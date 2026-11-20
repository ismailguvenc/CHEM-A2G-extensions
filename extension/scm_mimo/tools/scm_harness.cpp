// SPDX-License-Identifier: AGPL-3.0-or-later
// scm_harness: offline MIMO harness. nTx IQ streams -> CHEM -> nRx IQ streams.
//
//   scm_harness <run.json>
//
// pipeline "intermediate" (default): chem::Intermediate with a real
//   chem::emulator::Receiver (source node with nTx channels), worker thread,
//   UHD sc16 converters, and a capturing Transmitter (destination node with
//   nRx channels): the CHEM data path without UDP, including the core port
//   loop (with CHEM-EXT-MIMO: sum over source ports).
// The channel is either the scm_mimo extension ("ext": {...}, hook mode),
// or CHEM's own static per-port model ("chem": {"coeff": [[...]]} real
// matrix [src][dest], as set_mimo_matrix on the twin).
//
// run.json: {"fs":..., "fc":..., "frameSamples": 2000, "nTx": 2, "nRx": 2,
//   "in": ["tx0.cf32","tx1.cf32"], "out": ["rx0.cf32","rx1.cf32"],
//   "trajTx": [e,n,u] | "tx.csv", "trajRx": ..., "chem": {...}, "ext": {...}}
#include <boost/asio/io_service.hpp>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <fstream>
#include <iostream>
#include <mutex>
#include <uhd/convert.hpp>

#include "chem/aerpaw/node_characteristics.h"
#include "chem/channel/channel.h"
#include "chem/channel/intermediate.h"
#include "chem/common.h"
#include "chem/models/data_pool.hpp"
#include "chem/models/node_config.hpp"
#include "chem/models/signal.hpp"
#include "chem/models/usrp_channel.hpp"
#include "chem/node/receiver.h"
#include "chem/node/transmitter.h"
#include "scm_mimo.h"
#include "spdlog/spdlog.h"

using namespace chem;
using nlohmann::json;
using namespace std::chrono;

namespace {

NodeConfig makeCfg(const std::string& id, int nch, double fs, double fcHz, short txp, short rxp,
                   double srcPowerDbfs) {
    NodeConfig cfg;
    cfg.setId(id);
    cfg.setName(id);
    cfg.setNumChannels(static_cast<uint8_t>(nch));
    cfg.setIpAddress("127.0.0.1");
    cfg.setTxPortNumber(txp);
    cfg.setRxPortNumber(rxp);
    usrp::SampleRate sr;
    sr.setRxRate(fs);
    sr.setTxRate(fs);
    cfg.setSampleRate(sr);
    for (int c = 0; c < nch; ++c) {
        usrp::Channel ch;
        ch.setRxFreq(fcHz);
        ch.setTxFreq(fcHz);
        cfg.addChannel(ch);
    }
    cfg.setAntennaPattern("dipole");
    cfg.setInputFormat("fc32");
    cfg.setOutputFormat("fc32");
    chem::aerpaw::NodeCharacteristics chars{};
    chars.source_power_dbfs = srcPowerDbfs;
    cfg.setCharacteristics(chars);
    cfg.setDeviceGains(TX_GAIN, RX_GAIN);
    return cfg;
}

std::vector<fc> readCf32(const std::string& p) {
    std::ifstream f(p, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + p);
    f.seekg(0, std::ios::end);
    const size_t n = static_cast<size_t>(f.tellg()) / sizeof(fc);
    f.seekg(0);
    std::vector<fc> v(n);
    f.read(reinterpret_cast<char*>(v.data()), n * sizeof(fc));
    return v;
}

std::shared_ptr<a2g::Trajectory> loadTraj(const json& j) {
    if (j.is_string()) {
        auto t = std::make_shared<a2g::TableTrajectory>();
        if (!t->loadCsv(j.get<std::string>())) throw std::runtime_error("bad trajectory");
        return t;
    }
    return std::make_shared<a2g::StaticTrajectory>(
        a2g::Vec3{j[0].get<double>(), j[1].get<double>(), j[2].get<double>()});
}

class CaptureTransmitter : public emulator::Transmitter {
   public:
    CaptureTransmitter(boost::asio::io_service& io, const NodeConfig& c, iqPoolSPtr pool, int nch)
        : Transmitter(io, c, pool), m_pool(pool), m_out(nch) {
        uhd::convert::id_type id;
        id.input_format = "sc16_item32_le";
        id.num_inputs = 1;
        id.output_format = "fc32";
        id.num_outputs = 1;
        m_conv = uhd::convert::get_converter(id)();
        m_conv->set_scalar(1.0 / 32767.0);
    }
    void Add2Buff(std::unique_ptr<chem::Signal> sig) override {
        const size_t n = sig->getHeader()->number_of_samples;
        std::vector<std::vector<fc>> o(m_out.size(), std::vector<fc>(n));
        for (size_t c = 0; c < m_out.size(); ++c) {
            const void* in = sig->getData(c);
            void* op = o[c].data();
            m_conv->conv(in, op, n);
        }
        m_pool->release(sig->dataArray());
        std::lock_guard<std::mutex> g(m_mu);
        for (size_t c = 0; c < m_out.size(); ++c) m_out[c].insert(m_out[c].end(), o[c].begin(), o[c].end());
        ++m_count;
        m_cv.notify_all();
    }
    bool waitCount(size_t c, int ms) {
        std::unique_lock<std::mutex> l(m_mu);
        return m_cv.wait_for(l, milliseconds(ms), [&] { return m_count >= c; });
    }
    std::vector<std::vector<fc>> take() {
        std::lock_guard<std::mutex> g(m_mu);
        return m_out;
    }

   private:
    iqPoolSPtr m_pool;
    uhd::convert::converter::sptr m_conv;
    std::mutex m_mu;
    std::condition_variable m_cv;
    std::vector<std::vector<fc>> m_out;
    size_t m_count{0};
};

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: scm_harness run.json\n";
        return 2;
    }
    spdlog::set_level(spdlog::level::warn);
    json R = json::parse(std::ifstream(argv[1]));
    const double fs = R.at("fs").get<double>();
    const double fcHz = R.value("fc", 3.4e9);
    const size_t frame = R.value("frameSamples", 2000);
    const int nTx = R.value("nTx", 2), nRx = R.value("nRx", 2);
    const json C = R.value("chem", json::object());
    const json E = R.value("ext", json());
    const double inScale = std::pow(10.0, C.value("inputDbfs", -30.0) / 20.0);

    std::vector<std::vector<fc>> x;
    for (const auto& f : R.at("in")) x.push_back(readCf32(f.get<std::string>()));
    if (static_cast<int>(x.size()) != nTx) throw std::runtime_error("need nTx input files");
    size_t nS = x[0].size();
    for (auto& v : x) nS = std::min(nS, v.size());
    const size_t nFrames = nS / frame;
    auto trTx = loadTraj(R.at("trajTx"));
    auto trRx = loadTraj(R.at("trajRx"));

    NodeConfig ncTx = makeCfg("tx", nTx, fs, fcHz, 30000, 20000, C.value("inputDbfs", -30.0));
    NodeConfig ncRx = makeCfg("rx", nRx, fs, fcHz, 30001, 20001, C.value("inputDbfs", -30.0));

    boost::asio::io_service io;
    auto payload = std::make_shared<udpDataPool_t>(POOL_CAP * 4, "h_payload");
    auto iqPool = std::make_shared<iqPool_t>(256, frame * 4 * std::max(nTx, nRx) + HEADER_SIZE + 64, "h_iq");
    auto inter = std::make_shared<Intermediate>(fcHz, payload, std::vector<iqPoolSPtr>{iqPool},
                                                C.value("maxLatencyMs", 10000));
    auto rxNode = std::make_shared<emulator::Receiver>(io, ncTx, iqPool);
    auto capTx = std::make_shared<CaptureTransmitter>(io, ncRx, iqPool, nRx);
    inter->AttachDestination(capTx);
    inter->AttachSource(rxNode);
    const PropagationModel pl = C.value("pl", std::string("NONE")) == "NONE"
                                    ? PropagationModel::NONE
                                    : PropagationModelFromString(C.value("pl", std::string("FREE_SPACE")));
    inter->setPathLossMode(pl);
    inter->updatePathLoss(pl, -0.3f);
    auto it = inter->getChannelList().find({rxNode->getId(), capTx->getId()});
    if (it == inter->getChannelList().end()) {
        std::cerr << "channel not found\n";
        return 1;
    }
    Channel* ch = &it->second;
    ch->setShadowingSTD(0.0);
    ch->updateNoiseType(C.value("noise", false) ? NoiseType::AWGN : NoiseType::NONE);
    ch->setDopplerEnabled(false);
    // CHEM's static per-port model: real coefficient matrix [src][dest]
    if (C.contains("coeff")) {
        const auto& M = C["coeff"];
        for (int s = 0; s < nTx; ++s)
            for (int d = 0; d < nRx; ++d)
                ch->updateChCoeff(chId{static_cast<size_t>(d), static_cast<size_t>(s)},
                                  M[s][d].get<float>());
    }
    std::shared_ptr<ScmMimoExtension::LinkState> lst;
    if (E.is_object()) {
        scm::ScmConfig sc = scmParseConfig(E);
        sc.fc = fcHz;
        sc.fs = fs;
        lst = ScmMimoExtension::makeLinkState(sc, nTx, nRx, trTx, trRx, false);
        lst->link->setTime(0.0);
        ch->setActiveExtension("scm_mimo", false);
        if (!ScmMimoExtension::installHook(*ch, lst)) {
            std::cerr << "hook not compiled in\n";
            return 1;
        }
    }
    // geometry for CHEM's own path loss (fixed at the start; PL is not the subject here)
    {
        const auto a = trTx->at(0).pos, b = trRx->at(0).pos;
        ch->updateDistance(static_cast<float>(std::hypot(b.x - a.x, b.y - a.y)));
        ch->updateHeights(static_cast<float>(a.z), static_cast<float>(b.z));
    }

    uhd::convert::id_type id;
    id.input_format = "fc32";
    id.num_inputs = 1;
    id.output_format = "sc16_item32_le";
    id.num_outputs = 1;
    auto toSc16 = uhd::convert::get_converter(id)();
    toSc16->set_scalar(32767.0);

    std::vector<fc> in(frame);
    const auto wall0 = steady_clock::now();
    for (size_t k = 0; k < nFrames; ++k) {
        auto buf = iqPool->acquire();
        Header h{};
        h.size = static_cast<uint32_t>(frame * 4 * nTx);
        h.sample_rate = static_cast<uint64_t>(fs);
        h.number_of_samples = static_cast<uint32_t>(frame);
        h.total_frags = 1;
        h.number_of_channels = static_cast<uint8_t>(nTx);
        const int64_t now = duration_cast<nanoseconds>(steady_clock::now().time_since_epoch()).count();
        h.start = now;
        h.end = now + static_cast<int64_t>(frame / fs * 1e9);
        h.rt_tx_time = now;
        auto sig = std::make_unique<chem::Signal>(std::move(buf), h, fcHz, ncTx.getId());
        for (int c = 0; c < nTx; ++c) {
            for (size_t i = 0; i < frame; ++i) in[i] = x[c][k * frame + i] * static_cast<float>(inScale);
            const void* ip = in.data();
            void* op = sig->getData(c);
            toSc16->conv(ip, op, frame);
        }
        rxNode->Add2Buff(std::move(sig));
        if (!capTx->waitCount(k + 1, 5000)) {
            std::cerr << "timeout waiting for frame " << k << "\n";
            return 1;
        }
    }
    auto y = capTx->take();
    inter->DetachSource(rxNode);
    inter->DetachDestination(capTx);
    const double wall = duration<double>(steady_clock::now() - wall0).count();
    for (int c = 0; c < nRx; ++c) {
        std::ofstream fo(R.at("out")[c].get<std::string>(), std::ios::binary);
        fo.write(reinterpret_cast<const char*>(y[c].data()), y[c].size() * sizeof(fc));
    }
    std::cout << "{\"frames\":" << nFrames << ",\"samples\":" << y[0].size() << ",\"wall_s\":" << wall
              << ",\"stream_s\":" << nFrames * frame / fs
              << ",\"buffers_ext\":" << (lst ? lst->buffers : 0) << "}" << std::endl;
    return 0;
}
