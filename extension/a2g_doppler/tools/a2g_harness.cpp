// SPDX-License-Identifier: AGPL-3.0-or-later
// a2g_harness: offline harness that pushes an IQ file through CHEM code.
//
//   a2g_harness <run.json>
//
// pipeline "channel":      chem::Channel::processChannel() frame by frame
//                          (fc32, exact CHEM link processing).
// pipeline "intermediate": chem::Intermediate with a real chem::emulator::
//                          Receiver (source), worker thread, UHD sc16
//                          converters and a capturing Transmitter
//                          (destination) - the CHEM data path without UDP.
// Every channelUpdateMs of stream time the harness does what
// Coordinator::configChannels() does: horizontal distance, heights,
// elevation, azimuth and (if enabled) the single link Doppler from the
// trajectories. The a2g_doppler extension runs in "hook" mode (per-sample
// time-varying CIR) or "taps" mode (CIR pushed every tapUpdateMs).
#include <boost/asio/io_service.hpp>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <mutex>
#include <thread>
#include <uhd/convert.hpp>

#include "a2g_doppler.h"
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
#include "spdlog/spdlog.h"

using namespace chem;
using nlohmann::json;
using namespace std::chrono;
static const bool kDbg = std::getenv("A2G_DEBUG") != nullptr;
#define DBG(msg) do { if (kDbg) std::cerr << "[harness] " << msg << std::endl; } while (0)

namespace {

NodeConfig makeCfg(const std::string& id, double fs, double fcHz, short txp, short rxp,
                   double srcPowerDbfs) {
    NodeConfig cfg;
    cfg.setId(id);
    cfg.setName(id);
    cfg.setNumChannels(1);
    cfg.setIpAddress("127.0.0.1");
    cfg.setTxPortNumber(txp);
    cfg.setRxPortNumber(rxp);
    usrp::SampleRate sr;
    sr.setRxRate(fs);
    sr.setTxRate(fs);
    cfg.setSampleRate(sr);
    usrp::Channel ch;
    ch.setRxFreq(fcHz);
    ch.setTxFreq(fcHz);
    cfg.addChannel(ch);
    cfg.setAntennaPattern("dipole");
    cfg.setInputFormat("fc32");
    cfg.setOutputFormat("fc32");
    // RF characteristics with a reference source power, as on the twin.
    // (Without them CHEM's NaN sentinel is not detected under -ffast-math.)
    chem::aerpaw::NodeCharacteristics chars{};
    chars.source_power_dbfs = srcPowerDbfs;
    cfg.setCharacteristics(chars);
    cfg.setDeviceGains(TX_GAIN, RX_GAIN);  // explicit: NaN sentinel + -ffast-math
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
        if (!t->loadCsv(j.get<std::string>()))
            throw std::runtime_error("bad trajectory " + j.get<std::string>());
        return t;
    }
    return std::make_shared<a2g::StaticTrajectory>(
        a2g::Vec3{j[0].get<double>(), j[1].get<double>(), j[2].get<double>()});
}

// What Coordinator::configChannels() computes, from ENU trajectories.
struct Geo {
    double dist, dAlt, elev, az, doppler, h1, h2;
};
Geo coordinatorGeometry(const a2g::KinState& a, const a2g::KinState& b, double fcHz) {
    Geo g{};
    const double e = b.pos.x - a.pos.x, n = b.pos.y - a.pos.y, u = b.pos.z - a.pos.z;
    g.dist = std::sqrt(e * e + n * n);  // CHEM: horizontal (haversine) distance
    g.dAlt = u;
    g.elev = (g.dist > 0.01) ? std::atan(u / g.dist) * 180 / PI : 0.0;
    g.az = std::atan2(e, n) * 180 / PI;
    g.h1 = a.pos.z;
    g.h2 = b.pos.z;
    const double r = std::sqrt(e * e + n * n + u * u);
    if (r > 0.01) {
        const double rr = ((b.vel.x - a.vel.x) * e + (b.vel.y - a.vel.y) * n +
                           (b.vel.z - a.vel.z) * u) / r;
        g.doppler = -(rr / SPEED_OF_LIGHT) * fcHz;
    }
    return g;
}

PropagationModel plFrom(const std::string& s) {
    if (s == "NONE") return PropagationModel::NONE;
    return PropagationModelFromString(s);
}

// Destination that captures processed signals (no UDP).
class CaptureTransmitter : public emulator::Transmitter {
   public:
    CaptureTransmitter(boost::asio::io_service& io, const NodeConfig& c, iqPoolSPtr pool)
        : Transmitter(io, c, pool), m_pool(pool) {
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
        std::vector<fc> out(n);
        const void* in = sig->getData(0);
        void* o = out.data();
        m_conv->conv(in, o, n);
        m_pool->release(sig->dataArray());
        std::lock_guard<std::mutex> g(m_mu);
        m_out.insert(m_out.end(), out.begin(), out.end());
        ++m_count;
        m_cv.notify_all();
    }
    bool waitCount(size_t c, int ms) {
        std::unique_lock<std::mutex> l(m_mu);
        return m_cv.wait_for(l, milliseconds(ms), [&] { return m_count >= c; });
    }
    std::vector<fc> take() {
        std::lock_guard<std::mutex> g(m_mu);
        return m_out;
    }

   private:
    iqPoolSPtr m_pool;
    uhd::convert::converter::sptr m_conv;
    std::mutex m_mu;
    std::condition_variable m_cv;
    std::vector<fc> m_out;
    size_t m_count{0};
};

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: a2g_harness run.json\n";
        return 2;
    }
    spdlog::set_level(spdlog::level::warn);
    json R = json::parse(std::ifstream(argv[1]));
    const double fs = R.at("fs").get<double>();
    const double fcHz = R.value("fc", 3.4e9);
    const size_t frame = R.value("frameSamples", 2000);
    const std::string pipeline = R.value("pipeline", "channel");
    const json C = R.value("chem", json::object());
    const json E = R.value("ext", json());
    const double updS = C.value("channelUpdateMs", 100.0) * 1e-3;
    const double inScale = std::pow(10.0, C.value("inputDbfs", -30.0) / 20.0);

    auto x = readCf32(R.at("in").get<std::string>());
    const size_t nFrames = x.size() / frame;
    auto trTx = loadTraj(R.at("trajTx"));
    auto trRx = loadTraj(R.at("trajRx"));

    NodeConfig ncTx = makeCfg("tx", fs, fcHz, 30000, 20000, C.value("inputDbfs", -30.0));
    NodeConfig ncRx = makeCfg("rx", fs, fcHz, 30001, 20001, C.value("inputDbfs", -30.0));

    // ---------------- set up the CHEM link --------------------------------
    std::unique_ptr<Channel> standalone;
    std::shared_ptr<Intermediate> inter;
    boost::asio::io_service io;
    iqPoolSPtr iqPool;
    std::shared_ptr<CaptureTransmitter> capTx;
    std::shared_ptr<emulator::Receiver> rxNode;
    Channel* ch = nullptr;

    const PropagationModel pl = plFrom(C.value("pl", std::string("FREE_SPACE")));
    if (pipeline == "intermediate") {
        auto payload = std::make_shared<udpDataPool_t>(POOL_CAP * 4, "h_payload");
        iqPool = std::make_shared<iqPool_t>(256, frame * 4 + HEADER_SIZE + 64, "h_iq");
        DBG("make Intermediate");
        inter = std::make_shared<Intermediate>(fcHz, payload, std::vector<iqPoolSPtr>{iqPool},
                                               C.value("maxLatencyMs", 10000));
        // node "tx" is the radio that transmits (CHEM source = Receiver)
        DBG("make Receiver");
        rxNode = std::make_shared<emulator::Receiver>(io, ncTx, iqPool);
        DBG("make Transmitter");
        capTx = std::make_shared<CaptureTransmitter>(io, ncRx, iqPool);
        DBG("attach dest");
        inter->AttachDestination(capTx);
        DBG("attach src");
        inter->AttachSource(rxNode);
        DBG("attached; ids " << rxNode->getId() << " " << capTx->getId());
        inter->setPathLossMode(pl);
        inter->updatePathLoss(pl, C.value("groundCoeff", -0.3f));
        auto it = inter->getChannelList().find({rxNode->getId(), capTx->getId()});
        if (it == inter->getChannelList().end()) {
            std::cerr << "channel not found\n";
            return 1;
        }
        ch = &it->second;
    } else {
        standalone = std::make_unique<Channel>("tx", "rx", 1, 1);
        ch = standalone.get();
    }
    ch->setShadowingSTD(C.value("shadowStd", 0.0));
    ch->updateNoiseType(C.value("noise", false) ? NoiseType::AWGN : NoiseType::NONE);
    ch->update3gppScenario(C.value("scenario3gpp", std::string("RMa")));
    ch->setFrequencyOffsetHz(C.value("freqOffsetHz", 0.0));
    const bool coreDoppler = C.value("doppler", true);
    ch->setDopplerEnabled(coreDoppler);
    if (C.contains("cir")) {
        signal_v taps;
        for (const auto& t : C["cir"]) taps.emplace_back(t[0].get<float>(), t[1].get<float>());
        ch->updateChTaps({0, 0}, taps);
    }

    // ---------------- extension ---------------------------------------------
    std::shared_ptr<A2GDopplerExtension::LinkState> lst;
    std::string extMode;
    double tapUpdS = 1e-3;
    if (E.is_object()) {
        a2g::ModelConfig mc = a2gParseModelConfig(E);
        mc.fc = fcHz;
        mc.fs = fs;
        lst = A2GDopplerExtension::makeLinkState(mc, trTx, trRx, false);
        extMode = E.value("mode", std::string("hook"));
        tapUpdS = E.value("tapUpdateMs", 1.0) * 1e-3;
        ch->setDopplerEnabled(false);  // the extension owns the Doppler
        ch->setActiveExtension("a2g_doppler", false);
        if (extMode == "hook") {
            if (!A2GDopplerExtension::installHook(*ch, lst)) {
                std::cerr << "hook not compiled in\n";
                return 1;
            }
        }
    }

    // ---------------- state log ----------------------------------------------
    std::ofstream slog;
    if (R.contains("stateOut")) {
        slog.open(R["stateOut"].get<std::string>());
        slog << "t,dist_h,elev,chemDopplerHz,plDb";
        if (lst) {
            for (size_t i = 0; i < lst->cfg.paths.size(); ++i)
                slog << ",f" << i << ",tau" << i << "_us,amp" << i << "_db,sig" << i;
        }
        slog << "\n";
    }
    const double logEvery = R.value("logStateEveryMs", 10.0) * 1e-3;

    // ---------------- run ----------------------------------------------------
    std::vector<fc> y;
    y.reserve(nFrames * frame);
    signal_v in(std::max<size_t>(frame, MAX_NUM_SAMPLES)), out(std::max<size_t>(frame, MAX_NUM_SAMPLES));
    ChannelProcessParams p{};
    p.coeff = 1.0f;
    p.s_perBuff = frame;
    p.sample_ratio = 1.0f;
    p.sample_rate_hz = fs;
    p.propagationModel = pl;
    p.ground_coeff = C.value("groundCoeff", -0.3f);
    p.src_config = &ncTx;
    p.dest_config = &ncRx;
    p.freq = fcHz;

    uhd::convert::converter::sptr toSc16;
    if (pipeline == "intermediate") {
        uhd::convert::id_type id;
        id.input_format = "fc32";
        id.num_inputs = 1;
        id.output_format = "sc16_item32_le";
        id.num_outputs = 1;
        toSc16 = uhd::convert::get_converter(id)();
        toSc16->set_scalar(32767.0);
    }

    double nextUpd = 0, nextTap = 0, nextLog = 0;
    const auto wall0 = steady_clock::now();
    for (size_t k = 0; k < nFrames; ++k) {
        const double t = static_cast<double>(k * frame) / fs;
        if (t >= nextUpd) {
            const Geo g = coordinatorGeometry(trTx->at(t), trRx->at(t), fcHz);
            ch->updateDistance(static_cast<float>(g.dist));
            ch->updateAltitude(static_cast<float>(std::abs(g.dAlt)));
            ch->updateHeights(static_cast<float>(g.h1), static_cast<float>(g.h2));
            ch->updateElevation(static_cast<float>(g.elev));
            ch->updateAzimuth(static_cast<float>(g.az));
            if (ch->isDopplerEnabled()) ch->setDopplerHz(g.doppler);
            nextUpd += updS;
        }
        if (lst && extMode == "taps" && t >= nextTap) {
            std::lock_guard<std::mutex> gl(lst->mu);
            auto taps = lst->port(0, 0, fs).sampledCIR(t, 256);
            ch->updateChTaps({0, 0}, signal_v(taps.begin(), taps.end()));
            while (nextTap <= t) nextTap += tapUpdS;
        }
        if (slog.is_open() && t >= nextLog) {
            const Geo g = coordinatorGeometry(trTx->at(t), trRx->at(t), fcHz);
            slog << t << "," << g.dist << "," << g.elev << ","
                 << (ch->isDopplerEnabled() ? ch->getDopplerHz() : 0.0) << "," << ch->getPlDb();
            if (lst) {
                std::lock_guard<std::mutex> gl(lst->mu);
                auto snap = lst->port(0, 0, fs).snapshot(t);
                for (auto& ps : snap.paths)
                    slog << "," << ps.geoDopplerHz << "," << ps.excessDelayS * 1e6 << ","
                         << 20 * std::log10(std::max(1e-12, std::abs(ps.amp))) << "," << ps.sigmaHz;
            }
            slog << "\n";
            nextLog += logEvery;
        }
        for (size_t i = 0; i < frame; ++i) in[i] = x[k * frame + i] * static_cast<float>(inScale);

        if (pipeline == "intermediate") {
            auto buf = iqPool->acquire();
            Header h{};
            h.size = frame * 4;
            h.sample_rate = static_cast<uint64_t>(fs);
            h.number_of_samples = static_cast<uint32_t>(frame);
            h.total_frags = 1;
            h.number_of_channels = 1;
            const int64_t now = duration_cast<nanoseconds>(steady_clock::now().time_since_epoch()).count();
            h.start = now;
            h.end = now + static_cast<int64_t>(frame / fs * 1e9);
            h.rt_tx_time = now;
            auto sig = std::make_unique<chem::Signal>(std::move(buf), h, fcHz, ncTx.getId());
            const void* ip = in.data();
            void* op = sig->getData(0);
            toSc16->conv(ip, op, frame);
            if (k == 0) DBG("first Add2Buff");
            rxNode->Add2Buff(std::move(sig));
            if (k == 0) DBG("first frame queued");
            if (!capTx->waitCount(k + 1, 5000)) {
                std::cerr << "timeout waiting for frame " << k << "\n";
                return 1;
            }
        } else {
            ch->processChannel(in, out, p);
            y.insert(y.end(), out.begin(), out.begin() + frame);
        }
    }
    if (pipeline == "intermediate") {
        y = capTx->take();
        inter->DetachSource(rxNode);
        inter->DetachDestination(capTx);
    }
    const double wall = duration<double>(steady_clock::now() - wall0).count();
    std::ofstream fo(R.at("out").get<std::string>(), std::ios::binary);
    fo.write(reinterpret_cast<const char*>(y.data()), y.size() * sizeof(fc));
    std::cout << "{\"frames\":" << nFrames << ",\"samples\":" << y.size()
              << ",\"wall_s\":" << wall << ",\"stream_s\":" << nFrames * frame / fs
              << ",\"buffers_ext\":" << (lst ? lst->buffers : 0) << "}" << std::endl;
    return 0;
}
