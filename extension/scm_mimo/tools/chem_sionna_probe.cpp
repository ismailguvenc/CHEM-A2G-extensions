// SPDX-License-Identifier: AGPL-3.0-or-later
// chem_sionna_probe: runs CHEM's own "sionna" extension (SionnaExtension, the REST client that
// the deployed CHEM v1.2.4 also ships) inside a local CHEM build, against a local Sionna RT
// server (AERPAW-DT-SIONNA-EXTENSION, aerpaw Lake Wheeler scene). Two sandbox nodes (LW1, UAV)
// live in a real chem::Intermediate and a CHEM node map; the UAV walks along a position list.
// At each position the probe waits for the extension to compute paths and push the CIR, then
// sends an impulse through chem::Channel::processChannel() and logs what CHEM applies.
//
//   chem_sionna_probe run.json
// run.json: {"serverUrl": "http://127.0.0.1:8000", "sceneConfig": "aerpaw",
//   "sceneOffset": {"x":-3.7,"y":3.4,"z":0}, "referenceOrigin": {"lat":..,"lon":..,"alt":112},
//   "tx": {"lat":..,"lon":..,"alt":10}, "positions": "uav.csv" (t,lat,lon,alt AGL),
//   "dwellMs": 1500, "out": "probe.csv", "maxDepth": 3, "numSamples": 100000}
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <thread>

#include "chem/aerpaw/node_characteristics.h"
#include "chem/channel/channel.h"
#include "chem/channel/intermediate.h"
#include "chem/common.h"
#include "chem/models/data_pool.hpp"
#include "chem/models/node_config.hpp"
#include "chem/models/usrp_channel.hpp"
#include "chem/node/node_entity.h"
#include "sionna_client.h"
#include "spdlog/spdlog.h"

using namespace chem;
using nlohmann::json;

static NodeConfig makeCfg(const std::string& id, double fs, double fcHz) {
    NodeConfig cfg;
    cfg.setId(id);
    cfg.setName(id);
    cfg.setNumChannels(1);
    usrp::SampleRate sr;
    sr.setRxRate(fs);
    sr.setTxRate(fs);
    cfg.setSampleRate(sr);
    usrp::Channel ch;
    ch.setRxFreq(fcHz);
    ch.setTxFreq(fcHz);
    cfg.addChannel(ch);
    cfg.setAntennaPattern("isotropic");
    cfg.setInputFormat("fc32");
    cfg.setOutputFormat("fc32");
    chem::aerpaw::NodeCharacteristics chars{};
    chars.source_power_dbfs = -30.0;
    cfg.setCharacteristics(chars);
    cfg.setDeviceGains(0.0, 0.0);
    return cfg;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: chem_sionna_probe run.json\n";
        return 2;
    }
    spdlog::set_level(spdlog::level::warn);
    const json R = json::parse(std::ifstream(argv[1]));
    const double fcHz = R.value("fc", 3.4e9), fs = 1e6;  // sandbox nodes: the extension uses 1 MS/s
    auto payload = std::make_shared<udpDataPool_t>(POOL_CAP, "p_payload");
    auto iq = std::make_shared<iqPool_t>(8, 4096, "p_iq");
    std::map<double, std::shared_ptr<Intermediate>> inter;
    inter[fcHz] = std::make_shared<Intermediate>(fcHz, payload, std::vector<iqPoolSPtr>{iq}, 100);
    inter[fcHz]->addSandboxNode("LW1", 1);
    inter[fcHz]->addSandboxNode("UAV", 1);
    std::map<std::string, Node> nodes;
    auto cTx = makeCfg("LW1", fs, fcHz), cRx = makeCfg("UAV", fs, fcHz);
    nodes.emplace("LW1", Node(cTx));
    nodes.emplace("UAV", Node(cRx));
    const auto& tx = R.at("tx");
    nodes["LW1"].setLocation({tx["lat"].get<double>(), tx["lon"].get<double>(), tx["alt"].get<float>()});

    std::vector<std::array<double, 4>> pos;  // t, lat, lon, alt
    {
        std::ifstream f(R.at("positions").get<std::string>());
        std::string line;
        while (std::getline(f, line)) {
            if (line.empty() || !(std::isdigit(static_cast<unsigned char>(line[0])) || line[0] == '-')) continue;
            std::replace(line.begin(), line.end(), ',', ' ');
            std::istringstream ss(line);
            std::array<double, 4> p{};
            ss >> p[0] >> p[1] >> p[2] >> p[3];
            pos.push_back(p);
        }
    }
    nodes["UAV"].setLocation({pos[0][1], pos[0][2], static_cast<float>(pos[0][3])});

    SionnaExtension ext;
    ext.setContext(nodes, inter);
    json cfg = {{"serverUrl", R.value("serverUrl", std::string("http://127.0.0.1:8000"))},
                {"sceneConfig", R.value("sceneConfig", std::string("aerpaw"))},
                {"updateRateMs", R.value("updateRateMs", 300)},
                {"maxDepth", R.value("maxDepth", 3)},
                {"numSamples", R.value("numSamples", 100000)},
                {"referenceOrigin", R.at("referenceOrigin")},
                {"sceneOffset", R.at("sceneOffset")}};
    const auto t0 = std::chrono::steady_clock::now();
    if (!ext.onStart(cfg)) {
        std::cerr << "sionna extension did not start\n";
        return 1;
    }
    std::cerr << "scene " << ext.getSceneId() << " after "
              << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s\n";

    auto& chans = inter[fcHz]->getChannelList();
    auto it = chans.find({"LW1", "UAV"});
    if (it == chans.end()) {
        std::cerr << "no LW1->UAV channel\n";
        return 1;
    }
    Channel& ch = it->second;
    ch.updateNoiseType(NoiseType::NONE);
    ch.setShadowingSTD(0.0);
    ch.setDopplerEnabled(false);
    std::ofstream out(R.at("out").get<std::string>());
    out << "t,lat,lon,alt,active,taps,energy_db,peak_tap,h0_re,h0_im,h1_re,h1_im,h2_re,h2_im,pl_db,status_ms\n";
    const int dwell = R.value("dwellMs", 1500);
    const size_t n = 64;
    for (size_t k = 0; k < pos.size(); k += R.value("every", 1)) {
        const auto& p = pos[k];
        nodes["UAV"].setLocation({p[1], p[2], static_cast<float>(p[3])});
        std::this_thread::sleep_for(std::chrono::milliseconds(dwell));
        signal_v x(n, fc(0, 0)), y(n, fc(0, 0));
        x[0] = fc(0.0316f, 0.0f);  // -30 dBFS impulse = source power reference
        ChannelProcessParams prm{};
        prm.coeff = 1.0f;
        prm.s_perBuff = n;
        prm.sample_ratio = 1.0f;
        prm.sample_rate_hz = fs;
        prm.propagationModel = PropagationModel::NONE;
        prm.src_config = &cTx;
        prm.dest_config = &cRx;
        prm.freq = fcHz;
        ch.processChannel(x, y, prm);
        double e = 0, pk = 0;
        size_t ipk = 0;
        for (size_t i = 0; i < n; ++i) {
            const double v = std::norm(y[i]) / (0.0316 * 0.0316);
            e += v;
            if (v > pk) { pk = v; ipk = i; }
        }
        const auto st = ext.getStatus();
        out << p[0] << "," << p[1] << "," << p[2] << "," << p[3] << "," << ch.isExtensionActive() << ","
            << ch.getChTapCount({0, 0}) << "," << 10 * std::log10(std::max(e, 1e-30)) << "," << ipk;
        for (int i = 0; i < 3; ++i) out << "," << y[i].real() / 0.0316 << "," << y[i].imag() / 0.0316;
        out << "," << ch.getPlDb() << "," << st.value("avgCirResponseMs", -1.0) << "\n";
        out.flush();
        std::cerr << "t=" << p[0] << " taps=" << ch.getChTapCount({0, 0}) << " E=" << 10 * std::log10(std::max(e, 1e-30))
                  << " dB\n";
    }
    std::cout << ext.getStatus().dump() << std::endl;
    ext.onStop();
    return 0;
}
