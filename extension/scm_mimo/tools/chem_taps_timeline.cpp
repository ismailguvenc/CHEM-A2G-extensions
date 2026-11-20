// SPDX-License-Identifier: AGPL-3.0-or-later
// chem_taps_timeline: write a CIR timeline from the a2g_doppler or scm_mimo
// model, for the live-twin driver (pushes taps with the CHEM API set_cir
// per port pair) and for analysis.
//
//   chem_taps_timeline run.json
//
// run.json:
//   {"model": "a2g" | "scm", "config": {...extension JSON block...},
//    "fc": 3.4e9, "fs": 11.52e6, "rateHz": 1000, "durationS": 60,
//    "nTx": 2, "nRx": 2, "maxTaps": 8,
//    "trajTx": [e,n,u] | "tx.csv", "trajRx": [e,n,u] | "rx.csv",
//    "out": "timeline.bin", "infoOut": "timeline_info.csv"}
// out: one JSON header line, then float32 complex [K][nRx][nTx][maxTaps].
// The a2g model is SISO (nTx = nRx = 1).
#include <cmath>
#include <fstream>
#include <iostream>

#include "a2g/model.h"
#include "chem/nlohmann/json.hpp"
#include "scm/scm.h"

using nlohmann::json;

// parsers of the extensions, re-declared here without CHEM dependencies
namespace chem {
a2g::ModelConfig a2gParseModelConfig(const json& j, const a2g::ModelConfig& base);
scm::ScmConfig scmParseConfig(const json& j, const scm::ScmConfig& base);
}  // namespace chem

static std::shared_ptr<a2g::Trajectory> loadTraj(const json& j) {
    if (j.is_string()) {
        auto t = std::make_shared<a2g::TableTrajectory>();
        if (!t->loadCsv(j.get<std::string>())) throw std::runtime_error("bad trajectory");
        return t;
    }
    return std::make_shared<a2g::StaticTrajectory>(
        a2g::Vec3{j[0].get<double>(), j[1].get<double>(), j[2].get<double>()});
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: chem_taps_timeline run.json\n";
        return 2;
    }
    const json R = json::parse(std::ifstream(argv[1]));
    const std::string model = R.value("model", std::string("scm"));
    const double fc = R.value("fc", 3.4e9), fs = R.value("fs", 11.52e6);
    const double rate = R.value("rateHz", 1000.0), dur = R.value("durationS", 60.0);
    const int L = R.value("maxTaps", 8);
    auto tx = loadTraj(R.at("trajTx"));
    auto rx = loadTraj(R.at("trajRx"));
    const uint64_t K = static_cast<uint64_t>(std::llround(dur * rate));
    std::ofstream fo(R.at("out").get<std::string>(), std::ios::binary);
    std::ofstream info;
    if (R.contains("infoOut")) info.open(R["infoOut"].get<std::string>());
    std::vector<float> buf;

    if (model == "a2g") {
        a2g::ModelConfig c = chem::a2gParseModelConfig(R.value("config", json::object()), {});
        c.fc = fc;
        c.fs = fs;
        a2g::LinkChannel lc(c, tx, rx);
        json hdr = {{"model", "a2g"}, {"dopplerModel", a2g::toString(c.dopplerModel)}, {"fc", fc},
                    {"fs", fs}, {"rateHz", rate}, {"K", K}, {"nRx", 1}, {"nTx", 1}, {"L", L},
                    {"paths", c.paths.size()}};
        fo << hdr.dump() << "\n";
        if (info.is_open()) info << "t,alt,speed,f0,sigma0\n";
        for (uint64_t k = 0; k < K; ++k) {
            const double t = k / rate;
            auto taps = lc.sampledCIR(t, static_cast<size_t>(L));
            buf.assign(2 * L, 0.f);
            for (size_t l = 0; l < taps.size() && l < static_cast<size_t>(L); ++l) {
                buf[2 * l] = taps[l].real();
                buf[2 * l + 1] = taps[l].imag();
            }
            fo.write(reinterpret_cast<const char*>(buf.data()), buf.size() * sizeof(float));
            if (info.is_open() && k % static_cast<uint64_t>(std::max(1.0, rate / 100)) == 0) {
                const auto s = lc.snapshot(t);
                info << t << "," << s.uavAltM << "," << s.uavSpeed << ","
                     << (s.paths.empty() ? 0 : s.paths[0].geoDopplerHz) << ","
                     << (s.paths.empty() ? 0 : s.paths[0].sigmaHz) << "\n";
            }
        }
    } else {
        scm::ScmConfig c = chem::scmParseConfig(R.value("config", json::object()), {});
        c.fc = fc;
        c.fs = fs;
        c.updateS = 1.0 / rate;
        const int nTx = R.value("nTx", 2), nRx = R.value("nRx", 2);
        scm::ScmLink lk(c, nTx, nRx, tx, rx);
        lk.setTime(0.0);
        json hdr = {{"model", "scm"}, {"scenario", scm::toString(c.scenario)},
                    {"doppler", scm::toString(c.doppler)}, {"fc", fc}, {"fs", fs},
                    {"rateHz", rate}, {"K", K}, {"nRx", nRx}, {"nTx", nTx}, {"L", L},
                    {"groundSpacingWl", c.ground.spacingWl}, {"airSpacingWl", c.air.spacingWl}};
        fo << hdr.dump() << "\n";
        if (info.is_open()) info << "t,dist,elev,speed,fLos,sigma,fm,K_dB,drops,cap10,cap20,cond_dB\n";
        for (uint64_t k = 0; k < K; ++k) {
            const scm::Cir h = lk.cirAtPacket(k);
            buf.assign(static_cast<size_t>(nRx * nTx * L * 2), 0.f);
            for (int u = 0; u < nRx; ++u)
                for (int s = 0; s < nTx; ++s)
                    for (int l = 0; l < L && l < static_cast<int>(h[u][s].size()); ++l) {
                        const size_t o = static_cast<size_t>(((u * nTx + s) * L + l) * 2);
                        buf[o] = static_cast<float>(h[u][s][l].real());
                        buf[o + 1] = static_cast<float>(h[u][s][l].imag());
                    }
            fo.write(reinterpret_cast<const char*>(buf.data()), buf.size() * sizeof(float));
            if (info.is_open() && k % static_cast<uint64_t>(std::max(1.0, rate / 100)) == 0) {
                std::vector<std::vector<scm::cd>> H(nRx, std::vector<scm::cd>(nTx));
                for (int u = 0; u < nRx; ++u)
                    for (int s = 0; s < nTx; ++s)
                        for (const auto& v : h[u][s]) H[u][s] += v;
                const auto inf = lk.info();
                double cond = 0;
                if (nTx == 2) {
                    const auto sv = scm::singularValues2(H);
                    cond = 20 * std::log10(sv[0] / std::max(1e-12, sv[1]));
                }
                info << k / rate << "," << inf.dist << "," << inf.elevDeg << "," << inf.speed << ","
                     << inf.fLosHz << "," << inf.sigmaHz << "," << inf.fmHz << "," << inf.kDb << ","
                     << inf.nDrops << "," << scm::capacityBpsHz(H, 10.0) << ","
                     << scm::capacityBpsHz(H, 100.0) << "," << cond << "\n";
            }
        }
    }
    std::cout << "{\"packets\":" << K << "}" << std::endl;
    return 0;
}
