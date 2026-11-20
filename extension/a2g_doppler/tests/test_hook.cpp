// SPDX-License-Identifier: AGPL-3.0-or-later
// Integration tests: a2g_doppler through chem::Channel, the JSON parser and
// the live extension path (ExtensionRegistry context, node map, Intermediate).
#include <gtest/gtest.h>

#include <chrono>
#include <cmath>
#include <iostream>
#include <thread>

#include "a2g_doppler.h"
#include "chem/aerpaw/node_characteristics.h"
#include "chem/channel/channel.h"
#include "chem/channel/intermediate.h"
#include "chem/common.h"
#include "chem/models/data_pool.hpp"
#include "chem/models/node_config.hpp"
#include "chem/models/usrp_channel.hpp"
#include "chem/node/node_entity.h"

using namespace chem;
using nlohmann::json;

namespace {
NodeConfig makeCfg(const std::string& id, double fs, double fc) {
    NodeConfig cfg;
    cfg.setId(id);
    cfg.setName(id);
    cfg.setNumChannels(1);
    cfg.setIpAddress("127.0.0.1");
    usrp::SampleRate sr;
    sr.setRxRate(fs);
    sr.setTxRate(fs);
    cfg.setSampleRate(sr);
    usrp::Channel ch;
    ch.setRxFreq(fc);
    ch.setTxFreq(fc);
    cfg.addChannel(ch);
    cfg.setAntennaPattern("dipole");
    cfg.setInputFormat("fc32");
    cfg.setOutputFormat("fc32");
    chem::aerpaw::NodeCharacteristics chars{};
    chars.source_power_dbfs = -30.0;
    cfg.setCharacteristics(chars);
    cfg.setDeviceGains(TX_GAIN, RX_GAIN);
    return cfg;
}

double toneFreq(const signal_v& y, size_t n, double fs) {
    std::complex<double> acc = 0;
    for (size_t i = 1; i < n; ++i)
        acc += std::complex<double>(y[i].real(), y[i].imag()) *
               std::conj(std::complex<double>(y[i - 1].real(), y[i - 1].imag()));
    return std::arg(acc) * fs / (2 * M_PI);
}

ChannelProcessParams params(size_t n, double fs, double fc, const NodeConfig* a,
                            const NodeConfig* b) {
    ChannelProcessParams p{};
    p.coeff = 1.0f;
    p.s_perBuff = n;
    p.sample_ratio = 1.0f;
    p.sample_rate_hz = fs;
    p.src_channel_index = 0;
    p.dest_channel_index = 0;
    p.propagationModel = PropagationModel::NONE;
    p.ground_coeff = -0.3f;
    p.src_config = a;
    p.dest_config = b;
    p.freq = fc;
    return p;
}
}  // namespace

TEST(A2GParse, ReadsPathsAndDefaults) {
    json j = json::parse(R"({
        "K_dB": 9.0, "shape": "GAUSS", "updateMs": 1.0, "taps": "LOS_ONLY",
        "phaseUpdate": "START_HOLD",
        "paths": [{"type":"LOS"},
                  {"type":"GROUND","reflCoeff":-0.5},
                  {"type":"SCATTERER","point":[1,2,3],"gainDb":-10,"K_dB":3},
                  {"type":"FIXED","delayUs":2.0,"dopplerHz":-30}]})");
    auto c = a2gParseModelConfig(j);
    ASSERT_EQ(c.paths.size(), 4u);
    EXPECT_DOUBLE_EQ(c.paths[0].kDb, 9.0);
    EXPECT_TRUE(a2g::isSpecular(c.paths[1].kDb));  // LOS_ONLY: others specular
    EXPECT_DOUBLE_EQ(c.paths[1].reflCoeff, -0.5);
    EXPECT_DOUBLE_EQ(c.paths[2].kDb, 3.0);
    EXPECT_DOUBLE_EQ(c.paths[2].point.z, 3.0);
    EXPECT_NEAR(c.paths[3].delayS, 2e-6, 1e-15);
    EXPECT_DOUBLE_EQ(c.updateS, 1e-3);
    EXPECT_EQ(c.phaseUpdate, a2g::PhaseUpdate::START_HOLD);
}

#ifdef CHEM_EXT_TV_CIR_HOOK
TEST(A2GHook, ChannelAppliesGeometricDopplerOnce) {
    const double fs = 100e3, fc = 3.4e9;
    auto a = makeCfg("gt", fs, fc), b = makeCfg("uav", fs, fc);
    Channel ch("gt", "uav", 1, 1);
    ch.updateNoiseType(NoiseType::NONE);
    ch.setShadowingSTD(0.0);
    ch.setDopplerEnabled(false);  // the extension owns the Doppler
    auto gt = std::make_shared<a2g::StaticTrajectory>(a2g::Vec3{0, 0, 10});
    auto uavT = std::make_shared<a2g::TableTrajectory>();
    for (double t = -1; t < 50; t += 0.1) uavT->add(t, {300 - 10 * t, 0, 60}, {-10, 0, 0});
    a2g::ModelConfig mc;
    mc.fc = fc;
    mc.paths[0].kDb = a2g::kInf;
    mc.paths[0].shape = a2g::Shape::NONE;
    auto st = A2GDopplerExtension::makeLinkState(mc, gt, uavT, false);
    ASSERT_TRUE(A2GDopplerExtension::installHook(ch, st));
    EXPECT_TRUE(ch.hasTimeVaryingCIR());

    const size_t n = 10000;
    signal_v x(n, chem::fc(0.01f, 0.0f)), y(n);
    auto p = params(n, fs, fc, &a, &b);
    for (int k = 0; k < 5; ++k) ch.processChannel(x, y, p);
    const double f = toneFreq(y, n, fs);
    const double fExp = a2g::geometricDopplerHz(gt->at(0.45), uavT->at(0.45), fc);
    EXPECT_NEAR(f, fExp, 0.5);
    EXPECT_GT(fExp, 100.0);

    // Remove the hook: no Doppler (core Doppler is off)
    ch.clearTimeVaryingCIR();
    ch.processChannel(x, y, p);
    EXPECT_NEAR(toneFreq(y, n, fs), 0.0, 1e-3);
}
#endif  // CHEM_EXT_TV_CIR_HOOK

TEST(A2GExtension, LivePathAttachesToSandboxChannel) {
    const double fc = 3.4e9;
    auto payload = std::make_shared<udpDataPool_t>(POOL_CAP, "t_payload");
    auto iq = std::make_shared<iqPool_t>(8, 4096, "t_iq");
    std::map<double, std::shared_ptr<Intermediate>> inter;
    inter[fc] = std::make_shared<Intermediate>(fc, payload, std::vector<iqPoolSPtr>{iq}, 100);
    inter[fc]->addSandboxNode("gt", 1);
    inter[fc]->addSandboxNode("uav", 1);

    std::map<std::string, Node> nodes;
    nodes.emplace("gt", Node(makeCfg("gt", 1e6, fc)));
    nodes.emplace("uav", Node(makeCfg("uav", 1e6, fc)));
    nodes["gt"].setLocation({35.72750947, -78.69595819, 10.0f});
    nodes["uav"].setLocation({35.7290, -78.6960, 60.0f});
    nodes["uav"].setVelocity({0.0f, -5.0f, 0.0f});

    A2GDopplerExtension ext;
    ext.setContext(nodes, inter);
    json cfg = {{"enabled", true}, {"mode", "hook"}, {"pollMs", 5},
                {"links", json::array({{{"src", "gt"}, {"dest", "uav"}}})}};
    // enable the core Doppler first: the extension must turn it off
    inter[fc]->updateDopplerEnabled("gt", "uav", true);
    ASSERT_TRUE(ext.onStart(cfg));
    std::this_thread::sleep_for(std::chrono::milliseconds(150));
    auto& chans = inter[fc]->getChannelList();
    auto it = chans.find({"gt", "uav"});
    ASSERT_NE(it, chans.end());
#ifdef CHEM_EXT_TV_CIR_HOOK
    EXPECT_TRUE(it->second.hasTimeVaryingCIR());
#endif
    EXPECT_EQ(it->second.getActiveExtension(), "a2g_doppler");
    EXPECT_FALSE(it->second.isDopplerEnabled());
    // the reverse link is not selected
    auto it2 = chans.find({"uav", "gt"});
    ASSERT_NE(it2, chans.end());
    EXPECT_NE(it2->second.getActiveExtension(), "a2g_doppler");
    auto s = ext.getStatus();
    EXPECT_EQ(s["links"].size(), 1u);
    ext.onStop();
    EXPECT_FALSE(it->second.isExtensionActive());
#ifdef CHEM_EXT_TV_CIR_HOOK
    EXPECT_FALSE(it->second.hasTimeVaryingCIR());
#endif
}

TEST(A2GExtension, TapsModePushesCIR) {
    const double fc = 3.4e9;
    auto payload = std::make_shared<udpDataPool_t>(POOL_CAP, "t_payload2");
    auto iq = std::make_shared<iqPool_t>(8, 4096, "t_iq2");
    std::map<double, std::shared_ptr<Intermediate>> inter;
    inter[fc] = std::make_shared<Intermediate>(fc, payload, std::vector<iqPoolSPtr>{iq}, 100);
    inter[fc]->addSandboxNode("gt", 1);
    inter[fc]->addSandboxNode("uav", 1);
    std::map<std::string, Node> nodes;
    nodes.emplace("gt", Node(makeCfg("gt", 1e6, fc)));
    nodes.emplace("uav", Node(makeCfg("uav", 1e6, fc)));
    nodes["gt"].setLocation({35.72750947, -78.69595819, 10.0f});
    nodes["uav"].setLocation({35.7290, -78.6960, 60.0f});
    A2GDopplerExtension ext;
    ext.setContext(nodes, inter);
    json cfg = {{"mode", "taps"}, {"tapUpdateMs", 1.0}, {"sampleRateHz", 1e6},
                {"paths", json::array({{{"type", "LOS"}},
                                       {{"type", "FIXED"}, {"delayUs", 3.0}, {"gainDb", -6}}})}};
    ASSERT_TRUE(ext.onStart(cfg));
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    auto& ch = inter[fc]->getChannelList().find({"gt", "uav"})->second;
    EXPECT_EQ(ch.getChTapCount({0, 0}), 4u);  // taps at 0 and 3 samples
    EXPECT_GT(ext.getStatus()["tapPushes"].get<uint64_t>(), 10u);
    ext.onStop();
}

// CHEM note: unset device gains / RF characteristics use NaN sentinels, and
// CHEM's std::isfinite() checks are removed by -ffast-math, so
// processChannel() then returns NaN. With explicit gains and characteristics
// (as on the twin) the output is finite.
TEST(A2GHook, CoreOutputFiniteWithCharacteristics) {
    const double fs = 100e3, fcHz = 3.4e9;
    auto a = makeCfg("gt", fs, fcHz), b = makeCfg("uav", fs, fcHz);
    Channel ch("gt", "uav", 1, 1);
    ch.updateNoiseType(NoiseType::AWGN);
    ch.setShadowingSTD(0.0);
    ch.updateDistance(100.0f);
    const size_t n = 1000;
    signal_v x(n, chem::fc(0.0316f, 0.0f)), y(n);
    auto p = params(n, fs, fcHz, &a, &b);
    p.propagationModel = PropagationModel::FREE_SPACE;
    ch.processChannel(x, y, p);
    // NaN-safe under -ffast-math: a NaN fails both comparisons
    EXPECT_LT(std::abs(y[0]), 1e3f);
    EXPECT_GT(std::abs(y[0]), 0.0f);
    EXPECT_NEAR(ch.getPlDb(), 83.08, 0.1);
}

