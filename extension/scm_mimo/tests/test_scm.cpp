// SPDX-License-Identifier: AGPL-3.0-or-later
// Unit tests of the scm_mimo model library and the a2g Doppler presets.
#include <gtest/gtest.h>

#include <cmath>
#include <numeric>

#include "a2g/model.h"
#include "scm/scm.h"

using namespace scm;

namespace {
std::shared_ptr<a2g::Trajectory> fixedAt(double e, double n, double u) {
    return std::make_shared<a2g::StaticTrajectory>(a2g::Vec3{e, n, u});
}
std::shared_ptr<a2g::Trajectory> radialFlight(double speed, double alt) {
    auto t = std::make_shared<a2g::TableTrajectory>();
    for (int i = 0; i <= 400; ++i) {
        const double tt = i * 0.5;
        t->add(tt, {100.0 + speed * tt, 0.0, alt}, {speed, 0.0, 0.0});
    }
    return t;
}
ScmConfig base(Doppler d) {
    ScmConfig c;
    c.fc = 3.4e9;
    c.fs = 1.92e6;
    c.doppler = d;
    c.updateS = 0.5e-3;
    c.seed = 7;
    return c;
}
}  // namespace

TEST(ScmDrop, PowersDelaysAngles) {
    std::mt19937_64 rng(1);
    for (auto sc : {Scenario::A2G, Scenario::SMA, Scenario::UMA, Scenario::UMI}) {
        ScmConfig c;
        c.scenario = sc;
        for (int k = 0; k < 50; ++k) {
            const Drop d = drawDrop(c, Doppler::PROPOSED, rng);
            ASSERT_EQ(d.tau.size(), 6u);
            EXPECT_NEAR(std::accumulate(d.P.begin(), d.P.end(), 0.0), 1.0, 1e-12);
            EXPECT_DOUBLE_EQ(d.tau[0], 0.0);
            for (size_t n = 1; n < d.tau.size(); ++n) EXPECT_GE(d.tau[n], d.tau[n - 1]);
            for (const auto& q : d.q) {  // unit-rms Doppler quantiles
                double ms = 0;
                for (double v : q) ms += v * v;
                EXPECT_NEAR(std::sqrt(ms / q.size()), 1.0, 1e-9);
            }
        }
    }
}

TEST(ScmLink, UnitMeanPower) {
    // average |h|^2 over drops and time is 1 per port pair
    double p = 0;
    int cnt = 0;
    for (int seed = 0; seed < 40; ++seed) {
        auto c = base(Doppler::SCM);
        c.seed = 100 + seed;
        ScmLink L(c, 2, 2, fixedAt(0, 0, 10), radialFlight(5, 30));
        for (uint64_t k = 0; k < 2000; k += 50) {
            const auto H = L.narrowbandAtPacket(k);
            for (auto& r : H)
                for (auto& v : r) {
                    p += std::norm(v);
                    ++cnt;
                }
        }
    }
    EXPECT_NEAR(p / cnt, 1.0, 0.15);
}

TEST(ScmLink, LosFollowsGeometricDoppler) {
    auto c = base(Doppler::PROPOSED);
    c.kDb = 300;  // LOS only
    c.kFromDoppler = false;
    ScmLink L(c, 1, 1, fixedAt(0, 0, 10), radialFlight(5, 30));
    L.narrowbandAtPacket(1000);
    const auto inf = L.info();
    // UAV flies away radially at ~5 m/s: Doppler ~ -5 cos(elev)/lambda
    const double lam = a2g::kC0 / 3.4e9;
    const double expect = -5.0 * std::cos(inf.elevDeg * a2g::kPi / 180) / lam;
    EXPECT_NEAR(inf.fLosHz, expect, 0.5);
}

TEST(ScmLink, JakesHasNoLosAndTdlDHasLine) {
    ScmLink J(base(Doppler::JAKES), 2, 2, fixedAt(0, 0, 10), radialFlight(5, 30));
    J.narrowbandAtPacket(10);
    EXPECT_LT(J.info().kDb, -100);
    ScmLink D(base(Doppler::TDL_D), 2, 2, fixedAt(0, 0, 10), radialFlight(5, 30));
    D.narrowbandAtPacket(10);
    EXPECT_NEAR(D.info().kDb, 13.3, 1e-9);
    EXPECT_NEAR(D.info().fLosHz, D.info().fmHz * std::cos(a2g::kPi / 4), 1e-6);
}

TEST(ScmLink, StreamMatchesCirAtPacketBoundaries) {
    // impulse at a packet boundary: the output equals the CIR of that packet
    auto c = base(Doppler::SCM);
    c.fs = 23.04e6;
    c.scenario = Scenario::SMA;
    ScmLink A(c, 2, 2, fixedAt(0, 0, 10), radialFlight(5, 30));
    ScmLink B(c, 2, 2, fixedAt(0, 0, 10), radialFlight(5, 30));
    const size_t Lp = A.packetLen();
    const uint64_t k = 3;
    std::vector<cf> x(Lp * 6, cf(0, 0)), y(x.size());
    x[k * Lp] = cf(1, 0);
    const Cir h = B.cirAtPacket(k);
    for (size_t u = 0; u < 2; ++u)
        for (size_t s = 0; s < 2; ++s) {
            A.process(u, s, x.data(), y.data(), x.size());
            EXPECT_NEAR(y[k * Lp].real(), h[u][s][0].real(), 2e-5);
            EXPECT_NEAR(y[k * Lp].imag(), h[u][s][0].imag(), 2e-5);
            for (size_t l = 1; l < h[u][s].size() && l < 4; ++l) {
                // later taps see one sample of interpolation (tiny change)
                EXPECT_NEAR(std::abs(y[k * Lp + l]), std::abs(h[u][s][l]), 1e-2);
            }
        }
}

TEST(ScmLink, SpatialCorrelationDependsOnSpacing) {
    // wider ground-array spacing -> lower port correlation -> better conditioned
    auto meanCond = [](double spG) {
        double acc = 0;
        int cnt = 0;
        for (int seed = 0; seed < 60; ++seed) {
            auto c = base(Doppler::SCM);
            c.scenario = Scenario::SMA;
            c.ground.spacingWl = spG;
            c.seed = 500 + seed;
            ScmLink L(c, 2, 2, fixedAt(0, 0, 10), radialFlight(5, 30));
            const auto sv = singularValues2(L.narrowbandAtPacket(20));
            acc += 20 * std::log10(sv[0] / std::max(1e-12, sv[1]));
            ++cnt;
        }
        return acc / cnt;
    };
    EXPECT_GT(meanCond(0.5), meanCond(10.0));
}

TEST(ScmLink, CrossPolLosIsRankTwo) {
    auto c = base(Doppler::PROPOSED);
    c.kDb = 300;  // LOS only
    c.kFromDoppler = false;
    c.ground.xpol = true;
    c.air.xpol = true;
    ScmLink L(c, 2, 2, fixedAt(0, 0, 10), radialFlight(5, 30));
    const auto H = L.narrowbandAtPacket(5);
    EXPECT_NEAR(std::abs(H[0][1]), 0.0, 1e-9);
    EXPECT_NEAR(std::abs(H[1][0]), 0.0, 1e-9);
    EXPECT_NEAR(std::abs(H[0][0] + H[1][1]), 0.0, 1e-6);  // diag(1, -1) * common phase
    const auto sv = singularValues2(H);
    EXPECT_NEAR(sv[0] / sv[1], 1.0, 1e-6);
}

TEST(ScmHelpers, CapacityIdentity) {
    std::vector<std::vector<cd>> I{{1, 0}, {0, 1}};
    EXPECT_NEAR(capacityBpsHz(I, 100.0), 2 * std::log2(1 + 50.0), 1e-9);
    const auto sv = singularValues2(I);
    EXPECT_NEAR(sv[0], 1.0, 1e-12);
    EXPECT_NEAR(sv[1], 1.0, 1e-12);
}

TEST(A2gPresets, TdlTablesAndJakes) {
    a2g::ModelConfig c;
    c.paths = {a2g::PathSpec{}};
    c.paths[0].type = a2g::PathType::LOS;
    c.dopplerModel = a2g::DopplerModel::TDL_A;
    c.tdlDelaySpreadS = 100e-9;
    a2g::applyDopplerPreset(c);
    EXPECT_EQ(c.paths.size(), 23u);
    EXPECT_TRUE(c.normalizeTotal);
    EXPECT_NEAR(c.paths[22].delayS, 965.86e-9, 1e-12);
    c.dopplerModel = a2g::DopplerModel::TDL_D;
    a2g::applyDopplerPreset(c);
    EXPECT_EQ(c.paths.size(), 14u);
    EXPECT_TRUE(a2g::isSpecular(c.paths[0].kDb));
    EXPECT_TRUE(c.paths[0].useAoa);
    a2g::ModelConfig j;
    j.paths.resize(2);
    j.paths[1].type = a2g::PathType::GROUND;
    j.dopplerModel = a2g::DopplerModel::JAKES;
    a2g::applyDopplerPreset(j);
    for (const auto& p : j.paths) {
        EXPECT_FALSE(p.geoPhase);
        EXPECT_EQ(p.shape, a2g::Shape::JAKES);
        EXPECT_LE(p.kDb, -a2g::kSpecularDb);
    }
}

TEST(A2gPresets, TdlDLosLineAtPoint7Fm) {
    a2g::ModelConfig c;
    c.fs = 1.92e6;
    c.dopplerModel = a2g::DopplerModel::TDL_D;
    a2g::applyDopplerPreset(c);
    auto tx = fixedAt(0, 0, 10);
    auto rx = radialFlight(5, 30);
    a2g::LinkChannel lc(c, tx, rx);
    std::vector<a2g::cf> x(4000, a2g::cf(1, 0)), y(4000);
    lc.process(x.data(), y.data(), x.size());
    const auto& pk = lc.lastPacket();
    const double fm = 5.0 / (a2g::kC0 / 3.4e9);
    EXPECT_NEAR(pk.paths[0].dopplerHz, fm * std::cos(a2g::kPi / 4), 0.5);
    EXPECT_NEAR(pk.paths[1].dopplerHz, 0.0, 1e-9);
}
