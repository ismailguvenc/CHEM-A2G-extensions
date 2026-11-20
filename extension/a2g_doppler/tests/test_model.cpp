// SPDX-License-Identifier: AGPL-3.0-or-later
// Unit tests of the a2g_doppler model (no CHEM dependency).
#include <gtest/gtest.h>

#include <cmath>
#include <numeric>

#include "a2g/model.h"

using namespace a2g;

namespace {

constexpr double FC = 3.4e9;
const double LAM = kC0 / FC;

std::shared_ptr<Trajectory> line(Vec3 p0, Vec3 v, double T = 200.0) {
    auto t = std::make_shared<TableTrajectory>();
    for (double s = -1.0; s <= T; s += 0.1) t->add(s, p0 + s * v, v);
    return t;
}

// Unwrapped phase of a complex sequence
std::vector<double> uphase(const std::vector<cf>& y) {
    std::vector<double> p(y.size());
    double prev = 0, off = 0;
    for (size_t i = 0; i < y.size(); ++i) {
        double a = std::arg(cd(y[i].real(), y[i].imag()));
        if (i > 0) {
            double d = a - prev;
            if (d > kPi) off -= 2 * kPi;
            if (d < -kPi) off += 2 * kPi;
        }
        prev = a;
        p[i] = a + off;
    }
    return p;
}

double slopeHz(const std::vector<double>& ph, double fs, size_t a, size_t b) {
    // least-squares slope of phase -> frequency
    double n = b - a, sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (size_t i = a; i < b; ++i) {
        double x = (i - a) / fs;
        sx += x; sy += ph[i]; sxx += x * x; sxy += x * ph[i];
    }
    return (n * sxy - sx * sy) / (n * sxx - sx * sx) / (2 * kPi);
}

ModelConfig losOnly(double fs) {
    ModelConfig c;
    c.fc = FC;
    c.fs = fs;
    PathSpec los;
    los.type = PathType::LOS;
    los.kDb = kInf;
    los.shape = Shape::NONE;
    c.paths = {los};
    return c;
}

}  // namespace

TEST(Trajectory, HermiteExactForLinearMotion) {
    auto t = line({0, 0, 30}, {5, 1, 0});
    auto k = t->at(12.345);
    EXPECT_NEAR(k.pos.x, 5 * 12.345, 1e-9);
    EXPECT_NEAR(k.pos.y, 12.345, 1e-9);
    EXPECT_NEAR(k.vel.x, 5, 1e-9);
}

TEST(Trajectory, LiveTrajectoryEstimatesVelocity) {
    LiveTrajectory lt(0.0);
    for (int i = 0; i < 10; ++i) lt.push(i * 0.1, {i * 0.5, 0, 50}, {}, false);
    auto k = lt.at(0.55);
    EXPECT_NEAR(k.vel.x, 5.0, 1e-6);
    EXPECT_NEAR(k.pos.x, 2.75, 1e-6);
}

TEST(LosDoppler, ConstantVelocityToneMatchesGeometry) {
    const double fs = 20e3;
    auto gt = std::make_shared<StaticTrajectory>(Vec3{0, 0, 10});
    auto uav = line({200, 0, 30}, {-5, 0, 0});  // approaching
    LinkChannel lc(losOnly(fs), gt, uav);
    std::vector<cf> x(fs * 2, cf(1, 0)), y(x.size());
    lc.process(x.data(), y.data(), x.size());
    auto ph = uphase(y);
    const double f = slopeHz(ph, fs, 0, x.size());
    const double fExp = geometricDopplerHz(gt->at(1.0), uav->at(1.0), FC);
    EXPECT_GT(fExp, 50.0);  // approaching: positive Doppler (~ +56 Hz)
    EXPECT_NEAR(f, fExp, 0.01);
    for (auto v : y) EXPECT_NEAR(std::abs(v), 1.0, 1e-5);
}

TEST(LosDoppler, MidpointPhaseIsExactAtPacketBoundariesAndBounded) {
    // Overflight at 15 m/s, 30 m: Doppler rate up to ~ v^2/(lambda dmin)
    const double fs = 20e3;
    auto gt = std::make_shared<StaticTrajectory>(Vec3{0, 0, 10});
    auto uav = line({-150, 0, 30}, {15, 0, 0}, 30);
    ModelConfig c = losOnly(fs);
    c.updateS = 0.5e-3;
    LinkChannel lc(c, gt, uav);
    const size_t N = static_cast<size_t>(fs * 20);
    std::vector<cf> x(N, cf(1, 0)), y(N);
    lc.process(x.data(), y.data(), N);
    double maxErr = 0;
    for (size_t n = 0; n < N; n += 7) {
        const double t = n / fs;
        const double d = norm(uav->at(t).pos - gt->at(t).pos);
        const double exact = -2 * kPi * d / LAM;
        double e = std::arg(cd(y[n].real(), y[n].imag()) * std::polar(1.0, -exact));
        maxErr = std::max(maxErr, std::abs(e));
    }
    // bound pi/4 |fdot| Delta^2 with |fdot| <= v^2/(lambda dmin) = 127.6 Hz/s
    const double bound = kPi / 4 * 130.0 * c.updateS * c.updateS;
    EXPECT_LT(maxErr, bound + 2e-4);
}

TEST(LosDoppler, StartHoldAccumulatesError) {
    const double fs = 20e3;
    auto gt = std::make_shared<StaticTrajectory>(Vec3{0, 0, 10});
    auto uav = line({-150, 0, 30}, {15, 0, 0}, 30);
    ModelConfig c = losOnly(fs);
    c.updateS = 0.1;  // CHEM-like 100 ms hold
    c.phaseUpdate = PhaseUpdate::START_HOLD;
    LinkChannel lc(c, gt, uav);
    const size_t N = static_cast<size_t>(fs * 15);
    std::vector<cf> x(N, cf(1, 0)), y(N);
    lc.process(x.data(), y.data(), N);
    auto ph = uphase(y);
    // exact unwrapped phase difference
    const double d0 = norm(uav->at(0).pos - gt->at(0).pos);
    const double t = (N - 1) / fs;
    const double d1 = norm(uav->at(t).pos - gt->at(t).pos);
    const double exact = -2 * kPi * (d1 - d0) / LAM;
    const double err = (ph.back() - ph.front()) - exact;
    // predicted: -pi*Delta*(fD(t)-fD(0)) (report eq. 7), fD swings ~ +170 -> -170 Hz
    const double fD0 = geometricDopplerHz(gt->at(0), uav->at(0), FC);
    const double fD1 = geometricDopplerHz(gt->at(t), uav->at(t), FC);
    const double pred = -kPi * c.updateS * (fD1 - fD0);
    EXPECT_GT(std::abs(err), 10.0);       // large accumulated error
    EXPECT_NEAR(err, pred, 0.25 * std::abs(pred));
}

TEST(PerTap, ScattererHasOwnDopplerAndDelay) {
    const double fs = 1e6;
    auto gt = std::make_shared<StaticTrajectory>(Vec3{0, 0, 10});
    auto uav = line({300, 0, 50}, {-8, 3, 0});
    ModelConfig c = losOnly(fs);
    PathSpec sc;
    sc.type = PathType::SCATTERER;
    sc.point = {150, 400, 15};
    sc.gainDb = 0;
    sc.kDb = kInf;
    sc.shape = Shape::NONE;
    c.paths.push_back(sc);
    LinkChannel lc(c, gt, uav);
    // periodic impulses: output at kP + d gives the tap gains
    const size_t P = 2000, K = 500;
    std::vector<cf> x(P * K, cf(0, 0)), y(x.size());
    for (size_t k = 0; k < K; ++k) x[k * P] = cf(1, 0);
    lc.process(x.data(), y.data(), x.size());
    auto snap = lc.snapshot(0.5);
    const int d1 = static_cast<int>(std::lround(snap.paths[1].excessDelayS * fs));
    EXPECT_EQ(d1, 2);  // 553 m excess path = 1.84 us
    std::vector<cf> h0(K), h1(K);
    for (size_t k = 0; k < K; ++k) { h0[k] = y[k * P]; h1[k] = y[k * P + d1]; }
    const double fsTap = fs / P;
    const double f0 = slopeHz(uphase(h0), fsTap, 0, K);
    const double f1 = slopeHz(uphase(h1), fsTap, 0, K);
    EXPECT_NEAR(f0, snap.paths[0].geoDopplerHz, 0.05);
    EXPECT_NEAR(f1, snap.paths[1].geoDopplerHz, 0.05);
    EXPECT_GT(std::abs(f0 - f1), 5.0);  // distinct Dopplers per tap
}

TEST(PerTap, SharedLosModeUsesOneDoppler) {
    const double fs = 1e6;
    auto gt = std::make_shared<StaticTrajectory>(Vec3{0, 0, 10});
    auto uav = line({300, 0, 50}, {-8, 3, 0});
    ModelConfig c = losOnly(fs);
    c.dopplerMode = DopplerMode::SHARED_LOS;
    PathSpec sc;
    sc.type = PathType::SCATTERER;
    sc.point = {150, 400, 15};
    sc.kDb = kInf;
    sc.shape = Shape::NONE;
    c.paths.push_back(sc);
    LinkChannel lc(c, gt, uav);
    std::vector<cf> x(1000, cf(1, 0)), y(1000);
    lc.process(x.data(), y.data(), x.size());
    const auto& s = lc.lastPacket();
    EXPECT_NEAR(s.paths[0].dopplerHz, s.paths[1].dopplerHz, 1e-9);
}

TEST(GroundReflection, MatchesTwoRayClosedForm) {
    const double fs = 1e6;
    auto gt = std::make_shared<StaticTrajectory>(Vec3{0, 0, 12});
    auto uav = std::make_shared<StaticTrajectory>(Vec3{80, 0, 40});
    ModelConfig c = losOnly(fs);
    PathSpec g;
    g.type = PathType::GROUND;
    g.reflCoeff = -0.3;
    g.kDb = kInf;
    g.shape = Shape::NONE;
    c.paths.push_back(g);
    LinkChannel lc(c, gt, uav);
    std::vector<cf> x(4000, cf(1, 0)), y(4000);
    lc.process(x.data(), y.data(), x.size());
    const double pDb = 10 * std::log10(std::norm(cd(y[3000].real(), y[3000].imag())));
    EXPECT_NEAR(pDb, twoRayGainDb(80, 12, 40, FC, -0.3), 0.01);
    auto s = lc.snapshot(0);
    const double d1 = std::sqrt(80 * 80 + 28 * 28), d2 = std::sqrt(80 * 80 + 52 * 52);
    EXPECT_NEAR(s.paths[1].excessDelayS, (d2 - d1) / kC0, 1e-12);
}

TEST(Diffuse, GaussSosHasExactRmsWidth) {
    std::mt19937_64 rng(7);
    SosProcess g(Shape::GAUSS, 32, rng);
    double ms = 0;
    for (double f : g.unitFreqs()) ms += f * f;
    EXPECT_NEAR(std::sqrt(ms / 32), 1.0, 1e-9);
}

TEST(Diffuse, JakesRmsWidthIsFmOverSqrt2) {
    std::mt19937_64 rng(7);
    SosProcess g(Shape::JAKES, 64, rng);
    double ms = 0;
    for (double f : g.unitFreqs()) ms += f * f;
    EXPECT_NEAR(std::sqrt(ms / 64), 1.0 / std::sqrt(2.0), 1e-6);
}

TEST(Diffuse, EnsembleCorrelationMatchesGaussian) {
    // |R(tau)| of the diffuse part vs exp(-2 pi^2 sigma^2 tau^2)
    const double sig = 4.0, dt = 1e-3;
    const int lag = 50;  // 50 ms
    cd acc = 0;
    double p = 0;
    for (int s = 0; s < 400; ++s) {
        std::mt19937_64 rng(100 + s);
        SosProcess g(Shape::GAUSS, 32, rng);
        std::vector<cd> v(2000);
        for (auto& e : v) e = g.step(dt, sig, 0);
        for (size_t i = 0; i + lag < v.size(); ++i) {
            acc += v[i + lag] * std::conj(v[i]);
            p += std::norm(v[i]);
        }
    }
    const double R = std::abs(acc) / p;
    const double tau = lag * dt;
    EXPECT_NEAR(R, std::exp(-2 * kPi * kPi * sig * sig * tau * tau), 0.03);
}

TEST(Diffuse, RicianKFactorAndUnitPower) {
    const double fs = 2e3;
    auto gt = std::make_shared<StaticTrajectory>(Vec3{0, 0, 10});
    auto uav = std::make_shared<StaticTrajectory>(Vec3{200, 0, 70});
    ModelConfig c = losOnly(fs);
    c.paths[0].kDb = 8.5;
    c.paths[0].shape = Shape::GAUSS;
    c.paths[0].sigmaHz = 40.0;  // wide, for many independent samples
    c.nSinusoids = 64;
    double m2 = 0, m4 = 0;
    int n = 0;
    for (int s = 0; s < 20; ++s) {
        c.seed = 11 + s;
        LinkChannel lc(c, gt, uav);
        std::vector<cf> x(40000, cf(1, 0)), y(40000);
        lc.process(x.data(), y.data(), x.size());
        for (auto v : y) {
            const double a2 = std::norm(cd(v.real(), v.imag()));
            m2 += a2; m4 += a2 * a2; ++n;
        }
    }
    m2 /= n; m4 /= n;
    EXPECT_NEAR(m2, 1.0, 0.03);
    // moment-method K (Greenstein): K = sqrt(1-g)/(1-sqrt(1-g)), g = var/m2^2
    const double gm = (m4 - m2 * m2) / (m2 * m2);
    const double K = std::sqrt(1 - gm) / (1 - std::sqrt(1 - gm));
    EXPECT_NEAR(10 * std::log10(K), 8.5, 0.8);
}

TEST(Diffuse, WidthTableInterpolationAndScaling) {
    ModelConfig c = losOnly(1e3);
    auto gt = std::make_shared<StaticTrajectory>(Vec3{0, 0, 10});
    LinkChannel lc(c, gt, gt);
    const double sc = 3.4 / 3.51;
    EXPECT_NEAR(lc.sigmaAt(30, 5), 5.0 * sc, 1e-9);
    EXPECT_NEAR(lc.sigmaAt(60, 5), 4.25 * sc, 1e-9);
    EXPECT_NEAR(lc.sigmaAt(150, 5), 3.2 * sc, 1e-9);
    EXPECT_NEAR(lc.sigmaAt(10, 15), 5.0 * sc * 3, 1e-9);
}

TEST(Model, SeedReproducible) {
    auto gt = std::make_shared<StaticTrajectory>(Vec3{0, 0, 10});
    auto uav = line({200, 0, 50}, {-5, 0, 0});
    ModelConfig c = losOnly(5e3);
    c.paths[0].kDb = 3;
    c.paths[0].shape = Shape::GAUSS;
    std::vector<cf> x(5000, cf(1, 0)), y1(5000), y2(5000), y3(5000);
    LinkChannel a(c, gt, uav), b(c, gt, uav);
    a.process(x.data(), y1.data(), 5000);
    b.process(x.data(), y2.data(), 5000);
    c.seed = 99;
    LinkChannel d(c, gt, uav);
    d.process(x.data(), y3.data(), 5000);
    double e12 = 0, e13 = 0;
    for (int i = 0; i < 5000; ++i) { e12 += std::norm(y1[i] - y2[i]); e13 += std::norm(y1[i] - y3[i]); }
    EXPECT_EQ(e12, 0.0);
    EXPECT_GT(e13, 1.0);
}

TEST(Model, StreamSplitInvariance) {
    // processing in odd-sized chunks gives the same output as one call
    auto gt = std::make_shared<StaticTrajectory>(Vec3{0, 0, 10});
    auto uav = line({200, 0, 50}, {-5, 2, 0});
    ModelConfig c = losOnly(100e3);
    PathSpec f;
    f.type = PathType::FIXED;
    f.delayS = 37e-6;
    f.dopplerHz = 23.0;
    f.gainDb = -6;
    f.kDb = kInf;
    f.shape = Shape::NONE;
    c.paths.push_back(f);
    std::mt19937 r(1);
    std::normal_distribution<float> N;
    std::vector<cf> x(30000);
    for (auto& v : x) v = cf(N(r), N(r));
    std::vector<cf> y1(x.size()), y2(x.size());
    LinkChannel a(c, gt, uav), b(c, gt, uav);
    a.process(x.data(), y1.data(), x.size());
    size_t i = 0, k = 0;
    const size_t chunks[] = {1, 333, 1000, 4097, 77};
    while (i < x.size()) {
        size_t n = std::min(chunks[k++ % 5], x.size() - i);
        b.process(x.data() + i, y2.data() + i, n);
        i += n;
    }
    double e = 0, p = 0;
    for (size_t j = 0; j < x.size(); ++j) { e += std::norm(y1[j] - y2[j]); p += std::norm(y1[j]); }
    EXPECT_LT(e / p, 1e-10);
}

TEST(Model, FractionalDelayRealisesNonIntegerDelay) {
    const double fs = 1e6;
    auto gt = std::make_shared<StaticTrajectory>(Vec3{0, 0, 10});
    ModelConfig c = losOnly(fs);
    c.fractionalDelay = true;
    c.fracTaps = 32;
    c.paths[0].type = PathType::FIXED;
    c.paths[0].delayS = 2.5e-6;  // 2.5 samples
    LinkChannel lc(c, gt, gt);
    // slow complex tone (band-limited)
    const double f0 = 37e3;
    std::vector<cf> x(8000), y(8000);
    for (size_t n = 0; n < x.size(); ++n)
        x[n] = cf(std::cos(2 * kPi * f0 * n / fs), std::sin(2 * kPi * f0 * n / fs));
    lc.process(x.data(), y.data(), x.size());
    // realised delay = 2.5 + (fracTaps/2 - 1) samples
    const double D = 2.5 + (c.fracTaps / 2 - 1);
    double err = 0;
    for (size_t n = 200; n < 7800; ++n) {
        cd ref = std::polar(1.0, 2 * kPi * f0 * (n - D) / fs);
        err = std::max(err, std::abs(cd(y[n].real(), y[n].imag()) - ref));
    }
    EXPECT_LT(err, 0.02);
}
