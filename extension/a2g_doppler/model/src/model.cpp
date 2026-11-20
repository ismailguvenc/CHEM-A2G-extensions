// SPDX-License-Identifier: AGPL-3.0-or-later
// A2G per-tap Doppler channel model for CHEM (extension "a2g_doppler").
#include "a2g/model.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace a2g {

double norm(const Vec3& a) { return std::sqrt(a.x * a.x + a.y * a.y + a.z * a.z); }
double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

static std::string upper(std::string s) {
    for (auto& c : s) c = static_cast<char>(std::toupper(c));
    return s;
}

PathType pathTypeFromString(const std::string& s0) {
    const auto s = upper(s0);
    if (s == "LOS") return PathType::LOS;
    if (s == "GROUND" || s == "GROUND_REFLECTION" || s == "2RAY") return PathType::GROUND;
    if (s == "SCATTERER" || s == "POINT") return PathType::SCATTERER;
    if (s == "FIXED" || s == "TAP") return PathType::FIXED;
    throw std::invalid_argument("unknown path type: " + s0);
}
Shape shapeFromString(const std::string& s0) {
    const auto s = upper(s0);
    if (s == "NONE") return Shape::NONE;
    if (s == "GAUSS" || s == "GAUSSIAN") return Shape::GAUSS;
    if (s == "EXP" || s == "LAPLACE") return Shape::EXP;
    if (s == "JAKES" || s == "CLARKE") return Shape::JAKES;
    throw std::invalid_argument("unknown spectrum shape: " + s0);
}
WidthScaling widthScalingFromString(const std::string& s0) {
    const auto s = upper(s0);
    if (s == "NONE") return WidthScaling::NONE;
    if (s == "CARRIER") return WidthScaling::CARRIER;
    if (s == "CARRIER_SPEED") return WidthScaling::CARRIER_SPEED;
    throw std::invalid_argument("unknown width scaling: " + s0);
}
DopplerModel dopplerModelFromString(const std::string& s0) {
    const auto s = upper(s0);
    if (s == "PROPOSED" || s == "A2G" || s == "A2G_MEAS") return DopplerModel::PROPOSED;
    if (s == "JAKES" || s == "CLARKE") return DopplerModel::JAKES;
    if (s == "TDL_A" || s == "TDL-A") return DopplerModel::TDL_A;
    if (s == "TDL_D" || s == "TDL-D") return DopplerModel::TDL_D;
    throw std::invalid_argument("unknown Doppler model: " + s0);
}
std::string toString(DopplerModel m) {
    switch (m) {
        case DopplerModel::PROPOSED: return "PROPOSED";
        case DopplerModel::JAKES: return "JAKES";
        case DopplerModel::TDL_A: return "TDL_A";
        case DopplerModel::TDL_D: return "TDL_D";
    }
    return "?";
}

// TR 38.901 V17.0.0 Table 7.7.2-1 (TDL-A) and Table 7.7.2-4 (TDL-D).
// Transcribed for this extension; check against the specification before a
// release. TDL-D tap 1 = LOS -0.2 dB + Rayleigh -13.5 dB (K1 = 13.3 dB); the
// two parts are separate entries, so that only the specular part carries
// the LOS Doppler line and the Rayleigh part stays centred on 0 Hz.
const std::vector<TdlTap>& tdlTable(DopplerModel m) {
    constexpr double R = -kSpecularDb;  // Rayleigh
    constexpr double S = kSpecularDb;   // specular (LOS line)
    static const std::vector<TdlTap> A = {
        {0.0000, -13.4, R}, {0.3819, 0.0, R},   {0.4025, -2.2, R},  {0.5868, -4.0, R},
        {0.4610, -6.0, R},  {0.5375, -8.2, R},  {0.6708, -9.9, R},  {0.5750, -10.5, R},
        {0.7618, -7.5, R},  {1.5375, -15.9, R}, {1.8978, -6.6, R},  {2.2242, -16.7, R},
        {2.1718, -12.4, R}, {2.4942, -15.2, R}, {2.5119, -10.8, R}, {3.0582, -11.3, R},
        {4.0810, -12.7, R}, {4.4579, -16.2, R}, {4.5695, -18.3, R}, {4.7966, -18.9, R},
        {5.0066, -16.6, R}, {5.3043, -19.9, R}, {9.6586, -29.7, R}};
    static const std::vector<TdlTap> D = {
        {0.000, -0.2, S},       {0.000, -13.5, R}, {0.035, -18.8, R}, {0.612, -21.0, R}, {1.363, -22.8, R},
        {1.405, -17.9, R},      {1.804, -20.1, R}, {2.596, -21.9, R}, {1.775, -22.9, R},
        {4.042, -27.8, R},      {7.937, -23.6, R}, {9.424, -24.8, R}, {9.708, -30.0, R},
        {12.525, -27.7, R}};
    static const std::vector<TdlTap> none;
    if (m == DopplerModel::TDL_A) return A;
    if (m == DopplerModel::TDL_D) return D;
    return none;
}

void applyDopplerPreset(ModelConfig& c) {
    switch (c.dopplerModel) {
        case DopplerModel::PROPOSED:
            return;
        case DopplerModel::JAKES:
            // Jakes on every tap: delays and powers from the geometry, no
            // geometric Doppler line, Rayleigh fading with the classical
            // spectrum on each tap.
            for (auto& p : c.paths) {
                p.kDb = -kInf;
                p.shape = Shape::JAKES;
                p.geoPhase = false;
                if (p.type == PathType::FIXED) {
                    p.useAoa = false;
                    p.dopplerHz = 0.0;
                }
            }
            return;
        case DopplerModel::TDL_A:
        case DopplerModel::TDL_D: {
            c.paths.clear();
            for (const auto& t : tdlTable(c.dopplerModel)) {
                PathSpec p;
                p.type = PathType::FIXED;
                p.delayS = t.tauNorm * c.tdlDelaySpreadS;
                p.gainDb = t.powerDb;
                p.kDb = t.kDb;
                if (isSpecular(t.kDb)) {  // LOS part: line at fm cos(45 deg)
                    p.shape = Shape::NONE;
                    p.useAoa = true;
                    p.aoaDeg = 45.0;
                } else {
                    p.shape = Shape::JAKES;
                }
                c.paths.push_back(p);
            }
            c.normalizeTotal = true;
            return;
        }
    }
}

std::string toString(PathType t) {
    switch (t) {
        case PathType::LOS: return "LOS";
        case PathType::GROUND: return "GROUND";
        case PathType::SCATTERER: return "SCATTERER";
        case PathType::FIXED: return "FIXED";
    }
    return "?";
}
std::string toString(Shape s) {
    switch (s) {
        case Shape::NONE: return "NONE";
        case Shape::GAUSS: return "GAUSS";
        case Shape::EXP: return "EXP";
        case Shape::JAKES: return "JAKES";
    }
    return "?";
}

// ---------------------------------------------------------------------------
// Trajectories
// ---------------------------------------------------------------------------
static KinState hermiteEval(const std::vector<double>& ts,
                            const std::vector<Vec3>& ps,
                            const std::vector<Vec3>& vs, double t) {
    const size_t n = ts.size();
    if (n == 0) return {};
    if (n == 1 || t <= ts.front()) {
        return {ps.front() + (t - ts.front()) * vs.front(), vs.front()};
    }
    if (t >= ts.back()) {
        return {ps.back() + (t - ts.back()) * vs.back(), vs.back()};
    }
    const size_t i1 = static_cast<size_t>(
        std::upper_bound(ts.begin(), ts.end(), t) - ts.begin());
    const size_t i0 = i1 - 1;
    const double h = ts[i1] - ts[i0];
    if (h <= 0) return {ps[i1], vs[i1]};
    const double s = (t - ts[i0]) / h, s2 = s * s, s3 = s2 * s;
    const double h00 = 2 * s3 - 3 * s2 + 1, h10 = s3 - 2 * s2 + s,
                 h01 = -2 * s3 + 3 * s2, h11 = s3 - s2;
    const double d00 = 6 * s2 - 6 * s, d10 = 3 * s2 - 4 * s + 1,
                 d01 = -6 * s2 + 6 * s, d11 = 3 * s2 - 2 * s;
    KinState k;
    k.pos = h00 * ps[i0] + (h10 * h) * vs[i0] + h01 * ps[i1] + (h11 * h) * vs[i1];
    k.vel = (1.0 / h) * (d00 * ps[i0] + (d10 * h) * vs[i0] + d01 * ps[i1] +
                         (d11 * h) * vs[i1]);
    return k;
}

void TableTrajectory::add(double t, Vec3 p, Vec3 v) {
    if (!m_t.empty() && t <= m_t.back()) return;  // keep strictly increasing
    m_t.push_back(t);
    m_p.push_back(p);
    m_v.push_back(v);
}

bool TableTrajectory::loadCsv(const std::string& path) {
    std::ifstream f(path);
    if (!f) return false;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        const char c = line[0];
        if (!(std::isdigit(static_cast<unsigned char>(c)) || c == '-' || c == '.'))
            continue;  // header
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream ss(line);
        double t, e, n, u, ve = 0, vn = 0, vu = 0;
        if (!(ss >> t >> e >> n >> u)) continue;
        ss >> ve >> vn >> vu;
        add(t, {e, n, u}, {ve, vn, vu});
    }
    return !m_t.empty();
}

KinState TableTrajectory::at(double t) const { return hermiteEval(m_t, m_p, m_v, t); }

void LiveTrajectory::push(double t, Vec3 p, Vec3 v, bool hasVel) {
    std::lock_guard<std::mutex> g(m_mu);
    if (!m_t.empty() && t <= m_t.back()) return;
    if (!hasVel && !m_t.empty()) {
        const double dt = t - m_t.back();
        v = (1.0 / dt) * (p - m_p.back());
    }
    m_t.push_back(t);
    m_p.push_back(p);
    m_v.push_back(v);
    if (m_t.size() > m_keep) {
        m_t.erase(m_t.begin());
        m_p.erase(m_p.begin());
        m_v.erase(m_v.begin());
    }
}

KinState LiveTrajectory::at(double t) const {
    std::lock_guard<std::mutex> g(m_mu);
    return hermiteEval(m_t, m_p, m_v, t - m_delay);
}

bool LiveTrajectory::empty() const {
    std::lock_guard<std::mutex> g(m_mu);
    return m_t.empty();
}

// ---------------------------------------------------------------------------
// Diffuse sum-of-sinusoids process
// ---------------------------------------------------------------------------
static double normPpf(double p) {
    double lo = -9, hi = 9;
    for (int i = 0; i < 100; ++i) {
        const double mid = 0.5 * (lo + hi);
        const double c = 0.5 * (1.0 + std::erf(mid / std::sqrt(2.0)));
        (c < p ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}

SosProcess::SosProcess(Shape shape, int n, std::mt19937_64& rng) : m_shape(shape) {
    if (shape == Shape::NONE || n <= 0) return;
    std::uniform_real_distribution<double> U(0.0, 1.0);
    m_unit.resize(n);
    m_theta.resize(n);
    if (shape == Shape::JAKES) {
        // Angles of arrival equally spaced with a random rotation (Clarke /
        // Jakes isotropic scattering): f_i = fm cos(alpha_i), rms fm/sqrt(2).
        const double rot = U(rng);
        for (int i = 0; i < n; ++i) m_unit[i] = std::cos(2 * kPi * (i + rot) / n);
    } else {
        // Stratified quantiles of the target spectrum; rescaled so the rms
        // width is exactly 1 (report Sec. 4.5, "SoS with quantile frequencies").
        for (int i = 0; i < n; ++i) {
            const double p = (i + U(rng)) / n;
            if (shape == Shape::GAUSS) {
                m_unit[i] = normPpf(std::clamp(p, 1e-9, 1 - 1e-9));
            } else {  // EXP (two-sided exponential / Laplace)
                const double b = 1.0 / std::sqrt(2.0);
                const double q = p - 0.5;
                m_unit[i] = -b * (q < 0 ? -1.0 : 1.0) *
                            std::log(std::max(1e-12, 1 - 2 * std::abs(q)));
            }
        }
        double ms = 0, mean = 0;
        for (double f : m_unit) mean += f;
        mean /= n;
        for (double& f : m_unit) f -= mean;  // spectrum centred on 0
        for (double f : m_unit) ms += f * f;
        const double rms = std::sqrt(ms / n);
        for (double& f : m_unit) f /= rms;
    }
    for (int i = 0; i < n; ++i) m_theta[i] = 2 * kPi * U(rng);
}

cd SosProcess::value() const {
    if (m_unit.empty()) return {0, 0};
    cd s{0, 0};
    for (double th : m_theta) s += std::polar(1.0, th);
    return s / std::sqrt(static_cast<double>(m_theta.size()));
}

cd SosProcess::step(double dt, double sigmaHz, double fmHz) {
    if (m_unit.empty()) return {0, 0};
    const double scale = (m_shape == Shape::JAKES) ? fmHz : sigmaHz;
    for (size_t i = 0; i < m_unit.size(); ++i) {
        m_theta[i] = std::fmod(m_theta[i] + 2 * kPi * m_unit[i] * scale * dt, 2 * kPi);
    }
    return value();
}

// ---------------------------------------------------------------------------
// LinkChannel
// ---------------------------------------------------------------------------
LinkChannel::LinkChannel(const ModelConfig& cfg, std::shared_ptr<const Trajectory> tx,
                         std::shared_ptr<const Trajectory> rx)
    : m_cfg(cfg), m_tx(std::move(tx)), m_rx(std::move(rx)), m_rng(cfg.seed) {
    if (!(m_cfg.fs > 0)) throw std::invalid_argument("fs must be > 0");
    m_packetLen = std::max<uint64_t>(1, static_cast<uint64_t>(std::llround(m_cfg.updateS * m_cfg.fs)));
    std::mt19937_64 rng(cfg.seed);
    for (const auto& p : m_cfg.paths) {
        PathRt rt;
        rt.spec = p;
        const bool diffuse = p.shape != Shape::NONE && !isSpecular(p.kDb);
        if (diffuse) rt.sos = SosProcess(p.shape, m_cfg.nSinusoids, rng);
        m_paths.push_back(std::move(rt));
    }
    m_histLen = static_cast<size_t>(std::max(0, m_cfg.maxDelaySamples)) +
                static_cast<size_t>(m_cfg.fractionalDelay ? m_cfg.fracTaps : 0) + 1;
    m_hist.assign(m_histLen, cf(0, 0));
}

void LinkChannel::setTime(double t) {
    m_t0 = t;
    m_nTot = 0;
    m_nInPacket = m_packetLen;  // force a new packet at the next sample
}

double LinkChannel::sigmaAt(double altM, double speed) const {
    const auto& A = m_cfg.altTab;
    const auto& S = m_cfg.sigTab;
    double s = S.empty() ? 0.0 : S.front();
    if (A.size() == S.size() && A.size() >= 2) {
        if (altM <= A.front()) s = S.front();
        else if (altM >= A.back()) s = S.back();
        else {
            for (size_t i = 1; i < A.size(); ++i) {
                if (altM <= A[i]) {
                    const double w = (altM - A[i - 1]) / (A[i] - A[i - 1]);
                    s = S[i - 1] + w * (S[i] - S[i - 1]);
                    break;
                }
            }
        }
    }
    if (m_cfg.widthScaling != WidthScaling::NONE) s *= m_cfg.fc / m_cfg.widthRefFc;
    if (m_cfg.widthScaling == WidthScaling::CARRIER_SPEED)
        s *= std::max(speed, 0.0) / m_cfg.widthRefSpeed;
    return s;
}

double LinkChannel::pathLength(const PathSpec& p, const KinState& tx,
                               const KinState& rx) const {
    switch (p.type) {
        case PathType::LOS:
        case PathType::FIXED:
            return norm(tx.pos - rx.pos);
        case PathType::GROUND: {
            Vec3 img = rx.pos;
            img.z = 2 * m_cfg.groundZ - img.z;  // image of RX in the ground
            return norm(tx.pos - img);
        }
        case PathType::SCATTERER:
            return norm(tx.pos - p.point) + norm(p.point - rx.pos);
    }
    return 0;
}

void LinkChannel::uavState(const KinState& tx, const KinState& rx, double& alt,
                           double& speed) const {
    const KinState& u = (tx.pos.z >= rx.pos.z) ? tx : rx;
    alt = u.pos.z - m_cfg.groundZ;
    speed = norm(tx.vel - rx.vel);
}

double LinkChannel::fixedFreq(const PathSpec& p, double speed) const {
    if (p.useAoa) {
        const double fm = std::max(speed, m_cfg.speedFloor) / lambda(m_cfg.fc);
        return fm * std::cos(p.aoaDeg * kPi / 180.0);
    }
    return p.dopplerHz;
}

static double wrap(double x) {
    x = std::fmod(x, 2 * kPi);
    return x < 0 ? x + 2 * kPi : x;
}

void LinkChannel::startPacket(double tk) {
    const double D = static_cast<double>(m_packetLen) / m_cfg.fs;
    const double lam = lambda(m_cfg.fc);
    const KinState tx0 = m_tx->at(tk), rx0 = m_rx->at(tk);
    const KinState tx1 = m_tx->at(tk + D), rx1 = m_rx->at(tk + D);
    const double eps = 1e-3;
    const KinState txa = m_tx->at(tk + eps), rxa = m_rx->at(tk + eps);
    const KinState txb = m_tx->at(tk - eps), rxb = m_rx->at(tk - eps);
    double alt, speed;
    uavState(tx0, rx0, alt, speed);
    const double sigma = sigmaAt(alt, speed);
    const double fm = std::max(speed, m_cfg.speedFloor) / lam;
    const double dLos0 = norm(tx0.pos - rx0.pos);
    const double dLos1 = norm(tx1.pos - rx1.pos);
    const double kap = m_cfg.kappa;
    const bool firstPacket = !m_started;

    double fLos = (m_cfg.phaseUpdate == PhaseUpdate::MIDPOINT)
                      ? -kap * (dLos1 - dLos0) / (lam * D)
                      : -kap * (norm(txa.pos - rxa.pos) - norm(txb.pos - rxb.pos)) /
                            (2 * eps * lam);
    const double psiLosAnchor = wrap(-2 * kPi * kap * dLos0 / lam);

    m_last.t = tk;
    m_last.uavAltM = alt;
    m_last.uavSpeed = speed;
    m_last.paths.resize(m_paths.size());

    double pTot = 0;
    for (size_t i = 0; i < m_paths.size(); ++i) {
        auto& P = m_paths[i];
        const auto& s = P.spec;
        const double L0 = pathLength(s, tx0, rx0);
        const double L1 = pathLength(s, tx1, rx1);
        const double geo = -(pathLength(s, txa, rxa) - pathLength(s, txb, rxb)) / (2 * eps * lam);
        double amp = 1.0, exDelay = 0.0;
        switch (s.type) {
            case PathType::LOS: amp = 1.0; break;
            case PathType::GROUND: amp = s.reflCoeff * dLos0 / std::max(L0, 1e-3); break;
            case PathType::SCATTERER:
                amp = std::pow(10.0, s.gainDb / 20.0) * dLos0 / std::max(L0, 1e-3);
                break;
            case PathType::FIXED: amp = std::pow(10.0, s.gainDb / 20.0); break;
        }
        exDelay = (s.type == PathType::FIXED) ? s.delayS : (L0 - dLos0) / kC0;

        // --- frequency and phase -------------------------------------------
        const bool geometric = s.type != PathType::FIXED && s.geoPhase;
        if (m_cfg.dopplerMode == DopplerMode::SHARED_LOS) {
            // all paths rotate with the LOS phase (+ a fixed offset)
            if (firstPacket) P.psi = wrap(psiLosAnchor - 2 * kPi * (L0 - dLos0) / lam);
            else if (m_cfg.phaseUpdate == PhaseUpdate::MIDPOINT) {
                // keep the offset found at the first packet
            }
            P.freq = fLos;
        } else if (geometric) {
            P.freq = (m_cfg.phaseUpdate == PhaseUpdate::MIDPOINT)
                         ? -kap * (L1 - L0) / (lam * D)
                         : kap * geo;
            if (firstPacket || m_cfg.phaseUpdate == PhaseUpdate::MIDPOINT)
                P.psi = wrap(-2 * kPi * kap * L0 / lam);  // exact range phase
        } else {
            P.freq = fixedFreq(s, speed);
            if (firstPacket) P.psi = 0.0;
        }

        // --- diffuse part ----------------------------------------------------
        cd H{amp, 0.0};
        const double sp = (s.sigmaHz < 0) ? sigma : s.sigmaHz;
        if (!m_paths[i].sos.unitFreqs().empty()) {
            const cd g = firstPacket ? P.sos.value() : P.sos.step(D, sp, fm);
            const double K = kLinear(s.kDb);
            H = amp * (std::sqrt(K / (K + 1)) + std::sqrt(1 / (K + 1)) * g);
        }
        P.H = H;
        pTot += amp * amp;

        // --- delay ------------------------------------------------------------
        const double dly = std::max(0.0, exDelay * m_cfg.fs);
        if (m_cfg.fractionalDelay) {
            const int Lf = std::max(2, m_cfg.fracTaps);
            const int c = Lf / 2 - 1;
            const int di = static_cast<int>(std::floor(dly));
            const double mu = dly - di;
            P.dInt = std::min(di, m_cfg.maxDelaySamples);
            P.frac.assign(Lf, 0.f);
            double sum = 0;
            for (int m = 0; m < Lf; ++m) {
                const double x = m - c - mu;
                const double sinc = (std::abs(x) < 1e-12) ? 1.0 : std::sin(kPi * x) / (kPi * x);
                const double w = 0.5 - 0.5 * std::cos(2 * kPi * (m + 0.5 - mu) / Lf);  // Hann
                P.frac[m] = static_cast<float>(sinc * w);
                sum += sinc * w;
            }
            for (auto& v : P.frac) v = static_cast<float>(v / sum);
        } else {
            P.dInt = std::min(static_cast<int>(std::lround(dly)), m_cfg.maxDelaySamples);
            P.frac.clear();
        }

        auto& st = m_last.paths[i];
        st.lengthM = L0;
        st.excessDelayS = exDelay;
        st.dopplerHz = P.freq;
        st.geoDopplerHz = geometric ? kap * geo : P.freq;
        st.amp = amp;
        st.gain = H;
        st.phase = P.psi;
        st.sigmaHz = sp;
    }
    m_normScale = (m_cfg.normalizeTotal && pTot > 0) ? 1.0 / std::sqrt(pTot) : 1.0;

    // oscillator impairment (common to all paths), phase continuous
    if (m_cfg.osc.enable) {
        std::normal_distribution<double> N01(0, 1);
        if (tk >= m_oscNextJit) {
            m_oscJit = m_cfg.osc.jitterHzRms * N01(m_rng);
            m_oscNextJit = tk + m_cfg.osc.jitterHoldS;
        }
        // smooth random walk: 0.05 Hz/sqrt(s)
        m_oscRw += 0.05 * std::sqrt(D) * N01(m_rng);
    }
    m_tk = tk;
    m_started = true;
}

void LinkChannel::process(const cf* in, cf* out, size_t n) {
    if (n == 0) return;
    // extended input: history followed by the new samples
    thread_local std::vector<cf> ext;
    thread_local std::vector<cd> acc;
    ext.resize(m_histLen + n);
    std::copy(m_hist.begin(), m_hist.end(), ext.begin());
    std::copy(in, in + n, ext.begin() + m_histLen);
    acc.assign(n, cd(0, 0));

    size_t i = 0;
    while (i < n) {
        if (!m_started || m_nInPacket >= m_packetLen) {
            startPacket(m_t0 + static_cast<double>(m_nTot + i) / m_cfg.fs);
            m_nInPacket = 0;
        }
        const size_t seg = std::min<uint64_t>(n - i, m_packetLen - m_nInPacket);
        for (auto& P : m_paths) {
            const cd H = P.H * m_normScale;
            if (std::abs(H) == 0.0) {
                P.psi = wrap(P.psi + 2 * kPi * P.freq * seg / m_cfg.fs);
                continue;
            }
            cd ph = H * std::polar(1.0, P.psi);
            const cd stp = std::polar(1.0, 2 * kPi * P.freq / m_cfg.fs);
            const size_t base = m_histLen + i - static_cast<size_t>(P.dInt);
            if (P.frac.empty()) {
                for (size_t s = 0; s < seg; ++s) {
                    const cf x = ext[base + s];
                    acc[i + s] += ph * cd(x.real(), x.imag());
                    ph *= stp;
                }
            } else {
                const int Lf = static_cast<int>(P.frac.size());
                for (size_t s = 0; s < seg; ++s) {
                    cd xs(0, 0);
                    for (int m = 0; m < Lf; ++m) {
                        const cf x = ext[base + s - m];
                        xs += static_cast<double>(P.frac[m]) * cd(x.real(), x.imag());
                    }
                    acc[i + s] += ph * xs;
                    ph *= stp;
                }
            }
            P.psi = wrap(P.psi + 2 * kPi * P.freq * seg / m_cfg.fs);
        }
        if (m_cfg.osc.enable) {
            const double t = m_t0 + static_cast<double>(m_nTot + i) / m_cfg.fs;
            const double f = m_cfg.osc.driftHzPer30min / 1800.0 * t + m_oscRw + m_oscJit;
            const double dps = 2 * kPi * f / m_cfg.fs;
            for (size_t s = 0; s < seg; ++s) {
                acc[i + s] *= std::polar(1.0, m_oscPsi);
                m_oscPsi = wrap(m_oscPsi + dps);
            }
        }
        m_nInPacket += seg;
        i += seg;
    }
    m_nTot += n;
    for (size_t s = 0; s < n; ++s)
        out[s] = cf(static_cast<float>(acc[s].real()), static_cast<float>(acc[s].imag()));
    std::copy(ext.end() - static_cast<long>(m_histLen), ext.end(), m_hist.begin());
}

LinkSnapshot LinkChannel::snapshot(double t) const {
    LinkSnapshot snap;
    const double lam = lambda(m_cfg.fc);
    const KinState tx = m_tx->at(t), rx = m_rx->at(t);
    const double eps = 1e-3;
    const KinState txa = m_tx->at(t + eps), rxa = m_rx->at(t + eps);
    const KinState txb = m_tx->at(t - eps), rxb = m_rx->at(t - eps);
    double alt, speed;
    uavState(tx, rx, alt, speed);
    snap.t = t;
    snap.uavAltM = alt;
    snap.uavSpeed = speed;
    const double d0 = norm(tx.pos - rx.pos);
    for (const auto& P : m_paths) {
        const auto& s = P.spec;
        PathState st;
        st.lengthM = pathLength(s, tx, rx);
        st.excessDelayS = (s.type == PathType::FIXED) ? s.delayS : (st.lengthM - d0) / kC0;
        st.geoDopplerHz = (s.type == PathType::FIXED || !s.geoPhase)
                              ? fixedFreq(s, speed)
                              : -m_cfg.kappa * (pathLength(s, txa, rxa) - pathLength(s, txb, rxb)) /
                                    (2 * eps * lam);
        st.dopplerHz = st.geoDopplerHz;
        switch (s.type) {
            case PathType::LOS: st.amp = 1; break;
            case PathType::GROUND: st.amp = s.reflCoeff * d0 / st.lengthM; break;
            case PathType::SCATTERER: st.amp = std::pow(10.0, s.gainDb / 20) * d0 / st.lengthM; break;
            case PathType::FIXED: st.amp = std::pow(10.0, s.gainDb / 20); break;
        }
        st.gain = st.amp;
        st.phase = s.geoPhase ? wrap(-2 * kPi * m_cfg.kappa * st.lengthM / lam) : 0.0;
        st.sigmaHz = (s.sigmaHz < 0) ? sigmaAt(alt, speed) : s.sigmaHz;
        snap.paths.push_back(st);
    }
    return snap;
}

std::vector<cf> LinkChannel::sampledCIR(double t, size_t maxLen) {
    // Taps mode: CIR sampled at time t, held by CHEM until the next update.
    const double lam = lambda(m_cfg.fc);
    const LinkSnapshot snap = snapshot(t);
    const double dt = m_started ? std::max(0.0, t - m_tk) : 0.0;
    std::vector<cf> taps;
    double pTot = 0;
    for (size_t i = 0; i < m_paths.size(); ++i) {
        auto& P = m_paths[i];
        const auto& st = snap.paths[i];
        cd H{st.amp, 0};
        if (!P.sos.unitFreqs().empty()) {
            const cd g = m_started ? P.sos.step(dt, st.sigmaHz,
                                                std::max(snap.uavSpeed, m_cfg.speedFloor) / lam)
                                   : P.sos.value();
            const double K = kLinear(P.spec.kDb);
            H = st.amp * (std::sqrt(K / (K + 1)) + std::sqrt(1 / (K + 1)) * g);
        }
        double psi;
        if (P.spec.type == PathType::FIXED || !P.spec.geoPhase) {
            P.psi = wrap(P.psi + 2 * kPi * st.dopplerHz * dt);
            psi = P.psi;
        } else if (m_cfg.dopplerMode == DopplerMode::SHARED_LOS) {
            psi = snap.paths[0].phase;
        } else {
            psi = st.phase;
        }
        pTot += st.amp * st.amp;
        const size_t d = static_cast<size_t>(std::lround(std::max(0.0, st.excessDelayS * m_cfg.fs)));
        if (d >= maxLen) continue;
        if (taps.size() <= d) taps.resize(d + 1, cf(0, 0));
        const cd h = H * std::polar(1.0, psi);
        taps[d] += cf(static_cast<float>(h.real()), static_cast<float>(h.imag()));
    }
    if (m_cfg.normalizeTotal && pTot > 0)
        for (auto& v : taps) v *= static_cast<float>(1.0 / std::sqrt(pTot));
    m_tk = t;
    m_started = true;
    return taps;
}

double geometricDopplerHz(const KinState& tx, const KinState& rx, double fc) {
    const Vec3 r = rx.pos - tx.pos;
    const double d = norm(r);
    if (d <= 0) return 0;
    const Vec3 u = (1.0 / d) * r;
    const double rr = dot(rx.vel - tx.vel, u);  // range rate
    return -rr / (kC0 / fc);
}

double twoRayGainDb(double dH, double hTx, double hRx, double fc, double gamma) {
    const double d1 = std::sqrt(dH * dH + (hTx - hRx) * (hTx - hRx));
    const double d2 = std::sqrt(dH * dH + (hTx + hRx) * (hTx + hRx));
    const double k = 2 * kPi * fc / kC0;
    const cd e = 1.0 + gamma * (d1 / d2) * std::polar(1.0, -k * (d2 - d1));
    return 10 * std::log10(std::norm(e));
}

}  // namespace a2g
