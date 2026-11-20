// SPDX-License-Identifier: AGPL-3.0-or-later
// 3GPP TR 25.996 SCM for CHEM MIMO links (extension "scm_mimo").
#include "scm/scm.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace scm {

namespace {
constexpr double kPi = a2g::kPi;
constexpr double kDeg = kPi / 180.0;
constexpr int kM = 20;  // sub-paths per path (TR 25.996)

// TR 25.996 Table 5.2: sub-path offsets (deg) for 2 deg (BS) and 35 deg (MS)
// per-path angle spread.
const double kOffBs[10] = {0.0894, 0.2826, 0.4984, 0.7431, 1.0257,
                           1.3594, 1.7688, 2.2961, 3.0389, 4.3101};
const double kOffMs[10] = {1.5649, 4.9447, 8.7224, 13.0045, 17.9492,
                           23.7899, 30.9538, 40.1824, 53.1816, 75.4274};

std::string upper(std::string s) {
    for (auto& c : s) c = static_cast<char>(std::toupper(c));
    return s;
}
double wrap(double x) {
    x = std::fmod(x, 2 * kPi);
    return x < 0 ? x + 2 * kPi : x;
}
double normPpf(double p) {
    double lo = -9, hi = 9;
    for (int i = 0; i < 100; ++i) {
        const double mid = 0.5 * (lo + hi);
        const double c = 0.5 * (1.0 + std::erf(mid / std::sqrt(2.0)));
        (c < p ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}
}  // namespace

Scenario scenarioFromString(const std::string& s0) {
    const auto s = upper(s0);
    if (s == "A2G" || s == "SMA_A2G") return Scenario::A2G;
    if (s == "SMA" || s == "SUBURBAN") return Scenario::SMA;
    if (s == "UMA" || s == "URBAN") return Scenario::UMA;
    if (s == "UMI") return Scenario::UMI;
    if (s == "UMI_LOS") return Scenario::UMI_LOS;
    throw std::invalid_argument("unknown SCM scenario: " + s0);
}
Doppler dopplerFromString(const std::string& s0) {
    const auto s = upper(s0);
    if (s == "SCM" || s == "NATIVE") return Doppler::SCM;
    if (s == "PROPOSED" || s == "A2G") return Doppler::PROPOSED;
    if (s == "JAKES") return Doppler::JAKES;
    if (s == "TDL_D" || s == "TDL-D") return Doppler::TDL_D;
    throw std::invalid_argument("unknown SCM Doppler option: " + s0);
}
std::string toString(Scenario s) {
    switch (s) {
        case Scenario::A2G: return "A2G";
        case Scenario::SMA: return "SMA";
        case Scenario::UMA: return "UMA";
        case Scenario::UMI: return "UMI";
        case Scenario::UMI_LOS: return "UMI_LOS";
    }
    return "?";
}
std::string toString(Doppler d) {
    switch (d) {
        case Doppler::SCM: return "SCM";
        case Doppler::PROPOSED: return "PROPOSED";
        case Doppler::JAKES: return "JAKES";
        case Doppler::TDL_D: return "TDL_D";
    }
    return "?";
}

// ---------------------------------------------------------------------------
// Drop (TR 25.996 Sec. 5.3, steps 3-10)
// ---------------------------------------------------------------------------
Drop drawDrop(const ScmConfig& c, Doppler dop, std::mt19937_64& rng) {
    std::normal_distribution<double> N01(0.0, 1.0);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    const int N = std::max(1, c.nPaths);
    Drop d;
    d.tau.resize(N);
    d.P.resize(N);
    const bool umi = (c.scenario == Scenario::UMI || c.scenario == Scenario::UMI_LOS);
    double offBsScale = 1.0, aoaCoef = 0.2175;
    std::vector<double> aod(N);
    if (umi) {
        // delays U(0, 1.2 us); P'_n = 10^-(tau_n[us] + z_n/10), z ~ N(0, 3 dB)
        for (int n = 0; n < N; ++n) d.tau[n] = 1.2e-6 * U(rng);
        std::sort(d.tau.begin(), d.tau.end());
        const double t0 = d.tau[0];
        for (auto& t : d.tau) t -= t0;
        for (int n = 0; n < N; ++n)
            d.P[n] = std::pow(10.0, -(d.tau[n] * 1e6 + 3.0 * N01(rng) / 10.0));
        for (int n = 0; n < N; ++n) aod[n] = (-40.0 + 80.0 * U(rng)) * kDeg;
        offBsScale = 2.5;  // 5 deg per-path AS at the BS
        aoaCoef = 0.265;
        d.dsS = 0;
        d.asDeg = 0;
    } else {
        double muAS = 0.69, epsAS = 0.13, rAS = 1.2, muDS = -6.80, epsDS = 0.288, rDS = 1.4;
        if (c.scenario == Scenario::UMA) {
            muAS = 0.810; epsAS = 0.34; rAS = 1.3; muDS = -6.18; epsDS = 0.18; rDS = 1.7;
        }
        // correlated LSPs: rho(DS, AS) = 0.5
        const double w1 = N01(rng), w2 = N01(rng);
        const double a = w1, b = 0.5 * w1 + std::sqrt(0.75) * w2;
        const double sDS = std::pow(10.0, epsDS * a + muDS);
        const double sAS = std::pow(10.0, epsAS * b + muAS);
        d.dsS = sDS;
        d.asDeg = sAS;
        std::vector<double> tp(N);
        for (int n = 0; n < N; ++n) tp[n] = -rDS * sDS * std::log(std::max(1e-12, U(rng)));
        std::sort(tp.begin(), tp.end());
        for (int n = 0; n < N; ++n) d.tau[n] = tp[n] - tp[0];
        for (int n = 0; n < N; ++n)
            d.P[n] = std::exp((1.0 - rDS) * d.tau[n] / (rDS * sDS)) *
                     std::pow(10.0, -3.0 * N01(rng) / 10.0);
        // AoD offsets: N(0, (rAS sAS)^2), smallest |delta| to the shortest delay
        std::vector<double> dd(N);
        for (int n = 0; n < N; ++n) dd[n] = rAS * sAS * N01(rng);
        std::sort(dd.begin(), dd.end(), [](double x, double y) { return std::abs(x) < std::abs(y); });
        for (int n = 0; n < N; ++n) aod[n] = dd[n] * kDeg;
    }
    const double ps = std::accumulate(d.P.begin(), d.P.end(), 0.0);
    for (auto& p : d.P) p /= ps;
    d.dAod = aod;
    d.dAoa.resize(N);
    for (int n = 0; n < N; ++n) {
        const double sAoa = 104.12 * (1.0 - std::exp(-aoaCoef * std::abs(10.0 * std::log10(d.P[n]))));
        d.dAoa[n] = sAoa * N01(rng) * kDeg;
    }
    d.sAod.assign(N, std::vector<double>(kM));
    d.sAoa.assign(N, std::vector<double>(kM));
    d.Phi.assign(N, std::vector<double>(kM));
    d.PhiPol.assign(N, std::vector<std::array<double, 4>>(kM));
    d.q.assign(N, std::vector<double>(kM, 0.0));
    d.acc.assign(N, std::vector<double>(kM, 0.0));
    for (int n = 0; n < N; ++n) {
        std::vector<int> perm(kM);
        std::iota(perm.begin(), perm.end(), 0);
        std::shuffle(perm.begin(), perm.end(), rng);
        for (int m = 0; m < kM; ++m) {
            const double ob = kOffBs[m / 2] * (m % 2 ? -1.0 : 1.0) * offBsScale;
            const int pm = perm[m];
            const double om = kOffMs[pm / 2] * (pm % 2 ? -1.0 : 1.0);
            d.sAod[n][m] = ob * kDeg;
            d.sAoa[n][m] = om * kDeg;
            d.Phi[n][m] = 2 * kPi * U(rng);
            for (auto& p : d.PhiPol[n][m]) p = 2 * kPi * U(rng);
        }
        if (dop == Doppler::PROPOSED) {
            // stratified Gaussian quantiles, unit rms, random order
            std::vector<double> qq(kM);
            for (int m = 0; m < kM; ++m) qq[m] = normPpf(std::clamp((m + U(rng)) / kM, 1e-9, 1 - 1e-9));
            double mean = std::accumulate(qq.begin(), qq.end(), 0.0) / kM, ms = 0;
            for (auto& v : qq) v -= mean;
            for (double v : qq) ms += v * v;
            const double rms = std::sqrt(ms / kM);
            std::shuffle(qq.begin(), qq.end(), rng);
            for (int m = 0; m < kM; ++m) d.q[n][m] = qq[m] / rms;
        } else if (dop == Doppler::JAKES || dop == Doppler::TDL_D) {
            // isotropic scattering, independent of the spatial angles
            const double rot = U(rng);
            for (int m = 0; m < kM; ++m) d.q[n][m] = std::cos(2 * kPi * (m + rot) / kM);
            std::shuffle(d.q[n].begin(), d.q[n].end(), rng);
        }
    }
    return d;
}

// ---------------------------------------------------------------------------
// ScmLink
// ---------------------------------------------------------------------------
ScmLink::ScmLink(const ScmConfig& cfg, int nTx, int nRx, std::shared_ptr<const a2g::Trajectory> tx,
                 std::shared_ptr<const a2g::Trajectory> rx)
    : m_cfg(cfg), m_nTx(std::max(1, nTx)), m_nRx(std::max(1, nRx)), m_tx(std::move(tx)),
      m_rx(std::move(rx)), m_rng(cfg.seed) {
    if (!(m_cfg.fs > 0)) throw std::invalid_argument("fs must be > 0");
    m_packetLen = std::max<uint64_t>(1, static_cast<uint64_t>(std::llround(m_cfg.updateS * m_cfg.fs)));
    m_dt = static_cast<double>(m_packetLen) / m_cfg.fs;
    m_pairs.resize(static_cast<size_t>(m_nTx * m_nRx));
    for (auto& p : m_pairs) p.hist.assign(static_cast<size_t>(m_cfg.maxDelaySamples) + 1, cf(0, 0));
}

double ScmLink::sigmaAt(double alt, double speed) const {
    const auto& A = m_cfg.altTab;
    const auto& S = m_cfg.sigTab;
    double s = S.empty() ? 0.0 : S.front();
    if (A.size() == S.size() && A.size() >= 2) {
        if (alt <= A.front()) s = S.front();
        else if (alt >= A.back()) s = S.back();
        else
            for (size_t i = 1; i < A.size(); ++i)
                if (alt <= A[i]) {
                    s = S[i - 1] + (alt - A[i - 1]) / (A[i] - A[i - 1]) * (S[i] - S[i - 1]);
                    break;
                }
    }
    s *= m_cfg.fc / m_cfg.widthRefFc;
    if (m_cfg.speedScaling) s *= std::max(speed, 0.0) / m_cfg.widthRefSpeed;
    return s;
}

void ScmLink::setTime(double t0) {
    m_t0 = t0;
    m_pk.clear();
    m_nextK = 0;
    for (auto& p : m_pairs) {
        p.g = 0;
        std::fill(p.hist.begin(), p.hist.end(), cf(0, 0));
    }
}

double ScmLink::pairTime(size_t u, size_t s) const {
    return m_t0 + static_cast<double>(m_pairs[u * m_nTx + s].g) / m_cfg.fs;
}

void ScmLink::generateNext() {
    const uint64_t k = m_nextK++;
    const double t = m_t0 + static_cast<double>(k) * m_dt;
    const double lam = a2g::kC0 / m_cfg.fc;
    const double kw = 2 * kPi / lam;
    const a2g::KinState tx = m_tx->at(t), rx = m_rx->at(t);
    const a2g::KinState tx1 = m_tx->at(t + m_dt), rx1 = m_rx->at(t + m_dt);
    if (!m_started) m_groundIsTx = tx.pos.z <= rx.pos.z;
    const a2g::KinState& G = m_groundIsTx ? tx : rx;
    const a2g::KinState& A = m_groundIsTx ? rx : tx;
    const Vec3 rel = A.pos - G.pos;
    const double d0 = std::max(1e-3, a2g::norm(rel));
    const double d1 = std::max(1e-3, a2g::norm(tx1.pos - rx1.pos));
    const Vec3 uG2A = (1.0 / d0) * rel;          // ground -> air
    const Vec3 uA2G = (-1.0) * uG2A;             // air -> ground
    const double azG = std::atan2(uG2A.y, uG2A.x);
    const double azA = std::atan2(uA2G.y, uA2G.x);
    if (!m_axisSet) {
        m_axisG = m_cfg.ground.autoAxis ? azG + kPi / 2 : m_cfg.ground.axisAzDeg * kDeg;
        m_axisA = m_cfg.air.autoAxis ? azA + kPi / 2 : m_cfg.air.axisAzDeg * kDeg;
        m_axisSet = true;
    }
    const Vec3 v = A.vel - G.vel;
    const double speed = a2g::norm(v);
    const double fm = std::max(speed, m_cfg.speedFloor) / lam;
    const double alt = A.pos.z - m_cfg.groundZ;
    const double sigma = sigmaAt(alt, speed);
    const Doppler dop = m_cfg.doppler;

    // ---- LOS ray -------------------------------------------------------------
    const bool losScen = (m_cfg.scenario == Scenario::A2G || m_cfg.scenario == Scenario::UMI_LOS);
    const bool losOn = (dop == Doppler::TDL_D) || (losScen && dop != Doppler::JAKES);
    double kDb = m_cfg.kDb;
    if (m_cfg.kFromDoppler && dop == Doppler::TDL_D) kDb = 13.3;
    if (m_cfg.scenario == Scenario::UMI_LOS && dop != Doppler::TDL_D) kDb = 13.0 - 0.03 * d0;
    m_kLin = losOn ? std::pow(10.0, kDb / 10.0) : 0.0;
    double fLos = 0;
    if (dop == Doppler::TDL_D) {
        fLos = fm * std::cos(45.0 * kDeg);
        if (!m_started) m_psiLos = 0;
        else m_psiLos = wrap(m_psiLos + 2 * kPi * m_info.fLosHz * m_dt);
    } else {
        fLos = -(d1 - d0) / (lam * m_dt);       // midpoint rule (exact at boundaries)
        m_psiLos = wrap(-2 * kPi * d0 / lam);
    }
    const Vec3 axG{std::cos(m_axisG), std::sin(m_axisG), 0.0};
    const Vec3 axA{std::cos(m_axisA), std::sin(m_axisA), 0.0};
    // co-located cross-polarized pairs have no spatial phase between elements
    const double dG = m_cfg.ground.xpol ? 0.0 : m_cfg.ground.spacingWl * lam;
    const double dA = m_cfg.air.xpol ? 0.0 : m_cfg.air.spacingWl * lam;
    const bool xp = m_cfg.ground.xpol || m_cfg.air.xpol;
    // polarization index (0 = V, 1 = H) of a ground / air element
    auto polG = [&](int e) { return m_cfg.ground.xpol ? (e & 1) : 0; };
    auto polA = [&](int e) { return m_cfg.air.xpol ? (e & 1) : 0; };
    const double xpdA = std::pow(10.0, -m_cfg.xpdDb / 20.0);
    const int nG = m_groundIsTx ? m_nTx : m_nRx;
    const int nA = m_groundIsTx ? m_nRx : m_nTx;

    Packet P;
    P.k = k;
    P.t = t;
    P.psiLos = m_psiLos;
    P.fLos = fLos;
    P.aLos.assign(m_nRx, std::vector<cd>(m_nTx, cd(0, 0)));
    const double aL = std::sqrt(m_kLin / (m_kLin + 1.0));
    if (losOn) {
        const double pg = kw * dG * a2g::dot(uG2A, axG);
        const double pa = kw * dA * a2g::dot(uA2G, axA);
        for (int u = 0; u < m_nRx; ++u)
            for (int s = 0; s < m_nTx; ++s) {
                const int eG = m_groundIsTx ? s : u, eA = m_groundIsTx ? u : s;
                // LOS polarization matrix diag(1, -1) (TR 25.996 eq. 5.11)
                double pol = 1.0;
                if (xp) {
                    const int pr = polA(eA), pt = polG(eG);
                    pol = (pr != pt) ? 0.0 : (pr == 1 ? -1.0 : 1.0);
                }
                P.aLos[u][s] = aL * pol * std::polar(1.0, pg * eG + pa * eA);
            }
    }

    // ---- drops ------------------------------------------------------------------
    if (m_drops.empty() ||
        (m_drops.size() == 1 && a2g::norm(A.pos - m_drops.back().posStart) > m_cfg.dropDistanceM)) {
        Drop d = drawDrop(m_cfg, dop, m_rng);
        d.tStart = t;
        d.posStart = A.pos;
        m_drops.push_back(std::move(d));
        ++m_info.nDrops;
    }
    std::vector<double> w(m_drops.size(), 1.0);
    if (m_drops.size() == 2) {
        const double a = (m_cfg.crossfadeS > 0)
                             ? std::clamp((t - m_drops[1].tStart) / m_cfg.crossfadeS, 0.0, 1.0)
                             : 1.0;
        w[0] = std::cos(kPi / 2 * a);
        w[1] = std::sin(kPi / 2 * a);
        if (a >= 1.0) {
            m_drops.erase(m_drops.begin());
            w.erase(w.begin());
        }
    }
    int lMax = 0;
    for (const auto& d : m_drops)
        for (double tau : d.tau)
            lMax = std::max(lMax, std::min(m_cfg.maxDelaySamples, static_cast<int>(std::lround(tau * m_cfg.fs))));
    P.h.assign(m_nRx, std::vector<std::vector<cd>>(m_nTx, std::vector<cd>(lMax + 1, cd(0, 0))));
    const double aN = std::sqrt(1.0 / (m_kLin + 1.0));
    std::vector<cd> gPh(nG), aPh(nA);
    for (size_t di = 0; di < m_drops.size(); ++di) {
        Drop& d = m_drops[di];
        const int N = static_cast<int>(d.tau.size());
        for (int n = 0; n < N; ++n) {
            const int l = std::min(m_cfg.maxDelaySamples, static_cast<int>(std::lround(d.tau[n] * m_cfg.fs)));
            const double amp = w[di] * aN * std::sqrt(d.P[n] / kM);
            for (int m = 0; m < kM; ++m) {
                const double aod = azG + d.dAod[n] + d.sAod[n][m];
                const double aoa = azA + d.dAoa[n] + d.sAoa[n][m];
                const cd c = amp * std::polar(1.0, d.Phi[n][m] + d.acc[n][m]);
                const double pg = kw * dG * std::cos(aod - m_axisG);
                const double pa = kw * dA * std::cos(aoa - m_axisA);
                for (int e = 0; e < nG; ++e) gPh[e] = std::polar(1.0, pg * e);
                for (int e = 0; e < nA; ++e) aPh[e] = std::polar(1.0, pa * e);
                for (int u = 0; u < m_nRx; ++u)
                    for (int s = 0; s < m_nTx; ++s) {
                        const int eG = m_groundIsTx ? s : u, eA = m_groundIsTx ? u : s;
                        if (!xp) {
                            P.h[u][s][l] += c * gPh[eG] * aPh[eA];
                        } else {
                            const int pr = polA(eA), pt = polG(eG);
                            const int idx = pr * 2 + pt;  // vv, vh, hv, hh
                            const double g = (pr == pt) ? 1.0 : xpdA;
                            P.h[u][s][l] += c * g * std::polar(1.0, d.PhiPol[n][m][idx] - d.Phi[n][m]) *
                                            gPh[eG] * aPh[eA];
                        }
                    }
                // Doppler of this sub-path (held for the packet)
                double f = 0;
                switch (dop) {
                    case Doppler::SCM: {
                        const Vec3 ua{std::cos(aoa), std::sin(aoa), 0.0};
                        f = a2g::dot(v, ua) / lam;
                        break;
                    }
                    case Doppler::PROPOSED: f = sigma * d.q[n][m]; break;
                    case Doppler::JAKES:
                    case Doppler::TDL_D: f = fm * d.q[n][m]; break;
                }
                d.acc[n][m] = wrap(d.acc[n][m] + 2 * kPi * f * m_dt);
            }
        }
    }
    for (int l = 0; l <= lMax; ++l) {
        bool nz = false;
        for (int u = 0; u < m_nRx && !nz; ++u)
            for (int s = 0; s < m_nTx && !nz; ++s) nz = std::abs(P.h[u][s][l]) > 0;
        if (nz) P.nz.push_back(l);
    }
    m_info.t = t;
    m_info.dist = d0;
    m_info.elevDeg = std::asin(std::clamp(uG2A.z, -1.0, 1.0)) / kDeg;
    m_info.speed = speed;
    m_info.alt = alt;
    m_info.fLosHz = fLos;
    m_info.sigmaHz = sigma;
    m_info.fmHz = fm;
    m_info.kDb = losOn ? kDb : -300.0;
    m_info.groundIsTx = m_groundIsTx;
    m_started = true;
    m_pk.push_back(std::move(P));
}

const ScmLink::Packet& ScmLink::packet(uint64_t k) {
    if (!m_pk.empty() && k < m_pk.front().k)
        throw std::runtime_error("scm: packet already discarded");
    while (m_pk.empty() || m_pk.back().k < k) generateNext();
    return m_pk[static_cast<size_t>(k - m_pk.front().k)];
}

Cir ScmLink::cirAtPacket(uint64_t k) {
    const Packet& P = packet(k);
    Cir h = P.h;
    const cd rot = std::polar(1.0, P.psiLos);
    for (int u = 0; u < m_nRx; ++u)
        for (int s = 0; s < m_nTx; ++s) {
            if (h[u][s].empty()) h[u][s].assign(1, cd(0, 0));
            h[u][s][0] += P.aLos[u][s] * rot;
        }
    // stand-alone use: keep the timeline short
    while (m_pk.size() > 8 && m_pk.front().k + 4 < k) m_pk.pop_front();
    return h;
}

std::vector<std::vector<cd>> ScmLink::narrowbandAtPacket(uint64_t k) {
    const Cir h = cirAtPacket(k);
    std::vector<std::vector<cd>> H(m_nRx, std::vector<cd>(m_nTx, cd(0, 0)));
    for (int u = 0; u < m_nRx; ++u)
        for (int s = 0; s < m_nTx; ++s)
            for (const auto& v : h[u][s]) H[u][s] += v;
    return H;
}

void ScmLink::process(size_t u, size_t s, const cf* in, cf* out, size_t n) {
    if (n == 0) return;
    PairState& ps = m_pairs[u * m_nTx + s];
    const size_t hl = ps.hist.size();
    thread_local std::vector<cf> ext;
    ext.resize(hl + n);
    std::copy(ps.hist.begin(), ps.hist.end(), ext.begin());
    std::copy(in, in + n, ext.begin() + hl);
    const double twoPiOverFs = 2 * kPi / m_cfg.fs;
    size_t i = 0;
    while (i < n) {
        const uint64_t g = ps.g + i;
        const uint64_t k = g / m_packetLen;
        const uint64_t off = g % m_packetLen;
        const size_t seg = static_cast<size_t>(std::min<uint64_t>(n - i, m_packetLen - off));
        packet(k + 1);  // make sure both ends exist before taking references
        const Packet& P0 = packet(k);
        const Packet& P1 = packet(k + 1);
        // NLOS taps: linear interpolation between packet boundaries
        const auto& h0 = P0.h[u][s];
        const auto& h1 = P1.h[u][s];
        const size_t L = std::max(h0.size(), h1.size());
        thread_local std::vector<cd> hc, dh;
        hc.assign(L, cd(0, 0));
        dh.assign(L, cd(0, 0));
        const double inv = 1.0 / static_cast<double>(m_packetLen);
        std::vector<int> taps;
        for (size_t l = 0; l < L; ++l) {
            const cd a = l < h0.size() ? h0[l] : cd(0, 0);
            const cd b = l < h1.size() ? h1[l] : cd(0, 0);
            if (a == cd(0, 0) && b == cd(0, 0)) continue;
            dh[l] = (b - a) * inv;
            hc[l] = a + dh[l] * static_cast<double>(off);
            taps.push_back(static_cast<int>(l));
        }
        cd los = P0.aLos[u][s] * std::polar(1.0, P0.psiLos + twoPiOverFs * P0.fLos * static_cast<double>(off));
        const cd stp = std::polar(1.0, twoPiOverFs * P0.fLos);
        const bool losOn = std::abs(P0.aLos[u][s]) > 0;
        for (size_t q = 0; q < seg; ++q) {
            const size_t base = hl + i + q;
            const cf x0 = ext[base];
            cd acc = losOn ? los * cd(x0.real(), x0.imag()) : cd(0, 0);
            for (int l : taps) {
                const cf x = ext[base - static_cast<size_t>(l)];
                acc += hc[l] * cd(x.real(), x.imag());
                hc[l] += dh[l];
            }
            if (losOn) los *= stp;
            out[i + q] = cf(static_cast<float>(acc.real()), static_cast<float>(acc.imag()));
        }
        i += seg;
    }
    ps.g += n;
    std::copy(ext.end() - static_cast<long>(hl), ext.end(), ps.hist.begin());
    // discard packets that no pair needs any more
    uint64_t kMin = UINT64_MAX;
    for (const auto& p : m_pairs) kMin = std::min(kMin, p.g / m_packetLen);
    while (m_pk.size() > 2 && m_pk.front().k + 1 < kMin) m_pk.pop_front();
}

// ---------------------------------------------------------------------------
// analysis helpers
// ---------------------------------------------------------------------------
static cd det(std::vector<std::vector<cd>> A) {
    const size_t n = A.size();
    cd d(1, 0);
    for (size_t c = 0; c < n; ++c) {
        size_t p = c;
        for (size_t r = c + 1; r < n; ++r)
            if (std::abs(A[r][c]) > std::abs(A[p][c])) p = r;
        if (std::abs(A[p][c]) == 0) return cd(0, 0);
        if (p != c) {
            std::swap(A[p], A[c]);
            d = -d;
        }
        d *= A[c][c];
        for (size_t r = c + 1; r < n; ++r) {
            const cd f = A[r][c] / A[c][c];
            for (size_t k = c; k < n; ++k) A[r][k] -= f * A[c][k];
        }
    }
    return d;
}

double capacityBpsHz(const std::vector<std::vector<cd>>& H, double snr) {
    const size_t nr = H.size(), nt = H.empty() ? 0 : H[0].size();
    std::vector<std::vector<cd>> M(nt, std::vector<cd>(nt, cd(0, 0)));
    for (size_t i = 0; i < nt; ++i) {
        for (size_t j = 0; j < nt; ++j) {
            cd s(0, 0);
            for (size_t r = 0; r < nr; ++r) s += std::conj(H[r][i]) * H[r][j];
            M[i][j] = (i == j ? cd(1, 0) : cd(0, 0)) + snr / static_cast<double>(nt) * s;
        }
    }
    return std::log2(std::max(1e-300, std::abs(det(M))));
}

std::vector<double> singularValues2(const std::vector<std::vector<cd>>& H) {
    // eigenvalues of the 2x2 Hermitian H^H H
    const size_t nr = H.size();
    cd a(0, 0), b(0, 0), d(0, 0);
    for (size_t r = 0; r < nr; ++r) {
        a += std::norm(H[r][0]);
        d += std::norm(H[r][1]);
        b += std::conj(H[r][0]) * H[r][1];
    }
    const double tr = a.real() + d.real();
    const double dt = a.real() * d.real() - std::norm(b);
    const double disc = std::sqrt(std::max(0.0, tr * tr / 4 - dt));
    return {std::sqrt(std::max(0.0, tr / 2 + disc)), std::sqrt(std::max(0.0, tr / 2 - disc))};
}

}  // namespace scm
