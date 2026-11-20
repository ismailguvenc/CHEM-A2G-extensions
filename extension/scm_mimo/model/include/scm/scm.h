// SPDX-License-Identifier: AGPL-3.0-or-later
// 3GPP TR 25.996 Spatial Channel Model (SCM) for CHEM MIMO links
// (extension "scm_mimo", Nov 2026).
//
// One link = one ground node (fixed, "BS" in SCM terms) and one air node
// (UAV, "MS"). Each node has a uniform linear array (ULA). A drop draws
// N paths (clusters) with M = 20 sub-paths each (TR 25.996 Sec. 5.3):
//
//   h_{u,s,n}(t) = sqrt(P_n / M) sum_m exp(j(Phi_nm + k d_s sin(th_AoD,nm)
//                  + k d_u sin(th_AoA,nm) + phi_nm(t)))           (eq. 5.4)
//
// plus an optional LOS ray (K-factor split, eq. 5.5). The LOS ray uses the
// true 3-D direction from the trajectories, and its phase follows the path
// length exactly, as in the a2g_doppler model. Cluster angles are offsets
// from the current LOS direction, so the spatial structure turns with the
// geometry along a flight.
//
// The sub-path phase phi_nm(t) is the Doppler process. It is selectable, so
// the Doppler models of a2g_doppler and the MIMO model are coupled:
//   SCM     : f_nm = v . u_AoA,nm / lambda   (native SCM, eq. 5.4)
//   PROPOSED: f_nm = sigma_f(alt, v) q_nm    (narrow Gaussian diffuse
//             spectrum of the AERPAW letter; q_nm unit-rms quantiles);
//             LOS line at the geometric Doppler, K = 8.5 dB
//   JAKES   : f_nm = fm cos(alpha_nm), alpha uniform (isotropic), no LOS
//   TDL_D   : LOS line at fm cos(45 deg), K = 13.3 dB, Jakes diffuse
// The spatial signature (angles, array phases) is the same for all four.
//
// A new drop is drawn after dropDistanceM of travel and is cross-faded
// over crossfadeS (no phase jump).
//
// No CHEM dependency (uses the a2g model library for geometry).
#pragma once

#include <array>
#include <complex>
#include <cstdint>
#include <deque>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "a2g/model.h"

namespace scm {

using cf = std::complex<float>;
using cd = std::complex<double>;
using a2g::Vec3;

enum class Scenario { A2G, SMA, UMA, UMI, UMI_LOS };
enum class Doppler { SCM, PROPOSED, JAKES, TDL_D };
Scenario scenarioFromString(const std::string& s);
Doppler dopplerFromString(const std::string& s);
std::string toString(Scenario s);
std::string toString(Doppler d);

struct ArrayCfg {
    int n{2};               // elements
    double spacingWl{0.5};  // element spacing (wavelengths)
    bool autoAxis{true};    // axis perpendicular to the first LOS azimuth
    double axisAzDeg{0.0};  // ENU azimuth of the array axis (0 = east, CCW)
    // Cross-polarized pair (TR 25.996 Sec. 5.5.1 polarized option): n = 2
    // co-located elements, element 0 vertical, element 1 horizontal.
    bool xpol{false};
};

struct ScmConfig {
    double fc{3.4e9};
    double fs{1.92e6};
    Scenario scenario{Scenario::A2G};
    Doppler doppler{Doppler::PROPOSED};
    int nPaths{6};
    double xpdDb{8.0};  // NLOS cross-polarization discrimination (SCM default)
    double kDb{8.5};   // LOS K (A2G); TDL_D uses 13.3 dB, UMI_LOS 13 - 0.03 d
    bool kFromDoppler{true};  // TDL_D/JAKES set K; false: kDb is kept
    ArrayCfg ground{};
    ArrayCfg air{};
    double updateS{0.5e-3};   // coefficient packet; linear interpolation
    double dropDistanceM{50.0};
    double crossfadeS{1.0};
    int maxDelaySamples{256};
    // diffuse width table (letter, 3.51 GHz, 5 m/s), as a2g_doppler
    std::vector<double> altTab{30, 50, 70, 90, 110};
    std::vector<double> sigTab{5.0, 4.5, 4.0, 3.1, 3.2};
    double widthRefFc{3.51e9};
    double widthRefSpeed{5.0};
    bool speedScaling{true};
    double speedFloor{0.5};
    double groundZ{0.0};
    uint64_t seed{1234};
};

// One SCM drop (TR 25.996 steps 3-10). Angles in radians, relative to the
// LOS direction at each end.
struct Drop {
    double dsS{0}, asDeg{0};
    std::vector<double> tau, P;                 // [N]
    std::vector<double> dAod, dAoa;             // [N] path offsets (rad)
    std::vector<std::vector<double>> sAod, sAoa;// [N][M] sub-path offsets (rad)
    std::vector<std::vector<double>> Phi;       // [N][M] initial phases
    std::vector<std::vector<std::array<double, 4>>> PhiPol;  // [N][M] vv, vh, hv, hh phases
    std::vector<std::vector<double>> q;         // [N][M] unit Doppler (PROPOSED: rms 1; JAKES: cos)
    std::vector<std::vector<double>> acc;       // [N][M] Doppler phase accumulators
    double tStart{0};
    Vec3 posStart{};
};

Drop drawDrop(const ScmConfig& c, Doppler d, std::mt19937_64& rng);

// Dense CIR of all port pairs at one time: h[u][s][l] (u: RX port, s: TX port)
using Cir = std::vector<std::vector<std::vector<cd>>>;

struct LinkInfo {
    double t{0}, dist{0}, elevDeg{0}, speed{0}, alt{0};
    double fLosHz{0}, sigmaHz{0}, fmHz{0}, kDb{0};
    int nDrops{0};
    bool groundIsTx{true};
};

class ScmLink {
   public:
    // tx/rx: trajectories of the link source and destination nodes.
    ScmLink(const ScmConfig& cfg, int nTx, int nRx, std::shared_ptr<const a2g::Trajectory> tx,
            std::shared_ptr<const a2g::Trajectory> rx);

    // Stream processing for one port pair; y = sum_l h_{u,s}[l](t) x[n-l].
    // Every pair keeps its own sample counter and input history; all pairs
    // read the same coefficient timeline, so they stay coherent.
    void process(size_t u, size_t s, const cf* in, cf* out, size_t n);
    double pairTime(size_t u, size_t s) const;
    void setTime(double t0);  // re-anchor all pairs (clears the timeline)

    // CIR at packet boundary k (t = t0 + k*updateS), incl. the LOS ray.
    Cir cirAtPacket(uint64_t k);
    // Narrowband (sum over taps) channel matrix at packet k: H[u][s].
    std::vector<std::vector<cd>> narrowbandAtPacket(uint64_t k);
    LinkInfo info() const { return m_info; }
    const ScmConfig& config() const { return m_cfg; }
    int nTx() const { return m_nTx; }
    int nRx() const { return m_nRx; }
    uint64_t packetLen() const { return m_packetLen; }

   private:
    struct Packet {
        uint64_t k{0};
        double t{0};
        // NLOS dense taps h[u][s][l] at the packet start
        Cir h;
        // LOS: per-pair complex amplitude (array phases, power), phase at
        // packet start, frequency in the packet (Hz)
        std::vector<std::vector<cd>> aLos;
        double psiLos{0}, fLos{0};
        std::vector<int> nz;  // non-zero NLOS tap indices
    };
    const Packet& packet(uint64_t k);
    void generateNext();
    void steering(const Vec3& gPos, const Vec3& aPos, double& azG, double& azA) const;
    double sigmaAt(double alt, double speed) const;

    ScmConfig m_cfg;
    int m_nTx, m_nRx;
    std::shared_ptr<const a2g::Trajectory> m_tx, m_rx;
    std::mt19937_64 m_rng;
    uint64_t m_packetLen{1};
    double m_t0{0};
    double m_dt{0};
    std::deque<Packet> m_pk;      // timeline (k increasing)
    uint64_t m_nextK{0};
    std::vector<Drop> m_drops;    // active drops (1, or 2 while cross-fading)
    double m_kLin{0};
    double m_psiLos{0};
    bool m_started{false};
    bool m_groundIsTx{true};
    double m_axisG{0}, m_axisA{0};  // array axis azimuths (rad)
    bool m_axisSet{false};
    // per pair stream state
    struct PairState {
        uint64_t g{0};             // samples processed
        std::vector<cf> hist;      // input history
    };
    std::vector<PairState> m_pairs;  // index u*nTx+s
    LinkInfo m_info;
};

// Spectral efficiency helpers (analysis/tests)
double capacityBpsHz(const std::vector<std::vector<cd>>& H, double snrLin);
std::vector<double> singularValues2(const std::vector<std::vector<cd>>& H);  // 2x2 / Nx2 closed form

}  // namespace scm
