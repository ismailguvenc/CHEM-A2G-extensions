// SPDX-License-Identifier: AGPL-3.0-or-later
// A2G per-tap Doppler channel model for CHEM (extension "a2g_doppler").
//
// Model (report "AERPAW A2G Doppler Model for CHEM", v1, Nov 2026):
//   y[n] = sum_p H_p[k] * exp(j psi_p[n]) * x(n - tau_p[k] fs)
//   psi_p : continuous phase of path p; it follows the path length d_p(t):
//           psi_p(t) = -2 pi d_p(t) / lambda  (geometric Doppler, eq. (2))
//   H_p[k]: a_p (sqrt(K_p/(K_p+1)) + sqrt(1/(K_p+1)) g_p(t_k)), held for
//           one packet of Delta = updateS seconds (default 0.5 ms)
//   g_p   : independent unit-power narrow diffuse process per path
//           (Gaussian / exponential / Jakes spectrum, sum of sinusoids)
// Paths: LOS, specular ground reflection (image method), point scatterers
// (bistatic path length), and fixed taps (TDL style, constant Doppler).
//
// This file has no CHEM dependency, so the model can be unit-tested and
// used by offline tools on its own.
#pragma once

#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <memory>
#include <mutex>
#include <random>
#include <string>
#include <vector>

namespace a2g {

using cf = std::complex<float>;
using cd = std::complex<double>;

constexpr double kC0 = 299792458.0;
constexpr double kPi = 3.14159265358979323846;
// Finite sentinels: CHEM builds with -ffast-math, where inf/NaN tests are
// not reliable. K_dB >= kSpecularDb means 'specular only' (K = inf);
// K_dB <= -kSpecularDb means 'diffuse only' (K = 0, Rayleigh).
constexpr double kSpecularDb = 300.0;
constexpr double kInf = 1e300;  // any value >= kSpecularDb

struct Vec3 {
    double x{0}, y{0}, z{0};  // ENU metres (east, north, up)
};
inline Vec3 operator+(Vec3 a, Vec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3 operator-(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3 operator*(double s, Vec3 a) { return {s * a.x, s * a.y, s * a.z}; }
double norm(const Vec3& a);
double dot(const Vec3& a, const Vec3& b);

struct KinState {
    Vec3 pos;
    Vec3 vel;
};

// ---------------------------------------------------------------------------
// Trajectories
// ---------------------------------------------------------------------------
class Trajectory {
   public:
    virtual ~Trajectory() = default;
    virtual KinState at(double t) const = 0;
};

class StaticTrajectory : public Trajectory {
   public:
    explicit StaticTrajectory(Vec3 p) : m_p(p) {}
    KinState at(double) const override { return {m_p, {}}; }

   private:
    Vec3 m_p;
};

// Samples (t, pos, vel) with cubic Hermite interpolation (report Sec. 4.4).
// Outside the table the trajectory extrapolates with constant velocity.
class TableTrajectory : public Trajectory {
   public:
    void add(double t, Vec3 p, Vec3 v);
    // CSV columns: t,east,north,up,ve,vn,vu  (header line optional)
    bool loadCsv(const std::string& path);
    KinState at(double t) const override;
    size_t size() const { return m_t.size(); }

   private:
    std::vector<double> m_t;
    std::vector<Vec3> m_p, m_v;
};

// A trajectory shifted in time: at(t) = base(t - t0). Used for scripted
// flights (CSV in ENU) that start when the extension starts (v2).
class OffsetTrajectory : public Trajectory {
   public:
    OffsetTrajectory(std::shared_ptr<const Trajectory> base, double t0)
        : m_base(std::move(base)), m_t0(t0) {}
    KinState at(double t) const override { return m_base->at(t - m_t0); }

   private:
    std::shared_ptr<const Trajectory> m_base;
    double m_t0;
};

// Thread-safe trajectory fed by live position reports (CHEM node map).
// It evaluates with a fixed delay of one report interval so that it
// interpolates rather than extrapolates (report Sec. 4.4).
class LiveTrajectory : public Trajectory {
   public:
    explicit LiveTrajectory(double delayS = 0.0, size_t keep = 64)
        : m_delay(delayS), m_keep(keep) {}
    void push(double t, Vec3 p, Vec3 v, bool hasVel);
    KinState at(double t) const override;
    void setDelay(double d) { m_delay = d; }
    bool empty() const;

   private:
    mutable std::mutex m_mu;
    TableTrajectory m_tab;
    std::vector<double> m_t;
    std::vector<Vec3> m_p, m_v;
    double m_delay;
    size_t m_keep;
};

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------
enum class PathType { LOS, GROUND, SCATTERER, FIXED };
enum class Shape { NONE, GAUSS, EXP, JAKES };
enum class WidthScaling { NONE, CARRIER, CARRIER_SPEED };
enum class PhaseUpdate { MIDPOINT, START_HOLD };
// PER_TAP: each path has its own geometric Doppler (proposed).
// SHARED_LOS: all paths use the LOS Doppler (report base model, and the
//   behaviour of today's CHEM for multi-tap CIRs).
enum class DopplerMode { PER_TAP, SHARED_LOS };
// Doppler model preset (v2, Nov 2026). The preset rewrites the path list
// before the link is built (applyDopplerPreset):
//   PROPOSED: geometric phase per path + narrow diffuse part (letter model)
//   JAKES   : the configured paths keep their delays and powers, but each
//             tap is Rayleigh with a classical Jakes spectrum (fm = v/lambda)
//             and no geometric line ("Jakes on every tap")
//   TDL_A   : 3GPP TR 38.901 TDL-A (NLOS, 23 taps), Jakes on every tap
//   TDL_D   : 3GPP TR 38.901 TDL-D (LOS, 13 taps); tap 1 has a LOS line at
//             fm cos(45 deg) = 0.707 fm with K1 = 13.3 dB, all taps Jakes
// TDL delays are scaled by tdlDelaySpreadS (TR 38.901 Sec. 7.7.3).
enum class DopplerModel { PROPOSED, JAKES, TDL_A, TDL_D };

PathType pathTypeFromString(const std::string& s);
DopplerModel dopplerModelFromString(const std::string& s);
std::string toString(DopplerModel m);
Shape shapeFromString(const std::string& s);
WidthScaling widthScalingFromString(const std::string& s);
std::string toString(PathType t);
std::string toString(Shape s);

struct PathSpec {
    PathType type{PathType::LOS};
    Vec3 point{};           // SCATTERER: scatterer position (ENU m)
    double gainDb{0.0};     // SCATTERER/FIXED: power (dB) relative to LOS
    double reflCoeff{-0.3}; // GROUND: real reflection coefficient
    double delayS{0.0};     // FIXED: excess delay (s)
    double dopplerHz{0.0};  // FIXED: constant Doppler (Hz)
    bool useAoa{false};     // FIXED: f = fm cos(aoaDeg) instead of dopplerHz
    double aoaDeg{0.0};
    double kDb{kInf};       // Rician K (dB); >= kSpecularDb: specular only
    Shape shape{Shape::GAUSS};  // diffuse spectrum shape of this path
    double sigmaHz{-1.0};   // diffuse rms width; < 0: altitude table
    // false: a geometric path keeps its delay and amplitude, but its phase
    // does not follow the path length (no Doppler line). Used by JAKES.
    bool geoPhase{true};
};

struct OscImpairment {
    bool enable{false};
    double driftHzPer30min{125.0};
    double jitterHzRms{2.0};
    double jitterHoldS{0.1};
};

struct ModelConfig {
    double fc{3.4e9};
    double fs{1.92e6};
    double updateS{0.5e-3};
    PhaseUpdate phaseUpdate{PhaseUpdate::MIDPOINT};
    DopplerMode dopplerMode{DopplerMode::PER_TAP};
    double groundZ{0.0};
    std::vector<PathSpec> paths{PathSpec{}};
    // Diffuse width table (letter, 3.51 GHz, 5 m/s)
    std::vector<double> altTab{30, 50, 70, 90, 110};
    std::vector<double> sigTab{5.0, 4.5, 4.0, 3.1, 3.2};
    double widthRefFc{3.51e9};
    double widthRefSpeed{5.0};
    WidthScaling widthScaling{WidthScaling::CARRIER_SPEED};
    double kappa{1.0};       // LoS Doppler scale factor
    int nSinusoids{32};
    uint64_t seed{1234};
    bool fractionalDelay{false};
    int fracTaps{16};
    int maxDelaySamples{8192};
    bool normalizeTotal{false};  // false: LOS amplitude 1 (CHEM PL is LOS)
    double speedFloor{0.5};      // m/s; Jakes fm floor
    OscImpairment osc{};
    DopplerModel dopplerModel{DopplerModel::PROPOSED};
    double tdlDelaySpreadS{100e-9};  // TDL_A / TDL_D: rms delay spread
};

// 3GPP TR 38.901 Table 7.7.2-1 (TDL-A) and 7.7.2-4 (TDL-D): normalized delay,
// power (dB) and Rician K (dB) of each tap (K <= -kSpecularDb: Rayleigh).
struct TdlTap {
    double tauNorm;
    double powerDb;
    double kDb;
};
const std::vector<TdlTap>& tdlTable(DopplerModel m);
// Rewrites cfg.paths (and normalizeTotal for TDL) for cfg.dopplerModel.
// PROPOSED leaves the configuration as it is.
void applyDopplerPreset(ModelConfig& cfg);

// ---------------------------------------------------------------------------
// Diffuse process: sum of sinusoids with per-sinusoid phase accumulators,
// so that a time-varying width does not cause a discontinuity.
// ---------------------------------------------------------------------------
class SosProcess {
   public:
    SosProcess() = default;
    SosProcess(Shape shape, int n, std::mt19937_64& rng);
    // Advance by dt seconds with rms width sigmaHz (GAUSS/EXP) or max
    // Doppler fmHz (JAKES). Returns g at the new time.
    cd step(double dt, double sigmaHz, double fmHz);
    cd value() const;
    // Normalized spectral frequencies for the shape (unit rms width for
    // GAUSS/EXP, unit fm for JAKES).
    const std::vector<double>& unitFreqs() const { return m_unit; }

   private:
    Shape m_shape{Shape::NONE};
    std::vector<double> m_unit, m_theta;
};

// ---------------------------------------------------------------------------
// Per-path state at one packet boundary (for logs, taps mode and tests)
// ---------------------------------------------------------------------------
struct PathState {
    double lengthM{0};      // path length d_p
    double excessDelayS{0}; // (d_p - d_0)/c
    double dopplerHz{0};    // frequency used in the packet
    double geoDopplerHz{0}; // instantaneous -d'(t)/lambda
    double amp{0};          // specular amplitude a_p
    cd gain{0};             // H_p (amplitude incl. diffuse, no phase)
    double phase{0};        // psi_p at the packet start (rad)
    double sigmaHz{0};
};

struct LinkSnapshot {
    double t{0};
    double uavAltM{0};
    double uavSpeed{0};
    std::vector<PathState> paths;
};

// ---------------------------------------------------------------------------
// LinkChannel: time-varying TDL for one link (one TX port -> one RX port)
// ---------------------------------------------------------------------------
class LinkChannel {
   public:
    LinkChannel(const ModelConfig& cfg, std::shared_ptr<const Trajectory> tx,
                std::shared_ptr<const Trajectory> rx);

    // Stream processing. in/out may alias. Stream time advances by n/fs.
    void process(const cf* in, cf* out, size_t n);
    // Stream time of the next sample.
    double time() const { return m_t0 + double(m_nTot) / m_cfg.fs; }
    void setTime(double t);  // re-anchor (e.g. resync to wall clock)

    // Integer-sample CIR at time t (taps mode of the CHEM extension).
    std::vector<cf> sampledCIR(double t, size_t maxLen);
    // State at time t without advancing the stream.
    LinkSnapshot snapshot(double t) const;
    const LinkSnapshot& lastPacket() const { return m_last; }
    const ModelConfig& config() const { return m_cfg; }

    double sigmaAt(double altM, double speed) const;
    static double lambda(double fc) { return kC0 / fc; }

   private:
    struct PathRt {
        PathSpec spec;
        SosProcess sos;
        double psi{0};     // phase accumulator (rad), double precision
        double freq{0};    // Hz in current packet
        cd H{0};
        int dInt{0};
        std::vector<float> frac;  // fractional delay FIR (empty: integer)
        int fracOff{0};
    };
    void startPacket(double tk);
    double pathLength(const PathSpec& p, const KinState& tx,
                      const KinState& rx) const;
    double fixedFreq(const PathSpec& p, double speed) const;
    void uavState(const KinState& tx, const KinState& rx, double& alt,
                  double& speed) const;

    ModelConfig m_cfg;
    std::shared_ptr<const Trajectory> m_tx, m_rx;
    std::vector<PathRt> m_paths;
    std::vector<cf> m_hist;  // input history (maxDelay + frac)
    size_t m_histLen{0};
    double m_t0{0};
    uint64_t m_nTot{0};
    uint64_t m_nInPacket{0};
    uint64_t m_packetLen{1};
    bool m_started{false};
    double m_tk{0};
    double m_normScale{1.0};
    // oscillator impairment
    std::mt19937_64 m_rng;
    double m_oscPsi{0}, m_oscJit{0}, m_oscNextJit{0}, m_oscRw{0};
    LinkSnapshot m_last;
};

// Closed-form helpers used in tests and analysis
double geometricDopplerHz(const KinState& tx, const KinState& rx, double fc);
inline bool isSpecular(double kDb) { return kDb >= kSpecularDb; }
inline double kLinear(double kDb) { return kDb <= -kSpecularDb ? 0.0 : std::pow(10.0, kDb / 10.0); }
double twoRayGainDb(double dHoriz, double hTx, double hRx, double fc,
                    double gamma);  // |1 + G d1/d2 e^{-jk(d2-d1)}|^2 in dB

}  // namespace a2g
