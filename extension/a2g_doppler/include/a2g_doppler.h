// SPDX-License-Identifier: AGPL-3.0-or-later
// CHEM extension "a2g_doppler": per-tap A2G Doppler channel model.
//
// Two run modes:
//   "hook" (default): sample-accurate time-varying CIR. Needs the small
//          CHEM-EXT-HOOK in chem/channel/channel.{h,cpp} (see patch/).
//          Each tap gets its own continuous Doppler phase per sample; tap
//          gains and frequencies refresh every updateMs (default 0.5 ms).
//   "taps": no core change at all. A thread samples the CIR every
//          tapUpdateMs and pushes it with Intermediate::updateCIR(); CHEM
//          holds the taps until the next push (zero-order hold).
// In both modes the extension turns the CHEM link Doppler (one shift per
// link) off for the links it owns, so the shift is not applied twice.
#pragma once

#ifdef CHEM_EXT_A2G_DOPPLER

#include <atomic>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <tuple>

#include "a2g/model.h"
#include "chem/extensions/extension.h"

namespace chem {

// Parse the JSON block of config.json ("extensions": {"a2g_doppler": {...}})
// or of an EXT command into a model configuration (exposed for tests).
a2g::ModelConfig a2gParseModelConfig(const nlohmann::json& j,
                                     const a2g::ModelConfig& base = {});

class A2GDopplerExtension : public ChannelExtension {
   public:
    A2GDopplerExtension() = default;
    ~A2GDopplerExtension() override;

    std::string name() const override { return "a2g_doppler"; }
    bool onStart(const nlohmann::json& config) override;
    void onStop() override;
    bool isRunning() const override { return m_running.load(); }
    nlohmann::json handleCommand(const std::string& action,
                                 const nlohmann::json& params) override;
    nlohmann::json getStatus() const override;
    nlohmann::json getConfigSchema() const override;
    bool bypassesPathLoss() const override { return false; }

    // Install the model on one channel directly (used by tests and by the
    // offline harness; the live path calls this from the poll loop).
    // Returns false if the CHEM-EXT-HOOK is not compiled in.
    struct LinkState;
    static bool installHook(chem::Channel& ch, std::shared_ptr<LinkState> st);
    static std::shared_ptr<LinkState> makeLinkState(
        const a2g::ModelConfig& cfg, std::shared_ptr<const a2g::Trajectory> tx,
        std::shared_ptr<const a2g::Trajectory> rx, bool wallClockSync);

    struct LinkState {
        a2g::ModelConfig cfg;
        std::shared_ptr<const a2g::Trajectory> tx, rx;
        bool wallClockSync{true};
        double resyncS{0.05};
        std::mutex mu;
        // one time-varying channel per (dest port, src port)
        std::map<std::pair<size_t, size_t>, std::unique_ptr<a2g::LinkChannel>> ports;
        a2g::LinkChannel& port(size_t d, size_t s, double fs);
        uint64_t buffers{0};
    };

   private:
    void pollLoop();
    void syncNodes(double tNow);
    void attachLinks(double tNow);
    a2g::Vec3 toEnu(double lat, double lon, double alt) const;
    std::shared_ptr<a2g::LiveTrajectory> nodeTraj(const std::string& id);
    // v2: scripted trajectory (trajectoryFiles) if one is set for the node, else live
    std::shared_ptr<const a2g::Trajectory> linkTraj(const std::string& id);
    bool linkSelected(const std::string& src, const std::string& dest) const;

    std::atomic<bool> m_running{false};
    std::thread m_thread;
    mutable std::mutex m_mu;
    nlohmann::json m_rawCfg;
    a2g::ModelConfig m_model;
    std::string m_mode{"hook"};
    int m_pollMs{20};
    double m_tapUpdateMs{1.0};
    double m_tapsFs{1.92e6};
    size_t m_maxTaps{64};
    double m_trajDelayS{0.2};
    bool m_disableCoreDoppler{true};
    double m_lat0{35.72750947}, m_lon0{-78.69595819}, m_alt0{0.0};
    nlohmann::json m_links;  // empty: all links
    std::map<std::string, std::shared_ptr<a2g::LiveTrajectory>> m_nodes;
    std::map<std::string, std::shared_ptr<const a2g::Trajectory>> m_fileTraj;  // by name or id
    std::map<std::string, a2g::Vec3> m_staticNodes;
    // key: (freq, src, dest)
    std::map<std::tuple<double, std::string, std::string>, std::shared_ptr<LinkState>> m_linkStates;
    std::atomic<uint64_t> m_tapPushes{0};
};

}  // namespace chem

#endif  // CHEM_EXT_A2G_DOPPLER
