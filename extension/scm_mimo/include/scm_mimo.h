// SPDX-License-Identifier: AGPL-3.0-or-later
// CHEM extension "scm_mimo": 3GPP TR 25.996 SCM for MIMO links, with the
// Doppler options of a2g_doppler (SCM | PROPOSED | JAKES | TDL_D).
//
// Run modes:
//   "hook" (default): sample-accurate time-varying CIR per port pair through
//          the CHEM-EXT-HOOK (same hook as a2g_doppler). The CHEM core must
//          add the contributions of all source ports at each destination
//          port; the public CHEM snapshot (6b9e58f) keeps only the last
//          source port, so apply tools/apply_mimo_core.py (CHEM-EXT-MIMO).
//          The deployed CHEM v1.2.4 already adds the ports.
//   "taps": no core change. A thread pushes the CIR of every port pair with
//          Intermediate::updateCIR() every tapUpdateMs (zero-order hold).
// The extension owns the links it attaches to: it turns the CHEM link
// Doppler off and sets the real port coefficients to 1.
#pragma once

#ifdef CHEM_EXT_SCM_MIMO

#include <atomic>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <tuple>

#include "chem/extensions/extension.h"
#include "scm/scm.h"

namespace chem {

class Channel;

scm::ScmConfig scmParseConfig(const nlohmann::json& j, const scm::ScmConfig& base = {});

class ScmMimoExtension : public ChannelExtension {
   public:
    ScmMimoExtension() = default;
    ~ScmMimoExtension() override;

    std::string name() const override { return "scm_mimo"; }
    bool onStart(const nlohmann::json& config) override;
    void onStop() override;
    bool isRunning() const override { return m_running.load(); }
    nlohmann::json handleCommand(const std::string& action,
                                 const nlohmann::json& params) override;
    nlohmann::json getStatus() const override;
    nlohmann::json getConfigSchema() const override;
    bool bypassesPathLoss() const override { return false; }

    struct LinkState {
        std::unique_ptr<scm::ScmLink> link;
        bool wallClockSync{true};
        double resyncS{0.05};
        std::mutex mu;
        uint64_t buffers{0};
        // kept to rebuild the link when CHEM runs the link at another rate
        scm::ScmConfig cfg;
        int nTx{1}, nRx{1};
        std::shared_ptr<const a2g::Trajectory> tx, rx;
        uint64_t rebuilds{0};
    };
    static std::shared_ptr<LinkState> makeLinkState(
        const scm::ScmConfig& cfg, int nTx, int nRx,
        std::shared_ptr<const a2g::Trajectory> tx,
        std::shared_ptr<const a2g::Trajectory> rx, bool wallClockSync);
    // Install the per-port-pair hook (returns false without CHEM-EXT-HOOK).
    static bool installHook(chem::Channel& ch, std::shared_ptr<LinkState> st);

   private:
    void pollLoop();
    void syncNodes(double tNow);
    void attachLinks(double tNow);
    a2g::Vec3 toEnu(double lat, double lon, double alt) const;
    std::shared_ptr<a2g::LiveTrajectory> nodeTraj(const std::string& id);
    // v2: scripted trajectory (trajectoryFiles) if one is set for the node, else live
    std::shared_ptr<const a2g::Trajectory> linkTraj(const std::string& id);
    bool linkSelected(const std::string& src, const std::string& dest) const;
    int nodeChannels(const std::string& id) const;

    std::atomic<bool> m_running{false};
    std::thread m_thread;
    mutable std::mutex m_mu;
    nlohmann::json m_rawCfg;
    scm::ScmConfig m_cfg;
    std::string m_mode{"hook"};
    int m_pollMs{20};
    double m_tapUpdateMs{1.0};
    double m_tapsFs{1.92e6};
    double m_trajDelayS{0.2};
    double m_lat0{35.72750947}, m_lon0{-78.69595819}, m_alt0{0.0};
    nlohmann::json m_links;
    std::map<std::string, std::shared_ptr<a2g::LiveTrajectory>> m_nodes;
    std::map<std::string, std::shared_ptr<const a2g::Trajectory>> m_fileTraj;  // by name or id
    std::map<std::string, a2g::Vec3> m_staticNodes;
    std::map<std::tuple<double, std::string, std::string>, std::shared_ptr<LinkState>> m_linkStates;
    std::atomic<uint64_t> m_tapPushes{0};
};

}  // namespace chem

#endif  // CHEM_EXT_SCM_MIMO
