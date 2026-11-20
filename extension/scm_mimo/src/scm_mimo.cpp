// SPDX-License-Identifier: AGPL-3.0-or-later
// CHEM extension "scm_mimo": 3GPP TR 25.996 SCM for MIMO links.
#include "scm_mimo.h"

#ifdef CHEM_EXT_SCM_MIMO

#include <chrono>
#include <cmath>

#include "chem/channel/channel.h"
#include "chem/channel/intermediate.h"
#include "chem/common.h"
#include "spdlog/fmt/fmt.h"

namespace chem {

using nlohmann::json;

static double nowS() {
    using namespace std::chrono;
    return duration_cast<nanoseconds>(steady_clock::now().time_since_epoch()).count() * 1e-9;
}

static scm::ArrayCfg parseArray(const json& j, scm::ArrayCfg a) {
    if (!j.is_object()) return a;
    a.n = j.value("elements", a.n);
    a.spacingWl = j.value("spacingWl", a.spacingWl);
    a.xpol = j.value("xpol", a.xpol);
    if (j.contains("axisAzDeg")) {
        if (j["axisAzDeg"].is_string()) {
            a.autoAxis = true;  // "auto"
        } else {
            a.autoAxis = false;
            a.axisAzDeg = j["axisAzDeg"].get<double>();
        }
    }
    return a;
}

scm::ScmConfig scmParseConfig(const json& j, const scm::ScmConfig& base) {
    scm::ScmConfig c = base;
    if (!j.is_object()) return c;
    c.fc = j.value("fc", c.fc);
    c.fs = j.value("fs", c.fs);
    if (j.contains("scenario")) c.scenario = scm::scenarioFromString(j["scenario"].get<std::string>());
    if (j.contains("doppler")) c.doppler = scm::dopplerFromString(j["doppler"].get<std::string>());
    c.nPaths = j.value("nPaths", c.nPaths);
    c.xpdDb = j.value("xpd_dB", c.xpdDb);
    if (j.contains("K_dB")) {
        c.kDb = j["K_dB"].get<double>();
        c.kFromDoppler = false;
    }
    c.ground = parseArray(j.value("groundArray", json()), c.ground);
    c.air = parseArray(j.value("airArray", json()), c.air);
    if (j.contains("updateMs")) c.updateS = j["updateMs"].get<double>() * 1e-3;
    c.dropDistanceM = j.value("dropDistanceM", c.dropDistanceM);
    c.crossfadeS = j.value("crossfadeS", c.crossfadeS);
    c.maxDelaySamples = j.value("maxDelaySamples", c.maxDelaySamples);
    c.speedFloor = j.value("speedFloor", c.speedFloor);
    c.widthRefFc = j.value("widthRefCarrierHz", c.widthRefFc);
    c.widthRefSpeed = j.value("widthRefSpeed", c.widthRefSpeed);
    if (j.contains("widthScaling"))
        c.speedScaling = (j["widthScaling"].get<std::string>() == "carrier_speed");
    if (j.contains("diffuseWidthHz") && j["diffuseWidthHz"].is_array()) {
        c.altTab.clear();
        c.sigTab.clear();
        for (const auto& e : j["diffuseWidthHz"]) {
            c.altTab.push_back(e[0].get<double>());
            c.sigTab.push_back(e[1].get<double>());
        }
    }
    c.seed = j.value("seed", c.seed);
    return c;
}

std::shared_ptr<ScmMimoExtension::LinkState> ScmMimoExtension::makeLinkState(
    const scm::ScmConfig& cfg, int nTx, int nRx, std::shared_ptr<const a2g::Trajectory> tx,
    std::shared_ptr<const a2g::Trajectory> rx, bool wallClockSync) {
    auto st = std::make_shared<LinkState>();
    st->cfg = cfg;
    st->nTx = nTx;
    st->nRx = nRx;
    st->tx = tx;
    st->rx = rx;
    st->link = std::make_unique<scm::ScmLink>(cfg, nTx, nRx, std::move(tx), std::move(rx));
    st->wallClockSync = wallClockSync;
    if (wallClockSync) st->link->setTime(nowS());
    return st;
}

bool ScmMimoExtension::installHook(chem::Channel& ch, std::shared_ptr<LinkState> st) {
#ifdef CHEM_EXT_TV_CIR_HOOK
    ch.setTimeVaryingCIR([st](signal_v& buf, size_t n, double fs, size_t dPort, size_t sPort) {
        std::lock_guard<std::mutex> g(st->mu);
        if (std::abs(fs - st->link->config().fs) > 1e-6 && fs > 0) {
            // CHEM runs this link at another sample rate: rebuild the link at that rate
            scm::ScmConfig c = st->cfg;
            c.fs = fs;
            st->cfg = c;
            st->link = std::make_unique<scm::ScmLink>(c, st->nTx, st->nRx, st->tx, st->rx);
            if (st->wallClockSync) st->link->setTime(nowS());
            ++st->rebuilds;
        }
        auto& L = *st->link;
        if (dPort >= static_cast<size_t>(L.nRx()) || sPort >= static_cast<size_t>(L.nTx())) return;
        if (st->wallClockSync && dPort == 0 && sPort == 0) {
            const double tw = nowS();
            if (std::abs(L.pairTime(0, 0) - tw) > st->resyncS) L.setTime(tw);
        }
        L.process(dPort, sPort, buf.data(), buf.data(), n);
        ++st->buffers;
    });
    return true;
#else
    (void)ch;
    (void)st;
    return false;
#endif
}

ScmMimoExtension::~ScmMimoExtension() { onStop(); }

bool ScmMimoExtension::onStart(const json& config) {
    if (m_running.load()) onStop();
    {
        std::lock_guard<std::mutex> g(m_mu);
        m_rawCfg = config;
        m_cfg = scmParseConfig(config);
        m_mode = config.value("mode", std::string("hook"));
        m_pollMs = config.value("pollMs", 20);
        m_tapUpdateMs = config.value("tapUpdateMs", 1.0);
        m_tapsFs = config.value("sampleRateHz", 1.92e6);
        m_trajDelayS = config.value("trajectoryDelayMs", 200.0) * 1e-3;
        if (config.contains("referenceOrigin")) {
            const auto& o = config["referenceOrigin"];
            m_lat0 = o.value("lat", m_lat0);
            m_lon0 = o.value("lon", m_lon0);
            m_alt0 = o.value("alt", m_alt0);
        }
        if (config.contains("groundAltitudeM"))
            m_cfg.groundZ = config["groundAltitudeM"].get<double>() - m_alt0;
        m_links = config.value("links", json::array());
        m_staticNodes.clear();
        if (config.contains("nodes") && config["nodes"].is_object()) {
            for (auto& [id, v] : config["nodes"].items()) {
                if (v.contains("enu"))
                    m_staticNodes[id] = {v["enu"][0].get<double>(), v["enu"][1].get<double>(),
                                         v["enu"][2].get<double>()};
                else
                    m_staticNodes[id] = toEnu(v.value("lat", m_lat0), v.value("lon", m_lon0),
                                              v.value("alt", m_alt0));
            }
        }
#ifndef CHEM_EXT_TV_CIR_HOOK
        if (m_mode == "hook") {
            LOG_WARN("SCM_MIMO", "CHEM-EXT-HOOK not compiled in; using mode 'taps'");
            m_mode = "taps";
        }
#endif
        m_fileTraj.clear();
        if (config.contains("trajectoryFiles") && config["trajectoryFiles"].is_object()) {
            const double t0 = nowS();
            for (auto& [nm, path] : config["trajectoryFiles"].items()) {
                auto tab = std::make_shared<a2g::TableTrajectory>();
                if (!tab->loadCsv(path.get<std::string>())) {
                    LOG_WARN("SCM_MIMO", fmt::format("cannot read trajectory {}", path.get<std::string>()));
                    continue;
                }
                m_fileTraj[nm] = std::make_shared<a2g::OffsetTrajectory>(tab, t0);
            }
        }
        m_linkStates.clear();
        m_nodes.clear();
    }
    m_running.store(true);
    m_thread = std::thread(&ScmMimoExtension::pollLoop, this);
    LOG_INFO("SCM_MIMO", fmt::format("started: mode={} scenario={} doppler={} arrays {}x{}", m_mode,
                                     scm::toString(m_cfg.scenario), scm::toString(m_cfg.doppler),
                                     m_cfg.ground.n, m_cfg.air.n));
    return true;
}

void ScmMimoExtension::onStop() {
    if (!m_running.exchange(false)) return;
    if (m_thread.joinable()) m_thread.join();
    if (!m_intermediateMap) return;
    for (auto& [freq, inter] : *m_intermediateMap) {
        inter->withChannelsExclusive([&](auto& channels) {
            for (auto& [key, ch] : channels) {
                if (ch.getActiveExtension() != name()) continue;
#ifdef CHEM_EXT_TV_CIR_HOOK
                ch.clearTimeVaryingCIR();
#endif
                for (const auto& [id, c] : ch.getIndChannels()) ch.clearChTaps(id);
                ch.clearActiveExtension();
            }
        });
    }
    LOG_INFO("SCM_MIMO", "stopped");
}

a2g::Vec3 ScmMimoExtension::toEnu(double lat, double lon, double alt) const {
    constexpr double mPerDeg = 111319.5;  // as Coordinator::configChannels()
    const double latAvg = 0.5 * (lat + m_lat0) * a2g::kPi / 180.0;
    return {(lon - m_lon0) * mPerDeg * std::cos(latAvg), (lat - m_lat0) * mPerDeg, alt - m_alt0};
}

std::shared_ptr<a2g::LiveTrajectory> ScmMimoExtension::nodeTraj(const std::string& id) {
    auto it = m_nodes.find(id);
    if (it != m_nodes.end()) return it->second;
    auto t = std::make_shared<a2g::LiveTrajectory>(m_trajDelayS);
    m_nodes[id] = t;
    return t;
}

std::shared_ptr<const a2g::Trajectory> ScmMimoExtension::linkTraj(const std::string& id) {
    auto f = m_fileTraj.find(id);
    if (f == m_fileTraj.end() && m_nodeMap) {
        auto n = m_nodeMap->find(id);
        if (n != m_nodeMap->end()) f = m_fileTraj.find(n->second.getConfig().getName());
    }
    if (f != m_fileTraj.end()) return f->second;
    return nodeTraj(id);
}

bool ScmMimoExtension::linkSelected(const std::string& src, const std::string& dest) const {
    if (!m_links.is_array() || m_links.empty()) return true;
    for (const auto& l : m_links)
        if (l.value("src", std::string()) == src && l.value("dest", std::string()) == dest) return true;
    return false;
}

int ScmMimoExtension::nodeChannels(const std::string& id) const {
    if (!m_nodeMap) return 1;
    auto it = m_nodeMap->find(id);
    if (it == m_nodeMap->end()) return 1;
    return std::max<int>(1, it->second.getConfig().getNumChannels());
}

void ScmMimoExtension::syncNodes(double tNow) {
    if (!m_nodeMap) return;
    for (auto& [id, node] : *m_nodeMap) {
        auto tr = nodeTraj(id);
        const std::string& nm = node.getConfig().getName();
        auto st = m_staticNodes.find(id);
        if (st == m_staticNodes.end()) st = m_staticNodes.find(nm);
        if (st != m_staticNodes.end()) {
            tr->push(tNow, st->second, {}, true);
            continue;
        }
        auto loc = node.getLocation();
        if (!loc) continue;
        const a2g::Vec3 p = toEnu(loc->position.lat, loc->position.lon, loc->position.alt);
        const a2g::Vec3 v{loc->velocity.east, loc->velocity.north, loc->velocity.up};
        tr->push(tNow, p, v, loc->hasVelocity);
    }
}

void ScmMimoExtension::attachLinks(double tNow) {
    if (!m_intermediateMap) return;
    for (auto& [freq, inter] : *m_intermediateMap) {
        inter->withChannelsExclusive([&](auto& channels) {
            for (auto& [key, ch] : channels) {
                const auto& [src, dest] = key;
                if (!linkSelected(src, dest)) continue;
                auto tk = std::make_tuple(freq, src, dest);
                auto it = m_linkStates.find(tk);
                if (it == m_linkStates.end()) {
                    scm::ScmConfig c = m_cfg;
                    c.fc = freq;
                    if (m_mode == "taps") {
                        c.fs = m_tapsFs;
                        c.updateS = m_tapUpdateMs * 1e-3;
                    }
                    const int nTx = nodeChannels(src), nRx = nodeChannels(dest);
                    auto st = makeLinkState(c, nTx, nRx, linkTraj(src), linkTraj(dest), true);
                    it = m_linkStates.emplace(tk, st).first;
                    LOG_INFO("SCM_MIMO", fmt::format("link {} -> {} at {} MHz attached ({}x{}, {})", src, dest,
                                                     freq / 1e6, nRx, nTx, m_mode));
                }
                ch.setDopplerEnabled(false);
                if (ch.getActiveExtension() != name()) {
                    ch.setActiveExtension(name(), bypassesPathLoss());
                    for (const auto& [id, c] : ch.getIndChannels()) ch.updateChCoeff(id, 1.0f);
                    if (m_mode == "hook") installHook(ch, it->second);
                }
                if (m_mode == "taps") {
                    std::lock_guard<std::mutex> g(it->second->mu);
                    auto& L = *it->second->link;
                    const uint64_t k = static_cast<uint64_t>(
                        std::max(0.0, (tNow - L.pairTime(0, 0)) / L.config().updateS));
                    const scm::Cir h = L.cirAtPacket(k);
                    for (int u = 0; u < L.nRx(); ++u)
                        for (int s = 0; s < L.nTx(); ++s) {
                            signal_v t;
                            for (const auto& v : h[u][s])
                                t.emplace_back(static_cast<float>(v.real()), static_cast<float>(v.imag()));
                            // CHEM keys port CIRs as {dest port, src port}
                            ch.updateChTaps(chId{static_cast<size_t>(u), static_cast<size_t>(s)}, t);
                        }
                    m_tapPushes.fetch_add(1);
                }
            }
        });
    }
}

void ScmMimoExtension::pollLoop() {
    double nextSync = 0;
    while (m_running.load()) {
        const double t = nowS();
        {
            std::lock_guard<std::mutex> g(m_mu);
            if (t >= nextSync) {
                syncNodes(t);
                nextSync = t + m_pollMs * 1e-3;
            }
            attachLinks(t);
        }
        const double sleepMs = (m_mode == "taps") ? m_tapUpdateMs : m_pollMs;
        std::this_thread::sleep_for(std::chrono::microseconds(static_cast<int64_t>(sleepMs * 1e3)));
    }
}

json ScmMimoExtension::handleCommand(const std::string& action, const json& params) {
    json resp;
    if (action == "start" || action == "configure") {
        json cfg = m_rawCfg.is_object() ? m_rawCfg : json::object();
        for (auto& [k, v] : params.items()) cfg[k] = v;
        resp["status"] = onStart(cfg) ? "success" : "fail";
    } else if (action == "stop") {
        onStop();
        resp["status"] = "success";
    } else if (action == "status" || action == "snapshot") {
        resp = getStatus();
        resp["status"] = "success";
    } else {
        resp["status"] = "fail";
        resp["message"] = "Unknown action: " + action;
    }
    return resp;
}

json ScmMimoExtension::getStatus() const {
    std::lock_guard<std::mutex> g(m_mu);
    json s;
    s["running"] = m_running.load();
    s["mode"] = m_mode;
    s["scenario"] = scm::toString(m_cfg.scenario);
    s["doppler"] = scm::toString(m_cfg.doppler);
    s["tapPushes"] = m_tapPushes.load();
    json links = json::array();
    for (const auto& [k, st] : m_linkStates) {
        std::lock_guard<std::mutex> g2(st->mu);
        const auto inf = st->link->info();
        links.push_back({{"freqMHz", std::get<0>(k) / 1e6},
                         {"src", std::get<1>(k)},
                         {"dest", std::get<2>(k)},
                         {"ports", fmt::format("{}x{}", st->link->nRx(), st->link->nTx())},
                         {"buffers", st->buffers},
                         {"distM", inf.dist},
                         {"elevDeg", inf.elevDeg},
                         {"speed", inf.speed},
                         {"fLosHz", inf.fLosHz},
                         {"sigmaHz", inf.sigmaHz},
                         {"K_dB", inf.kDb},
                         {"drops", inf.nDrops}});
    }
    s["links"] = links;
    return s;
}

json ScmMimoExtension::getConfigSchema() const {
    return {
        {"mode", {{"type", "string"}, {"enum", {"hook", "taps"}}, {"default", "hook"}}},
        {"scenario", {{"type", "string"}, {"enum", {"A2G", "SMA", "UMA", "UMI", "UMI_LOS"}}, {"default", "A2G"}}},
        {"doppler", {{"type", "string"}, {"enum", {"SCM", "PROPOSED", "JAKES", "TDL_D"}}, {"default", "PROPOSED"}}},
        {"K_dB", {{"type", "number"}, {"default", 8.5}}},
        {"nPaths", {{"type", "integer"}, {"default", 6}}},
        {"groundArray", {{"type", "object"}, {"default", {{"elements", 2}, {"spacingWl", 0.5}, {"axisAzDeg", "auto"}}}}},
        {"airArray", {{"type", "object"}, {"default", {{"elements", 2}, {"spacingWl", 0.5}, {"axisAzDeg", "auto"}}}}},
        {"updateMs", {{"type", "number"}, {"default", 0.5}}},
        {"dropDistanceM", {{"type", "number"}, {"default", 50}}},
        {"crossfadeS", {{"type", "number"}, {"default", 1.0}}},
        {"tapUpdateMs", {{"type", "number"}, {"default", 1.0}}},
        {"seed", {{"type", "integer"}, {"default", 1234}}},
    };
}

}  // namespace chem

#endif  // CHEM_EXT_SCM_MIMO
