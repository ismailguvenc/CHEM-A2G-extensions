// SPDX-License-Identifier: AGPL-3.0-or-later
// CHEM extension "a2g_doppler": per-tap A2G Doppler channel model.
#include "a2g_doppler.h"

#ifdef CHEM_EXT_A2G_DOPPLER

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

static a2g::Vec3 vecFrom(const json& j) {
    if (j.is_array() && j.size() >= 3)
        return {j[0].get<double>(), j[1].get<double>(), j[2].get<double>()};
    return {};
}

a2g::ModelConfig a2gParseModelConfig(const json& j, const a2g::ModelConfig& base) {
    a2g::ModelConfig c = base;
    if (!j.is_object()) return c;
    c.fc = j.value("fc", c.fc);
    c.fs = j.value("fs", c.fs);
    if (j.contains("updateMs")) c.updateS = j["updateMs"].get<double>() * 1e-3;
    if (j.contains("phaseUpdate"))
        c.phaseUpdate = (j["phaseUpdate"].get<std::string>() == "START_HOLD")
                            ? a2g::PhaseUpdate::START_HOLD
                            : a2g::PhaseUpdate::MIDPOINT;
    if (j.contains("dopplerMode"))
        c.dopplerMode = (j["dopplerMode"].get<std::string>() == "SHARED_LOS")
                            ? a2g::DopplerMode::SHARED_LOS
                            : a2g::DopplerMode::PER_TAP;
    c.groundZ = j.value("groundZ", c.groundZ);
    c.kappa = j.value("kappa", c.kappa);
    c.nSinusoids = j.value("nSinusoids", c.nSinusoids);
    c.seed = j.value("seed", c.seed);
    c.fractionalDelay = j.value("fractionalDelay", c.fractionalDelay);
    c.fracTaps = j.value("fracTaps", c.fracTaps);
    c.maxDelaySamples = j.value("maxDelaySamples", c.maxDelaySamples);
    c.normalizeTotal = j.value("normalizeTotal", c.normalizeTotal);
    c.widthRefFc = j.value("widthRefCarrierHz", c.widthRefFc);
    c.widthRefSpeed = j.value("widthRefSpeed", c.widthRefSpeed);
    if (j.contains("widthScaling"))
        c.widthScaling = a2g::widthScalingFromString(j["widthScaling"].get<std::string>());
    if (j.contains("diffuseWidthHz") && j["diffuseWidthHz"].is_array()) {
        c.altTab.clear();
        c.sigTab.clear();
        for (const auto& e : j["diffuseWidthHz"]) {
            c.altTab.push_back(e[0].get<double>());
            c.sigTab.push_back(e[1].get<double>());
        }
    }
    const double kLos = j.value("K_dB", 8.5);
    const std::string shape = j.value("shape", std::string("GAUSS"));
    const std::string diffTaps = j.value("taps", std::string("ALL"));
    if (j.contains("oscImpairment")) {
        const auto& o = j["oscImpairment"];
        c.osc.enable = o.value("enable", false);
        c.osc.driftHzPer30min = o.value("driftHzPer30min", c.osc.driftHzPer30min);
        c.osc.jitterHzRms = o.value("jitterHzRms", c.osc.jitterHzRms);
        c.osc.jitterHoldS = o.value("jitterHoldMs", c.osc.jitterHoldS * 1e3) * 1e-3;
    }
    if (j.contains("paths") && j["paths"].is_array()) {
        c.paths.clear();
        for (const auto& p : j["paths"]) {
            a2g::PathSpec s;
            s.type = a2g::pathTypeFromString(p.value("type", std::string("LOS")));
            s.point = vecFrom(p.value("point", json::array()));
            s.gainDb = p.value("gainDb", 0.0);
            s.reflCoeff = p.value("reflCoeff", -0.3);
            s.delayS = p.value("delayUs", 0.0) * 1e-6;
            s.dopplerHz = p.value("dopplerHz", 0.0);
            if (p.contains("aoaDeg")) {
                s.aoaDeg = p["aoaDeg"].get<double>();
                s.useAoa = true;
            }
            // K: LOS defaults to K_dB; other paths default to specular
            // unless "taps":"ALL" asks for a diffuse part on every tap.
            double kd = (s.type == a2g::PathType::LOS) ? kLos
                        : (diffTaps == "ALL" ? kLos : a2g::kInf);
            if (p.contains("K_dB")) {
                // null or "inf": specular; "-inf": diffuse only (Rayleigh)
                if (p["K_dB"].is_null() || (p["K_dB"].is_string() && p["K_dB"] == "inf")) kd = a2g::kInf;
                else if (p["K_dB"].is_string() && p["K_dB"] == "-inf") kd = -a2g::kInf;
                else kd = p["K_dB"].get<double>();
            }
            s.kDb = kd;
            s.shape = a2g::shapeFromString(p.value("shape", shape));
            if (p.contains("sigmaHz")) s.sigmaHz = p["sigmaHz"].get<double>();
            c.paths.push_back(s);
        }
    } else {
        // v2: no path list -> the default paths get the documented defaults
        // (LOS K = K_dB, default 8.5 dB, with the configured diffuse shape)
        for (auto& s : c.paths) {
            if (s.type == a2g::PathType::LOS) s.kDb = kLos;
            s.shape = a2g::shapeFromString(shape);
        }
    }
    // v2: Doppler model preset (PROPOSED | JAKES | TDL_A | TDL_D)
    if (j.contains("dopplerModel"))
        c.dopplerModel = a2g::dopplerModelFromString(j["dopplerModel"].get<std::string>());
    if (j.contains("tdlDelaySpreadNs"))
        c.tdlDelaySpreadS = j["tdlDelaySpreadNs"].get<double>() * 1e-9;
    a2g::applyDopplerPreset(c);
    return c;
}

// ---------------------------------------------------------------------------
// LinkState
// ---------------------------------------------------------------------------
a2g::LinkChannel& A2GDopplerExtension::LinkState::port(size_t d, size_t s, double fs) {
    auto key = std::make_pair(d, s);
    auto it = ports.find(key);
    if (it == ports.end() || std::abs(it->second->config().fs - fs) > 1e-6) {
        a2g::ModelConfig c = cfg;
        c.fs = fs;
        c.seed = cfg.seed + 7919 * d + 104729 * s;  // independent per port pair
        auto lc = std::make_unique<a2g::LinkChannel>(c, tx, rx);
        if (wallClockSync) lc->setTime(nowS());
        it = ports.insert_or_assign(key, std::move(lc)).first;
    }
    return *it->second;
}

std::shared_ptr<A2GDopplerExtension::LinkState> A2GDopplerExtension::makeLinkState(
    const a2g::ModelConfig& cfg, std::shared_ptr<const a2g::Trajectory> tx,
    std::shared_ptr<const a2g::Trajectory> rx, bool wallClockSync) {
    auto st = std::make_shared<LinkState>();
    st->cfg = cfg;
    st->tx = std::move(tx);
    st->rx = std::move(rx);
    st->wallClockSync = wallClockSync;
    return st;
}

bool A2GDopplerExtension::installHook(chem::Channel& ch, std::shared_ptr<LinkState> st) {
#ifdef CHEM_EXT_TV_CIR_HOOK
    ch.setTimeVaryingCIR(
        [st](signal_v& buf, size_t n, double fs, size_t dPort, size_t sPort) {
            std::lock_guard<std::mutex> g(st->mu);
            auto& lc = st->port(dPort, sPort, fs);
            if (st->wallClockSync) {
                const double tw = nowS();
                if (std::abs(lc.time() - tw) > st->resyncS) lc.setTime(tw);
            }
            lc.process(buf.data(), buf.data(), n);
            ++st->buffers;
        });
    return true;
#else
    (void)ch;
    (void)st;
    return false;
#endif
}

// ---------------------------------------------------------------------------
// Extension life cycle
// ---------------------------------------------------------------------------
A2GDopplerExtension::~A2GDopplerExtension() { onStop(); }

bool A2GDopplerExtension::onStart(const json& config) {
    if (m_running.load()) onStop();
    {
        std::lock_guard<std::mutex> g(m_mu);
        m_rawCfg = config;
        m_model = a2gParseModelConfig(config);
        m_mode = config.value("mode", std::string("hook"));
        m_pollMs = config.value("pollMs", 20);
        m_tapUpdateMs = config.value("tapUpdateMs", 1.0);
        m_tapsFs = config.value("sampleRateHz", 1.92e6);
        m_maxTaps = config.value("maxTaps", 64);
        m_trajDelayS = config.value("trajectoryDelayMs", 200.0) * 1e-3;
        m_disableCoreDoppler = config.value("disableCoreDoppler", true);
        if (config.contains("referenceOrigin")) {
            const auto& o = config["referenceOrigin"];
            m_lat0 = o.value("lat", m_lat0);
            m_lon0 = o.value("lon", m_lon0);
            m_alt0 = o.value("alt", m_alt0);
        }
        if (config.contains("groundAltitudeM"))
            m_model.groundZ = config["groundAltitudeM"].get<double>() - m_alt0;
        m_links = config.value("links", json::array());
        m_staticNodes.clear();
        if (config.contains("nodes") && config["nodes"].is_object()) {
            for (auto& [id, v] : config["nodes"].items()) {
                if (v.contains("enu")) m_staticNodes[id] = vecFrom(v["enu"]);
                else m_staticNodes[id] = toEnu(v.value("lat", m_lat0), v.value("lon", m_lon0),
                                                v.value("alt", m_alt0));
            }
        }
#ifndef CHEM_EXT_TV_CIR_HOOK
        if (m_mode == "hook") {
            LOG_WARN("A2G_DOPPLER", "CHEM-EXT-HOOK not compiled in; using mode 'taps'");
            m_mode = "taps";
        }
#endif
        m_fileTraj.clear();
        if (config.contains("trajectoryFiles") && config["trajectoryFiles"].is_object()) {
            const double t0 = nowS();
            for (auto& [nm, path] : config["trajectoryFiles"].items()) {
                auto tab = std::make_shared<a2g::TableTrajectory>();
                if (!tab->loadCsv(path.get<std::string>())) {
                    LOG_WARN("A2G_DOPPLER", fmt::format("cannot read trajectory {}", path.get<std::string>()));
                    continue;
                }
                m_fileTraj[nm] = std::make_shared<a2g::OffsetTrajectory>(tab, t0);
            }
        }
        m_linkStates.clear();
        m_nodes.clear();
    }
    m_running.store(true);
    m_thread = std::thread(&A2GDopplerExtension::pollLoop, this);
    LOG_INFO("A2G_DOPPLER",
             fmt::format("started: mode={} paths={} updateMs={} K_LOS={} dB",
                         m_mode, m_model.paths.size(), m_model.updateS * 1e3,
                         m_model.paths.empty() ? 0.0 : m_model.paths[0].kDb));
    return true;
}

void A2GDopplerExtension::onStop() {
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
                ch.clearChTaps({0, 0});
                ch.clearActiveExtension();
            }
        });
    }
    LOG_INFO("A2G_DOPPLER", "stopped");
}

a2g::Vec3 A2GDopplerExtension::toEnu(double lat, double lon, double alt) const {
    // Same flat-earth conversion as Coordinator::configChannels()
    constexpr double mPerDeg = 111319.5;
    const double latAvg = 0.5 * (lat + m_lat0) * a2g::kPi / 180.0;
    return {(lon - m_lon0) * mPerDeg * std::cos(latAvg), (lat - m_lat0) * mPerDeg, alt - m_alt0};
}

std::shared_ptr<a2g::LiveTrajectory> A2GDopplerExtension::nodeTraj(const std::string& id) {
    auto it = m_nodes.find(id);
    if (it != m_nodes.end()) return it->second;
    auto t = std::make_shared<a2g::LiveTrajectory>(m_trajDelayS);
    m_nodes[id] = t;
    return t;
}

std::shared_ptr<const a2g::Trajectory> A2GDopplerExtension::linkTraj(const std::string& id) {
    auto f = m_fileTraj.find(id);
    if (f == m_fileTraj.end() && m_nodeMap) {
        auto n = m_nodeMap->find(id);
        if (n != m_nodeMap->end()) f = m_fileTraj.find(n->second.getConfig().getName());
    }
    if (f != m_fileTraj.end()) return f->second;
    return nodeTraj(id);
}

bool A2GDopplerExtension::linkSelected(const std::string& src, const std::string& dest) const {
    if (!m_links.is_array() || m_links.empty()) return true;
    for (const auto& l : m_links) {
        if (l.value("src", std::string()) == src && l.value("dest", std::string()) == dest)
            return true;
    }
    return false;
}

void A2GDopplerExtension::syncNodes(double tNow) {
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

void A2GDopplerExtension::attachLinks(double tNow) {
    if (!m_intermediateMap) return;
    for (auto& [freq, inter] : *m_intermediateMap) {
        inter->withChannelsExclusive([&](auto& channels) {
            for (auto& [key, ch] : channels) {
                const auto& [src, dest] = key;
                if (!linkSelected(src, dest)) continue;
                auto tk = std::make_tuple(freq, src, dest);
                auto it = m_linkStates.find(tk);
                if (it == m_linkStates.end()) {
                    a2g::ModelConfig c = m_model;
                    c.fc = freq;
                    auto st = makeLinkState(c, linkTraj(src), linkTraj(dest), true);
                    it = m_linkStates.emplace(tk, st).first;
                    LOG_INFO("A2G_DOPPLER", fmt::format("link {} -> {} at {} MHz attached ({})",
                                                        src, dest, freq / 1e6, m_mode));
                }
                if (m_disableCoreDoppler) ch.setDopplerEnabled(false);
                if (ch.getActiveExtension() != name()) {
                    ch.setActiveExtension(name(), bypassesPathLoss());
                    if (m_mode == "hook") installHook(ch, it->second);
                }
                if (m_mode == "taps") {
                    std::lock_guard<std::mutex> g(it->second->mu);
                    auto& lc = it->second->port(0, 0, m_tapsFs);
                    auto taps = lc.sampledCIR(tNow, m_maxTaps);
                    signal_v t(taps.begin(), taps.end());
                    ch.updateChTaps({0, 0}, t);
                    m_tapPushes.fetch_add(1);
                }
            }
        });
    }
}

void A2GDopplerExtension::pollLoop() {
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

json A2GDopplerExtension::handleCommand(const std::string& action, const json& params) {
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

json A2GDopplerExtension::getStatus() const {
    std::lock_guard<std::mutex> g(m_mu);
    json s;
    s["running"] = m_running.load();
    s["mode"] = m_mode;
    s["dopplerModel"] = a2g::toString(m_model.dopplerModel);
    s["tapPushes"] = m_tapPushes.load();
    json links = json::array();
    for (const auto& [k, st] : m_linkStates) {
        json l;
        l["freqMHz"] = std::get<0>(k) / 1e6;
        l["src"] = std::get<1>(k);
        l["dest"] = std::get<2>(k);
        std::lock_guard<std::mutex> g2(st->mu);
        l["buffers"] = st->buffers;
        if (!st->ports.empty()) {
            const auto& snap = st->ports.begin()->second->lastPacket();
            l["uavAltM"] = snap.uavAltM;
            l["uavSpeed"] = snap.uavSpeed;
            json ps = json::array();
            for (size_t i = 0; i < snap.paths.size(); ++i) {
                const auto& p = snap.paths[i];
                ps.push_back({{"type", a2g::toString(st->cfg.paths[i].type)},
                              {"dopplerHz", p.dopplerHz},
                              {"excessDelayUs", p.excessDelayS * 1e6},
                              {"ampDb", 20 * std::log10(std::max(1e-12, std::abs(p.amp)))},
                              {"sigmaHz", p.sigmaHz}});
            }
            l["paths"] = ps;
        }
        links.push_back(l);
    }
    s["links"] = links;
    return s;
}

json A2GDopplerExtension::getConfigSchema() const {
    return {
        {"mode", {{"type", "string"}, {"enum", {"hook", "taps"}}, {"default", "hook"}}},
        {"updateMs", {{"type", "number"}, {"default", 0.5}}},
        {"phaseUpdate", {{"type", "string"}, {"enum", {"MIDPOINT", "START_HOLD"}}, {"default", "MIDPOINT"}}},
        {"dopplerMode", {{"type", "string"}, {"enum", {"PER_TAP", "SHARED_LOS"}}, {"default", "PER_TAP"}}},
        {"dopplerModel", {{"type", "string"}, {"enum", {"PROPOSED", "JAKES", "TDL_A", "TDL_D"}}, {"default", "PROPOSED"}}},
        {"tdlDelaySpreadNs", {{"type", "number"}, {"default", 100.0}}},
        {"K_dB", {{"type", "number"}, {"default", 8.5}}},
        {"shape", {{"type", "string"}, {"enum", {"GAUSS", "EXP", "JAKES", "NONE"}}, {"default", "GAUSS"}}},
        {"taps", {{"type", "string"}, {"enum", {"ALL", "LOS_ONLY"}}, {"default", "ALL"}}},
        {"widthScaling", {{"type", "string"}, {"enum", {"none", "carrier", "carrier_speed"}}, {"default", "carrier_speed"}}},
        {"diffuseWidthHz", {{"type", "array"}, {"default", {{30, 5.0}, {50, 4.5}, {70, 4.0}, {90, 3.1}, {110, 3.2}}}}},
        {"paths", {{"type", "array"}, {"default", {{{"type", "LOS"}}, {{"type", "GROUND"}, {"reflCoeff", -0.3}}}}}},
        {"tapUpdateMs", {{"type", "number"}, {"default", 1.0}}},
        {"trajectoryDelayMs", {{"type", "number"}, {"default", 200}}},
        {"seed", {{"type", "integer"}, {"default", 1234}}},
    };
}

}  // namespace chem

#endif  // CHEM_EXT_A2G_DOPPLER
