# SPDX-License-Identifier: AGPL-3.0-or-later
"""Add the 'trajectoryFiles' option to a2g_doppler and scm_mimo (v2): a scripted ENU trajectory
(CSV t,east,north,up[,ve,vn,vu]) for a node, by node name or id, time 0 = extension start."""
import re

for hdr, src, cls in (('a2g_doppler/include/a2g_doppler.h', 'a2g_doppler/src/a2g_doppler.cpp', 'A2GDopplerExtension'),
                      ('scm_mimo/include/scm_mimo.h', 'scm_mimo/src/scm_mimo.cpp', 'ScmMimoExtension')):
    h = open(hdr, encoding='utf-8').read()
    if 'linkTraj' not in h:
        h = h.replace('    std::shared_ptr<a2g::LiveTrajectory> nodeTraj(const std::string& id);',
                      '    std::shared_ptr<a2g::LiveTrajectory> nodeTraj(const std::string& id);\n'
                      '    // v2: scripted trajectory (trajectoryFiles) if one is set for the node, else live\n'
                      '    std::shared_ptr<const a2g::Trajectory> linkTraj(const std::string& id);')
        h = h.replace('    std::map<std::string, std::shared_ptr<a2g::LiveTrajectory>> m_nodes;',
                      '    std::map<std::string, std::shared_ptr<a2g::LiveTrajectory>> m_nodes;\n'
                      '    std::map<std::string, std::shared_ptr<const a2g::Trajectory>> m_fileTraj;  // by name or id')
        open(hdr, 'w', encoding='utf-8', newline='\n').write(h)
    s = open(src, encoding='utf-8').read()
    if 'trajectoryFiles' not in s:
        # parse in onStart, after the static nodes block
        anchor = '        m_linkStates.clear();\n        m_nodes.clear();\n'
        assert s.count(anchor) == 1, src
        s = s.replace(anchor,
            '        m_fileTraj.clear();\n'
            '        if (config.contains("trajectoryFiles") && config["trajectoryFiles"].is_object()) {\n'
            '            const double t0 = nowS();\n'
            '            for (auto& [nm, path] : config["trajectoryFiles"].items()) {\n'
            '                auto tab = std::make_shared<a2g::TableTrajectory>();\n'
            '                if (!tab->loadCsv(path.get<std::string>())) {\n'
            '                    LOG_WARN("' + cls + '", fmt::format("cannot read trajectory {}", path.get<std::string>()));\n'
            '                    continue;\n'
            '                }\n'
            '                m_fileTraj[nm] = std::make_shared<a2g::OffsetTrajectory>(tab, t0);\n'
            '            }\n'
            '        }\n' + anchor)
        # helper after nodeTraj()
        m = re.search(r'std::shared_ptr<a2g::LiveTrajectory> ' + cls + r'::nodeTraj\(const std::string& id\) \{.*?\n\}\n', s, re.S)
        helper = ('\nstd::shared_ptr<const a2g::Trajectory> ' + cls + '::linkTraj(const std::string& id) {\n'
                  '    auto f = m_fileTraj.find(id);\n'
                  '    if (f == m_fileTraj.end() && m_nodeMap) {\n'
                  '        auto n = m_nodeMap->find(id);\n'
                  '        if (n != m_nodeMap->end()) f = m_fileTraj.find(n->second.getConfig().getName());\n'
                  '    }\n'
                  '    if (f != m_fileTraj.end()) return f->second;\n'
                  '    return nodeTraj(id);\n'
                  '}\n')
        s = s[:m.end()] + helper + s[m.end():]
        s = s.replace('nodeTraj(src), nodeTraj(dest), true)', 'linkTraj(src), linkTraj(dest), true)')
        open(src, 'w', encoding='utf-8', newline='\n').write(s)
    print(src, 'linkTraj(src)' in open(src).read())
