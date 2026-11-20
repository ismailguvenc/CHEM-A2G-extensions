#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Apply the marked CHEM-EXT edits for the scm_mimo extension.

Usage: python3 apply_mimo_core.py <CHEM root> [--no-sum]

Run a2g_doppler/tools/apply_core_hook.py first (scm_mimo uses the same
time-varying CIR hook and the a2g model library).

Edits (between '>>> CHEM-EXT' and '<<< CHEM-EXT' markers):
  1. CMakeLists.txt     option CHEM_EXT_SCM_MIMO + add_subdirectory
  2. coordinator.cpp    register ScmMimoExtension (under #ifdef)
  3. configs/config.json default "scm_mimo" block (enabled: false)
  4. CHEM-EXT-MIMO (skipped with --no-sum):
     channel.h          ChannelProcessParams::add_noise (default true)
     channel.cpp        no noise when add_noise is false
     intermediate.cpp   each destination port = sum over all source ports
                        (the public snapshot 6b9e58f keeps only the last
                        source port; the deployed CHEM v1.2.4 adds them)
Idempotent.
"""
import json
import os
import sys

root = sys.argv[1]
no_sum = '--no-sum' in sys.argv


def edit(rel, old, new, marker):
    p = os.path.join(root, rel)
    s = open(p, encoding='utf-8').read()
    if marker in s:
        print('skip (already applied):', rel, '|', marker)
        return
    if old not in s:
        raise SystemExit('anchor not found in %s:\n%s' % (rel, old))
    s = s.replace(old, new, 1)
    open(p, 'w', encoding='utf-8', newline='\n').write(s)
    print('edited:', rel, '|', marker)


edit('CMakeLists.txt',
     '# <<< CHEM-EXT a2g_doppler\n',
     '# <<< CHEM-EXT a2g_doppler\n'
     '# >>> CHEM-EXT scm_mimo: 3GPP SCM MIMO extension (needs a2g_doppler)\n'
     'option(CHEM_EXT_SCM_MIMO "Build the scm_mimo SCM MIMO extension" ON)\n'
     'if(CHEM_EXT_SCM_MIMO AND CHEM_EXT_A2G_DOPPLER AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/extensions/scm_mimo/CMakeLists.txt)\n'
     '  add_subdirectory(extensions/scm_mimo)\n'
     'endif()\n'
     '# <<< CHEM-EXT scm_mimo\n',
     'CHEM-EXT scm_mimo')

edit('src/chem/channel/coordinator.cpp',
     '// <<< CHEM-EXT a2g_doppler\n',
     '// <<< CHEM-EXT a2g_doppler\n'
     '// >>> CHEM-EXT scm_mimo\n'
     '#ifdef CHEM_EXT_SCM_MIMO\n'
     '#include "scm_mimo.h"\n'
     '#endif\n'
     '// <<< CHEM-EXT scm_mimo\n',
     '#include "scm_mimo.h"')
edit('src/chem/channel/coordinator.cpp',
     '        std::make_unique<chem::A2GDopplerExtension>());\n#endif\n',
     '        std::make_unique<chem::A2GDopplerExtension>());\n#endif\n'
     '#ifdef CHEM_EXT_SCM_MIMO\n'
     '    m_extensionRegistry.registerExtension(\n'
     '        std::make_unique<chem::ScmMimoExtension>());\n'
     '#endif\n',
     'ScmMimoExtension>()')

cp = os.path.join(root, 'configs/config.json')
txt = open(cp, encoding='utf-8').read()
if '"scm_mimo"' not in txt:
    block = {"enabled": False, "mode": "hook", "scenario": "A2G", "doppler": "PROPOSED",
             "K_dB": 8.5, "nPaths": 6, "updateMs": 0.5,
             "groundArray": {"elements": 2, "spacingWl": 0.5, "axisAzDeg": "auto"},
             "airArray": {"elements": 2, "spacingWl": 0.5, "axisAzDeg": "auto"},
             "dropDistanceM": 50, "crossfadeS": 1.0,
             "referenceOrigin": {"lat": 35.72750947, "lon": -78.69595819, "alt": 0.0},
             "links": [], "trajectoryDelayMs": 200, "seed": 1234}
    NL = chr(10)
    body = json.dumps(block, indent=4).replace(NL, NL + '        ')
    idx = txt.rstrip().rfind(NL + '    }')
    new_txt = txt[:idx] + ',' + NL + '        "scm_mimo": ' + body + txt[idx:]
    json.loads(new_txt)
    open(cp, 'w', encoding='utf-8', newline=NL).write(new_txt)
    print('edited: configs/config.json | scm_mimo')

if no_sum:
    sys.exit(0)

edit('include/chem/channel/channel.h',
     '    ChannelProcessTiming* timing{\n        nullptr};  // Optional timing collection (nullptr = disabled)\n};',
     '    ChannelProcessTiming* timing{\n        nullptr};  // Optional timing collection (nullptr = disabled)\n'
     '    // >>> CHEM-EXT-MIMO: false for the 2nd..Nth source port of a MIMO sum\n'
     '    bool add_noise{true};\n'
     '    // <<< CHEM-EXT-MIMO\n};',
     'CHEM-EXT-MIMO')
edit('src/chem/channel/channel.cpp',
     '    dsp::utils::fc32_mul_scalar_sum(dest, rx_scale, noise, dest,\n'
     '                                    params.s_perBuff);',
     '    // >>> CHEM-EXT-MIMO: the receiver noise is added once per destination port\n'
     '    if (!params.add_noise) {\n'
     '        dsp::utils::fc32_mul_scalar(dest, rx_scale, dest, params.s_perBuff);\n'
     '    } else\n'
     '    // <<< CHEM-EXT-MIMO\n'
     '    dsp::utils::fc32_mul_scalar_sum(dest, rx_scale, noise, dest,\n'
     '                                    params.s_perBuff);',
     'CHEM-EXT-MIMO')
edit('src/chem/channel/intermediate.cpp',
     '            local_params.coeff = coeff;\n'
     '            local_params.src_channel_index = c_src;\n'
     '            channel.processChannel(ctx.src_sig_arr[c_src], tl_dest_sig_arr[c_dest],\n'
     '                                   local_params);\n',
     '            local_params.coeff = coeff;\n'
     '            local_params.src_channel_index = c_src;\n'
     '            // >>> CHEM-EXT-MIMO: destination port = sum over source ports\n'
     '            // (receiver noise once). Without this, each source port\n'
     '            // overwrites the previous one and only the last one is kept.\n'
     '            thread_local signal_v tl_mimo_tmp;\n'
     '            if (c_src == mimo_first_src) {\n'
     '                local_params.add_noise = true;\n'
     '                channel.processChannel(ctx.src_sig_arr[c_src], tl_dest_sig_arr[c_dest],\n'
     '                                       local_params);\n'
     '            } else {\n'
     '                if (tl_mimo_tmp.size() < tl_dest_sig_arr[c_dest].size())\n'
     '                    tl_mimo_tmp.resize(tl_dest_sig_arr[c_dest].size(), fc(0.0f, 0.0f));\n'
     '                local_params.add_noise = false;\n'
     '                channel.processChannel(ctx.src_sig_arr[c_src], tl_mimo_tmp, local_params);\n'
     '                for (size_t i = 0; i < dest_per_buff; ++i)\n'
     '                    tl_dest_sig_arr[c_dest][i] += tl_mimo_tmp[i];\n'
     '                local_params.add_noise = true;\n'
     '            }\n'
     '            // <<< CHEM-EXT-MIMO\n',
     'CHEM-EXT-MIMO')
edit('src/chem/channel/intermediate.cpp',
     '        local_params.dest_channel_index = c_dest;\n'
     '        for (size_t c_src = 0; c_src < sig_ctx.src_n_ch; c_src++) {\n',
     '        local_params.dest_channel_index = c_dest;\n'
     '        // >>> CHEM-EXT-MIMO: first source port with a non-zero coefficient\n'
     '        size_t mimo_first_src = 0;\n'
     '        while (mimo_first_src + 1 < sig_ctx.src_n_ch &&\n'
     '               channel.getChCoeffs(c_dest, mimo_first_src) == 0.0f)\n'
     '            ++mimo_first_src;\n'
     '        // <<< CHEM-EXT-MIMO\n'
     '        for (size_t c_src = 0; c_src < sig_ctx.src_n_ch; c_src++) {\n',
     'mimo_first_src = 0')
print('done')
