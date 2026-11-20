#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Apply the minimal, clearly marked CHEM-EXT edits to a CHEM source tree.

Usage: python3 apply_core_hook.py <CHEM root> [--no-hook]

Edits (all between '>>> CHEM-EXT' and '<<< CHEM-EXT' markers):
  1. CMakeLists.txt          option CHEM_EXT_A2G_DOPPLER + add_subdirectory
  2. coordinator.cpp         register A2GDopplerExtension (under #ifdef)
  3. configs/config.json     default "a2g_doppler" block (enabled: false)
  4. channel.h / channel.cpp optional time-varying CIR hook (skipped with
                             --no-hook; the extension then runs in "taps" mode)
The script is idempotent: it does nothing if a marker is already present.
"""
import json
import os
import sys

root = sys.argv[1]
no_hook = '--no-hook' in sys.argv


def edit(rel, old, new, marker):
    p = os.path.join(root, rel)
    s = open(p, encoding='utf-8').read()
    if marker in s:
        print('skip (already applied):', rel)
        return
    if old not in s:
        raise SystemExit('anchor not found in %s:\n%s' % (rel, old))
    s = s.replace(old, new, 1)
    open(p, 'w', encoding='utf-8', newline='\n').write(s)
    print('edited:', rel)


# 1. top-level CMakeLists.txt -------------------------------------------------
edit('CMakeLists.txt',
     '# Format the project using the `clang-format` target',
     '# >>> CHEM-EXT a2g_doppler: per-tap A2G Doppler extension (own target/option)\n'
     'option(CHEM_EXT_A2G_DOPPLER "Build the a2g_doppler per-tap Doppler extension" ON)\n'
     'if(CHEM_EXT_A2G_DOPPLER AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/extensions/a2g_doppler/CMakeLists.txt)\n'
     '  add_subdirectory(extensions/a2g_doppler)\n'
     'endif()\n'
     '# <<< CHEM-EXT a2g_doppler\n\n'
     '# Format the project using the `clang-format` target',
     'CHEM-EXT a2g_doppler')

# 2. coordinator registration -------------------------------------------------
edit('src/chem/channel/coordinator.cpp',
     '#include "sionna_client.h"\n',
     '#include "sionna_client.h"\n'
     '// >>> CHEM-EXT a2g_doppler\n'
     '#ifdef CHEM_EXT_A2G_DOPPLER\n'
     '#include "a2g_doppler.h"\n'
     '#endif\n'
     '// <<< CHEM-EXT a2g_doppler\n',
     'CHEM-EXT a2g_doppler')
edit('src/chem/channel/coordinator.cpp',
     '        std::make_unique<chem::SandboxExtension>());\n}',
     '        std::make_unique<chem::SandboxExtension>());\n'
     '    // >>> CHEM-EXT a2g_doppler\n'
     '#ifdef CHEM_EXT_A2G_DOPPLER\n'
     '    m_extensionRegistry.registerExtension(\n'
     '        std::make_unique<chem::A2GDopplerExtension>());\n'
     '#endif\n'
     '    // <<< CHEM-EXT a2g_doppler\n}',
     'A2GDopplerExtension>()')

# 3. default config -------------------------------------------------------------
cp = os.path.join(root, 'configs/config.json')
txt = open(cp, encoding='utf-8').read()
if '"a2g_doppler"' not in txt:
    block = {
        "enabled": False,
        "mode": "hook",
        "updateMs": 0.5,
        "phaseUpdate": "MIDPOINT",
        "dopplerMode": "PER_TAP",
        "K_dB": 8.5,
        "shape": "GAUSS",
        "taps": "ALL",
        "diffuseWidthHz": [[30, 5.0], [50, 4.5], [70, 4.0], [90, 3.1], [110, 3.2]],
        "widthRefCarrierHz": 3.51e9, "widthRefSpeed": 5.0, "widthScaling": "carrier_speed",
        "paths": [{"type": "LOS"}, {"type": "GROUND", "reflCoeff": -0.3}],
        "referenceOrigin": {"lat": 35.72750947, "lon": -78.69595819, "alt": 0.0},
        "groundAltitudeM": 0.0,
        "trajectoryDelayMs": 200, "pollMs": 20, "tapUpdateMs": 1.0,
        "seed": 1234,
        "oscImpairment": {"enable": False}
    }
    NL = chr(10)
    body = json.dumps(block, indent=4).replace(NL, NL + '        ')
    # text insertion before the brace that closes "extensions" (4-space
    # indent), so the rest of the file keeps its format
    idx = txt.rstrip().rfind(NL + '    }')
    new_txt = txt[:idx] + ',' + NL + '        "a2g_doppler": ' + body + txt[idx:]
    json.loads(new_txt)  # must stay valid JSON
    open(cp, 'w', encoding='utf-8', newline=NL).write(new_txt)
    print('edited: configs/config.json')
else:
    print('skip (already applied): configs/config.json')

if no_hook:
    sys.exit(0)

# 4. time-varying CIR hook --------------------------------------------------------
edit('include/chem/channel/channel.h',
     '#include <cstdint>\n',
     '#include <cstdint>\n'
     '// >>> CHEM-EXT-HOOK (a2g_doppler)\n'
     '#include <functional>\n'
     '#include <memory>\n'
     '#define CHEM_EXT_TV_CIR_HOOK 1\n'
     '// <<< CHEM-EXT-HOOK\n',
     'CHEM-EXT-HOOK')
s = open(os.path.join(root, 'include/chem/channel/channel.h'), encoding='utf-8').read()
if 'setTimeVaryingCIR' not in s:
    s = s.replace(
        '    void setExtensionPathLossDb(float pl_db);\n    float getExtensionPathLossDb() const;\n',
        '    void setExtensionPathLossDb(float pl_db);\n    float getExtensionPathLossDb() const;\n\n'
        '    // >>> CHEM-EXT-HOOK (a2g_doppler): optional time-varying CIR stage.\n'
        '    // When set, it replaces the static CIR convolution of this link.\n'
        '    // Arguments: buffer (in/out), samples, sample rate (Hz), dest port,\n'
        '    // src port. The function must keep its own state between calls.\n'
        '    using TimeVaryingCIRFn =\n'
        '        std::function<void(signal_v&, size_t, double, size_t, size_t)>;\n'
        '    void setTimeVaryingCIR(TimeVaryingCIRFn fn);\n'
        '    void clearTimeVaryingCIR();\n'
        '    bool hasTimeVaryingCIR() const;\n'
        '    // <<< CHEM-EXT-HOOK\n', 1)
    s = s.replace(
        '    float _extensionPathLossDb{0.0f};\n};',
        '    float _extensionPathLossDb{0.0f};\n'
        '    // >>> CHEM-EXT-HOOK (a2g_doppler)\n'
        '    std::shared_ptr<TimeVaryingCIRFn> _tvCir;\n'
        '    // <<< CHEM-EXT-HOOK\n};', 1)
    open(os.path.join(root, 'include/chem/channel/channel.h'), 'w', encoding='utf-8', newline='\n').write(s)
    print('edited: include/chem/channel/channel.h (members)')

edit('src/chem/channel/channel.cpp',
     '    if (cir_active) {\n        if (_cirApplyMethod == dsp::channel::CIRApplyMethod::DIRECT_CONV) {',
     '    // >>> CHEM-EXT-HOOK (a2g_doppler): time-varying CIR replaces the static\n'
     '    // CIR stage when an extension installed one for this link.\n'
     '    if (auto tv = std::atomic_load(&_tvCir); tv && *tv) {\n'
     '        (*tv)(dest, params.s_perBuff, params.sample_rate_hz,\n'
     '              params.dest_channel_index, params.src_channel_index);\n'
     '    } else\n'
     '    // <<< CHEM-EXT-HOOK\n'
     '    if (cir_active) {\n        if (_cirApplyMethod == dsp::channel::CIRApplyMethod::DIRECT_CONV) {',
     'tv && *tv')
edit('src/chem/channel/channel.cpp',
     'float Channel::getExtensionPathLossDb() const {\n    return _extensionPathLossDb;\n}\n',
     'float Channel::getExtensionPathLossDb() const {\n    return _extensionPathLossDb;\n}\n\n'
     '// >>> CHEM-EXT-HOOK (a2g_doppler)\n'
     'void Channel::setTimeVaryingCIR(TimeVaryingCIRFn fn) {\n'
     '    std::atomic_store(&_tvCir, std::make_shared<TimeVaryingCIRFn>(std::move(fn)));\n'
     '}\n\n'
     'void Channel::clearTimeVaryingCIR() {\n'
     '    std::atomic_store(&_tvCir, std::shared_ptr<TimeVaryingCIRFn>());\n'
     '}\n\n'
     'bool Channel::hasTimeVaryingCIR() const {\n'
     '    auto tv = std::atomic_load(&_tvCir);\n'
     '    return tv && *tv;\n'
     '}\n'
     '// <<< CHEM-EXT-HOOK\n',
     'Channel::setTimeVaryingCIR')
print('done')
