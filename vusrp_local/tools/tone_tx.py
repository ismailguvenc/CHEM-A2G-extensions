#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Tone transmitter on a local V-USRP. usage: tone_tx.py <device name> <secs> <f_port0> [f_port1]
import os, sys, time, numpy as np, uhd
name, dur = sys.argv[1], float(sys.argv[2]); tones = [float(f) for f in sys.argv[3:]]
FS = 1e6; FC = 3.4e9; ch = list(range(len(tones)))
u = uhd.usrp.MultiUSRP(f"name={name},num_channels={len(ch)}")
for c in ch:
    u.set_tx_rate(FS, c); u.set_tx_freq(uhd.types.TuneRequest(FC), c); u.set_tx_gain(60, c)
sa = uhd.usrp.StreamArgs("fc32", "sc16"); sa.channels = ch; st = u.get_tx_stream(sa)
n = int(os.environ.get("TX_BUF", 10000)); t = np.arange(n) / FS
buf = np.vstack([0.3 * np.exp(2j * np.pi * f * t) for f in tones]).astype(np.complex64)
md = uhd.types.TXMetadata(); t0 = time.time(); k = 0
while time.time() - t0 < dur:
    st.send(buf, md); k += 1
md.end_of_burst = True; st.send(np.zeros((len(ch), 10), np.complex64), md); print("tx buffers", k, flush=True)
