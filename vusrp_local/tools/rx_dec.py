#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Capture on a local V-USRP at 1 MS/s; per port, mix by -f_k for each tone k and decimate by 100 -> 10 kHz.
# usage: rx_dec.py <device name> <secs> <out.npz> <nch> <f_tone0> [f_tone1]
import sys, time, numpy as np, uhd
name, dur, out, nch = sys.argv[1], float(sys.argv[2]), sys.argv[3], int(sys.argv[4]); tones = [float(f) for f in sys.argv[5:]]
FS = 1e6; D = 100; FC = 3.4e9; ch = list(range(nch))
u = uhd.usrp.MultiUSRP(f"name={name},num_channels={len(ch)}")
for c in ch:
    u.set_rx_rate(FS, c); u.set_rx_freq(uhd.types.TuneRequest(FC), c); u.set_rx_gain(40, c)
sa = uhd.usrp.StreamArgs("fc32", "sc16"); sa.channels = ch; st = u.get_rx_stream(sa)
N = int(FS * dur) // D * D; buf = np.zeros((nch, 50000), np.complex64); md = uhd.types.RXMetadata()
y = np.zeros((nch, len(tones), N // D), np.complex64); acc = np.zeros((nch, 0), np.complex64); got = 0; w = 0; errs = 0; pw = []
cmd = uhd.types.StreamCMD(uhd.types.StreamMode.start_cont); cmd.stream_now = True; st.issue_stream_cmd(cmd); wall0 = time.time()
while got < N:
    r = st.recv(buf, md, 2.0)
    if md.error_code != uhd.types.RXMetadataErrorCode.none: errs += 1
    if r == 0: continue
    r = min(r, N - got); acc = np.concatenate([acc, buf[:, :r]], axis=1); m = acc.shape[1] // D * D
    if m:
        n = np.arange(got - (acc.shape[1] - r), got - (acc.shape[1] - r) + m)
        for k, f in enumerate(tones):
            y[:, k, w:w + m // D] = (acc[:, :m] * np.exp(-2j * np.pi * f * n / FS)[None, :]).reshape(nch, -1, D).mean(axis=2)
        pw.append(float(np.mean(np.abs(acc[:, :m]) ** 2))); w += m // D; acc = acc[:, m:]
    got += r
st.issue_stream_cmd(uhd.types.StreamCMD(uhd.types.StreamMode.stop_cont))
np.savez(out, y=y[:, :, :w], fs=FS / D, wall0=wall0, errs=errs, pw=np.array(pw))
print(f"rx {w} samples at 10 kHz per port and tone, errors {errs}", flush=True)
