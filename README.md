# soemdsp-retrofusion

Five filters, none of them well-behaved. Feedback resonators and
chaotic waveshapers that would rather self-oscillate than sit quietly
and roll off frequencies -- the analog-filter grit that made old drum
machines and synths sound like they were fighting back.

## What's here

- **Yellowjacket Filter** — a feedback-modulated ellipse-oscillator
  filter through a one-pole stage. Grindy, easily produces
  square-wave-like output.
- **SuperLove Filter** — a trisaw-oscillator feedback resonator through
  a multi-pole ladder tap. Warm, bass-heavy, stably self-oscillating.
- **Chaotic Phase Locking Filter** — a feedback ellipse-waveshaper
  resonator producing phase-locked chaotic textures.
- **Resonator Filter** — a dual-phasor FM feedback resonator, 3 chaotic
  waveform-variation modes.
- **Human Filter** — a dual-phasor feedback network shaped by a
  bell/peak filter in the feedback path.

Each one compiled to a dependency-free WebAssembly module under
`native_modules/`.

All five are original designs, not based on any published third-party
filter (the RS-MET-lineage filters that used to be listed here --
RSMET, Ladder, and TB-303 -- now live in `soemdsp-rsmet`).
