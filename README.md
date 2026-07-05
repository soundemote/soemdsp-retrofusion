# soemdsp-retrofusion

Eight filters, none of them well-behaved. Feedback resonators, chaotic
waveshapers, and ladder circuits that would rather self-oscillate than
sit quietly and roll off frequencies -- the analog-filter grit that made
old drum machines and synths sound like they were fighting back.

## What's here

- **RSMET Filter** — a ladder filter preceded by a tanh soft clipper and
  noise injection stage, exponential frequency/resonance curves, 10 modes.
- **Ladder Filter** — the plain RS-MET-lineage ladder design underneath
  RSMET Filter, gain-compensated, resonant stages.
- **TB-303 Filter** — feedback highpass, resonance skewing, 15 output
  modes (LP/HP/BP at 6/12/18/24 dB/octave).
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

## Attribution

RSMET Filter, Ladder Filter, and TB-303 Filter are built in the lineage
of Robin Schmidt's (RS-MET) published filter designs, including
TeeBeeFilter / Open303. Credited here explicitly -- these three are not
original designs.
