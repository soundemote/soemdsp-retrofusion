# soemdsp-retrofusion

Six feedback/chaotic resonator filters, native WebAssembly, original
designs (no external reference implementation).

## Build target

`--target=wasm32 -O3 -nostdlib -fno-exceptions -fno-rtti`. Compiled
with clang++.

## Export convention

All 6 modules share the same shape:

```
int    {prefix}_create()
void   {prefix}_destroy(int handle)
double {prefix}_sample(int handle, double input, ...params, double sampleRate)
int    {prefix}_version()
```

`frequency`, `resonance`, `chaosAmount` are 0-1 normalized on every
module.

## Modules

| Module | Export prefix | `_sample` parameters after `input` | Modes | Description |
|---|---|---|---|---|
| Yellowjacket Filter | `soemdsp_yellowjacket_filter` | `frequency, resonance, chaosAmount, sampleRate` | — | Feedback-modulated ellipse-oscillator filter through a one-pole stage; `chaosAmount` drives filter cutoff scaling. |
| SuperLove Filter | `soemdsp_superlove_filter` | `frequency, resonance, chaosAmount, mode, sampleRate` | 0=LP18, 1=LP24, 2=HP6, 3=BP6 | Trisaw-oscillator feedback resonator through a multi-pole ladder tap. |
| Chaotic Phase Locking Filter | `soemdsp_chaotic_phase_locking_filter` | `frequency, resonance, chaosAmount, sampleRate` | — | Feedback ellipse-waveshaper resonator (no oscillator phasor) through a 12 dB lowpass and DC-blocking highpass. |
| Resonator Filter | `soemdsp_resonator_filter` | `frequency, resonance, chaosAmount, mode, sampleRate` | 0=Sinusoid, 1=Triangle, 2=Sawtooth | Dual-phasor FM feedback resonator through a one-pole lowpass and DC-blocking highpass. |
| Human Filter | `soemdsp_human_filter` | `frequency, resonance, chaosAmount, mode, sampleRate` | 0=BP6, 1=LP6, 2=LP12 | Dual-phasor feedback network shaped by a bell/peak filter in the feedback path, DC-blocking highpass on output. |
| Flower Child Rev3 Filter | `soemdsp_flower_child_rev3_filter` | `frequency, resonance, chaosAmount, sampleRate` (no `mode` param) | — | Ellipsoid waveshaper on a bipolar phasor, driven by 5 resonance-shaping curves, through two cascaded impulse-invariant-transform one-pole lowpass stages, feedback into its own phase modulation input. Extracted from `flower_child_filter`'s former mode 2 in `soemdsp-synthwave`; numerically verified byte-identical to the original mode 2 across multiple parameter sets and sample rates. |

Instance pool: `kMaxInstances` handles per module (see each `.cpp`).

Source: `native_modules/{yellowjacket_filter,superlove_filter,chaotic_phase_locking_filter,resonator_filter,human_filter,flower_child_rev3_filter}/`.
