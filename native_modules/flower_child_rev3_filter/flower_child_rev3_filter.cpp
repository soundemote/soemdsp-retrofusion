// soemdsp-native-module: flower_child_rev3_filter
// soemdsp-native-label: Flower Child Rev3 Filter
// soemdsp-native-target: flowerChildRev3Filter
// soemdsp-native-kind: filter
//
// Split out of flower_child_filter's mode 2 (Rev3): a different
// oscillator scheme from the filter's other modes -- an "ellipsoid"
// waveshaper on a bipolar phasor value, driven by 5 independent 4-node
// resonance-shaping curves (phase-mod amount, sine amplitude, sine-to-
// square morph, clip level, noise reduction), through two cascaded
// impulse-invariant-transform one-pole lowpass stages.
//
// One curve (clipLevelGraph) had a genuinely ambiguous node list in the
// original 2021 source (nodes added twice, with a shape/skew set-call
// that only lines up with the intended values under the reading that
// the later three-node redefinition (7 -> 7 -> 2) supersedes an
// earlier, apparently-superseded four-value block) -- reproduced here
// using that later, clearly-intentional definition.

namespace {

static const int kMaxInstances = 32;
static const double kPi     = 3.141592653589793238;
static const double kTwoPi  = 6.283185307179586476;
static const double kHalfPi = 1.5707963267948966192;

union DoubleBits {
  double d;
  unsigned long long u;
};

static double poly_sin_0_halfpi(double x) {
  const double x2 = x * x;
  return x * (1.0 + x2 * (-1.6666666666666667e-1 + x2 * (8.3333333333333329e-3 + x2 * (-1.9841269841269841e-4 + x2 * (2.7557319223985888e-6 + x2 * (-2.5052108385441720e-8 + x2 * 1.6059043836821614e-10))))));
}

// x must be in [0, pi]
static double dsp_sin_0_pi(double x) {
  if (x > kHalfPi) x = kPi - x;
  return poly_sin_0_halfpi(x);
}

static inline double dsp_floor(double x) {
  double xi = (double)(long long)x;
  return (x < xi) ? xi - 1.0 : xi;
}

// Full-range sin/cos via quadrant folding onto the [0, pi/2] polynomial.
static double dsp_sin(double x) {
  double wrapped = x - kTwoPi * dsp_floor(x / kTwoPi);
  double sign = 1.0;
  if (wrapped >= kPi) {
    wrapped -= kPi;
    sign = -1.0;
  }
  return sign * dsp_sin_0_pi(wrapped);
}

static double dsp_cos(double x) {
  return dsp_sin(x + kHalfPi);
}

// 2^f for f in [0,1), truncated Taylor series of e^(f*ln2) -- accurate to
// better than 1e-5 relative error, which is far more precision than a
// musical pitch-to-frequency conversion needs.
static double pow2_frac(double f) {
  const double c1 = 0.6931471805599453, c2 = 0.2402265069591007,
               c3 = 0.05550410866482158, c4 = 0.009618129107628477,
               c5 = 0.001333355814670365, c6 = 0.0001540353039338161;
  return 1.0 + f * (c1 + f * (c2 + f * (c3 + f * (c4 + f * (c5 + f * c6)))));
}

static double dsp_exp2(double x) {
  double xi = dsp_floor(x);
  double f = x - xi;
  double p = pow2_frac(f);
  long long n = (long long)xi;
  DoubleBits bits;
  bits.d = p;
  long long expBits = (long long)((bits.u >> 52) & 0x7FF);
  expBits += n;
  if (expBits < 1) expBits = 1;
  if (expBits > 2046) expBits = 2046;
  bits.u = (bits.u & ~(0x7FFULL << 52)) | ((unsigned long long)expBits << 52);
  return bits.d;
}

static inline double dsp_sqrt(double x) {
  if (x <= 0.0) return 0.0;
  double guess = x;
  for (int i = 0; i < 24; i++) guess = 0.5 * (guess + x / guess);
  return guess;
}

static inline double clampd(double v, double lo, double hi) {
  return v < lo ? lo : (v > hi ? hi : v);
}

static inline double jmap01(double v, double outMin, double outMax) {
  return outMin + (outMax - outMin) * v;
}

static inline double jmapGeneral(double v, double srcMin, double srcMax, double dstMin, double dstMax) {
  return dstMin + (dstMax - dstMin) * (v - srcMin) / (srcMax - srcMin);
}

static inline double dsp_exp(double x) {
  double clamped = x < -40.0 ? -40.0 : (x > 40.0 ? 40.0 : x);
  return dsp_exp2(clamped * 1.4426950408889634);
}

static double dsp_ln(double x) {
  if (x <= 0.0) return -700.0;
  DoubleBits bits;
  bits.d = x;
  long long expBits = (long long)((bits.u >> 52) & 0x7FF);
  int e = (int)(expBits - 1023);
  bits.u = (bits.u & ~(0x7FFULL << 52)) | (1023ULL << 52);
  double m = bits.d;
  double t = (m - 1.0) / (m + 1.0);
  double t2 = t * t;
  double series = t * (1.0 + t2 * (1.0 / 3.0 + t2 * (1.0 / 5.0 + t2 * (1.0 / 7.0 + t2 * (1.0 / 9.0)))));
  return (double)e * 0.6931471805599453 + 2.0 * series;
}

// rsOnePoleFilter LOWPASS_IIT (matched-Z-transform one-pole): a1=exp(-w),
// b0=1-a1, y[n]=b0*x[n]+a1*y[n-1]. Exact, per soemdsp/filter/OnePoleFilter.hpp.
static inline double onePoleIitCoefficient(double cutoffHz, double sampleRate) {
  double w = clampd(kTwoPi * cutoffHz / sampleRate, 1e-9, kPi * 0.98);
  return dsp_exp(-w);
}

static inline double onePoleIitStep(double* y1, double input, double a1) {
  double b0 = 1.0 - a1;
  *y1 = b0 * input + a1 * (*y1);
  return *y1;
}

static inline double rationalCurve(double p, double skew) {
  return ((1.0 + skew) * p) / (1.0 - skew + 2.0 * skew * p);
}

// Generic N-node soemdsp::utility::Graph evaluator (shape 1=RATIONAL,
// 2=EXPONENTIAL, else linear).
struct GraphNode {
  double x, y, skew;
  int shape;
};

static double evalGraph(const GraphNode* nodes, int count, double x) {
  if (count <= 0) return 0.0;
  if (x < nodes[0].x) return nodes[0].y;
  int i = -1;
  for (int k = 0; k < count; k++) {
    if (nodes[k].x > x) { i = k; break; }
  }
  if (i < 0) return nodes[count - 1].y;
  if (i == 0) return nodes[0].y;
  const GraphNode& n1 = nodes[i - 1];
  const GraphNode& n2 = nodes[i];
  if (n2.x - n1.x < 1e-9) return 0.5 * (n1.y + n2.y);
  double p = (x - n1.x) / (n2.x - n1.x);
  if (n2.shape == 1) return n1.y + (n2.y - n1.y) * rationalCurve(p, n2.skew);
  if (n2.shape == 2) {
    double c = 0.5 * (n2.skew + 1.0);
    double a = 2.0 * dsp_ln((1.0 - c) / c);
    return n1.y + (n2.y - n1.y) * (1.0 - dsp_exp(p * a)) / (1.0 - dsp_exp(a));
  }
  return n1.y + (n2.y - n1.y) * p;
}

static inline double pitchToFreq(double pitch) {
  return 440.0 * dsp_exp2((pitch - 69.0) / 12.0);
}

// waveshape::ellipse -- phase is unipolar [0,1). A/B_sin/B_cos/C per the
// original signature; Flower Child always calls it with A=0, B_sin=0,
// B_cos=1, C=ellipseC.
static inline double waveEllipse(double phase, double ellipseC) {
  double sinX = dsp_sin(phase * kTwoPi);
  double cosX = dsp_cos(phase * kTwoPi);
  double sqrtVal = dsp_sqrt(cosX * cosX + (ellipseC * sinX) * (ellipseC * sinX));
  if (sqrtVal < 1e-12) sqrtVal = 1e-12;
  return cosX / sqrtVal;
}

// Simple xorshift32 PRNG for the (nearly-unfiltered in the original --
// its bandpass corners sit at ~0.125Hz and ~20.6kHz) noise contribution.
static inline double nextNoiseBipolar(unsigned int* state) {
  unsigned int x = *state;
  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 5;
  *state = x;
  return ((double)x / 4294967295.0) * 2.0 - 1.0;
}

struct FlowerChildRev3State {
  bool active;
  double phase;
  unsigned int rngState;
  double feedback;
  double lpf1Y1, lpf2Y1;  // one-pole IIT stages
};

static FlowerChildRev3State gPool[kMaxInstances];

}  // namespace

extern "C" int soemdsp_flower_child_rev3_filter_create() {
  for (int i = 0; i < kMaxInstances; i++) {
    if (!gPool[i].active) {
      FlowerChildRev3State& s = gPool[i];
      s.phase = 0.0;
      s.rngState = 0x9E3779B9u + (unsigned int)(i + 1) * 2654435761u;
      s.feedback = 0.0;
      s.lpf1Y1 = 0.0;
      s.lpf2Y1 = 0.0;
      s.active = true;
      return i + 1;
    }
  }
  return 0;
}

extern "C" void soemdsp_flower_child_rev3_filter_destroy(int handle) {
  if (handle < 1 || handle > kMaxInstances) return;
  gPool[handle - 1].active = false;
}

extern "C" double soemdsp_flower_child_rev3_filter_sample(
  int handle,
  double input,
  double frequency,   // 0..1 normalized slider, matches original
  double resonance,   // 0..1
  double chaosAmount,  // 0..1
  double sampleRate
) {
  if (handle < 1 || handle > kMaxInstances) return 0.0;
  FlowerChildRev3State& s = gPool[handle - 1];

  const double safeRate = sampleRate < 1.0 ? 44100.0 : sampleRate;
  const double freqNorm = clampd(frequency, 0.0, 1.0);
  const double reso = clampd(resonance, 0.0, 1.0);
  const double chaos = clampd(chaosAmount, 0.0, 1.0);

  const double masterPitch = jmap01(freqNorm, -120.0, 105.0);
  const double masterFrequency = pitchToFreq(masterPitch);
  const double fmAmount = pitchToFreq(-48.377);
  const double lpf1Cutoff = pitchToFreq(jmapGeneral(masterPitch, -120.0, 120.0, 90.0, 180.0));
  const double lpf2Cutoff = pitchToFreq(jmapGeneral(masterPitch, -120.0, 120.0, 80.0, 130.0));
  const double lpf1A = onePoleIitCoefficient(lpf1Cutoff, safeRate);
  const double lpf2A = onePoleIitCoefficient(lpf2Cutoff, safeRate);

  const GraphNode phaseModGraph[4] = {
    {0, 0.0, 0, 0}, {0.5, -0.017446, 0.9, 1}, {0.6, -0.017575, 0.0, 1}, {1.0, -0.0147, 0.6, 1},
  };
  const GraphNode sineAmpGraph[4] = {
    {0, 4.44777, 0, 0}, {0.5, 8.6687, 0.9, 1}, {0.6, 8.6687, 0.0, 1}, {1.0, 2.0, 0.6, 1},
  };
  const GraphNode sineToSquareGraph[4] = {
    {0, 0.6792, 0, 0}, {0.5, 0.9552, 0.9, 1}, {0.6, 0.9552, 0.0, 1}, {1.0, 0.001, 0.6, 1},
  };
  const GraphNode clipLevelGraph[3] = { {0.0, 7.0, 0, 0}, {0.7, 7.0, 0.0, 1}, {1.0, 2.0, 0.6, 1} };
  const GraphNode noiseGraph[3] = { {0.0, 0.0, 0, 0}, {0.8, 0.1, 0, 0}, {1.0, 1.0, 0.0, 1} };

  const double pmAmount = evalGraph(phaseModGraph, 4, reso);
  const double sineAmp = evalGraph(sineAmpGraph, 4, reso);
  const double sineToSquare = evalGraph(sineToSquareGraph, 4, reso);
  const double clipLevelRaw = evalGraph(clipLevelGraph, 3, reso);
  const double clipLevel = sineAmp < clipLevelRaw ? sineAmp : clipLevelRaw;
  const double noiseReduction = evalGraph(noiseGraph, 3, reso);
  const double chaosAmount4x = chaos * 4.0;

  double in = s.feedback + clampd(-1.0 * input, -clipLevel, clipLevel);
  const double f = masterFrequency * in * fmAmount;
  const double noiseTerm = masterFrequency * nextNoiseBipolar(&s.rngState) * chaosAmount4x * noiseReduction;

  const double incAmt = (f + noiseTerm) / safeRate;
  s.phase = s.phase + incAmt;
  s.phase = s.phase - dsp_floor(s.phase);
  const double bipolarPhasor = 2.0 * s.phase - 1.0;
  const double phasorOut = bipolarPhasor + pmAmount * s.feedback;

  const double ellipseOut = sineAmp * waveEllipse(phasorOut, sineToSquare);

  double feedback = onePoleIitStep(&s.lpf1Y1, ellipseOut, lpf1A);
  feedback = onePoleIitStep(&s.lpf2Y1, feedback, lpf2A);
  s.feedback = feedback;

  return feedback * 0.15;
}

extern "C" int soemdsp_flower_child_rev3_filter_version() {
  return 1;
}
