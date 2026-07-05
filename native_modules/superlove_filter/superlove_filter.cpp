// soemdsp-native-module: superlove_filter
// soemdsp-native-label: SuperLove Filter
// soemdsp-native-target: superloveFilter
// soemdsp-native-kind: filter
//
// A trisaw-oscillator feedback resonator through a multi-pole ladder
// filter tap. Four modes ported from the original SuperLove_LP18,
// SuperLove_LP24, SuperLove_HP6, and SuperLove_BP6 -- a warm, bass-heavy,
// stably self-oscillating design (Japan's answer to the Flower Child
// filter wars, per the original help text).

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

static double dsp_sin_0_pi(double x) {
  if (x > kHalfPi) x = kPi - x;
  return poly_sin_0_halfpi(x);
}

static double dsp_cos_0_pi(double x) {
  double y = kHalfPi - x;
  if (y < 0.0) return -poly_sin_0_halfpi(-y);
  return poly_sin_0_halfpi(y);
}

static double dsp_tan_neg_halfquarter(double x) {
  const double ax = -x;
  const double s = poly_sin_0_halfpi(ax);
  const double c = poly_sin_0_halfpi(kHalfPi - ax);
  return (c == 0.0) ? -1e15 : -(s / c);
}

static inline double dsp_floor(double x) {
  double xi = (double)(long long)x;
  return (x < xi) ? xi - 1.0 : xi;
}

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

static inline double clampd(double v, double lo, double hi) {
  return v < lo ? lo : (v > hi ? hi : v);
}

static inline double jmap01(double v, double outMin, double outMax) {
  return outMin + (outMax - outMin) * v;
}

static inline double pitchToFreq(double pitch) {
  return 440.0 * dsp_exp2((pitch - 69.0) / 12.0);
}

// soemdsp's standalone curve() tension function (distinct from
// curve::Rational). 0->0, 1->1, tension 0 is linear.
static inline double curveShape(double v, double tension) {
  double denom = 2.0 * tension * v - tension - 1.0;
  if (denom == 0.0) return v;
  return (tension * v - v) / denom;
}

static inline double rationalCurve(double p, double skew) {
  return ((1.0 + skew) * p) / (1.0 - skew + 2.0 * skew * p);
}

// Generic N-node soemdsp::utility::Graph evaluator (shape 1=RATIONAL,
// 2=EXPONENTIAL, else linear), matching Graph::getValue exactly: clamps
// below the first node and at/beyond the last node.
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

// waveshape::trisaw(phaseCycles, morph) -- exact, see waveshapes.cpp.
static double waveTrisaw(double phaseCycles, double morph) {
  double phaseRad = phaseCycles * kTwoPi;
  phaseRad = phaseRad - kTwoPi * dsp_floor(phaseRad / kTwoPi);
  double morphRad = morph * kTwoPi;
  double sourceMin, sourceMax, targetMin, targetRange;
  if (phaseRad > morphRad) {
    sourceMin = morphRad; sourceMax = kTwoPi; targetMin = 1.0; targetRange = -1.0;
  } else {
    sourceMin = 0.0; sourceMax = morphRad; targetMin = 0.0; targetRange = 1.0;
  }
  double sourceRange = sourceMax - sourceMin;
  double uni;
  if (sourceMin == sourceMax) uni = sourceMin;
  else uni = targetMin + (targetRange * (phaseRad - sourceMin)) / sourceRange;
  return 2.0 * uni - 1.0;
}

// General N-pole ladder cascade tap (resonance fixed at 0, matching the
// original SuperLove filters which never call setResonance on their
// ladder instances).
static double ladderTapStep(double y[5], double input, double a, int mode, int stages) {
  double c[5] = {0, 0, 0, 0, 0};
  if (mode == 1) {
    c[stages] = 1.0;
  } else if (mode == 2) {
    static const double hp[4][5] = {
      {1.0, -1.0, 0.0, 0.0, 0.0},
      {1.0, -2.0, 1.0, 0.0, 0.0},
      {1.0, -3.0, 3.0, -1.0, 0.0},
      {1.0, -4.0, 6.0, -4.0, 1.0},
    };
    for (int i = 0; i <= stages; i++) c[i] = hp[stages - 1][i];
  } else if (mode == 3) {
    static const double bp[4][5] = {
      {0.0, 2.0, -2.0, 0.0, 0.0},
      {0.0, 2.0, -2.0, 0.0, 0.0},
      {0.0, 0.0, 3.0, -3.0, 0.0},
      {0.0, 0.0, 4.0, -8.0, 4.0},
    };
    for (int i = 0; i < 5; i++) c[i] = bp[stages - 1][i];
  }
  double y0 = input;
  y0 = y0 / (1.0 + y0 * y0);
  y[1] = y0 + a * (y0 - y[1]);
  y[2] = y[1] + a * (y[1] - y[2]);
  y[3] = y[2] + a * (y[2] - y[3]);
  y[4] = y[3] + a * (y[3] - y[4]);
  y[0] = y0;
  return c[0] * y[0] + c[1] * y[1] + c[2] * y[2] + c[3] * y[3] + c[4] * y[4];
}

static inline double ladderCoefficient(double cutoffHz, double sampleRate) {
  double rawWc = kTwoPi * cutoffHz / sampleRate;
  double wc = clampd(rawWc, 1e-9, kPi * 0.98);
  double s = dsp_sin_0_pi(wc);
  double c = dsp_cos_0_pi(wc);
  double t = dsp_tan_neg_halfquarter(0.25 * (wc - kPi));
  double denom = s - c * t;
  if (denom > -1e-12 && denom < 1e-12) denom = (denom >= 0.0) ? 1e-12 : -1e-12;
  return t / denom;
}

static inline double nextNoiseBipolar(unsigned int* state) {
  unsigned int x = *state;
  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 5;
  *state = x;
  return ((double)x / 4294967295.0) * 2.0 - 1.0;
}

struct SuperLoveState {
  bool active;
  double feedbackSignal;
  double filterY[5];
  double dcY[5];
  unsigned int rngState;
};

static SuperLoveState gPool[kMaxInstances];

}  // namespace

extern "C" int soemdsp_superlove_filter_create() {
  for (int i = 0; i < kMaxInstances; i++) {
    if (!gPool[i].active) {
      SuperLoveState& s = gPool[i];
      s.feedbackSignal = 0.0;
      for (int j = 0; j < 5; j++) { s.filterY[j] = 0.0; s.dcY[j] = 0.0; }
      s.rngState = 0x85EBCA6Bu + (unsigned int)(i + 1) * 0xC2B2AE35u;
      s.active = true;
      return i + 1;
    }
  }
  return 0;
}

extern "C" void soemdsp_superlove_filter_destroy(int handle) {
  if (handle < 1 || handle > kMaxInstances) return;
  gPool[handle - 1].active = false;
}

extern "C" double soemdsp_superlove_filter_sample(
  int handle,
  double input,
  double frequency,  // 0..1 normalized
  double resonance,  // 0..1 (decay/density depending on mode)
  double chaosAmount,  // 0..1 (noise/shape depending on mode)
  int mode,          // 0=LP18, 1=LP24, 2=HP6, 3=BP6
  double sampleRate
) {
  if (handle < 1 || handle > kMaxInstances) return 0.0;
  SuperLoveState& s = gPool[handle - 1];

  const double safeRate = sampleRate < 1.0 ? 44100.0 : sampleRate;
  const double freqNorm = clampd(frequency, 0.0, 1.0);
  const double reso = clampd(resonance, 0.0, 1.0);
  const double chaos = clampd(chaosAmount, 0.0, 1.0);
  const int safeMode = mode < 0 ? 0 : (mode > 3 ? 3 : mode);

  if (safeMode <= 1) {
    // LP18 / LP24
    const GraphNode resonanceGraph[2] = { {0, 0, 0, 0}, {1, -2.7175, -0.85, 2} };
    const GraphNode noiseGraph[3] = { {0, 0.00, 0, 0}, {0.75, 0.05, -0.7, 2}, {1, 0.10, 0.6, 2} };

    const double cutoffHz = clampd(pitchToFreq(jmap01(freqNorm, -12.0, 135.0)), 0.0, 0.5 * safeRate);
    const double mod = evalGraph(resonanceGraph, 2, reso);
    const double noiseAmp = evalGraph(noiseGraph, 3, chaos);
    const double shape = chaos;

    s.feedbackSignal = mod * s.feedbackSignal + input;
    // noiseFilter is a near-unfiltered bandpass in the original (corners at
    // ~0.125Hz and ~20.6kHz), so raw white noise stands in for it here --
    // same documented approximation used in flower_child_filter.
    const double pm = nextNoiseBipolar(&s.rngState) * noiseAmp;
    double oscValue = -waveTrisaw(s.feedbackSignal + 0.25725 + pm, shape);

    const double a = ladderCoefficient(cutoffHz, safeRate);
    const int ladderMode = 1;
    const int stages = safeMode == 0 ? 3 : 4;
    s.feedbackSignal = ladderTapStep(s.filterY, oscValue, a, ladderMode, stages);

    const double dcCutoff = safeMode == 0 ? 10.0 : 5.0;
    const int dcStages = safeMode == 0 ? 3 : 1;
    const double dcA = ladderCoefficient(dcCutoff, safeRate);
    const double dcOut = ladderTapStep(s.dcY, s.feedbackSignal, dcA, 2, dcStages);

    return dcOut * 1.02;
  } else if (safeMode == 2) {
    // HP6
    const GraphNode resonanceGraph[2] = { {0, -0.2, 0, 0}, {1, 1.3, -0.85, 2} };
    const double mod = evalGraph(resonanceGraph, 2, reso);
    const double shape = 1.0 - chaos;

    s.feedbackSignal = mod * s.feedbackSignal + input;
    double oscValue = -waveTrisaw(s.feedbackSignal + 0.75, shape);

    const double lpCutoff = safeRate * 0.5;
    const double lpA = ladderCoefficient(lpCutoff, safeRate);
    double fb = ladderTapStep(s.filterY, oscValue * 0.1, lpA, 1, 1);

    const double cutoffHz = clampd(pitchToFreq(jmap01(freqNorm, -12.0, 135.0)), 0.0, 0.5 * safeRate);
    const double hpA = ladderCoefficient(cutoffHz, safeRate);
    fb = ladderTapStep(s.dcY, fb, hpA, 2, 1);
    fb *= 10.0;

    s.feedbackSignal = fb;
    return -fb * 0.31;
  } else {
    // BP6
    const GraphNode resonanceGraph[2] = { {0, -0.2, 0, 0}, {1, 1.3, -0.85, 2} };
    const double mod = evalGraph(resonanceGraph, 2, reso);
    const double shape = 1.0 - chaos;

    s.feedbackSignal = mod * s.feedbackSignal + input;
    double oscValue = -waveTrisaw(s.feedbackSignal + 0.75, shape);

    const double cutoffHz = clampd(pitchToFreq(jmap01(freqNorm, -12.0, 135.0)), 0.0, 0.5 * safeRate);
    const double a = ladderCoefficient(cutoffHz, safeRate);
    double fb = ladderTapStep(s.filterY, oscValue * 0.1, a, 3, 1);
    fb *= 10.0;

    s.feedbackSignal = fb;
    return fb;
  }
}

extern "C" int soemdsp_superlove_filter_version() {
  return 1;
}
