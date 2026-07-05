// soemdsp-native-module: yellowjacket_filter
// soemdsp-native-label: Yellowjacket Filter
// soemdsp-native-target: yellowjacketFilter
// soemdsp-native-kind: filter
//
// A feedback-modulated ellipse-oscillator filter, one-pole (6dB/octave)
// output stage, with a resonance-vs-frequency shaping curve controlling
// both the oscillator's waveshape and its feedback gain. Ported from the
// original Yellowjacket_BP -- grindy, heavily overdriven, easily produces
// square-wave-like output at most resonance settings.

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

static inline double pitchToFreq(double pitch) {
  return 440.0 * dsp_exp2((pitch - 69.0) / 12.0);
}

// Exact soemdsp::utility::Graph::getValue, 2-node EXPONENTIAL curve.
static double evalExponentialGraph2(double x, double y0, double y1, double skew) {
  if (x <= 0.0) return y0;
  if (x >= 1.0) return y1;
  const double c = 0.5 * (skew + 1.0);
  const double a = 2.0 * dsp_ln((1.0 - c) / c);
  const double expA = dsp_exp(a);
  const double expIA = dsp_exp(x * a);
  return y0 + (y1 - y0) * (1.0 - expIA) / (1.0 - expA);
}

// Exact soemdsp::curve::Rational::get(p).
static inline double rationalCurve(double p, double skew) {
  return ((1.0 + skew) * p) / (1.0 - skew + 2.0 * skew * p);
}

// 3-node Graph, RATIONAL shape at node2, matching Yellowjacket's
// resVfreqGraph: node0/node1 track resonance (flat segment between them),
// node2 is a FIXED constant (never moved after construction).
static double evalResVFreqGraph(double x, double n0y, double breakpoint, double n2yFixed, double skew) {
  if (x < 0.0) return n0y;
  if (x >= 1.0) return n2yFixed;
  if (x < breakpoint) return n0y;
  double p = (x - breakpoint) / (1.0 - breakpoint);
  return n0y + (n2yFixed - n0y) * rationalCurve(p, skew);
}

// waveshape::ellipse, full signature (A, B_sin, B_cos, C).
static double waveEllipseFull(double phase, double A, double bSin, double bCos, double C) {
  double sinX = dsp_sin(phase * kTwoPi);
  double cosX = dsp_cos(phase * kTwoPi);
  double apc = A + cosX;
  double sqrtVal = dsp_sqrt(apc * apc + (C * sinX) * (C * sinX));
  if (sqrtVal < 1e-12) sqrtVal = 1e-12;
  return (apc * bCos + (C * sinX) * bSin) / sqrtVal;
}

static inline double onePoleCoefficient(double cutoffHz, double sampleRate) {
  double rawWc = kTwoPi * cutoffHz / sampleRate;
  double wc = clampd(rawWc, 1e-9, kPi * 0.98);
  double s = dsp_sin_0_pi(wc);
  double c = dsp_cos_0_pi(wc);
  double t = dsp_tan_neg_halfquarter(0.25 * (wc - kPi));
  double denom = s - c * t;
  if (denom > -1e-12 && denom < 1e-12) denom = (denom >= 0.0) ? 1e-12 : -1e-12;
  return t / denom;
}

static inline double onePoleStep(double* y1, double input, double a) {
  double y0 = input;
  y0 = y0 / (1.0 + y0 * y0);
  *y1 = y0 + a * (y0 - *y1);
  return *y1;
}

struct YellowjacketState {
  bool active;
  double phase;
  double filterY1;
  double oscSelfMod;
  double lastOutValue;
};

static YellowjacketState gPool[kMaxInstances];

}  // namespace

extern "C" int soemdsp_yellowjacket_filter_create() {
  for (int i = 0; i < kMaxInstances; i++) {
    if (!gPool[i].active) {
      YellowjacketState& s = gPool[i];
      s.phase = 0.0;
      s.filterY1 = 0.0;
      s.oscSelfMod = 0.0;
      s.lastOutValue = 0.0;
      s.active = true;
      return i + 1;
    }
  }
  return 0;
}

extern "C" void soemdsp_yellowjacket_filter_destroy(int handle) {
  if (handle < 1 || handle > kMaxInstances) return;
  gPool[handle - 1].active = false;
}

extern "C" double soemdsp_yellowjacket_filter_sample(
  int handle,
  double input,
  double frequency,  // 0..1 normalized
  double resonance,  // 0..1
  double chaosAmount,  // 0..1 (drives filter cutoff scaling, matching original)
  double sampleRate
) {
  if (handle < 1 || handle > kMaxInstances) return 0.0;
  YellowjacketState& s = gPool[handle - 1];

  const double safeRate = sampleRate < 1.0 ? 44100.0 : sampleRate;
  const double freqNorm = clampd(frequency, 0.0, 1.0);
  const double reso = clampd(resonance, 0.0, 1.0);
  const double chaos = clampd(chaosAmount, 0.0, 1.0);

  double maxPitch, resDropPoint;
  if (safeRate <= 44100.0) { maxPitch = 87.7; resDropPoint = 0.77; }
  else if (safeRate <= 88200.0) { maxPitch = 96.0; resDropPoint = 0.82; }
  else if (safeRate <= 132300.0) { maxPitch = 96.0; resDropPoint = 0.83; }
  else if (safeRate <= 176400.0) { maxPitch = 96.0; resDropPoint = 0.86; }
  else if (safeRate <= 220500.0) { maxPitch = 96.0; resDropPoint = 0.89; }
  else if (safeRate <= 264600.0) { maxPitch = 96.0; resDropPoint = 0.90; }
  else { maxPitch = 96.0; resDropPoint = 0.95; }

  const double pitch = jmap01(freqNorm, -156.0, 96.0);
  const double frequencyHz = pitchToFreq(pitch < maxPitch ? pitch : maxPitch);
  const double cutoffHz = frequencyHz * jmap01(chaos, 4.56415, 0.972007);

  const double newResNormalized = evalResVFreqGraph(freqNorm, reso, resDropPoint, 0.2, 0.57);
  const double ellipseC = evalExponentialGraph2(newResNormalized, 7.6024, 0.00001, 0.99);
  const double feedbackGain = evalExponentialGraph2(newResNormalized, 20.0, -0.0429102, 0.99);

  const double a = onePoleCoefficient(cutoffHz, safeRate);

  double inputSignal = clampd(input * 4.0, -7.0, 7.0);
  inputSignal = s.oscSelfMod + 1.04025 * inputSignal + s.lastOutValue;

  const double incAmt = (frequencyHz * 1.9400625 * inputSignal) / safeRate;
  s.phase = s.phase + incAmt;
  s.phase = s.phase - dsp_floor(s.phase);

  double oscValue = waveEllipseFull(s.phase, 0.0, -0.71286768918541499, 0.70129855105756955, ellipseC);
  oscValue *= 0.635417;

  inputSignal = onePoleStep(&s.filterY1, oscValue, a);

  s.oscSelfMod = inputSignal * 20.0;

  const double out = 1.3892758936011171 * oscValue;
  s.lastOutValue = out * 0.5 * feedbackGain;

  return out;
}

extern "C" int soemdsp_yellowjacket_filter_version() {
  return 1;
}
