// soemdsp-native-module: rsmet_filter
// soemdsp-native-label: RSMET Filter
// soemdsp-native-target: rsmetFilter
// soemdsp-native-kind: filter
//
// A general-purpose ladder filter (the same design as this repo's
// ladder_filter module) preceded by a tanh soft clipper and noise
// injection stage, matching the original RsmetFilter family (10 modes:
// LP6/12/18/24, HP6/12/18/24, BP6, BP12). Frequency and resonance are
// mapped through exact exponential easing curves
// (soemdsp::utility::Graph 2-node EXPONENTIAL shape) rather than being
// linear, matching the original design.

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

// Natural log for x > 0, via IEEE754 exponent extraction + an atanh-series
// on the mantissa (m in [1,2), t=(m-1)/(m+1) in [0, 1/3], fast convergence).
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

static inline double dsp_tanh(double x) {
  return 1.0 - 2.0 / (dsp_exp(2.0 * x) + 1.0);
}

static inline double clampd(double v, double lo, double hi) {
  return v < lo ? lo : (v > hi ? hi : v);
}

// Exact soemdsp::utility::Graph::getValue for a 2-node EXPONENTIAL curve
// (nodes at x=0 and x=1). Clamps outside [0,1] like Graph::getValue does.
static double evalExponentialGraph2(double x, double y0, double y1, double skew) {
  if (x <= 0.0) return y0;
  if (x >= 1.0) return y1;
  const double c = 0.5 * (skew + 1.0);
  const double a = 2.0 * dsp_ln((1.0 - c) / c);
  const double expA = dsp_exp(a);
  const double expIA = dsp_exp(x * a);
  return y0 + (y1 - y0) * (1.0 - expIA) / (1.0 - expA);
}

static void compute_mix(int mode, int stages, double c[5], double* s_out) {
  for (int i = 0; i < 5; i++) c[i] = 0.0;
  if (mode == 0) {
    c[0] = 1.0;
    *s_out = 0.125;
  } else if (mode == 1) {
    c[stages] = 1.0;
    *s_out = stages * 0.25;
  } else if (mode == 2) {
    static const double hp[4][5] = {
      {1.0, -1.0,  0.0,  0.0, 0.0},
      {1.0, -2.0,  1.0,  0.0, 0.0},
      {1.0, -3.0,  3.0, -1.0, 0.0},
      {1.0, -4.0,  6.0, -4.0, 1.0},
    };
    for (int i = 0; i <= stages; i++) c[i] = hp[stages - 1][i];
    *s_out = stages * 0.25;
  } else {
    static const double bp[4][5] = {
      {0.0, 2.0, -2.0,  0.0, 0.0},
      {0.0, 2.0, -2.0,  0.0, 0.0},
      {0.0, 0.0,  3.0, -3.0, 0.0},
      {0.0, 0.0,  4.0, -8.0, 4.0},
    };
    for (int i = 0; i < 5; i++) c[i] = bp[stages - 1][i];
    *s_out = 0.125;
  }
}

static inline double safe(double x) {
  return (x * 0.0 == 0.0) ? x : 0.0;
}

static inline double nextNoiseBipolar(unsigned int* state) {
  unsigned int x = *state;
  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 5;
  *state = x;
  return ((double)x / 4294967295.0) * 2.0 - 1.0;
}

// rsmet_filter's 10 modes, each mapping to a (ladderMode, stages) pair on
// the shared ladder core -- matches the original Rsmet_LP6..Rsmet_BP12
// subclasses exactly (each was just a thin wrapper picking one ladder tap).
static void modeToLadder(int rsmetMode, int* ladderMode, int* stages) {
  static const int table[10][2] = {
    {1, 1}, {1, 2}, {1, 3}, {1, 4},  // LP6, LP12, LP18, LP24
    {2, 1}, {2, 2}, {2, 3}, {2, 4},  // HP6, HP12, HP18, HP24
    {3, 1}, {3, 4},                  // BP6 (BP_6_6), BP12 (BP_12_12)
  };
  int idx = rsmetMode < 0 ? 0 : (rsmetMode > 9 ? 9 : rsmetMode);
  *ladderMode = table[idx][0];
  *stages = table[idx][1];
}

struct RsmetState {
  double y[5];
  bool active;
  unsigned int rngState;
};

static RsmetState gPool[kMaxInstances];

}  // namespace

extern "C" int soemdsp_rsmet_filter_create() {
  for (int i = 0; i < kMaxInstances; i++) {
    if (!gPool[i].active) {
      RsmetState& s = gPool[i];
      for (int j = 0; j < 5; j++) s.y[j] = 0.0;
      s.rngState = 0x2545F491u + (unsigned int)(i + 1) * 2246822519u;
      s.active = true;
      return i + 1;
    }
  }
  return 0;
}

extern "C" void soemdsp_rsmet_filter_destroy(int handle) {
  if (handle < 1 || handle > kMaxInstances) return;
  gPool[handle - 1].active = false;
}

extern "C" double soemdsp_rsmet_filter_sample(
  int handle,
  double input,
  double frequency,  // 0..1 normalized
  double resonance,  // 0..1 normalized
  double chaosAmount,  // 0..1
  int mode,          // 0..9, see modeToLadder
  double sampleRate
) {
  if (handle < 1 || handle > kMaxInstances) return 0.0;
  RsmetState& s = gPool[handle - 1];

  const double safeRate = sampleRate < 1.0 ? 44100.0 : sampleRate;
  const double freqNorm = clampd(frequency, 0.0, 1.0);
  const double resoNorm = clampd(resonance, 0.0, 1.0);
  const double chaos = clampd(chaosAmount, 0.0, 1.0);

  const double cutoffHz = clampd(evalExponentialGraph2(freqNorm, 3.0, 20000.0, -0.95), 0.000001, safeRate * 0.49);
  const double feedback = clampd(evalExponentialGraph2(resoNorm, 0.0, 1.0, 0.5), 0.0, 0.999);

  int ladderMode, stages;
  modeToLadder(mode, &ladderMode, &stages);

  const double wc = clampd((2.0 * kPi * cutoffHz) / safeRate, 1e-9, kPi * 0.98);
  const double sine = dsp_sin_0_pi(wc);
  const double cosine = dsp_cos_0_pi(wc);
  const double tangent = dsp_tan_neg_halfquarter(0.25 * (wc - kPi));
  double a = sine - cosine * tangent;
  a = (a > -1e-12 && a < 1e-12) ? (a >= 0.0 ? 1e-12 : -1e-12) : a;
  a = tangent / a;

  double c[5];
  double mixS;
  compute_mix(ladderMode, stages, c, &mixS);

  const double b = 1.0 + a;
  const double denom = 1.0 + a * a + 2.0 * a * cosine;
  const double safeDenom = denom < 1e-12 ? 1e-12 : denom;
  const double g2 = (b * b) / safeDenom;
  const double g2sq = g2 * g2 < 1e-12 ? 1e-12 : g2 * g2;
  const double k = feedback / g2sq;
  const double g = 1.0 + mixS * k;

  double inputSignal = dsp_tanh(input * 2.0);
  if (chaos > 0.0) {
    inputSignal += nextNoiseBipolar(&s.rngState) * chaos;
  }

  const double safeIn = safe(g * inputSignal - k * s.y[4]);
  double y0 = safeIn / (1.0 + safeIn * safeIn);
  const double ny1 = safe(y0 + a * (y0 - s.y[1]));
  const double ny2 = safe(ny1 + a * (ny1 - s.y[2]));
  const double ny3 = safe(ny2 + a * (ny2 - s.y[3]));
  const double ny4 = safe(ny3 + a * (ny3 - s.y[4]));
  s.y[0] = safe(y0);
  s.y[1] = ny1;
  s.y[2] = ny2;
  s.y[3] = ny3;
  s.y[4] = ny4;

  const double out = c[0] * s.y[0] + c[1] * s.y[1] + c[2] * s.y[2] + c[3] * s.y[3] + c[4] * s.y[4];
  return safe(out) * 0.41;
}

extern "C" int soemdsp_rsmet_filter_version() {
  return 1;
}
