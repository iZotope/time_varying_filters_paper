/*
Copyright (c) 2014, Aaron S Wishnick
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.
* Neither the name of the Aaron S Wishnick, iZotope, Inc., nor the
  names of its contributors may be used to endorse or promote products
  derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include <iostream>
#include <cmath>
#include <array>
#include <vector>
#include <iomanip>
#include <string>
#include <complex>
#include <functional>
#include <map>
#include <random>
#include <Eigen/Dense>
#include <memory>
#include "sndfile.h"
using namespace std;
using namespace Eigen;

static array<double, 5> lowpass_coefficients(double w0, double q) {
  array<double, 5> coeffs;

  const double cosw0 = cos(w0);
  const double alpha = sin(w0) / (2 * q);

  const double a0 = 1 + alpha;

  coeffs[0] = 0.5 * (1 - cosw0);
  coeffs[1] = 1 - cosw0;
  coeffs[2] = 0.5 * (1 - cosw0);
  coeffs[3] = -2 * cosw0;
  coeffs[4] = 1 - alpha;

  for (auto& x : coeffs) {
    x /= a0;
  }

  return coeffs;
}

static array<double, 5> peaking_coefficients(double w0, double q, double gain) {
  array<double, 5> coeffs;

  const double A = pow(10.0, gain / 40);
  const double alpha = sin(w0) / (2 * q);

  coeffs[0] = 1 + alpha * A;
  coeffs[1] = -2 * cos(w0);
  coeffs[2] = 1 - alpha * A;
  double a0 = 1 + alpha / A;
  coeffs[3] = -2 * cos(w0);
  coeffs[4] = 1 - alpha / A;

  for (auto& x : coeffs) {
    x /= a0;
  }

  return coeffs;
}

static vector<double> generate_sine_wave(int length, double cycles_per_sample) {
  vector<double> signal(length);

  if (cycles_per_sample == 0) {
    signal.assign(length, 1.0);
    return signal;
  }

  for (int i = 0; i < length; ++i) {
    signal[i] = sin(2 * M_PI * i * cycles_per_sample);
  }
  return signal;
}

// Smooth a signal with a linear-phase Hann filter.
static vector<double> smooth_hann(const vector<double>& x, int radius) {
  int n = x.size();
  int m = 2 * radius + 1;

  vector<double> y(n + m - 1);

  vector<double> h(m);
  double norm = 0;
  for (int i = 0; i < m; ++i) {
    h[i] = 0.5 - 0.5 * cos((i + 0.5) * 2.0 * M_PI / m);
    norm += h[i];
  }
  for (int i = 0; i < m; ++i) {
    h[i] /= norm;
  }

  for (int i = 0; i < n + m - 1; ++i) {
    double sum = 0;
    for (int j = 0; j < m; ++j) {
      int idx = radius + i - j;
      if (idx < 0) {
        break;
      }
      sum += h[j] * x[idx];
    }
    y[i] = sum;
  }

  return y;
}

// Generate a Hann window.
static vector<double> generate_hann_window(int n) {
  vector<double> w(n);
  for (int i = 0; i < n; ++i) {
    w[i] = 0.5 - 0.5 * cos(2.0 * M_PI * i / (n - 1));
  }
  return w;
}

// Compute the gain of a biquad at frequency w0.
static double biquad_gain(array<double, 5> coeffs, double w0) {
  const double b0 = coeffs[0], b1 = coeffs[1], b2 = coeffs[2];
  const double a1 = coeffs[3], a2 = coeffs[4];

  auto z = polar(1.0, w0);
  auto z2 = z * z;

  auto h = (b0 + b1 / z + b2 / z2) / (1.0 + a1 / z + a2 / z2);
  return abs(h);
}

// High-quality reference signal. Assumes that the filter input is a sinusoid.
// Smooths the gain of the filter coefficients, and applies it to the signal.
static vector<double> reference(vector<double> x,
                                const vector<array<double, 5>>& coeffs,
                                double w0, int radius) {
  const int n = x.size();

  vector<double> gains(n);
  transform(begin(coeffs), end(coeffs), begin(gains),
            [w0](array<double, 5> coeffs) { return biquad_gain(coeffs, w0); });

  gains = smooth_hann(gains, radius);

  transform(begin(x), end(x), begin(gains), begin(x),
            [](double sig, double gain) { return sig * gain; });

  return x;
}

// Low-quality anchor. Applies the gain to the signal, with no smoothing, and
// inserts an impulse whenever coefficients change.
static vector<double>
anchor(vector<double> x, const vector<array<double, 5>>& coeffs, double w0) {
  const int n = x.size();

  vector<double> gains(n);
  transform(begin(coeffs), end(coeffs), begin(gains),
            [w0](array<double, 5> coeffs) { return biquad_gain(coeffs, w0); });

  double max_gain = accumulate(begin(gains), end(gains), 0.0,
                               [](double a, double b) { return max(a, b); });
  double transient_gain = 3 * max_gain;

  for (int i = 0; i < n; ++i) {
    if (i != 0 && gains[i] != gains[i - 1]) {
      // If there was a change in gain, insert a click.
      x[i] += transient_gain;
    } else {
      x[i] *= gains[i];
    }
  }

  return x;
}

// Vanilla Direct Form II implementation.
static vector<double> time_varying_df2(vector<double> x,
                                       const vector<array<double, 5>>& coeffs) {
  double w1 = 0, w2 = 0;
  const int n = x.size();
  for (int i = 0; i < n; ++i) {
    auto cur_coeffs = coeffs[i];

    double w = x[i] - cur_coeffs[3] * w1 - cur_coeffs[4] * w2;
    double y = cur_coeffs[0] * w + cur_coeffs[1] * w1 + cur_coeffs[2] * w2;
    x[i] = y;

    w2 = w1;
    w1 = w;
  }

  return x;
}

// Vanilla SVF implementation.
static vector<double> time_varying_svf(vector<double> x,
                                       const vector<array<double, 5>>& coeffs) {
  double s1 = 0, s2 = 0;
  const int n = x.size();

  for (int i = 0; i < n; ++i) {
    const double b0 = coeffs[i][0], b1 = coeffs[i][1], b2 = coeffs[i][2];
    const double a1 = coeffs[i][3], a2 = coeffs[i][4];

    const complex<double> neg_sqrt = sqrt(complex<double>(-1 - a1 - a2));
    const complex<double> pos_sqrt = sqrt(complex<double>(-1 + a1 - a2));

    const double g = (neg_sqrt / pos_sqrt).real();
    const double R = (a2 - 1) / (neg_sqrt * pos_sqrt).real();
    const double TwoR = 2 * R;
    const double g2 = g * g;
    const double chp = (b0 - b1 + b2) / (1 - a1 + a2);
    const double cbp = -2 * (b0 - b2) / (neg_sqrt * pos_sqrt).real();
    const double clp = (b0 + b1 + b2) / (1 + a1 + a2);

    double u1 = (x[i] - (g + TwoR) * s1 - s2) / (1 + g2 + TwoR * g);
    double u2 = g * u1 + s1;
    double u3 = g * u2 + s2;

    s1 += 2.0 * g * u1;
    s2 += 2.0 * g * u2;

    x[i] = chp * u1 + cbp * u2 + clp * u3;
  }

  return x;
}

// SVF, with the stabilization method from "Stability of Recursive Time-Varying
// Digital Filters by State Vector Transformation."
static vector<double>
svf_rabenstein_czarnach(vector<double> x,
                        const vector<array<double, 5>>& coeffs) {
  Matrix<double, 2, 1> s = Matrix<double, 2, 1>::Zero();
  const int n = x.size();

  double R_prev;

  for (int i = 0; i < n; ++i) {
    const double b0 = coeffs[i][0], b1 = coeffs[i][1], b2 = coeffs[i][2];
    const double a1 = coeffs[i][3], a2 = coeffs[i][4];

    const complex<double> neg_sqrt = sqrt(complex<double>(-1 - a1 - a2));
    const complex<double> pos_sqrt = sqrt(complex<double>(-1 + a1 - a2));

    const double g = (neg_sqrt / pos_sqrt).real();
    const double R = (a2 - 1) / (neg_sqrt * pos_sqrt).real();
    const double g2 = g * g;

    const double chp = (b0 - b1 + b2) / (1 - a1 + a2);
    const double cbp = -2 * (b0 - b2) / (neg_sqrt * pos_sqrt).real();
    const double clp = (b0 + b1 + b2) / (1 + a1 + a2);

    Matrix<double, 2, 2> A;
    A(0, 0) = -2 * R;
    A(0, 1) = -1;
    A(1, 0) = 1;
    A(1, 1) = 0;

    Matrix<double, 2, 1> B;
    B(0, 0) = 1;
    B(1, 0) = 0;

    Matrix<double, 2, 2> H =
        (Matrix<double, 2, 2>::Identity() - g * A).inverse();

    Matrix<double, 2, 2> Pa = Matrix<double, 2, 2>::Identity() + 2 * g * A * H;
    Matrix<double, 2, 1> Pb = 2 * g2 * A * H * B + 2 * g * B;
    Matrix<double, 2, 1> Pc;
    Pc(0, 0) = cbp - 2 * R * chp;
    Pc(1, 0) = clp - chp;
    double Pd = chp;

    Matrix<double, 2, 2> TT;
    double Rn = i != 0 ? R_prev : R;
    TT(0, 0) = (sqrt(complex<double>(Rn * Rn - 1)) /
                sqrt(complex<double>(R * R - 1))).real();
    TT(0, 1) = -(complex<double>(0, 1) * (R - Rn) /
                 sqrt(complex<double>(R * R - 1))).real();
    TT(1, 0) = 0;
    TT(1, 1) = 1;

    s = TT * s;

    // Compute output.
    double y = Pc.transpose() * s + Pd * x[i];

    // State update
    s = Pa * s + Pb * x[i];

    R_prev = R;

    x[i] = y;
  }

  return x;
}

template <int N>
static Matrix<double, N, N>
ComputeRabensteinK(Matrix<double, N, N> A1, Matrix<double, N, 1> b1,
                   Matrix<double, N, N> A2, Matrix<double, N, 1> b2) {
  Matrix<double, N, N> K = b1 * b2.transpose();

  auto Ae1 = A1;
  auto Ae2 = A2;
  int i = 0;
  while (i++ < 1000000 && (Ae1.template lpNorm<Infinity>() > 1e-12 ||
                           Ae2.template lpNorm<Infinity>() > 1e-12)) {
    Matrix<double, N, 1> f1 = Ae1 * b1;
    Matrix<double, N, 1> f2 = Ae2 * b2;

    K += f1 * f2.transpose();

    Ae1 *= A1;
    Ae2 *= A2;
  }

  return K;
}

// Direct Form 2 Transposed
static vector<double> tdf2(vector<double> x,
                           const vector<array<double, 5>>& coeffs) {
  const int n = x.size();
  Matrix<double, 2, 1> s = Matrix<double, 2, 1>::Zero();

  for (int i = 0; i < n; ++i) {
    const double b0 = coeffs[i][0], b1 = coeffs[i][1], b2 = coeffs[i][2];
    const double a1 = coeffs[i][3], a2 = coeffs[i][4];

    Matrix<double, 2, 2> A;
    A << -a1, 1.0, -a2, 0.0;

    Matrix<double, 2, 1> b;
    b << b1 - b0* a1, b2 - b0* a2;

    Matrix<double, 2, 1> c;
    c << 1.0, 0.0;

    double d = b2;

    double y = c.transpose() * s + d * x[i];
    s = A * s + b * x[i];

    x[i] = y;
  }

  return x;
}

static vector<double>
tdf2_rabenstein_czarnach(vector<double> x,
                         const vector<array<double, 5>>& coeffs) {
  const int n = x.size();
  Matrix<double, 2, 1> s = Matrix<double, 2, 1>::Zero();

  for (int i = 0; i < n; ++i) {
    const double b0 = coeffs[i][0], b1 = coeffs[i][1], b2 = coeffs[i][2];
    const double a1 = coeffs[i][3], a2 = coeffs[i][4];

    Matrix<double, 2, 2> A;
    A << -a1, 1.0, -a2, 0.0;

    Matrix<double, 2, 1> b;
    b << b1 - b0* a1, b2 - b0* a2;

    Matrix<double, 2, 1> c;
    c << 1.0, 0.0;

    double d = b2;

    double r = sqrt(a2);
    double t = acos(a1 / (-2 * r));

    auto prev_coeffs = i != 0 ? coeffs[i - 1] : coeffs[i];
    double rp = sqrt(prev_coeffs[4]);
    double tp = acos(prev_coeffs[3] / (-2 * rp));

    Matrix<double, 2, 2> TT;
    TT << 1, 0, r* sin(t - tp) / sin(tp), r * sin(t) / (rp * sin(tp));
    s = TT * s;

    double y = c.transpose() * s + d * x[i];
    s = A * s + b * x[i];

    x[i] = y;
  }

  return x;
}

// SVF, with transient minimization as "Minimization of Transient Signals in
// Recursive Time-Varying Digital Filters".
static vector<double> svf_rabenstein(vector<double> x,
                                     const vector<array<double, 5>>& coeffs) {
  Matrix<double, 2, 1> s = Matrix<double, 2, 1>::Zero();
  const int n = x.size();

  Matrix<double, 2, 2> Pa_prev;
  Matrix<double, 2, 1> Pb_prev;

  for (int i = 0; i < n; ++i) {
    const double b0 = coeffs[i][0], b1 = coeffs[i][1], b2 = coeffs[i][2];
    const double a1 = coeffs[i][3], a2 = coeffs[i][4];

    const complex<double> neg_sqrt = sqrt(complex<double>(-1 - a1 - a2));
    const complex<double> pos_sqrt = sqrt(complex<double>(-1 + a1 - a2));

    const double g = (neg_sqrt / pos_sqrt).real();
    const double R = (a2 - 1) / (neg_sqrt * pos_sqrt).real();
    const double g2 = g * g;

    const double chp = (b0 - b1 + b2) / (1 - a1 + a2);
    const double cbp = -2 * (b0 - b2) / (neg_sqrt * pos_sqrt).real();
    const double clp = (b0 + b1 + b2) / (1 + a1 + a2);

    Matrix<double, 2, 2> A;
    A(0, 0) = -2 * R;
    A(0, 1) = -1;
    A(1, 0) = 1;
    A(1, 1) = 0;

    Matrix<double, 2, 1> B;
    B(0, 0) = 1;
    B(1, 0) = 0;

    Matrix<double, 2, 2> H =
        (Matrix<double, 2, 2>::Identity() - g * A).inverse();

    Matrix<double, 2, 2> Pa = Matrix<double, 2, 2>::Identity() + 2 * g * A * H;
    Matrix<double, 2, 1> Pb = 2 * g2 * A * H * B + 2 * g * B;
    Matrix<double, 2, 1> Pc;
    Pc(0, 0) = cbp - 2 * R * chp;
    Pc(1, 0) = clp - chp;
    double Pd = chp;

    if (i != 0 && coeffs[i] != coeffs[i - 1]) {
      auto Knn = ComputeRabensteinK(Pa_prev, Pb_prev, Pa_prev, Pb_prev);
      auto Knp = ComputeRabensteinK(Pa_prev, Pb_prev, Pa, Pb);

      Pa_prev = Pa;
      Pb_prev = Pb;

      Matrix<double, 2, 2> T = Knp.transpose() * Knn.inverse();
      Pa = Pa * T;
      Pc = (Pc.transpose() * T).transpose();
    } else {
      Pa_prev = Pa;
      Pb_prev = Pb;
    }

    // Compute output.
    double y = Pc.transpose() * s + Pd * x[i];

    // State update
    s = Pa * s + Pb * x[i];

    x[i] = y;
  }

  return x;
}

// Truncated output switching, as in "Suppression of Transients in Time-Varying
// Recursive Filters for Audio Signals".
static vector<double> truncated_output_switching_filter(
    vector<double> x, const vector<array<double, 5>> coeffs, int Na) {
  const int n = x.size();
  vector<double> ys(n);
  for (int i = 0; i < n; ++i) {
    auto cur_coeffs = coeffs[i];
    double w1 = 0, w2 = 0;
    double y = 0;

    int start = max(i - Na + 1, 0);
    for (int j = start; j <= i; ++j) {
      double w = x[j] - cur_coeffs[3] * w1 - cur_coeffs[4] * w2;
      y = cur_coeffs[0] * w + cur_coeffs[1] * w1 + cur_coeffs[2] * w2;

      w2 = w1;
      w1 = w;
    }

    ys[i] = y;
  }

  return ys;
}

// Output switching, equivalent to Zetterberg-Zhang.
static vector<double>
output_switching_filter(vector<double> x,
                        const vector<array<double, 5>> coeffs) {
  return truncated_output_switching_filter(x, coeffs, x.size() + 1);
  const int n = x.size();
  vector<double> ys(n);
  for (int i = 0; i < n; ++i) {
    auto cur_coeffs = coeffs[i];
    double w1 = 0, w2 = 0;
    double y = 0;

    for (int j = 0; j <= i; ++j) {
      double w = x[j] - cur_coeffs[3] * w1 - cur_coeffs[4] * w2;
      y = cur_coeffs[0] * w + cur_coeffs[1] * w1 + cur_coeffs[2] * w2;

      w2 = w1;
      w1 = w;
    }

    ys[i] = y;
  }

  return ys;
}

// Gold/Rader or coupled/normal form.
static vector<double> gold_rader(vector<double> x,
                                 const vector<array<double, 5>> coeffs) {
  const int n = x.size();
  double u1 = 0, u2 = 0;

  for (int i = 0; i < n; ++i) {
    const double a1 = coeffs[i][3], a2 = coeffs[i][4];
    double alpha = -a1 / 2;
    double beta = -0.5 * sqrt(4 * a2 - a1 * a1);

    const double b0 = coeffs[i][0], b1 = coeffs[i][1], b2 = coeffs[i][2];
    double alpha2 = alpha * alpha, beta2 = beta * beta;
    double k1 = b0 - b2 / (alpha2 + beta2);
    double k2 = (b2 * alpha + (b1 + b0 * alpha) * (alpha2 + beta2)) /
                (beta * (alpha2 + beta2));
    double k3 = b2 / (alpha2 + beta2);

    double y = k1 * u1 + k2 * u2 + k3 * x[i];

    double u1_new = alpha * u1 - beta * u2 + x[i];
    double u2_new = beta * u1 + alpha * u2;
    u1 = u1_new;
    u2 = u2_new;
    x[i] = y;
  }

  return x;
}

static void write_wav(const char* path, double sampling_rate, int channel_count,
                      const vector<double>& samples) {
  SF_INFO sfinfo;
  sfinfo.samplerate = static_cast<int>(sampling_rate + 0.5);
  sfinfo.channels = channel_count;
  sfinfo.format = SF_FORMAT_WAV | SF_FORMAT_FLOAT;
  unique_ptr<SNDFILE, int (*)(SNDFILE*)> file(sf_open(path, SFM_WRITE, &sfinfo),
                                              &sf_close);

  if (!file) {
    cerr << "Error opening \"" << path << "\" for writing. "
         << sf_strerror(nullptr) << endl;
    return;
  }

  // Make the data multichannel.
  vector<float> multichannel_samples;
  multichannel_samples.reserve(channel_count * samples.size());
  for (auto x : samples) {
    for (int i = 0; i < channel_count; ++i) {
      multichannel_samples.push_back(x);
    }
  }

  sf_write_float(file.get(), multichannel_samples.data(),
                 multichannel_samples.size());
}

int main(int argc, char* argv[]) {
  if (argc != 4) {
    cerr << "Usage: " << argv[0]
         << " [filter type] [sine_hz|dc] [output_dir]\n";
    return -1;
  }

  double cycles_per_second;
  bool dc_input;
  try {
    if (strcmp(argv[2], "dc") == 0) {
      dc_input = true;
      cycles_per_second = 100;
    } else {
      dc_input = false;
      cycles_per_second = stod(argv[2]);
    }
  }
  catch (...) {
    cerr << "Invalid stimulus \"" << argv[2] << "\".\n";
    return -2;
  }

  string output_dir(argv[3]);

  const double sampling_rate = 48000;
  const double cycles_per_sample = cycles_per_second / sampling_rate;
  const double q = 6.0;
  const double cutoff_scale = 0.2;
  const double start_cutoff =
      max(cycles_per_second, 100.0) * (1 - cutoff_scale);
  const double end_cutoff = max(cycles_per_second, 100.0) * (1 + cutoff_scale);
  const int constant_length = 48000;
  const int end_constant_length = sampling_rate;
  const int start_vary_pos = constant_length;
  const int max_ramp_length = 200;
  const int length = constant_length + max_ramp_length + end_constant_length;

  auto ramp_fn = [=](int i, int ramp_length, double start_val, double end_val) {
    if (i < start_vary_pos) {
      return start_val;
    }
    if (i < start_vary_pos + ramp_length) {
      return start_val + static_cast<double>(i - start_vary_pos) / ramp_length *
                             (end_val - start_val);
    }
    return end_val;
  };

  const double start_w0 = start_cutoff * 2 * M_PI / sampling_rate;
  const double end_w0 = end_cutoff * 2 * M_PI / sampling_rate;
  auto vary_freq_fn = [=](int i, int ramp_length = 0) {
    return ramp_fn(i, ramp_length, start_w0, end_w0);
  };

  const double start_gain = -4, end_gain = 4;
  const double w0 = cycles_per_second * 2 * M_PI / sampling_rate;
  auto vary_gain_fn = [=](int i, int ramp_length = 0) {
    return ramp_fn(i, ramp_length, start_gain, end_gain);
  };

  const double start_q = 0.6, end_q = 4;
  auto vary_q_fn = [=](int i, int ramp_length = 0) {
    return ramp_fn(i, ramp_length, start_q, end_q);
  };

  auto input_signal =
      generate_sine_wave(length, dc_input ? 0 : cycles_per_sample);

  map<string, function<array<double, 5>(int i)>> filter_types = {
      {"lowpass_freq",
       [=](int i) { return lowpass_coefficients(vary_freq_fn(i), q); }},
      {"lowpass_q",
       [=](int i) { return lowpass_coefficients(w0, vary_q_fn(i)); }},
      {"peaking_freq", [=](int i) {
        return peaking_coefficients(vary_freq_fn(i), q, end_gain);
      }},
      {"peaking_gain",
       [=](int i) { return peaking_coefficients(w0, q, vary_gain_fn(i)); }},
      {"peaking_q", [=](int i) {
        return peaking_coefficients(end_w0, vary_q_fn(i), end_gain);
      }},
  };

  auto filter_iter = filter_types.find(argv[1]);
  if (filter_iter == end(filter_types)) {
    cerr << "Invalid filter type \"" << argv[1] << "\". Types:\n";
    for (const auto& p : filter_types) {
      cerr << "\t" << p.first << endl;
    }
    return -1;
  }
  auto coeffs_fn = filter_iter->second;

  vector<array<double, 5>> coeffs(length);
  for (int i = 0; i < length; ++i) {
    coeffs[i] = coeffs_fn(i);
  }

  typedef function<vector<double>(vector<double>,
                                  const vector<array<double, 5>>&)> filter_fn;
  map<string, filter_fn> methods = {
      {"df2", &time_varying_df2},
      {"gr", &gold_rader},
      {"svf", &time_varying_svf},
      {"svf_r", &svf_rabenstein},
      {"svf_rc", &svf_rabenstein_czarnach},
      {"tdf2", &tdf2},
      {"tdf2_rc", &tdf2_rabenstein_czarnach},
      {"zz", &output_switching_filter},
      {"reference",
       [=](vector<double> x, const vector<array<double, 5>>& coeffs) {
        return reference(std::move(x), coeffs, w0,
                         static_cast<int>(sampling_rate * 0.01));
      }},
      {
       "anchor", [=](vector<double> x, const vector<array<double, 5>>& coeffs) {
                   return anchor(std::move(x), coeffs, w0);
                 },
      }};

  // Reference and anchor aren't needed for testing DC input.
  if (dc_input) {
    methods.erase("reference");
    methods.erase("anchor");
  }

  // Process with each filter structure.
  map<string, vector<double>> outputs;
  transform(begin(methods), end(methods), inserter(outputs, outputs.end()),
            [&](pair<string, filter_fn> cur_method) {
    return make_pair(cur_method.first, cur_method.second(input_signal, coeffs));
  });

  // Normalize peak levels.
  double max_level =
      accumulate(begin(outputs), end(outputs), 0.0,
                 [](double cur, const pair<string, vector<double>>& p) {
        double new_max_level =
            accumulate(begin(p.second), end(p.second), 0.0,
                       [](double cur, double x) { return max(cur, abs(x)); });
        return max(cur, new_max_level);
      });
  for (auto& p : outputs) {
    transform(begin(p.second), end(p.second), begin(p.second),
              [max_level](double x) { return x / max_level; });
  }

  // Apply fade-in and fade-out, so the signal starting and ending doesn't
  // distract from the filter transition.
  // This only gets done for sinusoidal input, since it's used for the
  // subjective tests. For the objective tests, the transition doesn't matter
  // because there's no listening component.
  if (!dc_input) {
    auto fade_samples = static_cast<int>(sampling_rate * 0.05 + 0.5);
    auto hann_window = generate_hann_window(2 * fade_samples);
    hann_window.resize(fade_samples);
    for (auto& p : outputs) {
      auto& xs = p.second;
      transform(begin(xs), begin(xs) + fade_samples, begin(hann_window),
                begin(xs), std::multiplies<double>());
      transform(end(xs) - fade_samples, end(xs), begin(hann_window),
                end(xs) - fade_samples,
                [](double x, double w) { return x * (1.0 - w); });
    }
  }

  // Write the output.
  const int channel_count = 2;
  for (auto& output : outputs) {
    auto output_path = output_dir + "/" + output.first + ".wav";
    write_wav(output_path.c_str(), sampling_rate, channel_count, output.second);
  }
}
