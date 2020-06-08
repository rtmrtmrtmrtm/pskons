#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <algorithm>
#include <ctype.h>
#include <complex>
#include "fft.h"
#include "util.h"
#include "demod.h"
#include "taps.h"

//
// "PSK31: A New Radio-Teletype Mode", Peter Martinez, G3PLX,
// QEX July/Aug 1999.
//
// Also Moe Wheatley's PSKcore documentation and filter taps.
//
// by Robert Morris, AB1HL.
//

//
// adjustable parameters.
//

double max_hz = 3000;

double coarse_blocks = 3;
double coarse_weight = 0.5;
double coarse_top_n = 15;

double sync_buckets_coarse = 6;
double sync_buckets_fine = 12;
double sync_weight_coarse = 0.08;
double sync_weight_fine = 0.1;
double sync_thresh = 0.60;

double hz_fine_thresh = 0.65; // 0.57; // switch to fine hz adjustment
double coarse_hz_n = 4;
double fine_hz_n = 3;
double hz_fine = 0.07;
double hz_weight_coarse = 0.03;
double hz_weight_fine = 0.5;

double infant_time = 0.2; // grace time for new Signals
double keep_time = 10.0;   // grace time for Signals that were once good

// adjust_hz_phase()
double hz_phase_thresh = 0.85; // s->hz_phase_q_
double hz_phase_max = 0.25; // hz
double hz_phase_weight = 0.8;

// weights for various weighted averages of q
double hz_q_weight = 0.01;
double phase_q_weight = 0.1;
double sync_q_weight = 0.3;
double cull_q_weight = 0.3;
double good_q_weight = 0.25;

double cull_q = 0.6;      // drop Signals weaker than this
double good_q = 0.8;      // keep for at least keep_time

Demod::Demod(int rate)
{
  rate_ = rate;
  block_ = round(rate_ / 31.25);
  assert(fabs((rate_/31.25) - block_) < 0.00001);

  clock_ = 0;

  // recent signal strength, used by coarse(), to average.
  recent_.resize(4096);
  for(ulong i = 0; i < recent_.size(); i++)
    recent_[i] = 0;

  // simple one-symbol-long pulse shape.
  shape_ = raised_cosine(block_);

  // more serious raised-cosine low-pass filter.
  assert(rate == 8000); // for taps
  taps_ = taps();

  sigx_ = 0;
  sign_ = 0;
}

Demod::~Demod()
{
  for(ulong i = 0; i < signals_.size(); i++){
    delete signals_[i];
  }
}

void
Demod::clear_watches()
{
  watch_.clear();
}

//
// ui.cc wants to know about decodes near this hz.
//
void
Demod::watch(double hz, watch_cb_t cb, void *arg)
{
  Watch w;
  w.hz = hz;
  w.cb = cb;
  w.arg = arg;
  watch_.push_back(w);
}

//
// look for candidate signals -- the FFT buckets
// with highest amplitude. FFT with buckets narrower
// than 31.25, and sum them up.
//
void
Demod::coarse(const std::vector<double> &samples, int i0)
{
  double bin_hz = 31.25 / coarse_blocks;

  ulong nn = block_ * coarse_blocks;
  assert(nn <= samples.size());

  int start = (i0 + block_/2) - nn / 2;
  if(start < 0)
    start = 0;

  assert(start >= 0 && start + nn <= samples.size());

  std::vector<std::complex<double>> bins = one_fft(samples, start, nn);

  for(ulong i = 0; i < bins.size() && i*bin_hz < max_hz; i++){
    double a = std::abs(bins[i]);
    recent_[i] = (1-coarse_weight)*recent_[i] + coarse_weight*a;
  }

  //
  // sum up adjacent 31.25-hz worth of the small recent_ bins.
  //

  struct Bin {
    double hz;
    double str;
  };
  std::vector<Bin> r;

  ulong cb = round(coarse_blocks);

  for(ulong i = 0; i < bins.size() - cb && i*bin_hz < max_hz; i++){
    if(i * bin_hz < 150)
      continue;
    double sum = 0;
    for(ulong j = i; j < i + cb; j++){
      sum += recent_[j];
    }
    Bin b;
    b.str = sum;
    if((cb % 2) == 0){
      b.hz = (i + cb/2) * bin_hz - bin_hz/2;
    } else {
      b.hz = (i + cb/2) * bin_hz;
    }
    r.push_back(b);
  }

  std::sort(r.begin(), r.end(),
            [ ](const Bin &a, const Bin &b) -> bool
            { return a.str > b.str; } );

  double max_excursion = 31.25 / (2 * coarse_blocks);

  // ensure strongest are in Signals.
  mu_.lock();
  for(int i = 0; i < coarse_top_n; i++){
    double hz = r[i].hz;
    int got = 0;
    for(ulong j = 0; j < signals_.size(); j++){
      if(fabs(hz - signals_[j]->hz0_) < max_excursion ||
         fabs(hz - signals_[j]->hz_) < max_excursion){
        got = 1;
      }
    }
    if(got == 0){
      Signal *s = new Signal(hz);
      s->start_ = clock_;
      signals_.push_back(s);
    }
  }
  mu_.unlock();
}

//
// get rid of lowest-quality Signals.
//
void
Demod::cull()
{
  
  std::vector<Signal*> ns;

  mu_.lock();

  sigx_ += signals_.size();
  sign_ += 1;

  for(ulong i = 0; i < signals_.size(); i++){
    Signal *s = signals_[i];

    int keep = 0;

    if(clock_ - s->start_ < infant_time){
      keep = 1;
    } else if(clock_ - s->last_good_ < keep_time){
      keep = 1;
    } else if(s->cull_q_ > cull_q){
      keep = 1;
    } else {
      keep = 0;
    }

    if(s->good_q_ > good_q){
      s->last_good_ = clock_;
    }
    
    if(keep){
      ns.push_back(s);
    } else {
      delete s;
    }
  }
  signals_ = ns;

  mu_.unlock();
}

void
Demod::got(const std::vector<double> &samples)
{
  // append at most one block at a time,
  // so we can have each Signal inspect
  // the input a block at a time in the
  // middle of a 5-block window.
  for(ulong i = 0; i < samples.size(); ){
    ulong n = samples.size() - i;
    if(n > block_)
      n = block_;

    samples_.insert(samples_.end(), samples.begin() + i, samples.begin() + i + n);

    // has an entire block accumulated since
    // we last processed one? the last processed
    // signal (if any) started at samples_[block_].
    // we'd like to process next the symbol at samples_[2*block_].

    while(samples_.size() >= 5*block_){
      coarse(samples_, block_*2);
      coarse(samples_, block_*2 + block_/2);

      tick_all();
      samples_.erase(samples_.begin(), samples_.begin() + block_);

      clock_ += block_ / (double) rate_;

      cull();
    }

    i += n;
  }
}

void
Demod::tick_all()
{
  for(ulong i = 0; i < signals_.size(); i++){
    tick(signals_[i]);
  }
}

//
// basically does a one-bit FFT for various sub-symbol
// offsets, averages the power seen at each offset,
// and sets s->sync_ to the offset with the most power.
// should be called with off=2*block_ when this symbol
// and the previous have opposite phases (i.e. for "0"),
// in order that an incorrect offset will see some
// samples with the opposite (canceling) phase.
// the filter array idea for different offsets is
// from Moe Wheatley's PSKCore.
//
void
Demod::adjust_sync(const std::vector<double> &samples, Signal *s, int off)
{
  // initial phase of reference tone.
  double tone_phase = 0.0;

  // reference phase advances by this much each sample.
  double tone_inc = 2 * M_PI / (rate_ / s->hz_);

  // multiply by cos and sin of reference tone.
  // basically mixing down to s->hz_ and
  // converting to I/Q.
  // the goal is to sum these up over proposed
  // symbol times to see if the sum vector is long.
  // this is essentially extracting one FFT bucket,
  // from an FFT of length block_.
  std::vector<double> x;
  std::vector<double> y;
  for(ulong i = 0; i < samples.size(); i++){
    double xx = cos(tone_phase) * samples[i];
    double yy = sin(tone_phase) * samples[i];
    x.push_back(xx);
    y.push_back(yy);
    tone_phase += tone_inc;
  }

  // maintain buckets for both coarse and fine so that
  // switching doesn't have a discontinuity.

  for(int ii = 0; ii < sync_buckets_coarse; ii++){
    int i0 = off - block_ + ii * (block_ / sync_buckets_coarse);

    double sum = 0;
    for(int cnt = 0; cnt < 2; i0 += block_, cnt++){
      double xsum = 0;
      double ysum = 0;

      // look at each of this symbol's samples,
      // and sum up vectors.

#if 1
      int j0 = i0 + block_/2 - taps_.size()/2;
      for(int j = 0; j < taps_.size(); j++){
        if(j0+j >= 0 && j0+j < x.size()){
          xsum += x[j0+j] * taps_[j];
          ysum += y[j0+j] * taps_[j];
        }
      }
#else
      for(int j = 0; j < block_; j++){
        xsum += x[i0+j] * shape_[j];
        ysum += y[i0+j] * shape_[j];
      }
#endif

      double amp = sqrt(xsum*xsum + ysum*ysum);
      sum += amp;
    }

    s->coarse_buckets_[ii] = (1.0 - sync_weight_coarse) *
      s->coarse_buckets_[ii] + sync_weight_coarse * sum;
  }

  for(int ii = 0; ii < sync_buckets_fine; ii++){
    int i0 = off - block_ + ii * (block_ / sync_buckets_fine);

    double sum = 0;
    for(int cnt = 0; cnt < 2; i0 += block_, cnt++){
      double xsum = 0;
      double ysum = 0;

      // look at each of this symbol's samples,
      // and sum up vectors.

      for(int j = 0; j < block_; j++){
        xsum += x[i0+j] * shape_[j];
        ysum += y[i0+j] * shape_[j];
      }

      double amp = sqrt(xsum*xsum + ysum*ysum);
      sum += amp;
    }

    s->fine_buckets_[ii] = (1.0 - sync_weight_fine) *
      s->fine_buckets_[ii] + sync_weight_fine * sum;
  }

  if(s->sync_q_ < sync_thresh){
    int mxi = -1;
    for(int i = 0; i < sync_buckets_coarse; i++){
      if(mxi < 0 || s->coarse_buckets_[i] > s->coarse_buckets_[mxi]){
        mxi = i;
      }
    }
    s->sync_ = mxi * (block_ / sync_buckets_coarse);
  } else {
    int mxi = -1;
    for(int i = 0; i < sync_buckets_fine; i++){
      if(mxi < 0 || s->fine_buckets_[i] > s->fine_buckets_[mxi]){
        mxi = i;
      }
    }
    s->sync_ = mxi * (block_ / sync_buckets_fine);
  }

  if(s->sync_ > (int)block_/2)
    s->sync_ -= block_;
  if(s->sync_ < -(int)(block_/2))
    s->sync_ += block_;

  assert(s->sync_ >= -(int)block_/2 && s->sync_ <= (int)block_/2);
}

//
// mix a signal at hz down to 0 hz, with I/Q to preserve negative
// frequencies.
//
std::vector<std::complex<double>>
Demod::mix(const std::vector<double> a, double hz, int rate)
{
  std::vector<std::complex<double>> out(a.size());

  double phase = 0;
  for(int i = 0; i < a.size(); i++){
    double x = cos(phase) * a[i];
    double y = sin(phase) * a[i];
    out[i] = std::complex<double>(x, y);
    phase += 2 * M_PI / (rate / hz);
  }

  return out;
}

//
// demodulate the symbol that starts at samples[off],
// at the indicated hz. compares with phase of the
// previous symbol.
// q is quality, 0..1, reflects phase difference.
//
void
Demod::demod_bit(const std::vector<double> &samples, uint off,
                 double hz, int &bit, double &phase_diff, double &q)
{
  assert(off >= 0);
  assert(off >= block_);
  assert(off - block_ >= 0);
  assert(off + block_ <= samples.size());

  // correlate against both a sin wave and a cos
  // wave, together they will indicate what the phase is.
  // then atan2 to find angle. we're really mixing down
  // to an I/Q baseband, which yields a sum and difference,
  // and we need to average away the sum signal.

  // calculate phase of previous symbol and of this symbol.
  // relative to reference tone.

  std::vector<std::complex<double>> mixed = mix(samples, hz, rate_);

  std::vector<double> theta;
  
  for(int si = 0; si < 2; si++){
    double xsum = 0, ysum = 0;
    int i0 = off - block_ + si*block_ + block_/2 - taps_.size()/2;
    for(int i = 0; i < taps_.size(); i++){
      if(i0+i >= 0 && i0+i < mixed.size()){
        double xx = mixed[i0+i].real();
        double yy = mixed[i0+i].imag();
        xsum += xx * taps_[i];
        ysum += yy * taps_[i];
      }
    }
    theta.push_back(atan2(ysum, xsum));
  }

  // are they the same-ish phase, or different?
  // the range is -pi .. pi
  double d = theta[1] - theta[0];
  if(d < 0){
    d = 0 - d;
  }
  if(d > M_PI){
    d = (2 * M_PI) - d;
  }

  if(d >= M_PI/2){
    // different phase
    phase_diff = d;
    bit = 0;
    q = (d - M_PI / 2) / (M_PI / 2);
  } else {
    // same phase
    phase_diff = d;
    bit = 1;
    q = 1.0 - d / (M_PI / 2);
  }
}

//
// use phase change from symbol to symbol to adjust s->hz_.
//
// XXX similar to demod_bit()...
//
bool
Demod::adjust_hz_phase(const std::vector<double> &samples, Signal *s)
{
  if(s->phase_q_ < hz_phase_thresh)
    return false;
  
  std::vector<std::complex<double>> mixed = mix(samples, s->hz_, rate_);

  std::vector<double> theta;
  
#if 1
  int off = 2*block_ + s->sync_;
  for(int si = 0; si < 2; si++){
    double xsum = 0, ysum = 0;
    int i0 = off - block_ + si*block_ + block_/2 - taps_.size()/2;
    for(int i = 0; i < taps_.size(); i++){
      if(i0+i >= 0 && i0+i < mixed.size()){
        double xx = mixed[i0+i].real();
        double yy = mixed[i0+i].imag();
        xsum += xx * taps_[i];
        ysum += yy * taps_[i];
      }
    }
    theta.push_back(atan2(ysum, xsum));
  }
#else
  for(int si = 0; si < 2; si++){
    double xsum = 0, ysum = 0;
    int i0 = 2*block_ + s->sync_ - block_ + si*block_;
    for(int i = 0; i < block_; i++){
      xsum += mixed[i0+i].real();
      ysum += mixed[i0+i].imag();
    }
    theta.push_back(atan2(ysum, xsum));
  }
#endif

  // atan2() yields -pi .. +pi
  // adjust phase difference range to 0 .. pi
  double d = theta[1] - theta[0];
  if(d < 0){
    d = 0 - d;
  }
  if(d > M_PI){
    d = (2 * M_PI) - d;
  }

  double th0 = theta[0];
  double th1 = theta[1];
  if(d >= M_PI / 2){
    // different phase, probably a "0".
    // so simulate same phase.
    th1 += M_PI;
  }
  while(th1 > th0 + M_PI)
    th1 -= 2 * M_PI;
  while(th1 < th0 - M_PI)
    th1 += 2 * M_PI;

  // convert from radians per symbol to cycles per second.
  double dhz = ((th1 - th0) * 32) / (2 * M_PI);

  if(fabs(dhz) < hz_phase_max){
    s->hz_ = (1 - hz_phase_weight) * s->hz_ + (s->hz_ - dhz) * hz_phase_weight;
    return true;
  } else {
    return false;
  }
}

//
// AFC
//
double
Demod::adjust_hz_fine(const std::vector<double> &samples, Signal *s)
{
  int hz_n;
  double candidate[13]; // try a few near s->hz_

  double inc = hz_fine;
  hz_n = fine_hz_n;
  candidate[0] = s->hz_;
#if 0
  candidate[1] = s->hz_ - inc;
  candidate[2] = s->hz_ + inc;
  candidate[3] = s->hz_ - 2*inc;
  candidate[4] = s->hz_ + 2*inc;
  candidate[5] = s->hz_ - 4*inc;
  candidate[6] = s->hz_ + 4*inc;
  candidate[7] = s->hz_ - 8*inc;
  candidate[8] = s->hz_ + 8*inc;
  candidate[9] = s->hz_ - 16*inc;
  candidate[10] = s->hz_ + 16*inc;
  candidate[11] = s->hz_ - 32*inc;
  candidate[12] = s->hz_ + 32*inc;
#else
  candidate[1] = s->hz_ - inc;
  candidate[2] = s->hz_ + inc;
  candidate[3] = s->hz_ - 2*inc;
  candidate[4] = s->hz_ + 2*inc;
  candidate[5] = s->hz_ - 3*inc;
  candidate[6] = s->hz_ + 3*inc;
  candidate[7] = s->hz_ - 4*inc;
  candidate[8] = s->hz_ + 4*inc;
  candidate[9] = s->hz_ - 5*inc;
  candidate[10] = s->hz_ + 5*inc;
  candidate[11] = s->hz_ - 6*inc;
  candidate[12] = s->hz_ + 6*inc;
#endif
  
  double vv[13];
  assert(hz_n <= sizeof(vv) / sizeof(vv[0]));

  for(int i = 0; i < hz_n; i++){
    double dummy_q;
    int dummy_bit;
    double ph;
    demod_bit(samples, 2*block_ + s->sync_, candidate[i],
              dummy_bit, ph, dummy_q);
    vv[i] = ph;
  }

  double weight = hz_weight_fine;
  
  int pi = -1;
  for(int i = 0; i < hz_n; i++){
    if(pi < 0 || vv[i] > vv[pi]){
      pi = i;
    }
  }
  double fhz = (1 - weight) * s->hz_ + candidate[pi] * weight;

  double max_excursion = 31.25 / (2 * coarse_blocks);

  if(fhz > s->hz0_ + max_excursion)
    fhz = s->hz0_ + max_excursion;
  if(fhz < s->hz0_ - max_excursion)
    fhz = s->hz0_ - max_excursion;

  return fhz;
}

double
Demod::adjust_hz_coarse(const std::vector<double> &samples, Signal *s)
{
  double bin_hz = 31.25 / coarse_blocks;
  double inc = bin_hz / coarse_hz_n;

  double hzv[(int)coarse_hz_n];
  for(int i = 0; i < coarse_hz_n; i++){
    hzv[i] = s->hz0_ - bin_hz/2 + inc/2 + i*inc;
    double dummy_q;
    int dummy_bit;
    double ph;
    demod_bit(samples, 2*block_ + s->sync_, hzv[i],
              dummy_bit, ph, dummy_q);
    s->hz_buckets_[i] = (1 - hz_weight_coarse) * s->hz_buckets_[i] +
      hz_weight_coarse * ph;
  }

  int mxi = -1;
  for(int i = 0; i < coarse_hz_n; i++){
    if(mxi < 0 || s->hz_buckets_[i] > s->hz_buckets_[mxi]){
      mxi = i;
    }
  }

  return hzv[mxi];
}

void
Demod::adjust_hz(const std::vector<double> &samples, Signal *s)
{
  if(s->hz_q_ > hz_fine_thresh){
    s->hz_ = adjust_hz_fine(samples, s);
  } else {
    s->hz_ = adjust_hz_coarse(samples, s);
  }
}

//
// process a new symbol's-worth of samples for
// a single signal that we've been decoding.
// the new samples are in samples_[2*block_ .. 3*block_],
// plus sync_.
//
void
Demod::tick(Signal *s)
{
  assert(samples_.size() >= 5*block_);

  if(s->started_ == 0){
    s->started_ = 1;
    adjust_sync(samples_, s, 1*block_);
    adjust_sync(samples_, s, 2*block_);
    adjust_sync(samples_, s, 3*block_);
  }

  int bit;
  double dummy_diff;
  double q;
  demod_bit(samples_, 2*block_ + s->sync_, s->hz_, bit, dummy_diff, q);

  s->hz_q_ = (1 - hz_q_weight) * s->hz_q_ + hz_q_weight * q;
  s->phase_q_ = (1 - phase_q_weight) * s->phase_q_ + phase_q_weight * q;
  s->sync_q_ = (1 - sync_q_weight) * s->sync_q_ + sync_q_weight * q;
  s->cull_q_ = (1 - cull_q_weight) * s->cull_q_ + cull_q_weight * q;
  s->good_q_ = (1 - good_q_weight) * s->good_q_ + good_q_weight * q;
  
  got_bit(s, bit);

  if(adjust_hz_phase(samples_, s) == false){
    if(bit == 0){
      adjust_hz(samples_, s);
    }
  }
  if(bit == 0){
    // and adjust symbol boundary -- s->sync_.
    // estimates how much symbol overlaps with a symbol
    // of a different phase. so phases must be different.
    adjust_sync(samples_, s, 2*block_);
  }
}

struct varithing varicode[] = {
  { "00100", " " },
  { "00101111111100", ""}, // back space XXX
  { "001111100", "\n" },
  { "001110100", ""}, // line feed
  { "0011111111100", "!" },
  { "0010101111100", "\"" },
  { "0011111010100", "#" },
  { "0011101101100", "$" },
  { "00101101010100", "%" },
  { "00101011101100", "&" },
  { "0010111111100", "'" },
  { "001111101100", "(" },
  { "001111011100", ")" },
  { "0010110111100", "*" },
  { "0011101111100", "+" },
  { "00111010100", "," },
  { "0011010100", "-" },
  { "00101011100", "." },
  { "0011010111100", "/" },
  { "001011011100", "0" },
  { "001011110100", "1" },
  { "001110110100", "2" },
  { "001111111100", "3" },
  { "0010111011100", "4" },
  { "0010101101100", "5" },
  { "0010110101100", "6" },
  { "0011010110100", "7" },
  { "0011010101100", "8" },
  { "0011011011100", "9" },
  { "001111010100", ":" },
  { "0011011110100", ";" },
  { "0011110110100", "<" },
  { "00101010100", "=" },
  { "0011101011100", ">" },
  { "00101010111100", "?" },
  { "00101011110100", "@" },
  { "00111110100", "A" },
  { "001110101100", "B" },
  { "001010110100", "C" },
  { "001011010100", "D" },
  { "00111011100", "E" },
  { "001101101100", "F" },
  { "001111110100", "G" },
  { "0010101010100", "H" },
  { "00111111100", "I" },
  { "0011111110100", "J" },
  { "0010111110100", "K" },
  { "001101011100", "L" },
  { "001011101100", "M" },
  { "001101110100", "N" },
  { "001010101100", "O" },
  { "001101010100", "P" },
  { "0011101110100", "Q" },
  { "001010111100", "R" },
  { "00110111100", "S" },
  { "00110110100", "T" },
  { "0010101011100", "U" },
  { "0011011010100", "V" },
  { "0010101110100", "W" },
  { "0010111010100", "X" },
  { "0010111101100", "Y" },
  { "00101010110100", "Z" },
  { "0011111011100", "[" },
  { "0011110111100", "\\" },
  { "0011111101100", "]" },
  { "00101011111100", "^" },
  { "0010110110100", "_" },
  { "00101101111100", "`" },
  { "00101100", "a" },
  { "00101111100", "b" },
  { "0010111100", "c" },
  { "0010110100", "d" },
  { "001100", "e" },
  { "0011110100", "f" },
  { "00101101100", "g" },
  { "0010101100", "h" },
  { "00110100", "i" },
  { "0011110101100", "j" },
  { "001011111100", "k" },
  { "001101100", "l" },
  { "0011101100", "m" },
  { "00111100", "n" },
  { "0011100", "o" },
  { "0011111100", "p" },
  { "0011011111100", "q" },
  { "001010100", "r" },
  { "001011100", "s" },
  { "0010100", "t" },
  { "0011011100", "u" },
  { "00111101100", "v" },
  { "00110101100", "w" },
  { "001101111100", "x" },
  { "00101110100", "y" },
  { "0011101010100", "z" },
  { "00101011011100", "{" },
  { "0011011101100", "|" },
  { "00101011010100", "}" },
  { "00101101011100", "~" },
  { 0, 0 },
};

//
// drop leading bits from bits_.
// ok < 0 means they were junk.
// ok > 0 means they were good bits.
// ok == 0 means unknown.
//
// caller must hold mu_
//
void
Demod::drop(Signal *s, unsigned int n, int ok)
{
  assert(n <= s->bits_.size());
  s->bits_.erase(s->bits_.begin(), s->bits_.begin() + n);

  if(ok < 0){
    emit(s, "~"); // keep emitting junk so it's clearly low quality
  }
}

void
Demod::got_bit(Signal *s, int bit)
{
  mu_.lock();
  
  s->bits_.push_back(bit);

  // if we see two zeros, then a legal varicode character,
  // consume them.

  // if we other than two zeros, skip until we see two zeros.

  // if we see two zeros, but then not a valid varicode
  // character, and then two zeros, skip to the second two zeros.

  while(s->bits_.size() >= 3){
    if(s->bits_[0] != 0){
      // 1xx
      // keep looking for 001
      drop(s, 1, -1);
    } else if(s->bits_[1] != 0){
      // 01x
      // keep looking for 001
      drop(s, 2, -1);
    } else if(s->bits_[2] == 0){
      // 000 -- idle.
      drop(s, 1, 0);
    } else {
      // we have 00 and at least one more bit.
      // do the bits after 00 form a complete character, including trailing 00?
      // do the bits after 00 start a valid character?
      int prefix = -1;
      int match = -1;
      for(ulong i = 0; varicode[i].c; i++){
        for(uint j = 0; ; j++){
          if(varicode[i].b[j] == '\0'){
            match = i;
            break;
          }
          if(j >= s->bits_.size()){
            prefix = i;
            break;
          }
          if(varicode[i].b[j] == '0'){
            if(s->bits_[j] != 0){
              break;
            }
          } else {
            if(s->bits_[j] != 1){
              break;
            }
          }
        }
      }

      if(match >= 0){
        const char *txt = varicode[match].c;
        emit(s, txt);
        // delete leading 00, and varicode character, but not trailing 00.
        drop(s, strlen(varicode[match].b) - 2, 1);
      } else if(prefix >= 0){
        // don't consume anything, since bits_
        // holds a prefix for some character.
        break;
      } else {
        // no match, no prefix.
        if(0 && fabs(s->hz_ - 1200) < 10){
          printf(" <");
          for(int i = 0; i < s->bits_.size(); i++){
            printf("%c", s->bits_[i] ? '1' : '0');
          }
          printf("> ");
        }
        drop(s, 1, -1);
      }
    }
  }

  mu_.unlock();
}

//
// caller must hold mu_
//
void
Demod::emit(Signal *s, const char *txt)
{
  s->text_ += txt;
  
  double max_excursion = 31.25 / (2 * coarse_blocks);
  for(int i = 0; i < watch_.size(); i++){
    double whz = watch_[i].hz;
    double diff = fabs(s->hz_ - whz);
    if(diff <= max_excursion){
      // is this the closest Signal to the watch hz?
      int closest = 1;
      for(ulong i = 0; i < signals_.size(); i++){
        int close = fabs(signals_[i]->hz_ - whz) < max_excursion;
        int better = signals_[i]->hz_q_ > s->hz_q_;
        if(close && better){
          closest = 0;
          break;
        }
      }
      if(closest){
        watch_[i].cb(watch_[i].arg, whz, s->hz_, txt);
      }
    }
  }
}

//
// return a copy of signals_ for ui.cc to display.
//
std::vector<Signal>
Demod::signals()
{
  mu_.lock();
  std::vector<Signal> v;
  for(int i = 0; i < signals_.size(); i++){
    // the locking is not great, so copy selected fields,
    // rather than all.
    Signal *s = signals_[i];
    Signal ss(s->hz_);
    ss.hz_ = s->hz_;
    ss.hz0_ = s->hz0_;
    ss.sync_ = s->sync_;
    ss.text_ = s->text_;
    ss.start_ = s->start_;
    ss.hz_q_ = s->hz_q_;
    ss.phase_q_ = s->phase_q_;
    ss.sync_q_ = s->sync_q_;
    ss.cull_q_ = s->cull_q_;
    ss.good_q_ = s->good_q_;
    ss.last_good_ = s->last_good_;
    v.push_back(ss);
  }
  mu_.unlock();
  return v;
}

//
// for benchmark() in ui.cc.
//

double dummy_var = 0.0;

struct opt_var vars [] =
  {
   { "hz_weight_coarse", &hz_weight_coarse, { 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.1, 0.2, 0.3, 0.4, -9999 } },
   { "hz_weight_fine", &hz_weight_fine, { 0.1, 0.3, 0.5, 0.7, 0.9, 1.0, -9999 } },
   { "hz_fine_thresh", &hz_fine_thresh, { 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75,  -9999 } },
   { "hz_q_weight", &hz_q_weight, { 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.02, 0.05, 0.1, -9999 } },
   { "hz_phase_weight", &hz_phase_weight, { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, -9999 } },
   { "sync_weight_coarse", &sync_weight_coarse, { 0.05, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.3, -9999 } },
   { "coarse_hz_n", &coarse_hz_n, { 3, 4, 5, 6, 7, 8, 10, 12, -9999 } },
   { "fine_hz_n", &fine_hz_n, { 3, 5, 7, 9, -9999 } },
   { "hz_fine", &hz_fine, { 0.01, 0.015, 0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.4, 0.5, -9999 } },
   { "hz_phase_thresh", &hz_phase_thresh, { 0.6, 0.7, 0.8, 0.9, 0.95, -9999 } },
   { "hz_phase_max", &hz_phase_max, { 0.05, 0.1, 0.20, 0.30, 0.40, 0.5, 1, -9999 } },
   { "coarse_blocks", &coarse_blocks, { 1, 2, 3, 4, -9999 } },
   { "coarse_weight", &coarse_weight, { 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, -9999 } },
   { "sync_weight_fine", &sync_weight_fine, { 0.025, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, -9999 } },
   { "sync_thresh", &sync_thresh, { 0.55, 0.575, 0.6, 0.625, 0.65, 0.7, 0.75, 0.8, -9999 } },
   { "sync_buckets_coarse", &sync_buckets_coarse, { 4, 5, 6, 7, 8, 11, 15, -9999 } },
   { "sync_buckets_fine", &sync_buckets_fine, { 9, 11, 12, 13, 15, 20, 22, 23, 24, -9999 } },
   { "sync_q_weight", &sync_q_weight, { 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, -9999 } },
   { "phase_q_weight", &phase_q_weight, { 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, -9999 } },
   { "coarse_top_n", &coarse_top_n, { 10, 12, 14, 16, 18, 20, 25, 30, -9999 } },
   { "good_q_weight", &good_q_weight, { 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, -9999 } },
   { "good_q", &good_q, { 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, -9999 } },
   { "keep_time", &keep_time, { 2, 3, 4, 5, 6, 8, 10, 12, -9999 } },
   { "infant_time", &infant_time, { 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 1.0, 2, -9999 } },
   { "cull_q", &cull_q, { 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, -9999 } },
   { "cull_q_weight", &cull_q_weight, { 0.1, 0.2, 0.3, 0.5, 0.65, 0.8, -9999 } },
   { "defaults", &dummy_var, { 0, 0, 0, -9999 } },
   { 0, 0, { } },
};
