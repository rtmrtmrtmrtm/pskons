#include "mod.h"
#include "util.h"
#include <math.h>
#include <assert.h>


//
// streaming filter: you give it input samples a
// bit at a time, it hands you output a bit at a time.
//
class Filter {
private:

  std::vector<double> taps_;
  int n_; // taps_.size()

  std::vector<double> buf_;

public:

  Filter(const std::vector<double> & taps) {
    taps_ = taps;
    n_ = taps_.size();
  }

  bool go1(double in, double &out) {
    buf_.push_back(in);
    if(buf_.size() >= n_){
      assert(buf_.size() == n_);
      double x = 0;
      for(int i = 0; i < n_; i++){
        x += buf_[i] * taps_[n_ - i - 1];
      }
      buf_.erase(buf_.begin());
      out = x;
      return true;
    } else {
      // no output yet.
      return false;
    }
  }

  std::vector<double> go(const std::vector<double> &samples) {
    std::vector<double> out;
    for(int i = 0; i < samples.size(); i++){
      double x;
      bool valid = go1(samples[i], x);
      if(valid)
        out.push_back(x);
    }
    return out;
  }
};

//
// includes raised-cosine low-pass transmit filter.
//
std::vector<double>
bits2psk(const std::vector<int> &bits, double hz, int rate, int starting, int ending)
{
  static double last_phase = 1;

  //
  // for each sample time, the phase, 1 or -1.
  //
  std::vector<double> phases;

  int bitno = 0;

  if(starting){
    //
    // start with zero-amplitude signal to start w/o click.
    //
    while(phases.size() < (bitno+1) * (rate / 31.25)){
      phases.push_back(0);
    }
    bitno++;

    //
    // now a dummy symbol as reference for phase changes.
    //
    while(phases.size() < (bitno+1) * (rate / 31.25)){
      phases.push_back(last_phase);
    }
    bitno++;
  }
  
  for(int i = 0; i < bits.size(); i++){
    // 1-bit leave phase unchanged.
    // 0-bit swaps the phases.
    double this_phase = bits[i] ? last_phase : -last_phase;

    while(phases.size() < (bitno+1) * (rate / 31.25)){
      phases.push_back(this_phase);
    }
    bitno++;

    last_phase = this_phase;
  }

  if(ending){
    //
    // and fade to nothing at end to avoid a click.
    //
    while(phases.size() < (bitno+1) * (rate / 31.25)){
      phases.push_back(0);
    }
    bitno++;
  }

  //
  // transmit shaping filter on the phases.
  //
  static Filter *ff = 0;
  if(ff == 0){
    ff = new Filter(raised_cosine(round(rate / 31.25)));
  }
  std::vector<double> shaped_phases = ff->go(phases);

  std::vector<double> samples;
  static double theta = 0.0;
  for(int i = 0; i < shaped_phases.size(); i++){
    double x = shaped_phases[i] * cos(theta);
    samples.push_back(x);
    theta += 2 * M_PI / (rate / hz);
    if(theta >= 2 * M_PI)
      theta -= 2 * M_PI;
  }

  return samples;
}
