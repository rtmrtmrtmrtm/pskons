#ifndef DEMOD_H
#define DEMOD_H

#include <vector>
#include <string>
#include <mutex>

//
// an ongoing transmission that we're tracking.
//
class Signal {
 public:
  double hz_;   // last known hz
  double hz0_;  // original coarse hz
  int sync_;    // symbols start this many samples before/after nominal.
  std::vector<int> bits_;
  std::string text_; // all text
  double start_; // UNIX time when this Signal first created.
  bool pin_;

  // quality based on phase differences.
  double hz_q_;
  double phase_q_;
  double sync_q_;
  double cull_q_;
  double good_q_;
  double last_good_;

  double coarse_buckets_[32];
  double fine_buckets_[32];

  double hz_buckets_[32];

  int started_;
  
  Signal(double hz) {
    hz_ = hz;
    hz0_ = hz;
    sync_ = 0;
    hz_q_ = 0.5;
    phase_q_ = 0.5;
    sync_q_ = 0.5;
    cull_q_ = 0.5;
    good_q_ = 0.5;
    pin_ = false;
    last_good_ = 0;
    for(int i = 0; i < 32; i++){
      coarse_buckets_[i] = 0;
      fine_buckets_[i] = 0;
      hz_buckets_[i] = 0;
    }
    started_ = 0;
  }
};

typedef void (*watch_cb_t)(void*, double, double, const char *);
struct Watch {
  double hz;
  watch_cb_t cb;
  void *arg;
};

class Demod {
 private:

  // protects signals_ from ui thread's call to signals().
  std::mutex mu_;
  
  int rate_;
  int block_;
  double clock_; // seconds, driven by sample arrivals

  // recent average signal strength, coarse, 31.25/2 per bin.
  std::vector<double> recent_;

  // samples_ is a window of five symbols (block) of samples.
  std::vector<double> samples_;

  std::vector<Signal*> signals_;

  // integral of signals_.size()
  long long sigx_;
  int sign_;

  // raised_cosine(block_)
  std::vector<double> shape_;

  // big low-pass filter
  std::vector<double> taps_;

  std::vector<Watch> watch_;
    
 public:

  Demod(int rate);
  ~Demod();

  void got(const std::vector<double> &);
  void coarse(const std::vector<double> &samples, int i0);
  void tick_all();
  void tick(Signal *);
  void one_adjust_sync(Signal *, int);
  void adjust_sync(const std::vector<double> &, Signal *, int);
  void adjust_hz(const std::vector<double> &, Signal *);
  double adjust_hz_coarse(const std::vector<double> &, Signal *);
  double adjust_hz_fine(const std::vector<double> &, Signal *);
  bool adjust_hz_phase(const std::vector<double> &, Signal *s);
  void demod_bit(const std::vector<double> &, unsigned int off, double hz, int &bit, double &phase_diff, double &q);
  void got_bit(Signal *, int);
  void drop(Signal *, unsigned int, int);
  void emit(Signal *, const char *);
  void cull();
  void init_signal(Signal *);
  double demod3(const std::vector<std::complex<double>> &, uint off,
                double hz, int bit1, int bit2);
  std::vector<std::complex<double>> mix(const std::vector<double> a, double hz, int rate);

  std::vector<Signal> signals();

  void watch(double, watch_cb_t, void*);
  void clear_watches();

  int rate() { return rate_; }

  // average number of signals
  double nsig() { return sigx_ / (double) sign_; }
};

struct opt_var {
  const char *name;
  double *v;
  double vals[20];
};
extern struct opt_var vars[];

struct varithing {
  const char *b;
  const char *c;
};
extern struct varithing varicode[];

#endif
