#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <stdio.h>
#include <mutex>
#include <map>
#include <list>
#include <thread>
#include <algorithm>
#include <liquid/liquid.h>
#include <complex>
#include "util.h"
#include "demod.h"
#include "snd.h"
#include <string>
#include <math.h>
#include <assert.h>
#include <sys/ioctl.h>
#include <termios.h>
#include "sndfile.h"
#include "bench.h"
#include "mod.h"

void
usage()
{
  fprintf(stderr, "Usage: pskons [-only hz] -file files...\n");
  fprintf(stderr, "       pskons -opt\n");
  fprintf(stderr, "       pskons -bench\n");
  fprintf(stderr, "       pskons -card X [-out Y]\n");
  fprintf(stderr, "       pskons -cardfile xxx.wav\n");
  fprintf(stderr, "       pskons -levels X\n");
  fprintf(stderr, "       pskons -gen file hz text\n");
  snd_list();
  exit(1);
}

volatile int nsignals;

class Line {
public:
  double hz_;
  std::string text_;
  double last_; // UNIX time we last saw q_ > 0.75
  double q_; // from demod.cc, based on phase differences
};

std::vector<Line> lines;
std::mutex lines_mu;

struct layout {
  int rows;
  int cols;
  int pan_rows; // all signals
  int rx_rows;  // the one received signal
  int tx_rows;  // my transmission
};

void
get_layout(layout &lay)
{
  int rows = 24;
  int cols = 80;

  struct winsize ws;
  if(ioctl(1, TIOCGWINSZ, &ws) == 0){
    rows = ws.ws_row;
    cols = ws.ws_col;
  }

  lay.rows = rows;
  lay.cols = cols;
  lay.tx_rows = rows / 4;
  lay.rx_rows = rows / 4;
  lay.pan_rows = rows - lay.tx_rows - lay.rx_rows - 1;
}

//
// squeeze out repeated ~
//
std::string
simplify(std::string s)
{
  // squeeze out repeated ~
  std::string ss;
  for(int i = 0; i < s.size(); i++){
    if(s[i] == '~' && i > 0 && s[i-1] == '~'){
      // nothing
    } else {
      ss.push_back(s[i]);
    }
  }

  return ss;
}

//
// prints exactly n lines, the last n lines of s.
// does not print the final newline!
//
void
print_n(const char *prefix, std::string s, int n, int cols)
{
  std::vector<std::string> rxlines;
  std::string line;
  for(int i = 0; i < s.size(); i++){
    if(s[i] == '\n' || s[i] == '\r' || line.size() > cols - 2){
      if(line.size() > 0){
        rxlines.push_back(line);
        line = "";
      }
    } else {
      line.push_back(s[i]);
    }
  }
  int i0 = rxlines.size() - n + 1;
  if(line.size() == 0)
    i0 -= 1;
  if(i0 < 0)
    i0 = 0;

  for(int i = 0; i < n; i++){
    printf("%s", prefix);
    if(i0 < rxlines.size()){
      printf("%s", rxlines[i0].c_str());
    } else if(i0 == rxlines.size() && line.size() > 0){
      printf("%s", line.c_str());
    }
    i0++;
    if(i+1 < n)
      printf("\n");
  }
}

//
// the signal we're talking to.
//
double rx_hz = 1500;
std::mutex rx_buf_mu;
FILE *qso_fp = 0;
std::string rx_buf;

// buffer of characters I want to transmit.
std::mutex tx_buf_mu;
std::string tx_buf;

volatile int transmitting = 0;

void
draw_screen()
{
  layout lay;

  get_layout(lay);
  
  // clear screen
  printf("\033[H");  // home
  printf("\033[2J"); // clear

  printf("%4d ", nsignals);

  if(transmitting){
    printf("TX ");
  } else {
    printf("RX ");
  }
  if(rx_hz > 0){
    printf("%.0f", rx_hz);
  }
  printf("\n");
  
  for(int i = 0; i < lay.pan_rows; i++){
    lines_mu.lock();
    if(i >= lines.size())
      break;
    Line ll = lines[i];
    lines_mu.unlock();
    
    printf("%c %4.0f %.1f ",
           'A' + i,
           ll.hz_,
           ll.q_);
    
    std::string lt = simplify(ll.text_);
    
    int j0 = lt.size() - lay.cols + 16;
    if(j0 < 0)
      j0 = 0;
    for(int j = j0; j < lt.size(); j++){
      int c = lt[j];
      if(c == '\n' || c == '\r'){
        printf(" ");
      } else {
        printf("%c", c);
      }
    }
    printf("\n");
  }

  // the selected signal.
  rx_buf_mu.lock();
  std::string s = rx_buf;
  rx_buf_mu.unlock();
  s = simplify(s);
  print_n("- ", s, lay.rx_rows, lay.cols-1);
  printf("\n");

  // my transmitted text, if any.
  tx_buf_mu.lock();
  std::string tmp = tx_buf;
  tx_buf_mu.unlock();
  print_n("> ", tmp, lay.tx_rows, lay.cols-1);
  
  fflush(stdout);
}

void
kb_loop(Demod *dm)
{
  int state = -1;
  while(1){
    unsigned char c;
    if(read(0, &c, 1) != 1){
      exit(1);
    }
    if(state < 0 && c == '\001'){
      // control-A
      state = '\001';
    } else if(c == '\030'){
      // control-X
      exit(0);
    } else if(state == '\001'){
      // control-A <letter>
      // chooses a signal to receive.

      int i = c - 'a';
      lines_mu.lock();
      if(i >= 0 && i < lines.size()){
        rx_hz = lines[i].hz_;
        dm->clear_watches();
        rx_buf_mu.lock();
        rx_buf.clear();
        rx_buf = lines[i].text_;
        if(qso_fp){
          time_t clock;
          time(&clock);
          struct tm tmx;
          gmtime_r(&clock, &tmx);
          fprintf(qso_fp, "\n\n----- %02d/%02d/%04d %02d:%02d\n\n%s",
                  tmx.tm_mon+1,
                  tmx.tm_mday,
                  tmx.tm_year+1900,
                  tmx.tm_hour,
                  tmx.tm_min,
                  rx_buf.c_str());
        }
        rx_buf_mu.unlock();
        dm->watch(rx_hz,
                  [](void *, double, double, const char *txt) {
                    rx_buf_mu.lock();
                    rx_buf += txt;
                    if(qso_fp)
                      fprintf(qso_fp, "%s", txt);
                    rx_buf_mu.unlock();
                  },
                  0);
      }
      lines_mu.unlock();

      state = -1;
    } else if(c == '\010' || c == '\177'){
      tx_buf_mu.lock();
      if(tx_buf.size() > 0 && transmitting){
        tx_buf.pop_back();
      }
      tx_buf_mu.unlock();
    } else {
      tx_buf_mu.lock();
      tx_buf.push_back((char)c);
      tx_buf_mu.unlock();
    }
  }
}

void
tx_loop(SoundOut *sout)
{
  int tx_i = 0; // how far we've gotten in tx_buf[].
  double session_start; // of this continuous transmission
  double session_samples;

  while(1){
    unsigned int c;
    bool c_valid = false;
    
    tx_buf_mu.lock();
    if(tx_buf.size() > tx_i){
      c = tx_buf[tx_i];
      c_valid = true;
      tx_i++;
    }
    tx_buf_mu.unlock();

    if(c_valid || transmitting){
      std::vector<int> bits;
      int starting = 0;
      int ending = 0;

      if(transmitting == 0){
        transmitting = 1;
        starting = 1;
        session_start = now();
        session_samples = 0;
        // preamble: transmission starts with zeros (reversals).
        for(int i = 0; i < 32; i++){
          bits.push_back(0);
        }
      }

      if(c_valid){
        const char *bbb = 0;
        for(int j = 0; varicode[j].c; j++){
          if(varicode[j].c[0] == c){
            bbb = varicode[j].b;
            break;
          }
        }
        if(bbb){
          // each varicode[j].b starts and ends with 00;
          // skip the starting 00.
          for(int i = 2; bbb[i]; i++){
            bits.push_back(bbb[i] == '0' ? 0 : 1);
          }
        } else {
          printf(" <oops %c> ", c); fflush(stdout);
        }
      }
            
      if(bits.size() == 0){
        // user stopped typing.
        // postamble: send trailer of 32 1s, then stop sending.
        transmitting = 0;
        ending = 1;
        for(int i = 0; i < 32; i++){
          bits.push_back(1);
        }
      }

      if(bits.size() > 0){
        if(sout == 0){
          fprintf(stderr, "no output card\n");
          exit(1);
        }
        int rate = sout->rate();

        std::vector<double> samples = bits2psk(bits, rx_hz, rate, starting, ending);

        std::vector<short int> shorts;
        for(int i = 0; i < samples.size(); i++){
          short int x = samples[i] * 0.44 * 32767;
          shorts.push_back(x);
        }

        // emit on sound card.
        sout->write(shorts);

        session_samples += shorts.size();
        double et = session_start + session_samples/(double)rate;

        while(now() < et - 0.2){
          usleep((et - now() - 0.2) * 1000000);
        }

      }
    } else {
      usleep(500 * 1000);
    }
  }
}

void
rx_loop(SoundIn *sin, Demod *dm)
{
  // if we're reducing 48000 to 8000, then some fancy
  // footwork in case we read a number of bytes that's
  // not divisible by six.
  std::vector<double> excess;

#if 0
  // record input audio.
  rename("record.wav", "record.wav.old");
  SF_INFO sf;
  sf.channels = 1;
  sf.samplerate = sin->rate();
  sf.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
  SNDFILE *f = sf_open("record.wav", SFM_WRITE, &sf);
  assert(f);
#endif

  while(1){
    double xxx;
    std::vector<double> v = sin->get(1024, xxx);

#if 0
    sf_write_double(f, v.data(), v.size());
    sf_write_sync(f);
#endif
    
    if(excess.size() > 0){
      v.insert(v.begin(), excess.begin(), excess.end());
      excess.clear();
    }
    
    if(v.size() == 0)
      usleep(50 * 1000);

    if(dm->rate() != sin->rate()){
      int factor = sin->rate() / dm->rate();
      int n = v.size() % factor;
      if(n != 0){
        // move last n from v to excess.
        excess.insert(excess.begin(), v.end()-n, v.end());
        v.erase(v.end()-n, v.end());
      }
      v = thin(v, factor);
    }
    
    if(transmitting){
      // suppress receive audio.
      std::vector<double> vv(v.size());
      for(int i = 0; i < vv.size(); i++)
        vv[i] = 0;
      dm->got(vv);
    } else {
      dm->got(v);
    }
  }
}
  

struct termios save_tt;

void
screen_main(SoundIn *sin, SoundOut *sout)
{
  struct termios tt;
  if(tcgetattr(0, &tt) != 0){
    fprintf(stderr, "tcgetattr failed\n");
    exit(1);
  }
  save_tt = tt;
  tt.c_lflag &= ~ICANON;
  if(tcsetattr(0, TCSANOW, &tt) != 0){
    fprintf(stderr, "tcsetattr failed\n");
    exit(1);
  }

  qso_fp = fopen("qso-trace.txt", "a");
  if(qso_fp)
    setbuf(qso_fp, 0);
  
  sin->start();

  if(sout){
    sout->start();
  }

  int demodrate;
  if((sin->rate() % 8000) == 0 && sin->rate() > 8000){
    // reduce rate to 8000
    demodrate = 8000;
  } else {
    demodrate = sin->rate();
  }
  
  Demod *dm = new Demod(demodrate);
  
  std::thread in_th( [ sin, dm ] () {
                       rx_loop(sin, dm);
                     } );

  std::thread kb_th( [ dm ] () { kb_loop(dm); } );

  std::thread tx_th( [ sout ] () { tx_loop(sout); } );

  while(1){
    layout lay;
    get_layout(lay);

    lines_mu.lock();
    if(lay.pan_rows > lines.size()){
      for(int i = lines.size(); i < lay.pan_rows; i++){
        Line x;
        x.last_ = 0;
        x.hz_ = -1;
        x.q_ = 0;
        lines.push_back(x);
      }
    }

    while(lay.pan_rows < lines.size()){
      lines.pop_back();
    }
    lines_mu.unlock();

    std::vector<Signal> sv = dm->signals();
    nsignals = sv.size();

    std::vector<bool> found(sv.size());

    // for each existing display line, look for best match in sv[].
    double tolerance = 10.4; // hz
    lines_mu.lock();
    for(int i = 0; i < lines.size(); i++){
      if(lines[i].hz_ < 0.1)
        continue;
      int best = -1;
      for(int j = 0; j < sv.size(); j++){
        if(fabs(lines[i].hz_ - sv[j].hz_) < tolerance){
          found[j] = true;
          if(best < 0 || sv[j].hz_q_ > sv[best].hz_q_){
            best = j;
          }
        }
      }
      if(best >= 0){
        lines[i].text_ = sv[best].text_;
        lines[i].q_ = sv[best].hz_q_;
        if(sv[best].hz_q_ > 0.75){ // XXX
          lines[i].last_ = now();
          lines[i].hz_ = sv[best].hz_;
        }
      }
    }
    lines_mu.unlock();

    // try to fit new signals onto the screen,
    // replacing the weakest.
    for(int i = 0; i < sv.size(); i++){
      if(sv[i].hz_ < 0.1)
        continue;
      if(found[i] == false){
        lines_mu.lock();
        int worst = -1;
        for(int j = 0; j < lines.size(); j++){
          if(sv[i].hz_q_ > lines[j].q_ && now() - lines[j].last_ > 10 &&
             (rx_hz <= 0.0 || fabs(rx_hz - lines[j].hz_) > tolerance)){
            if(worst < 0 || lines[j].q_ < lines[worst].q_){
              worst = j;
            }
          }
        }
        if(worst >= 0){
          lines[worst].hz_ = sv[i].hz_;
          lines[worst].q_ = sv[i].hz_q_;
          lines[worst].text_ = sv[i].text_;
          lines[worst].last_ = 0; // now();
          found[i] = true;
        }
        lines_mu.unlock();
      }
    }

    draw_screen();

    usleep(200*1000);
  }
}

void
genfile(const char *file, double hz, const char *text)
{
  int rate = 8000;
  std::vector<int> bits;

  // preamble: transmission starts with zeros (reversals).
  for(int i = 0; i < 32; i++){
    bits.push_back(0);
  }

  for(int i = 0; text[i]; i++){
    int c = text[i];
    const char *bbb = 0;
    for(int j = 0; varicode[j].c; j++){
      if(varicode[j].c[0] == c){
        bbb = varicode[j].b;
        break;
      }
    }
    if(bbb){
      // each varicode[j].b starts and ends with 00;
      // skip the starting 00.
      for(int i = 2; bbb[i]; i++){
        bits.push_back(bbb[i] == '0' ? 0 : 1);
      }
    } else {
      printf("oops '%c'\n", c);
    }
  }

  // postamble: trailer of 1s.
  for(int i = 0; i < 32; i++){
    bits.push_back(1);
  }

  std::vector<double> samples = bits2psk(bits, hz, rate, 1, 1);

  if(0){
    writetxt(samples, "a1.txt");
    // add some multi-path
    int spread = 0.5 * 256;
    spread += 0;
    std::vector<double> s1;
    for(int i = 0; i < samples.size(); i++){
      double x = samples[i];
      if(i >= spread)
        x += 1.0 * samples[i-spread];
      s1.push_back(x);
    }
    samples = s1;
    writetxt(samples, "a2.txt");
  }

  writewav(samples, file, rate);
  writetxt(samples, "gen.txt");
}

int
main(int argc, char *argv[])
{
  double only = -1;
  int incard = -1;
  int outcard = -1;
  const char *cardfile = 0;

  if(argc < 2)
    usage();

  int ai = 1;
  while(ai < argc){
    if(strcmp(argv[ai], "-opt") == 0 && ai+1 == argc){
      optimize();
      ai++;
    } else if(strcmp(argv[ai], "-bench") == 0 && ai+1 == argc){
      double nsig;
      benchmark(1, nsig);
      ai++;
    } else if(strcmp(argv[ai], "-only") == 0 && ai+1 < argc){
      ai++;
      only = atof(argv[ai]);
      ai++;
    } else if(strcmp(argv[ai], "-file") == 0){
      ai++;
      for(; ai < argc; ai++){
        int filerate;
        std::vector<double> v = readwav(argv[ai], filerate);
        int demodrate;
        if((filerate % 8000) == 0 && filerate > 8000){
          // reduce to 8000
          v = thin(v, filerate / 8000);
          demodrate = 8000;
        } else {
          demodrate = filerate;
        }
        Demod *dm = new Demod(demodrate);
        if(only >= 0){
          dm->watch(only,
                    [ ](void *arg, double, double, const char *txt) {
                      printf("%s", txt);
                      fflush(stdout);
                    },
                    0);
        }
        dm->got(v);
      }
    } else if(strcmp(argv[ai], "-card") == 0 && ai+1 < argc){
      ai++;
      incard = atoi(argv[ai]);
      ai++;
    } else if(strcmp(argv[ai], "-cardfile") == 0 && ai+1 < argc){
      ai++;
      cardfile = argv[ai];
      ai++;
    } else if(strcmp(argv[ai], "-out") == 0 && ai+1 < argc){
      ai++;
      outcard = atoi(argv[ai]);
      ai++;
    } else if(strcmp(argv[ai], "-levels") == 0 && ai+1 < argc){
      ai++;
      incard = atoi(argv[ai]);
      ai++;
      levels(incard);
      exit(0);
    } else if(strcmp(argv[ai], "-gen") == 0 && ai+3 < argc){
      ai++;
      genfile(argv[ai], atof(argv[ai+1]), argv[ai+2]);
      ai += 3;
    } else {
      usage();
    }
  }

  if(incard >= 0){
    CardSoundIn *sin = new CardSoundIn(incard);
    SoundOut *sout = 0;
    if(outcard >= 0){
      sout = new SoundOut(outcard);
    }
    screen_main(sin, sout);
  }

  if(cardfile){
    FileSoundIn *sin = new FileSoundIn(cardfile);
    screen_main(sin, 0);
  }
}
