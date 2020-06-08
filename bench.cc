#include "bench.h"
#include "util.h"
#include "demod.h"
#include <vector>
#include <string>
#include <mutex>
#include <thread>
#include <unistd.h>
#include <stdlib.h>
#include <sys/types.h>
#ifndef __linux__
#include <sys/sysctl.h>
#endif
#include <assert.h>
#include <string.h>

struct {
  const char *wavname;
  double hz;
  const char *wanted;
} files[] = {
#if 1
             { "g.wav", 1104, "g1104.txt" },
             { "g.wav", 1730, "g1730.txt" },
             { "g.wav", 1557, "g1557.txt" },
             { "g.wav", 1283, "g1283.txt" },
             { "g.wav", 1200, "g1200.txt" },
             { "g.wav", 1182, "g1182.txt" },
             { "a.wav", 1608, "1608.txt" },
             { "a.wav", 1165, "a1165.txt" },
             { "a.wav", 912, "a912.txt" },
             { "a.wav", 1480, "a1480.txt" },
             { "b.wav", 1136, "b1136.txt" },
             { "b.wav", 1647, "b1647.txt" },
             { "b.wav", 1940, "b1940.txt" },
             { "d.wav", 1658, "d1658.txt" },
             { "d.wav", 352, "d352.txt" },
             { "d.wav", 1866, "d1866.txt" },
             { "e.wav", 1675, "e1675.txt" },
             // { "e.wav", 986, "e986.txt" },
             { "e.wav", 983, "e986x.txt" },
             { "f.wav", 974, "f974.txt" },
             { "f.wav", 1156, "f1156.txt" },
             { "f.wav", 1329, "f1329.txt" },
             { "f.wav", 1514, "f1514.txt" },
#endif
#if 0
             { "f.wav", 974, "f974.txt" },
             { "f.wav", 1156, "f1156.txt" },
             { "f.wav", 1329, "f1329.txt" },
             { "f.wav", 1514, "f1514.txt" },
#endif
             { 0, 0, 0 }
};


//
// hz[] is an array of hz's to monitor.
// wanted[i] is the name of the file containing what we expect to see at hz[i].
//
std::vector<double>
benchmark1(const char *wavname, std::vector<double> hz, std::vector<const char *> wanted,
           double &nsig)
{
  assert(hz.size() == wanted.size());
  
  int filerate;
  std::vector<double> v = readwav(wavname, filerate);
  int demodrate;
  if((filerate % 8000) == 0 && filerate > 8000){
    // reduce to 8000
    v = thin(v, filerate / 8000);
    demodrate = 8000;
  } else {
    demodrate = filerate;
  }

  std::vector<std::string> out(hz.size());
  
  Demod *dm = new Demod(demodrate);
  for(int i = 0; i < hz.size(); i++){
    dm->watch(hz[i],
              [ ](void*arg, double, double, const char *txt) {
                std::string *ooo = (std::string*) arg;
                (*ooo) += txt;
              },
              (void*)&(out[i]));
  }

  dm->got(v);

  nsig = dm->nsig();
  
  delete dm;

  // unique number to name temporary files.
  static std::mutex mu;
  static int bseq = 0;
  int seq;
  mu.lock();
  seq = bseq;
  bseq++;
  mu.unlock();

  char odi[256];
  char odo[256];

  std::vector<double> scores;

  for(int i = 0; i < hz.size(); i++){
    sprintf(odi, "/tmp/odi-%d-%d", getpid(), seq);
    sprintf(odo, "/tmp/odo-%d-%d", getpid(), seq);
  
    unlink(odi);
    FILE *fp = fopen(odi, "w");
    assert(fp);
    fprintf(fp, "%s", out[i].c_str());
    fclose(fp);
  
    unlink(odo);
    char cmd[512];
    sprintf(cmd, "../rtty-diff.pl %s %s > %s",
            wanted[i], odi, odo);
    system(cmd);
    
    double score = 0;
    {
      FILE *fp = fopen(odo, "r");
      assert(fp);
      int nscan = fscanf(fp, "%lf", &score);
      fclose(fp);
      assert(nscan == 1);
    }

    scores.push_back(score);
    
    unlink(odo);
    unlink(odi);
  }

  return scores;
}

int
get_ncpus()
{
  int count = 0;
  
#ifdef __linux__
  FILE *fp = popen("nproc", "r");
  if(fp){
    fscanf(fp, "%d", &count);
    pclose(fp);
    if(count > 0)
      return count;
  }
#else
  // mac, FreeBSD
  size_t count_len = sizeof(count);
  sysctlbyname("hw.ncpu", &count, &count_len, NULL, 0);
  if(count > 0)
    return count;
#endif
  
  return 1;
}
  

double
benchmark(int verbose, double &nsig)
{
  int ncpus = get_ncpus();

  ncpus /= 2;
  if(ncpus < 1){
    ncpus = 1;
  }

  std::vector<std::thread*> threads;
  
  // XXX one per files[]
  double scv[100];
  double nsv[100]; // average number of signals
  bool done[100];
  bool printed[100];

  int entries = 0;
  for(int i = 0; files[i].wavname; i++){
    done[i] = false;
    printed[i] = false;
    entries++;
  }
  assert(entries < sizeof(scv) / sizeof(scv[0])); // XXX

  int fi;
  for(fi = 0; fi < entries; ){
    if(threads.size() >= ncpus){
      threads[0]->join();
      delete threads[0];
      threads.erase(threads.begin());
    }

    for(int i = 0; i < entries; i++){
      if(done[i] && printed[i] == false){
        printed[i] = true;
        if(verbose)
          printf("%s %.1f %s : %.3f %.0f\n", files[i].wavname, files[i].hz, files[i].wanted, scv[i], nsv[i]);
      }
    }

    std::vector<double> hzv;
    std::vector<const char *> wantedv;

    int fi1 = fi;
    for( ; fi1 < entries && strcmp(files[fi1].wavname, files[fi].wavname) == 0; fi1++){
      hzv.push_back(files[fi1].hz);
      wantedv.push_back(files[fi1].wanted);
    }

    int n = fi1 - fi;
    std::thread *th = new std::thread( [ fi, n, hzv, wantedv, &scv, &nsv, &done ] () {
                                 
                                 double nsig;
                                 std::vector<double> scores;
                                 scores = benchmark1(files[fi].wavname, hzv, wantedv, nsig);

                                 assert(scores.size() == hzv.size());
                                 assert(scores.size() == n);

                                 for(int i = 0; i < n; i++){
                                   scv[fi+i] = scores[i];
                                   nsv[fi+i] = nsig;
                                   done[fi+i] = true;
                                 }
                               } );

    fi += n;

    threads.push_back(th);
  }

  while(threads.size() > 0){
    threads[0]->join();
    delete threads[0];
    threads.erase(threads.begin());

    for(int i = 0; i < entries; i++){
      if(done[i] && printed[i] == false){
        printed[i] = true;
        if(verbose)
          printf("%s %.1f %s : %.3f %.0f\n", files[i].wavname, files[i].hz, files[i].wanted, scv[i], nsv[i]);
      }
    }
  }

  double sum = 0;
  double nsigsum = 0;
  for(int i = 0; i < entries; i++){
    sum += scv[i];
    nsigsum += nsv[i];
  }
  
  if(verbose)
    printf("%.3f %.0f\n", sum / entries, nsigsum / entries);
  nsig = nsigsum / entries;
  return sum / entries;
}

void
optimize()
{
  for(int i = 0; vars[i].name; i++){
    double def = *vars[i].v;
    for(int j = 0; vars[i].vals[j] > -9998; j++){
      printf("%s %.2f : ", vars[i].name, vars[i].vals[j]);
      fflush(stdout);
      *vars[i].v = vars[i].vals[j];

      double nsig = 0;
      double score = benchmark(0, nsig);

      printf("%.3f %.0f\n", score, nsig);
    }
    *vars[i].v = def;
  }
}
