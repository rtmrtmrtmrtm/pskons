#ifndef FFT_H
#define FFT_H

#include <vector>
#include <complex>

std::vector<std::complex<double>> one_fft(const std::vector<double> &samples,
                                          int i0, int block);

std::vector<double> hilbert_shift(const std::vector<double> &x,
                                  double hz0, double hz1, int rate);


#endif
