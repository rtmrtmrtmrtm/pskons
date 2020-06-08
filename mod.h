#ifndef MOD_H
#define MOD_H

#include <vector>

std::vector<double> bits2psk(const std::vector<int> &bits, double hz,
                             int rate, int starting, int ending);
 
#endif
