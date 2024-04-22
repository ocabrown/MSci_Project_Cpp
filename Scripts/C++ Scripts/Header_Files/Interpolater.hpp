#ifndef Interpolater_hpp
#define Interpolater_hpp

#include <vector>

int findNearestNeighbourIndex(double value, std::vector<double> &x);
std::vector<double> interp1(std::vector<double> &x, std::vector<double> &y, std::vector<double> &x_new);

#endif
