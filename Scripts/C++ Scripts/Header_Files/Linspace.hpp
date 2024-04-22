#ifndef Linspace_hpp
#define Linspace_hpp

#include <vector>

std::vector<double> linspace(double start_in, double end_in, int num_in);
std::vector<double> linspace(double start_in, double end_in, int num_in, int max_stencil);
std::vector<double> geomspace(double start, double stop, int num, bool endpoint);

#endif
