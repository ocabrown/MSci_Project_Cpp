#ifndef Combined_Opacities_Function_hpp
#define Combined_Opacities_Function_hpp

#include <vector>

double Planck_Deriv(double l, double T);
double Q(double x);
std::vector<double> dust_opa(std::vector<double> ls, double a);
double Dust_Ross(double T, double a);
double Combined_Opacity(double a_d, double d_to_g, double P, double T);

#endif
