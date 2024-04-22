#ifndef Max_Acc_Rate_Calc_hpp
#define Max_Acc_Rate_Calc_hpp

#include <vector>

std::vector< std::vector<double> > Sigma_Calc(double Sigma_outer, double q, double H_r, double nu, int num_sig);
std::pair< std::vector< std::vector<double> >, std::vector<bool> > Max_Acc_Rate(double Sigma_outer, double Ms, double Mp, double Tp, double Rp, double r_in, double r_out, double alpha, int num_sig, bool Const_Sigma, bool Crida_Comparison, bool L09);

#endif
