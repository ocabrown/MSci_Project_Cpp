#ifndef Evolution_hpp
#define Evolution_hpp

#include <vector>
#include <array>

std::pair< std::vector< std::vector<double> >, std::vector< std::vector< std::vector<double> > > > EvolverIsoO(double mc, double me, double rp, double pd, double Td, double kd, double L_guess_ini, double acc_l, int n_i, double acc_m, std::array<double, 2> p_ana_c_ini, std::vector<double> m_grid, double t_tot, int n_t, double Sigma_in, double Ms, double Rp, double r_in, double r_out, double alpha, int num_sig, bool L09);

#endif
