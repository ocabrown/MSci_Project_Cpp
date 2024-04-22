#ifndef Convergence_Tester_hpp
#define Convergence_Tester_hpp

#include <vector>
#include <array>

std::vector< std::vector<double> > ConvergenceIsoT(double mc, double me, double rp, double pd, double Td, std::vector<double> num_p, double acc_m, std::array<double, 2> p_ana_c_ini);
std::vector< std::vector<double> > ConvergenceIsoD(double mc, double me, double rp, double pd, double Td, double kd, double lp, std::vector<double> num_p);
std::vector< std::vector<double> > ConvergenceIsoO(double mc, double me, double rp, double pd, double Td, double kd, double lp, std::vector<double> num_p, double acc_m, std::array<double, 2> p_ana_c_ini);
std::vector< std::vector<double> > ConvergenceFreO(double mc, double me, double rp, double pd, double Td, double kd, double lp, double met, double acc_k, std::vector<double> num_p, double acc_m, std::array<double, 2> p_ana_c_ini);
std::vector< std::vector<double> > ConvergenceComO(double mc, double me, double rp, double pd, double Td, double kd, double lp, double ad, double d_to_g, double acc_k, std::vector<double> num_p, double acc_m, std::array<double, 2> p_ana_c_ini);
std::vector< std::vector<double> > ConvergenceSigma(double Sigma_outer, double Ms, double Mp, double Tp, double Rp, double r_in, double r_out, double alpha, std::vector<double> num_sig);
std::vector< std::vector<double> > ConvergenceMdotMax(double Sigma_outer, double Ms, double Mp, double Tp, double Rp, double r_in, double r_out, double alpha, std::vector<double> num_sig);
std::vector< std::vector<double> > ConvergenceIsoO_LF(double mc, double me, double rp, double pd, double Td, double kd, double L_guess, double acc_l, double n, std::vector<double> m_grid_ini);
std::vector< std::vector<double> > ConvergenceFreO_LF(double mc, double me, double rp, double pd, double Td, double kd, double L_guess, double met, double acc_k, double acc_l, double n, std::vector<double> m_grid_ini);
std::vector< std::vector<double> > ConvergenceComO_LF(double mc, double me, double rp, double pd, double Td, double kd, double L_guess, double ad, double d_to_g, double acc_k, double acc_l, double n, std::vector<double> m_grid_ini);

#endif
