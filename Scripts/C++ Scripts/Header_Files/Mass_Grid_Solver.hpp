#ifndef Mass_Grid_Solver_hpp
#define Mass_Grid_Solver_hpp

#include <vector>
#include <array>

std::pair< std::vector<double>, std::array<double, 2> > Mass_Grid_Solver_IsoT(double mc, double me, double rp, double pd, double Td, int num_p, double acc_m, std::array<double, 2> p_ana_c_ini);
// No IsoD m_grid solver
std::pair< std::vector<double>, std::array<double, 2> > Mass_Grid_Solver_IsoO(double mc, double me, double rp, double pd, double Td, double kd, double lp, int num_p, double acc_m, std::array<double, 2> p_ana_c_ini);
std::pair< std::vector<double>, std::array<double, 2> > Mass_Grid_Solver_FreO(double mc, double me, double rp, double pd, double Td, double kd, double lp, double met, double acc_k, int num_p, double acc_m, std::array<double, 2> p_ana_c_ini);
std::pair< std::vector<double>, std::array<double, 2> > Mass_Grid_Solver_ComO(double mc, double me, double rp, double pd, double Td, double kd, double lp, double ad, double d_to_g, double acc_k, int num_p, double acc_m, std::array<double, 2> p_ana_c_ini);

#endif
