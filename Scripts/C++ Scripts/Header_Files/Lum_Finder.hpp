#ifndef Lum_Finder_hpp
#define Lum_Finder_hpp

#include <vector>

std::pair< std::vector< std::vector<double> >, std::vector< std::vector<double> > > L_Finder_IsoO(double mc, double me, double rp, double pd, double Td, double kd, double L_guess, double acc_l, int n, std::vector<double> m_grid_ini);
std::pair< std::vector< std::vector<double> >, std::vector< std::vector<double> > > L_Finder_FreO(double mc, double me, double rp, double pd, double Td, double kd, double L_guess, double met, double acc_k, double acc_l, int n, std::vector<double> m_grid_ini);
std::pair< std::vector< std::vector<double> >, std::vector< std::vector<double> > > L_Finder_ComO(double mc, double me, double rp, double pd, double Td, double kd, double L_guess, double ad, double d_to_g, double acc_k, double acc_l, int n, std::vector<double> m_grid_ini);

#endif
