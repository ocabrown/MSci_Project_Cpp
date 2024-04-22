#include "./../Header_Files/Evolution.hpp"

#include "./../Header_Files/Mass_Grid_Solver.hpp"
#include "./../Header_Files/Lum_Finder.hpp"
#include "./../Header_Files/Max_Acc_Rate_Calc.hpp"

#include <vector>
#include <array>
#include <iostream>





std::pair< std::vector< std::vector<double> >, std::vector< std::vector< std::vector<double> > > > EvolverIsoO(double mc, double me, double rp, double pd, double Td, double kd, double L_guess_ini, double acc_l, int n_i, double acc_m, std::array<double, 2> p_ana_c_ini, std::vector<double> m_grid, double t_tot, int n_t, double Sigma_in, double Ms, double Rp, double r_in, double r_out, double alpha, int num_sig, bool L09)
{
    const double G = 6.67430e-8;                // Gravitational constant
    const double Ma_E = 5.972e27;               // Earth mass [g]
    const double L_Sun = 3.846e33;              // Solar luminosity [erg/s]
    const double pc = 3.2;                      // Core density [g/cm^3]
    const double pi = 3.1415926535;             // Pi
    const double R = pow(3. * mc * Ma_E /(4 * pi * pc), 1./3.);    // Core Radius [cm]
    
//    auto returned_values = Mass_Grid_Solver_IsoO(mc, me, rp, pd, Td, kd, L_guess_ini, n_i, acc_m, p_ana_c_ini);
//    m_grid = returned_values.first;             // mass grid
//    p_ana_c_ini = returned_values.second;       // updated density fit parameters
    
    std::pair< std::vector< std::vector<double> >, std::vector< std::vector<double> > > IsoO_LF_Data = L_Finder_IsoO(mc, me, rp, pd, Td, kd, L_guess_ini, acc_l, n_i, m_grid);
    double L0 = IsoO_LF_Data.first[0][0];
    std::vector< std::vector<double> > vari = IsoO_LF_Data.second;
    
    n_t += 1;
    
    std::vector<double> empty(n_t, 0.);
    std::vector< std::vector<double> > empty_vari = vari;
    for (int i = 0; i < empty_vari.size(); i++)
    {
        std::fill(empty_vari[i].begin(), empty_vari[i].end(), 0.);
    }
    
    std::vector< std::vector< std::vector<double> > > varis(n_t, empty_vari);
    
    t_tot *= (365.*24.*60.*60.);                // Convert input total time in years into s
    double dt = t_tot / (n_t - 1.);
    
    std::vector<double> M, L, t, M_acc, M_acc_max, M_acc_true;
    M = empty;
    L = empty;
    t = empty;
    M_acc = empty;
    M_acc_max = empty;
    M_acc_true = empty;
    
    M[0] = (mc + me) * Ma_E;
    L[0] = L0;
    varis[0] = vari;
    double tn = 0;
    t[0] = tn;
    
    double Ln, Mdot_max;
    std::vector< std::vector<double> > Max_Acc_Rate_Vals;
    std::vector<double> Sigma, R_sig, extras;
    std::pair< std::vector< std::vector<double> >, std::vector< std::vector<double> > > IsoO_LF_Data_n;
    
    for (int n = 0; n < n_t - 1; n++)
    {
        Ln = L[n] * L_Sun;
//        M[n+1] = M[n] + (Ln * R * dt) / (G * M[n]);
//        M_acc[n] = (M[n+1] - M[n]) / dt;
        M_acc[n] = (Ln * R) / (G * M[n]);
        
        if (L09)
        {
            Max_Acc_Rate_Vals = Max_Acc_Rate(Sigma_in, Ms, M[n], Td, Rp, r_in, r_out, alpha, num_sig, false, false, true).first;
        }
        else
        {
            Max_Acc_Rate_Vals = Max_Acc_Rate(Sigma_in, Ms, M[n], Td, Rp, r_in, r_out, alpha, num_sig, false, false, false).first;
        }
        Mdot_max = Max_Acc_Rate_Vals[0][0];
        Sigma = Max_Acc_Rate_Vals[1];
        R_sig = Max_Acc_Rate_Vals[2];
        extras = Max_Acc_Rate_Vals[3];
        
        M_acc_max[n] = Mdot_max;
        double factor = 365. * 24. * 60. * 60. / Ma_E;
        std::cout << M[n]/Ma_E << ", " << M_acc[n]*factor << ", " << M_acc_max[n]*factor << std::endl;
        
        if (M_acc[n] > M_acc_max[n])
        {
            std::cout << "Max Reached!" << std::endl;
            M_acc_true[n] = M_acc_max[n];
            M[n+1] = M[n] + (dt * M_acc_max[n]);
        }
        else
        {
            M_acc_true[n] = M_acc[n];
            M[n+1] = M[n] + (dt * M_acc[n]);
        }
        
        me += (M[n+1] - M[n]) / Ma_E;
        
        for (double m_element : m_grid)
        {
            m_element *= M[n+1] / M[n];
        }
        
        IsoO_LF_Data_n = L_Finder_IsoO(mc, me, rp, pd, Td, kd, L[n], acc_l, n_i, m_grid);
        L[n+1] = IsoO_LF_Data_n.first[0][0];
        varis[n+1] = IsoO_LF_Data_n.second;
        
        tn += dt;
        t[n+1] = tn;
    }
    
    for (double t_element : t)
    {
        t_element /= (365*24*60*60);            // Convert time in s into years for plotting
    }
    
    if (M_acc[M_acc.size()-1] > M_acc_max[M_acc_max.size()-1])
    {
        std::cout << "Max Reached!" << std::endl;
        M_acc_true[M_acc_true.size()-1] = M_acc_max[M_acc_max.size()-1];
        M[M.size()-1] = M[M.size()-2] + (dt * M_acc_max[M_acc_max.size()-1]);
    }
    else
    {
        M_acc_true[M_acc_true.size()-1] = M_acc[M_acc.size()-1];
        M[M.size()-1] = M[M.size()-2] + (dt * M_acc[M_acc.size()-1]);
    }
    
    std::pair< std::vector< std::vector<double> >, std::vector< std::vector< std::vector<double> > > > IsoO_Evo_Data = {{L, M, t, M_acc, M_acc_max, M_acc_true}, varis};
    
    return IsoO_Evo_Data;
}
