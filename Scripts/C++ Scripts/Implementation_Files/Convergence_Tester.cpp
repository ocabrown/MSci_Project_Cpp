#include "./../Header_Files/Convergence_Tester.hpp"

#include "./../Header_Files/Max_Acc_Rate_Calc.hpp"
#include "./../Header_Files/IsoT_Model.hpp"
#include "./../Header_Files/IsoD_Model.hpp"
#include "./../Header_Files/IsoO_Model.hpp"
#include "./../Header_Files/FreO_Model.hpp"
#include "./../Header_Files/ComO_Model.hpp"
#include "./../Header_Files/Mass_Grid_Solver.hpp"
#include "./../Header_Files/Lum_Finder.hpp"

#include <vector>
#include <array>





std::vector< std::vector<double> > ConvergenceIsoT(double mc, double me, double rp, double pd, double Td, std::vector<double> num_p, double acc_m, std::array<double, 2> p_ana_c_ini)
{
    std::vector<double> L2norms_p;
    std::vector<double> L2norms_dP_dr;
    std::vector<double> min_dP_drs;
    
    for (int i = 0; i < num_p.size(); i++)
    {
        auto returned_values = Mass_Grid_Solver_IsoT(mc, me, rp, pd, Td, (int)num_p[i], acc_m, p_ana_c_ini);
        std::vector<double> m_grid = returned_values.first;
        std::array<double, 2> p_ana_c = returned_values.second;

        Iso_Temperature_Giant_Planet IsoT(mc, me, rp, pd, Td, (int)num_p[i], true, m_grid, p_ana_c);
        IsoT.Planet_Formation();
        std::vector<double> IsoT_L2norm_data = IsoT.L2norm();
        
        L2norms_p.push_back(IsoT_L2norm_data[0]);
        L2norms_dP_dr.push_back(IsoT_L2norm_data[1]);
        min_dP_drs.push_back(IsoT_L2norm_data[2]);
    }
    
    std::vector< std::vector<double> > L2norms_vars = {L2norms_p, L2norms_dP_dr, min_dP_drs, num_p};
    
    return L2norms_vars;
}



std::vector< std::vector<double> > ConvergenceIsoD(double mc, double me, double rp, double pd, double Td, double kd, double lp, std::vector<double> num_p)
{
    std::vector<double> L2norms_p;
    std::vector<double> L2norms_dP_dr;
    std::vector<double> min_dP_drs;
    std::vector<double> L2norms_dP_dr_lower;
    std::vector<double> L2norms_dP_dr_upper;
    
    for (int i = 0; i < num_p.size(); i++)
    {
        Iso_Density_Giant_Planet IsoD(mc, me, rp, pd, Td, kd, lp, (int)num_p[i]);
        IsoD.Planet_Formation();
        std::vector<double> IsoD_L2norm_data = IsoD.L2norm();
        
        L2norms_p.push_back(IsoD_L2norm_data[0]);
        L2norms_dP_dr.push_back(IsoD_L2norm_data[1]);
        min_dP_drs.push_back(IsoD_L2norm_data[2]);
        L2norms_dP_dr_lower.push_back(IsoD_L2norm_data[3]);
        L2norms_dP_dr_upper.push_back(IsoD_L2norm_data[4]);
    }
    
    std::vector< std::vector<double> > L2norms_vars = {L2norms_p, L2norms_dP_dr, min_dP_drs, L2norms_dP_dr_lower, L2norms_dP_dr_upper, num_p};
    
    return L2norms_vars;
}



std::vector< std::vector<double> > ConvergenceIsoO(double mc, double me, double rp, double pd, double Td, double kd, double lp, std::vector<double> num_p, double acc_m, std::array<double, 2> p_ana_c_ini)
{
    std::vector<double> r0s;
    std::vector<double> L2norms_dP_dr;
    std::vector<double> min_dP_drs;
    
    num_p.push_back(10.*num_p.back());
    
    for (int i = 0; i < num_p.size(); i++)
    {
        auto returned_values = Mass_Grid_Solver_IsoO(mc, me, rp, pd, Td, kd, lp, (int)num_p[i], acc_m, p_ana_c_ini);
        std::vector<double> m_grid = returned_values.first;
        std::array<double, 2> p_ana_c = returned_values.second;

        Iso_Opacity_Giant_Planet IsoO(mc, me, rp, pd, Td, kd, lp, (int)num_p[i], true, m_grid, p_ana_c);
        std::vector< std::vector<double> > IsoO_Sim_Data = IsoO.Planet_Formation("Both");
        r0s.push_back(IsoO_Sim_Data[1][1]);
        if (i == num_p.size()-1)
        {
            continue;
        }
        std::vector<double> IsoO_L2norm_data = IsoO.L2norm();
        
        L2norms_dP_dr.push_back(IsoO_L2norm_data[0]);
        min_dP_drs.push_back(IsoO_L2norm_data[1]);
    }
    
    num_p.pop_back();
    
    std::vector< std::vector<double> > L2norms_vars = {r0s, L2norms_dP_dr, min_dP_drs, num_p};
    
    return L2norms_vars;
}



std::vector< std::vector<double> > ConvergenceFreO(double mc, double me, double rp, double pd, double Td, double kd, double lp, double met, double acc_k, std::vector<double> num_p, double acc_m, std::array<double, 2> p_ana_c_ini)
{
    std::vector<double> r0s;
    
    num_p.push_back(10.*num_p.back());
    
    for (int i = 0; i < num_p.size(); i++)
    {
        auto returned_values = Mass_Grid_Solver_FreO(mc, me, rp, pd, Td, kd, lp, met, acc_k, (int)num_p[i], acc_m, p_ana_c_ini);
        std::vector<double> m_grid = returned_values.first;
        std::array<double, 2> p_ana_c = returned_values.second;

//        std::vector<double> m_grid;
//        Fre_Opacity_Giant_Planet FreO(mc, me, rp, pd, Td, kd, lp, met, acc_k, (int)num_p[i], false, m_grid, p_ana_c_ini);
        Fre_Opacity_Giant_Planet FreO(mc, me, rp, pd, Td, kd, lp, met, acc_k, (int)num_p[i], true, m_grid, p_ana_c);
        std::vector< std::vector<double> > FreO_Sim_Data = FreO.Planet_Formation();
        r0s.push_back(FreO_Sim_Data[1][1]);
    }
    
    num_p.pop_back();
    
    std::vector< std::vector<double> > L2norms_vars = {r0s, num_p};
    
    return L2norms_vars;
}



std::vector< std::vector<double> > ConvergenceComO(double mc, double me, double rp, double pd, double Td, double kd, double lp, double ad, double d_to_g, double acc_k, std::vector<double> num_p, double acc_m, std::array<double, 2> p_ana_c_ini)
{
    std::vector<double> r0s;
    
    num_p.push_back(10.*num_p.back());
    
    for (int i = 0; i < num_p.size(); i++)
    {
        auto returned_values = Mass_Grid_Solver_ComO(mc, me, rp, pd, Td, kd, lp, ad, d_to_g, acc_k, (int)num_p[i], acc_m, p_ana_c_ini);
        std::vector<double> m_grid = returned_values.first;
        std::array<double, 2> p_ana_c = returned_values.second;

//        std::vector<double> m_grid;
//        Com_Opacity_Giant_Planet ComO(mc, me, rp, pd, Td, kd, lp, ad, d_to_g, acc_k, (int)num_p[i], false, m_grid, p_ana_c_ini);
        Com_Opacity_Giant_Planet ComO(mc, me, rp, pd, Td, kd, lp, ad, d_to_g, acc_k, (int)num_p[i], true, m_grid, p_ana_c);
        std::vector< std::vector<double> > ComO_Sim_Data = ComO.Planet_Formation();
        r0s.push_back(ComO_Sim_Data[1][1]);
    }
    
    num_p.pop_back();
    
    std::vector< std::vector<double> > L2norms_vars = {r0s, num_p};
    
    return L2norms_vars;
}



std::vector< std::vector<double> > ConvergenceSigma(double Sigma_outer, double Ms, double Mp, double Tp, double Rp, double r_in, double r_out, double alpha, std::vector<double> num_sig)
{
    std::vector<double> S0s;
    std::vector<double> Sigma;
    
    num_sig.push_back(10.*num_sig.back());
    
    for (int i = 0; i < num_sig.size(); i++)
    {
        int num;
        if (i == 0)
        {
            num = (int)num_sig[i];
        }
        else
        {
            num = (int)num_sig[i]+1;
        }
        Sigma = Max_Acc_Rate(Sigma_outer, Ms, Mp, Tp, Rp, r_in, r_out, alpha, num, false, false, false).first[1];
        S0s.push_back(Sigma[0]);
    }
    
    num_sig.pop_back();
    
    std::vector< std::vector<double> > L2norms_vars = {S0s, num_sig};
    
    return L2norms_vars;
}



std::vector< std::vector<double> > ConvergenceMdotMax(double Sigma_outer, double Ms, double Mp, double Tp, double Rp, double r_in, double r_out, double alpha, std::vector<double> num_sig)
{
    std::vector<double> Mdot_rel_diffs;
    double Mdot_max, Mdot_max_ana, Mdot_rel_diff;
    
    for (int i = 0; i < num_sig.size(); i++)
    {
        int num;
        if (i == 0)
        {
            num = (int)num_sig[i];
        }
        else
        {
            num = (int)num_sig[i]+1;
        }
        std::vector< std::vector<double> > Mdots = Max_Acc_Rate(Sigma_outer, Ms, Mp, Tp, Rp, r_in, r_out, alpha, num, true, false, false).first;
        Mdot_max = Mdots[0][0];
        Mdot_max_ana = Mdots[1][0];
        Mdot_rel_diff = abs((Mdot_max - Mdot_max_ana) / Mdot_max_ana);
        Mdot_rel_diffs.push_back(Mdot_rel_diff);
    }
    
    std::vector< std::vector<double> > L2norms_vars = {Mdot_rel_diffs, num_sig};
    
    return L2norms_vars;
}



std::vector< std::vector<double> > ConvergenceIsoO_LF(double mc, double me, double rp, double pd, double Td, double kd, double L_guess, double acc_l, double n, std::vector<double> m_grid_ini)
{
    const double Ma_E = 5.972e27;               // Earth mass [g]
    const double pc = 3.2;                      // Core density [g/cm^3]
    const double pi = 3.1415926535;             // Pi
    
    double rc = pow(3. * mc * Ma_E /(4 * pi * pc), 1./3.);
    
    std::pair< std::vector< std::vector<double> >, std::vector< std::vector<double> > > IsoO_LF_Data = L_Finder_IsoO(mc, me, rp, pd, Td, kd, L_guess, acc_l, n, m_grid_ini);
    
    std::vector<double> Ls = IsoO_LF_Data.first[1];
    std::vector<double> rs = IsoO_LF_Data.first[3];
    std::vector<double> n_is = IsoO_LF_Data.first[4];
    
    double L_mean_sum = 0.;
    double sum_num = 4.;
    for (int i = (int)Ls.size()-(int)sum_num; i < (int)Ls.size(); i++)
    {
        L_mean_sum += Ls[i];
    }
    double L_mean = L_mean_sum/sum_num;
    
    std::vector<double> L_rel_diff, r_rel_diff;
    for (int i = 0; i < Ls.size(); i++)
    {
        L_rel_diff.push_back(abs((Ls[i] - L_mean) / L_mean));
        r_rel_diff.push_back(abs((rs[i] - rc) / rc));
    }
    
    std::vector< std::vector<double> > L2norms_vars = {L_rel_diff, r_rel_diff, n_is};
    
    return L2norms_vars;
}



std::vector< std::vector<double> > ConvergenceFreO_LF(double mc, double me, double rp, double pd, double Td, double kd, double L_guess, double met, double acc_k, double acc_l, double n, std::vector<double> m_grid_ini)
{
    const double Ma_E = 5.972e27;               // Earth mass [g]
    const double pc = 3.2;                      // Core density [g/cm^3]
    const double pi = 3.1415926535;             // Pi
    
    double rc = pow(3. * mc * Ma_E /(4 * pi * pc), 1./3.);
    
    std::pair< std::vector< std::vector<double> >, std::vector< std::vector<double> > > FreO_LF_Data = L_Finder_FreO(mc, me, rp, pd, Td, kd, L_guess, met, acc_k, acc_l, n, m_grid_ini);
    
    std::vector<double> Ls = FreO_LF_Data.first[1];
    std::vector<double> rs = FreO_LF_Data.first[3];
    std::vector<double> n_is = FreO_LF_Data.first[4];
    
    double L_mean_sum = 0.;
    double sum_num = 4.;
    for (int i = (int)Ls.size()-(int)sum_num; i < (int)Ls.size(); i++)
    {
        L_mean_sum += Ls[i];
    }
    double L_mean = L_mean_sum/sum_num;
    
    std::vector<double> L_rel_diff, r_rel_diff;
    for (int i = 0; i < Ls.size(); i++)
    {
        L_rel_diff.push_back(abs((Ls[i] - L_mean) / L_mean));
        r_rel_diff.push_back(abs((rs[i] - rc) / rc));
    }
    
    std::vector< std::vector<double> > L2norms_vars = {L_rel_diff, r_rel_diff, n_is};
    
    return L2norms_vars;
}



std::vector< std::vector<double> > ConvergenceComO_LF(double mc, double me, double rp, double pd, double Td, double kd, double L_guess, double ad, double d_to_g, double acc_k, double acc_l, double n, std::vector<double> m_grid_ini)
{
    const double Ma_E = 5.972e27;               // Earth mass [g]
    const double pc = 3.2;                      // Core density [g/cm^3]
    const double pi = 3.1415926535;             // Pi
    
    double rc = pow(3. * mc * Ma_E /(4 * pi * pc), 1./3.);
    
    std::pair< std::vector< std::vector<double> >, std::vector< std::vector<double> > > ComO_LF_Data = L_Finder_ComO(mc, me, rp, pd, Td, kd, L_guess, ad, d_to_g, acc_k, acc_l, n, m_grid_ini);
    
    std::vector<double> Ls = ComO_LF_Data.first[1];
    std::vector<double> rs = ComO_LF_Data.first[3];
    std::vector<double> n_is = ComO_LF_Data.first[4];
    
    double L_mean_sum = 0.;
    double sum_num = 4.;
    for (int i = (int)Ls.size()-(int)sum_num; i < (int)Ls.size(); i++)
    {
        L_mean_sum += Ls[i];
    }
    double L_mean = L_mean_sum/sum_num;
    
    std::vector<double> L_rel_diff, r_rel_diff;
    for (int i = 0; i < Ls.size(); i++)
    {
        L_rel_diff.push_back(abs((Ls[i] - L_mean) / L_mean));
        r_rel_diff.push_back(abs((rs[i] - rc) / rc));
    }
    
    std::vector< std::vector<double> > L2norms_vars = {L_rel_diff, r_rel_diff, n_is};
    
    return L2norms_vars;
}
