#include "./../Header_Files/Mass_Grid_Solver.hpp"

#include "./../Header_Files/IsoT_Model.hpp"
#include "./../Header_Files/IsoO_Model.hpp"
#include "./../Header_Files/FreO_Model.hpp"
#include "./../Header_Files/ComO_Model.hpp"

#include <vector>
#include <array>





std::pair< std::vector<double>, std::array<double, 2> > Mass_Grid_Solver_IsoT(double mc, double me, double rp, double pd, double Td, int num_p, double acc_m, std::array<double, 2> p_ana_c_ini)
{
    double e_m0 = 1.;
    double e_m1 = 1.;
    
    std::vector<double> m_grid;
    std::array<double, 2> p_ana_c_guess = p_ana_c_ini;
    std::array<double, 2> p_ana_c_guess_new;
    int iter = 0;
    
    while(e_m0 > acc_m or e_m1 > acc_m)
    {
        iter += 1;
        Iso_Temperature_Giant_Planet IsoT(mc, me, rp, pd, Td, num_p, false, m_grid, p_ana_c_guess);
        IsoT.Planet_Formation();
        p_ana_c_guess_new = IsoT.DensityFit();
        
        e_m0 = abs(p_ana_c_guess_new[0]/p_ana_c_guess[0] - 1);
        e_m1 = abs(p_ana_c_guess_new[1]/p_ana_c_guess[1] - 1);
        
        p_ana_c_guess = p_ana_c_guess_new;
        m_grid = IsoT.getVariables()[0];
    }
    
    //Iso_Temperature_Giant_Planet IsoT(mc, me, rp, pd, Td, num_p, false, m_grid, p_ana_c_guess_new);
    //m_grid = IsoT.getVariables()[0];
    
    std::pair< std::vector<double>, std::array<double, 2> > grid_vars = {m_grid, p_ana_c_guess_new};
    
    return grid_vars;
}



std::pair< std::vector<double>, std::array<double, 2> > Mass_Grid_Solver_IsoO(double mc, double me, double rp, double pd, double Td, double kd, double lp, int num_p, double acc_m, std::array<double, 2> p_ana_c_ini)
{
    double e_m0 = 1.;
    double e_m1 = 1.;
    
    std::vector<double> m_grid;
    std::array<double, 2> p_ana_c_guess = p_ana_c_ini;
    std::array<double, 2> p_ana_c_guess_new;
    int iter = 0;
    
    while(e_m0 > acc_m or e_m1 > acc_m)
    {
        iter += 1;
        Iso_Opacity_Giant_Planet IsoO(mc, me, rp, pd, Td, kd, lp, num_p, false, m_grid, p_ana_c_guess);
        IsoO.Planet_Formation("Both");
        p_ana_c_guess_new = IsoO.DensityFit();
        
        e_m0 = abs(p_ana_c_guess_new[0]/p_ana_c_guess[0] - 1);
        e_m1 = abs(p_ana_c_guess_new[1]/p_ana_c_guess[1] - 1);
        
        p_ana_c_guess = p_ana_c_guess_new;
        m_grid = IsoO.getVariables()[0];
    }
    
    //Iso_Opacity_Giant_Planet IsoO(mc, me, rp, pd, Td, kd, lp, num_p, false, m_grid, p_ana_c_guess);
    //m_grid = IsoO.getVariables()[0];
    
    std::pair< std::vector<double>, std::array<double, 2> > grid_vars = {m_grid, p_ana_c_guess_new};
    
    return grid_vars;
}



std::pair< std::vector<double>, std::array<double, 2> > Mass_Grid_Solver_FreO(double mc, double me, double rp, double pd, double Td, double kd, double lp, double met, double acc_k, int num_p, double acc_m, std::array<double, 2> p_ana_c_ini)
{
    double e_m0 = 1.;
    double e_m1 = 1.;
    
    std::vector<double> m_grid;
    std::array<double, 2> p_ana_c_guess = p_ana_c_ini;
    std::array<double, 2> p_ana_c_guess_new;
    int iter = 0;
    
    while(e_m0 > acc_m or e_m1 > acc_m)
    {
        iter += 1;
        Fre_Opacity_Giant_Planet FreO(mc, me, rp, pd, Td, kd, lp, met, acc_k, num_p, false, m_grid, p_ana_c_guess);
        FreO.Planet_Formation();
        p_ana_c_guess_new = FreO.DensityFit();
        
        e_m0 = abs(p_ana_c_guess_new[0]/p_ana_c_guess[0] - 1);
        e_m1 = abs(p_ana_c_guess_new[1]/p_ana_c_guess[1] - 1);
        
        p_ana_c_guess = p_ana_c_guess_new;
        m_grid = FreO.getVariables()[0];
    }
    
    //Fre_Opacity_Giant_Planet FreO(mc, me, rp, pd, Td, kd, lp, met, acc_k, num_p, false, m_grid, p_ana_c_guess);
    //m_grid = FreO.getVariables()[0];
    
    std::pair< std::vector<double>, std::array<double, 2> > grid_vars = {m_grid, p_ana_c_guess_new};
    
    return grid_vars;
}



std::pair< std::vector<double>, std::array<double, 2> > Mass_Grid_Solver_ComO(double mc, double me, double rp, double pd, double Td, double kd, double lp, double ad, double d_to_g, double acc_k, int num_p, double acc_m, std::array<double, 2> p_ana_c_ini)
{
    double e_m0 = 1.;
    double e_m1 = 1.;
    
    std::vector<double> m_grid;
    std::array<double, 2> p_ana_c_guess = p_ana_c_ini;
    std::array<double, 2> p_ana_c_guess_new;
    int iter = 0;
    
    while(e_m0 > acc_m or e_m1 > acc_m)
    {
        iter += 1;
        Com_Opacity_Giant_Planet ComO(mc, me, rp, pd, Td, kd, lp, ad, d_to_g, acc_k, num_p, false, m_grid, p_ana_c_guess);
        ComO.Planet_Formation();
        p_ana_c_guess_new = ComO.DensityFit();
        
        e_m0 = abs(p_ana_c_guess_new[0]/p_ana_c_guess[0] - 1);
        e_m1 = abs(p_ana_c_guess_new[1]/p_ana_c_guess[1] - 1);
        
        p_ana_c_guess = p_ana_c_guess_new;
        m_grid = ComO.getVariables()[0];
    }
    
    //Com_Opacity_Giant_Planet ComO(mc, me, rp, pd, Td, kd, lp, ad, d_to_g, acc_k, num_p, false, m_grid, p_ana_c_guess);
    //m_grid = ComO.getVariables()[0];
    
    std::pair< std::vector<double>, std::array<double, 2> > grid_vars = {m_grid, p_ana_c_guess_new};
    
    return grid_vars;
}
