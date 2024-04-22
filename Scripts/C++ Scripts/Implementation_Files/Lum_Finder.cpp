#include "./../Header_Files/Lum_Finder.hpp"

#include "./../Header_Files/IsoO_Model.hpp"
#include "./../Header_Files/FreO_Model.hpp"
#include "./../Header_Files/ComO_Model.hpp"

#include <vector>
#include <array>
#include <string>





std::pair< std::vector< std::vector<double> >, std::vector< std::vector<double> > > L_Finder_IsoO(double mc, double me, double rp, double pd, double Td, double kd, double L_guess, double acc_l, int n, std::vector<double> m_grid_ini)
{
    const double Ma_E = 5.972e27;               // Earth mass [g]
    const double pc = 3.2;                      // Core density [g/cm^3]
    const double pi = 3.1415926535;             // Pi
    
    std::vector<double> rc = {pow(3. * mc * Ma_E /(4 * pi * pc), 1./3.)};
    
    std::vector<double> L = {L_guess};
    double L_step = 1./2.;
    double e_r = 1.;
    int n_i = 0;
    
    std::array<double, 2> p_ana_c;
    std::vector<double> n_is, Ls, rs;
    std::vector<std::string> Plus_Minus;
    std::vector< std::vector<double> > var;
    double r_inner;
    
    while (abs(e_r) > acc_l)
    {
        n_i += 1;
        n_is.push_back(n_i);
        
        Iso_Opacity_Giant_Planet IsoO(mc, me, rp, pd, Td, kd, L[0], n, true, m_grid_ini, p_ana_c);
        var = IsoO.Planet_Formation("Both");
        r_inner = var[1][1];
        
        Ls.push_back(L[0]);
        rs.push_back(r_inner);
        
        e_r = r_inner / rc[0] - 1.;
        
        if (abs(e_r) >  acc_l)
        {
            if (e_r > 0.)
            {
                Plus_Minus.push_back("+");
                if (n_i > 1)
                {
                    if (Plus_Minus[Plus_Minus.size()-2] == "+")
                    {
                        L_step *= 2.;
                    }
                }
                L[0] *= (1. + L_step);
            }
            else if (e_r < 0.)
            {
                Plus_Minus.push_back("-");
                if (n_i > 1)
                {
                    if (Plus_Minus[Plus_Minus.size()-2] == "-")
                    {
                        L_step *= 2.;
                    }
                }
                L[0] *= (1. - L_step);
            }
        }
        
        else if (abs(e_r) <= acc_l)
        {
            L[0] = L[0];
        }
        
        L_step /= 2.;
    }
    
    std::pair< std::vector< std::vector<double> >, std::vector< std::vector<double> > > result = {{L, Ls, rc, rs, n_is}, var};
    
    return result;
}



std::pair< std::vector< std::vector<double> >, std::vector< std::vector<double> > > L_Finder_FreO(double mc, double me, double rp, double pd, double Td, double kd, double L_guess, double met, double acc_k, double acc_l, int n, std::vector<double> m_grid_ini)
{
    const double Ma_E = 5.972e27;               // Earth mass [g]
    const double pc = 3.2;                      // Core density [g/cm^3]
    const double pi = 3.1415926535;             // Pi
    
    std::vector<double> rc = {pow(3. * mc * Ma_E /(4 * pi * pc), 1./3.)};
    
    std::vector<double> L = {L_guess};
    double L_step = 1./2.;
    double e_r = 1.;
    int n_i = 0;
    
    std::array<double, 2> p_ana_c;
    std::vector<double> n_is, Ls, rs;
    std::vector<std::string> Plus_Minus;
    std::vector< std::vector<double> > var;
    std::vector<bool> breaks;
    double r_inner;
    
    while (abs(e_r) > acc_l)
    {
        n_i += 1;
        n_is.push_back(n_i);
        Fre_Opacity_Giant_Planet FreO(mc, me, rp, pd, Td, kd, L[0], met, acc_k, n, true, m_grid_ini, p_ana_c);
        var = FreO.Planet_Formation();
        r_inner = var[1][1];
        
        Ls.push_back(L[0]);
        rs.push_back(r_inner);
        
        breaks.push_back(FreO.getInfo().second);
        
        if (breaks[breaks.size()-1])
        {
            Plus_Minus.push_back("-");
            n_is.pop_back();
            Ls.pop_back();
            rs.pop_back();
            L[0] *= (1. - L_step);
            continue;
        }
        
        e_r = r_inner / rc[0] - 1.;
        
        if (abs(e_r) >  acc_l)
        {
            if (e_r > 0.)
            {
                Plus_Minus.push_back("+");
                if (n_i > 1)
                {
                    if (Plus_Minus[Plus_Minus.size()-2] == "+")
                    {
                        L_step *= 2.;
                    }
                }
                L[0] *= (1. + L_step);
            }
            else if (e_r < 0.)
            {
                Plus_Minus.push_back("-");
                if (n_i > 1)
                {
                    if (Plus_Minus[Plus_Minus.size()-2] == "-")
                    {
                        L_step *= 2.;
                    }
                }
                L[0] *= (1. - L_step);
            }
        }
        
        else if (abs(e_r) <= acc_l)
        {
            L[0] = L[0];
        }
        
        L_step /= 2.;
    }
    
    std::pair< std::vector< std::vector<double> >, std::vector< std::vector<double> > > result = {{L, Ls, rc, rs, n_is}, var};
    
    return result;
}



std::pair< std::vector< std::vector<double> >, std::vector< std::vector<double> > > L_Finder_ComO(double mc, double me, double rp, double pd, double Td, double kd, double L_guess, double ad, double d_to_g, double acc_k, double acc_l, int n, std::vector<double> m_grid_ini)
{
    const double Ma_E = 5.972e27;               // Earth mass [g]
    const double pc = 3.2;                      // Core density [g/cm^3]
    const double pi = 3.1415926535;             // Pi
    
    std::vector<double> rc = {pow(3. * mc * Ma_E /(4 * pi * pc), 1./3.)};
    
    std::vector<double> L = {L_guess};
    double L_step = 1./2.;
    double e_r = 1.;
    int n_i = 0;
    
    std::array<double, 2> p_ana_c;
    std::vector<double> n_is, Ls, rs;
    std::vector<std::string> Plus_Minus;
    std::vector< std::vector<double> > var;
    std::vector<bool> breaks;
    double r_inner;
    
    while (abs(e_r) > acc_l)
    {
        n_i += 1;
        n_is.push_back(n_i);
        Com_Opacity_Giant_Planet ComO(mc, me, rp, pd, Td, kd, L[0], ad, d_to_g, acc_k, n, true, m_grid_ini, p_ana_c);
        var = ComO.Planet_Formation();
        r_inner = var[1][1];
        
        Ls.push_back(L[0]);
        rs.push_back(r_inner);
        
        breaks.push_back(ComO.getInfo().second);
        
        if (breaks[breaks.size()-1])
        {
            Plus_Minus.push_back("-");
            n_is.pop_back();
            Ls.pop_back();
            rs.pop_back();
            L[0] *= (1. - L_step);
            continue;
        }
        
        e_r = r_inner / rc[0] - 1.;
        
        if (abs(e_r) >  acc_l)
        {
            if (e_r > 0.)
            {
                Plus_Minus.push_back("+");
                if (n_i > 1)
                {
                    if (Plus_Minus[Plus_Minus.size()-2] == "+")
                    {
                        L_step *= 2.;
                    }
                }
                L[0] *= (1. + L_step);
            }
            else if (e_r < 0.)
            {
                Plus_Minus.push_back("-");
                if (n_i > 1)
                {
                    if (Plus_Minus[Plus_Minus.size()-2] == "-")
                    {
                        L_step *= 2.;
                    }
                }
                L[0] *= (1. - L_step);
            }
        }
        
        else if (abs(e_r) <= acc_l)
        {
            L[0] = L[0];
        }
        
        L_step /= 2.;
    }
    
    std::pair< std::vector< std::vector<double> >, std::vector< std::vector<double> > > result = {{L, Ls, rc, rs, n_is}, var};
    
    return result;
}
