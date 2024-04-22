#include "./../Header_Files/IsoT_Model.hpp"

#include "./../Header_Files/Linspace.hpp"
#include "./../Header_Files/Linear_Fit.hpp"

#include <vector>
#include <array>
#include <iostream>
#include <boost/math/tools/minima.hpp>      //https://www.boost.org/doc/libs/1_70_0/libs/math/doc/html/math_toolkit/brent_minima.html
using boost::math::tools::brent_find_minima;





Iso_Temperature_Giant_Planet::Iso_Temperature_Giant_Planet(double mc, double me, double rp, double pd, double Td, int num_p, bool grid_given, std::vector<double> m_grid, std::array<double, 2> p_ana_c) : mc(mc), me(me), rp(rp), pd(pd), Td(Td), num_p(num_p), grid_given(grid_given), m_grid(m_grid), p_ana_c(p_ana_c)                                       // Custom constructor
{
    core_mass_max = mc * Ma_E;
    env_mass_max = me * Ma_E;
    Mp = core_mass_max + env_mass_max;
    R_H = (5.2*AU) * pow(Mp/(3.*Ma_S), 1./3.);
    cs2 = R * Td;
//    rad_hi = (G*Mp/((cs2/1.)+((G*Mp)/(0.25*R_H))));
//    std::cout << rad_hi/R_J << std::endl;
    rad_hi = rp * R_J;
    den_hi = pd;
    tem_hi = Td;
    pre_hi = pd * R * Td;
    
    // Making empty arrays for other variables with boundary conditions
    std::vector<double> empty(num_p + 2 * num_g, 0.0);
    r = empty;
    p = empty;
    T = empty;
    P = empty;
    s = empty;
    p_ana = empty;
    dP_dr_num = empty;
    dP_dr_ana = empty;
    
    // Putting initial values in each array
    r[hi] = rad_hi;
    p[hi] = den_hi;
    p_ana[hi] = den_hi;
    T[hi] = tem_hi;
    P[hi] = pre_hi;
    
    if (grid_given)
    {
        m = m_grid;
    }
    
    else
    {
        // Calculating envelope mass array so r has Uniform logarithmic spacing
        env_mass_min = 1.e20;          // [g]
        
        // Lambda for integrating density function of r to give m
        auto m_r_func = [this, p_ana_c](double r)
        {
            double a = p_ana_c[0];
            double b = p_ana_c[1];
            
            double integral = -(4.*pi*(pow(a*r,2.) + 2.*a*r + 2.) * exp(-(a*r + b))) / pow(a,3.);
            double C = env_mass_max - (-(4.*pi*(pow(a*rad_hi,2.) + 2.*a*rad_hi + 2.) * exp(-(a*rad_hi + b))) / pow(a,3.));
            return integral + C;
        };
        
        // Lambda to calculate minimum mass
        auto m_min = [this, m_r_func](double r)
        {
            return env_mass_min - m_r_func(r);
        };
        
        // Finding the two radii where the mass flips from +ve to -ve
        std::vector<double> rad_lo_guess = linspace(1.e8, 1.e12, 10000);
        std::vector<double> m_mins;
        for (int x = 0; x < rad_lo_guess.size(); x++)
        {
            double rad_lo_x = rad_lo_guess[x];
            double m_min_x = m_min(rad_lo_x);
            m_mins.push_back(m_min_x);
        }
        
        std::vector<double> rad_lo_lims;
        for (int x = 0; x < m_mins.size() - 1; x++)
        {
            if (m_mins[x] / m_mins[x+1] < 0.)
            {
                rad_lo_lims = {rad_lo_guess[x], rad_lo_guess[x+1]};
                break;
            }
        }
        
        // Brents method to find the zero point between the two radii
        double rad_lo = brent_find_minima(m_min, rad_lo_lims[0], rad_lo_lims[1], double_bits).first;
        
        // Creating m array from value found such that inner radius gives 0 mass
        std::vector<double> r_lin_geom = geomspace(rad_lo, rad_hi, num_p, true);
        
        for (int i = 0; i < r_lin_geom.size(); i++)
        {
            m.push_back(m_r_func(r_lin_geom[i]));
        }
        
        for (int i = 0; i < num_g; i++)
        {
            m.insert(m.begin(), 0.);
            m.push_back(0.);
        }
    }
    
//    // Uniform logarithmic spacing for me (envelope mass)
//    env_mass_min = 1.e20;
//    std::vector<double> m1 = geomspace(env_mass_min, 0.9*env_mass_max, static_cast<int>(0.35*num_p), true);
//    std::vector<double> m2 = geomspace((0.9*env_mass_max), env_mass_max, static_cast<int>((0.65*num_p)+1), true);
//    m2.erase(m2.begin());
//    m1.insert(m1.end(), m2.begin(), m2.end());
//    m = m1;
//    //std::vector<double> m = geomspace(env_mass_min, env_mass_max, num_p, true);
//    for (int i = 0; i < num_g; i++)
//    {
//        m.insert(m.begin(), 0.);
//        m.push_back(0.);
//    }
}


Iso_Temperature_Giant_Planet::~Iso_Temperature_Giant_Planet(){}       // Destructor


std::vector< std::vector<double> > Iso_Temperature_Giant_Planet::getVariables(void)
{
    return {m, r, p, T, P, s, p_ana, dP_dr_num, dP_dr_ana};
}


void Iso_Temperature_Giant_Planet::printVariable(std::vector<double> variable)
{
    for (int i = hi + 1; i >= lo - 1; i--)
    {
        std::cout << variable[i] << ", ";
    };
    std::cout << std::endl;
}


std::vector< std::vector<double> > Iso_Temperature_Giant_Planet::Planet_Formation(void)
{
    for (int i = hi; i > lo; i--)
    {
        r[i-1] = pow((pow(r[i],3.) + 3.*(m[i-1]-m[i])/(4.*pi*p[i])), 1./3.);
        T[i-1] = T[i];
        P[i-1] = P[i] - G*(m[i-1]+core_mass_max)*(m[i-1]-m[i])/(4.*pi*pow(r[i-1],4.));
        p[i-1] = P[i-1]/(R*T[i-1]);
        s[i-1] = T[i-1]/pow(P[i-1], (gamma-1.)/gamma);
        
        p_ana[i-1] = den_hi * exp(G*(m[i-1]+core_mass_max)/(R*tem_hi) * ((1./r[i-1]) - (1./rad_hi)));
        dP_dr_num[i-1] = (P[i-1]-P[i])/(r[i-1]-r[i]);
        dP_dr_ana[i-1] = -G*(m[i-1]+core_mass_max)*p[i-1]/pow(r[i-1], 2.);
    }
    
    variables = {m, r, p, T, P, s, p_ana, dP_dr_num, dP_dr_ana};
    
    return variables;
}


std::vector<double> Iso_Temperature_Giant_Planet::L2norm(void)
{
    std::vector<double> m_(variables[0].begin()+lo+num_p/10, variables[0].begin()+hi-num_p/10);
    std::vector<double> r_(variables[1].begin()+lo+num_p/10, variables[1].begin()+hi-num_p/10);
    std::vector<double> p_(variables[2].begin()+lo+num_p/10, variables[2].begin()+hi-num_p/10);
    std::vector<double> T_(variables[3].begin()+lo+num_p/10, variables[3].begin()+hi-num_p/10);
    std::vector<double> P_(variables[4].begin()+lo+num_p/10, variables[4].begin()+hi-num_p/10);
    std::vector<double> s_(variables[5].begin()+lo+num_p/10, variables[5].begin()+hi-num_p/10);
    std::vector<double> p_ana_(variables[6].begin()+lo+num_p/10, variables[6].begin()+hi-num_p/10);
    std::vector<double> dP_dr_num_(variables[7].begin()+lo+num_p/10, variables[7].begin()+hi-num_p/10);
    std::vector<double> dP_dr_ana_(variables[8].begin()+lo+num_p/10, variables[8].begin()+hi-num_p/10);
    
    auto L2norm_calc = [](std::vector<double> vec1, std::vector<double> vec2)
    {
        std::vector<double> result;
        for (int i = 0; i < vec1.size(); i++)
        {
            result.push_back(pow(((vec1[i] - vec2[i]) / vec1[i]), 2.));
        }
        return result;
    };
    
    auto sum = [](std::vector<double> vec)
    {
        double summation = 0.;
        for (int i = 0; i < vec.size(); i++)
        {
            summation += vec[i];
        }
        return summation;
    };
    
    double min_dP_dr = *min_element(dP_dr_num_.begin(),dP_dr_num_.end());
    
    std::vector<double> L2norms_p = L2norm_calc(p_ana_, p_);
    std::vector<double> L2norms_dP_dr = L2norm_calc(dP_dr_ana_, dP_dr_num_);
    double L2norm_p = sqrt(sum(L2norms_p)/p_.size());
    double L2norm_dP_dr = sqrt(sum(L2norms_dP_dr)/dP_dr_num_.size());
    
    L2norm_vars = {L2norm_p, L2norm_dP_dr, min_dP_dr};
    
    return L2norm_vars;
}


std::array<double, 2> Iso_Temperature_Giant_Planet::DensityFit(void)
{
    std::vector<double> r_fit(variables[1].begin()+lo, variables[1].begin()+hi);
    std::vector<double> p_fit(variables[2].begin()+lo, variables[2].begin()+hi);
    std::vector<double> p_fit_log;
    for (int i = 0; i < p_fit.size(); i++)
    {
        p_fit_log.push_back(log(p_fit[i]));
    }
    
    std::array<double, 2> popt = lin_fit(r_fit, p_fit_log);
    popt[0] *= -1.;
    popt[1] *= -1.;
    
    return popt;
}
