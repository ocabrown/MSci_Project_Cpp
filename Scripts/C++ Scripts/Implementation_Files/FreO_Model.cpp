#include "./../Header_Files/FreO_Model.hpp"

#include "./../Header_Files/Linspace.hpp"
#include "./../Header_Files/Linear_Fit.hpp"
#include "./../Header_Files/Sign.hpp"
#include "./../Header_Files/Freedman_Opacities_Function.hpp"

#include <vector>
#include <array>
#include <iostream>
#include <boost/math/tools/minima.hpp>      //https://www.boost.org/doc/libs/1_70_0/libs/math/doc/html/math_toolkit/brent_minima.html
using boost::math::tools::brent_find_minima;





Fre_Opacity_Giant_Planet::Fre_Opacity_Giant_Planet(double mc, double me, double rp, double pd, double Td, double kd, double lp, double met, double acc_k, int num_p, bool grid_given, std::vector<double> m_grid, std::array<double, 2> p_ana_c) : mc(mc), me(me), rp(rp), pd(pd), Td(Td), kd(kd), lp(lp), met(met), acc_k(acc_k), num_p(num_p), grid_given(grid_given), m_grid(m_grid), p_ana_c(p_ana_c)                                     // Custom constructor
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
    opa_hi = kd;
    pre_hi = pd * R * Td;
    lum = lp * L_Sun;
    Z = met;
    
    // Making empty arrays for other variables with boundary conditions
    std::vector<double> empty(num_p + 2 * num_g, 0.0);
    r = empty;
    p = empty;
    T = empty;
    k = empty;
    P = empty;
    s = empty;
    
    // Putting initial values in each array
    r[hi] = rad_hi;
    p[hi] = den_hi;
    T[hi] = tem_hi;
    k[hi] = opa_hi;
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
}


Fre_Opacity_Giant_Planet::~Fre_Opacity_Giant_Planet(){}       // Destructor


std::vector< std::vector<double> > Fre_Opacity_Giant_Planet::getVariables(void)
{
    return {m, r, p, T, k, P, s};
}


std::pair< std::pair< std::vector<int>, std::vector< std::pair<int,int> > >, bool > Fre_Opacity_Giant_Planet::getInfo(void)
{
    return {{rad_or_conv, boundary_i}, broken};
}


void Fre_Opacity_Giant_Planet::printVariable(std::vector<double> variable)
{
    for (int i = hi + 1; i >= lo - 1; i--)
    {
        std::cout << variable[i] << ", ";
    };
    std::cout << std::endl;
}


std::vector< std::vector<double> > Fre_Opacity_Giant_Planet::Planet_Formation(void)
{
    for (int i = hi; i > lo; i--)
    {
        r[i-1] = pow((pow(r[i],3.) + 3.*(m[i-1]-m[i])/(4.*pi*p[i])), 1./3.);
        if ((pow(r[i],3.) + 3.*(m[i-1]-m[i])/(4.*pi*p[i])) <= 0.)
        {
            broken = true;
            std::cout << "Broke" << std::endl;
            break;
        }
//        if (r[i-1] < 1.64686048e9)
//        {
//            r[i-1] = 1.64686048e9;
//            broken = true;
//            break
//        }
        
        P[i-1] = P[i] - G*(m[i-1]+core_mass_max)*(m[i-1]-m[i])/(4.*pi*pow(r[i-1],4.));
        
        // Assuming adiabatic:
        T[i-1] = T[i] * pow(P[i-1]/P[i], (gamma-1.)/gamma);
        k[i-1] = Freedman_Opacity(Z, P[i-1], T[i-1]);
        
        double A_adia = T[i] * pow(P[i-1]/P[i], (gamma-1.)/gamma);
        
        double Grad_Rad = 3.*k[i-1]*lum*P[i]/(16.*pi*a*c*G*(m[i-1]+core_mass_max)*pow(T[i],4));
        
        if (abs(Grad_Rad) < abs(Grad_Adia))     // Radiative
        {
            rad_or_conv.push_back(1);
            double B = T[i];
            std::vector<double> A_guess = linspace(B+1.e-10, A_adia, 10);
            
            auto g_func = [this, i, B](double A)
            {
                double C = 3.*(m[i-1]-m[i])*lum/(16.*pow(pi,2.)*a*c*pow(r[i-1],4.));
                double F = Freedman_Opacity(Z, P[i-1], A);
                return (pow(A,4.) - pow(B,4.))/(C*F) + 1.;
            };
            
            std::vector<double> gs;
            
            for (int x = 0; x < A_guess.size(); x++)
            {
                double A = A_guess[x];
                double g = g_func(A);
                gs.push_back(g);
            }
            
            double A_exact;
            bool all_same_sign = true;
            for (double element : gs)
            {
                if (sign_val(element) != sign_val(gs[0]))
                {
                    all_same_sign = false;
                    break;
                }
            }
            
            if (all_same_sign)
            {
                A_exact = A_adia;
            }
            else
            {
                std::vector<double> A_lims;
                for (int y = 0; y < gs.size()-1; y++)
                {
                    if (sign_val(gs[y]) != sign_val(gs[y+1]))
                    {
                        A_lims = {A_guess[y], A_guess[y+1]};
                    }
                }
                A_exact = brent_find_minima(g_func, A_lims[0], A_lims[1], double_bits).first;
            }
            T[i-1] = A_exact;
        }
        
//        if (abs(Grad_Rad) < abs(Grad_Adia))     // Radiative
//        {
//            double kR = k[-i];
//            double TR = T[i];
//            double kL = kR;                        // Initial guess
//            double TL = TR;                        // Initial guess
//            std::vector<double> kLs;
//            std::vector<double> TLs;
//            double k_avr = sqrt(kR * kL);
//            double e_k = 1.;
//            double e_T = 1.;
//            int n_i = 0;
//            std::vector<int> n_is;
//
//            while ((e_k > acc_k) or (e_T > acc_k))
//            {
//                n_i += 1;
//                n_is.push_back(n_i);
//
//                double kL_old = kL;
//                double TL_old = TL;
//                TL = pow(pow(TR,4.) - 3.*k_avr*lum*(m[i-1]-m[i])/(16.*a*c*pow(pi,2.)*pow(r[i-1],4.)), 1./4.);
//                kL = Freedman_Opacity(Z, P[i-1], TL);
//
//                kLs.push_back(kL);
//                TLs.push_back(TL);
//
//                e_k = abs(kL/kL_old - 1.);
//                e_T = abs(TL/TL_old - 1.);
//
//                k_avr = sqrt(kR * kL);
//            }
//
//            k[i-1] = kL;
//            T[i-1] = TL;
//
//            T[i-1] = pow(pow(T[i],4.) - 3.*k[i-1]*lum*(m[i-1]-m[i])/(16.*a*c*pow(pi,2.)*pow(r[i-1],4.)), 1./4.);
//            p[i-1] = P[i-1]/(R*T[i-1]);
//        }
        
        else                                    // Convective
        {
            rad_or_conv.push_back(0);
            T[i-1] = T[i] * pow(P[i-1]/P[i], (gamma-1.)/gamma);
        }
        
        k[i-1] = Freedman_Opacity(Z, P[i-1], T[i-1]);
        p[i-1] = P[i-1]/(R*T[i-1]);
        s[i-1] = T[i-1]/pow(P[i-1], (gamma-1.)/gamma);
        
        if (i > 2)
        {
            if (rad_or_conv[rad_or_conv.size()-1] - rad_or_conv[rad_or_conv.size()-2] == 1)
            {
                boundary_i.push_back({i,1});
            }
            else if (rad_or_conv[rad_or_conv.size()-1] - rad_or_conv[rad_or_conv.size()-2] == -1)
            {
                boundary_i.push_back({i,0});
            }
        }
    }
    
    variables = {m, r, p, T, k, P, s};
    
    return variables;
}


std::array<double, 2> Fre_Opacity_Giant_Planet::DensityFit(void)
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
