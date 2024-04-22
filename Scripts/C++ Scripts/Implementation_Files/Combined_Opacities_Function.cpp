#include "./../Header_Files/Combined_Opacities_Function.hpp"

#include "./../Header_Files/Linspace.hpp"
#include "./../Header_Files/Freedman_Opacities_Function.hpp"

#include <vector>





//-----------------------------------------------------------------------------
//
//
// Mean opacity - Rosseland
//
//
//-----------------------------------------------------------------------------

const double h = 6.626196e-27;
const double c = 2.998e10;
const double k = 1.380649e-16;

double planck_a = 2. * h * pow(c,2.);
double planck_b = h * c / k;

//-----------------------------------------------------------------------------
// Planck B(T) = (2hc^2/l^5) / (e^(hc/lkT) - 1)
// dB(T)/dT = (2h^2c^3/l^6kT^2) * (e^(hc/lkT)/((e^(hc/lkT)-1)^2))
// wav: wavelength in cm
// T: temperature in K
// returns: Energy dens in erg cm^-3 s^-1 sr^-1 K^-1

double Planck_Deriv(double l, double T)
{
    double dB_dT = (2.*pow(h,2.)*pow(c,3.)/(pow(l,6.)*k*pow(T,2.))) * (1. / ((exp(h*c/(l*k*T))-1.) * (1.-exp(-h*c/(l*k*T)))));
    return dB_dT;
}



//-----------------------------------------------------------------------------
// Dust effective absorption function
// x: wavelength in cm
// returns: Exctinction efficiency [1], see e.g. Mordasini 2014b (Analytical dust model paper) and refs therein

double Q(double x)
{
    if (x < 1.e-3)
    {
        return 2.;
    }
    else if (x < 0.5)
    {
        return 2. + 4. * x;
    }
    else if (x > 8./3.)
    {
        return 0.3 * pow(x,-1.);
    }
    else
    {
        return 0.8 * pow(x,-2.);
    }
}



//-----------------------------------------------------------------------------
// Wrapper function
// l: wavelength in cm
// a: dust size  in cm
// output: opacity in cm^2 per gram of dust at all given x

const double pi = 3.1415926535;

std::vector<double> dust_opa(std::vector<double> ls, double a)
{
    std::vector<double> xs;
    for (int i = 0; i < ls.size(); i++)
    {
        xs.push_back(ls[i] / (2. * pi * a));
    }
    std::vector<double> d_o;
    for (double x : xs)
    {
        d_o.push_back(Q(x)/a);
    }
    return d_o;
}



////-----------------------------------------------------------------------------
//// Dust effective absorption function
//// x: 2 pi a / l
//// returns: Exctinction efficiency [1], see e.g. Mordasini 2014b eqn 54 (Analytical dust model paper) and refs therein
//
//double Q(double x)
//{
//    if (x < 0.375)
//    {
//        return 0.3 * x;
//    }
//    else if (0.375 <= x < 2.188)
//    {
//        return 0.8 * pow(x,2.);
//    }
//    else if (2.188 <= x < 1000)
//    {
//        return 2. + 4. / x;
//    }
//    else
//    {
//        return 2;
//    }
//}
//
//
//
////-----------------------------------------------------------------------------
//// Wrapper function
//// l: wavelength in cm
//// a: dust size  in cm
//// output: opacity in cm^2/g of dust at all given x
//
//std::vector<double> dust_opa(std::vector<double> l, double a)
//{
//    std::vector<double> xs;
//    for (int i = 0; i < ls.size(); i++)
//    {
//        xs.push_back(2. * pi * a / l[i]);
//    }
//    std::vector<double> d_o;
//    for (double x : xs)
//    {
//        d_o.push_back(Q(x)/a);
//    }
//    return d_o;
//}



//-----------------------------------------------------------------------------
// On the fly Rosseland mean opacity for dust of a given size a and temperature field T
// T: temperature [K]
// a: particle size [cm]
// returns: Rosseland mean opacity, cm^2/g

double Dust_Ross(double T, double a)
{
    std::vector<double> lin_arr = linspace(-4., 4., 1000);
    std::vector<double> ls;
    for (int i = 0; i < lin_arr.size(); i++)
    {
        ls.push_back(pow(10., lin_arr[i])*1.e-4);               // mum to cm
    }
    ls.insert(ls.begin(), 0.);
    std::vector<double> dls;
    for (int i = 1; i < ls.size(); i++)
    {
        dls.push_back(ls[i]-ls[i-1]);
    }
    ls.erase(ls.begin());
    
    std::vector<double> k_dust = dust_opa(ls, a);
    double norm = 0.;
    double dust_norm = 0.;
    for (int i = 0; i < ls.size(); i++)
    {
        norm += Planck_Deriv(ls[i],T) * dls[i];
        dust_norm += (1/k_dust[i]) * Planck_Deriv(ls[i],T) * dls[i];
    }
    
    if (T < 1500.)
    {
        return norm / dust_norm;
    }
    else
    {
        return 0.;                                              // Due to evaporation
    }
}



//-----------------------------------------------------------------------------
//
//
// Combined Rosseland mean opacity for dust and gas under a Rosseland addition rule assumption
//
//
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// P: gas pressure [dyne/cm^2], ref value: 1e6=1bar
// T: temperature [K]
// a: particle size [cm]
// dust_to_gas: mass ratio of dust to gas, sets metallicity for gas opacities as well
// returns: cm^2/g

double Combined_Opacity(double ad, double d_to_g, double P, double T)
{
    double gas_opacity = Freedman_Opacity(P, T, log10(d_to_g/1.e-2));
    double dust_opacity = Dust_Ross(T,ad) * d_to_g;
    double opacity = gas_opacity + dust_opacity;
    
    return opacity;
}
