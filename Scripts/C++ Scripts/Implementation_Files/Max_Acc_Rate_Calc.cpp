#include "./../Header_Files/Max_Acc_Rate_Calc.hpp"

#include "./../Header_Files/Linspace.hpp"
#include "./../Header_Files/Sign.hpp"

#include <vector>
#include <string>
#include <iostream>





std::vector< std::vector<double> > Sigma_Calc(double Sigma_outer, double q, double H_r, double nu, int num_sig)
{
    const double G = 1.;                        // Gravitational constant
    const double AU = 1.495978707e10;           // Astronomical unit [cm]
    
    const double Rp = 1.;
    const double Ms = 1.;
    
    const double Omega_p = sqrt(G * Ms / pow(Rp, 3.));
    const double R_H = Rp * pow(q/3., 1./3.);
    
    std::vector<double> r = linspace(0.1, 3., int(num_sig));
    std::vector<double> Sigma(num_sig, 0.);
    
    Sigma[Sigma.size()-1] = Sigma_outer;
    
    double r_i, delta_i, Omega_i, a_dd_i, tg_i, t_f_nu_Omega_i, H_r_term_i, t_f_nu_Omega_r_i, Sigma_log;
    
    for (int i = 0; i < r.size(); i++)
    {
        r_i = r[r.size()-i];
        delta_i = r_i - Rp;
        
        Omega_i = sqrt(G * Ms / pow(r_i, 3.));
        a_dd_i = ((1. / 8.) * pow(abs(delta_i/R_H), -1.2)) + (200. * pow(abs(delta_i/R_H), -10.));
        
        tg_i = 0.35 * pow(q,2.) * pow(Rp,5.) * pow(Omega_p,2.) * (double)sign_val(delta_i) * r_i * pow(1./delta_i,4.);
        t_f_nu_Omega_i = 3. * nu * Omega_i / 4.;
        H_r_term_i = pow(H_r,2.) * R_H * Rp * pow(Omega_p,2.) * r_i * a_dd_i;
        t_f_nu_Omega_r_i = 3. * nu * Omega_i * r_i / 2.;
        
        Sigma_log = log(Sigma[Sigma.size()-i]) + (((tg_i - t_f_nu_Omega_i) / (H_r_term_i + t_f_nu_Omega_r_i)) * (r[r.size()-(i+1)] - r_i));
        Sigma[Sigma.size()-(i+1)] = exp(Sigma_log);
    }
    
    for (int i = 0; i < r.size(); i++)
    {
        r[i] *= 5.2 * AU;
        Sigma[i] *= 800;
    }
    
    std::vector< std::vector<double> > Sigma_Vals = {Sigma, r};
    
    return Sigma_Vals;
}



std::pair< std::vector< std::vector<double> >, std::vector<bool>> Max_Acc_Rate(double Sigma_outer, double Ms, double Mp, double Tp, double Rp, double r_in, double r_out, double alpha, int num_sig, bool Const_Sigma, bool Crida_Comparison, bool L09)
{
    const double G = 6.67430e-8;                // Gravitational constant
    const double mu = 2.3;                      // Relative mass
    const double R = 8.314e7/mu;                // Gas constant [erg/mol/K]
    const double Ma_S = 1.989e33;               // Sun mass [g]
    const double AU = 1.495978707e13;           // Astronomical unit [cm]
    
    Ms *= Ma_S;
    Rp *= AU;
    
    if (L09)
    {
        double cs2 = R * Tp;                           // Sound speed [cm/s]
        double R_H = (5.2 * AU) * pow(Mp/(3.*Ma_S), 1./3.);
        r_in = 1.e-10;
        r_out = (G * Mp / ((cs2 / 1.) + ((G * Mp) / (0.25 * R_H)))) / R_H;
    }
    
    double H_r, nu;
    
    if (Crida_Comparison)
    {
        H_r = 0.05;
        nu = 3.e14;
        Mp = 1.e-3 * Ms;
        Sigma_outer = 1. / sqrt(3.);
//        std::cout << "H/r = " << H_r << " and \u03BD = " << nu/1.e14 << std::endl;
    }
    else
    {
        double cs_p = sqrt(R * Tp);
        double H_p = cs_p * sqrt(pow(Rp,3.) / (G * Ms));
        H_r = H_p / Rp;
        nu = alpha * cs_p * H_p;
//        std::cout << "H/r = " << H_r << " and \u03BD = " << nu/1.e14 << std::endl;
    }
    
    double Omega_p = sqrt(G * Ms / pow(Rp, 3.));
    double R_H = Rp * pow(Mp/(3 * Ms), 1./3.);
    
    std::vector<double> r_d;
    
    if (Mp == 0.)
    {
        r_d = linspace(0.5*AU, 3.*Rp, num_sig);
    }
    else
    {
        std::vector<double> r_d1 = linspace(0.5*AU, Rp-(r_out*R_H), int(num_sig/5)+1);
        r_d1.pop_back();
        std::vector<double> r_d2 = linspace(Rp-(r_out*R_H), Rp-(r_in*R_H), int(num_sig/5));
        std::vector<double> r_d3 = linspace(Rp-(r_in*R_H), Rp+(r_in*R_H), int(num_sig/5)+2);
        r_d3.erase(r_d3.begin());
        r_d3.pop_back();
        std::vector<double> r_d4 = linspace(Rp+(r_in*R_H), Rp+(r_out*R_H), int(num_sig/5));
        std::vector<double> r_d5 = linspace(Rp+(r_out*R_H), 3*Rp, int(num_sig/5)+1);
        r_d5.erase(r_d5.begin());
        r_d1.insert(r_d1.end(), r_d2.begin(), r_d2.end());
        r_d1.insert(r_d1.end(), r_d3.begin(), r_d3.end());
        r_d1.insert(r_d1.end(), r_d4.begin(), r_d4.end());
        r_d1.insert(r_d1.end(), r_d5.begin(), r_d5.end());
        r_d = r_d1;
    }
    
    std::vector<double> Mdot_maxs;
    std::vector<double> Sigma(r_d.size(), 0.);
    std::vector<double> Mdot_max_ana;
    
    if (Const_Sigma)
    {
        for (int i = 0; i < Sigma.size(); i++)
        {
            Sigma[i] = Sigma_outer;
        }
        double Mdot_max_ana_in = 2. * Sigma_outer * Omega_p * pow(Rp,2.) * (sqrt(1. - (r_in * R_H / Rp)) - sqrt(1. - (r_out * R_H / Rp)));
        double Mdot_max_ana_out = 2. * Sigma_outer * Omega_p * pow(Rp,2.) * (sqrt(1. + (r_out * R_H / Rp)) - sqrt(1. + (r_in * R_H / Rp)));
        Mdot_max_ana.push_back(Mdot_max_ana_in + Mdot_max_ana_out);
    }
    else
    {
        Sigma[Sigma.size()-1] = Sigma_outer;
        double rd_i, delta_i, Omega_i, a_dd_i, tg_i, t_f_nu_Omega_i, H_r_term_i, t_f_nu_Omega_r_i, Sigma_log;
        for (int i = 1; i < r_d.size(); i++)
        {
            rd_i = r_d[r_d.size()-i];
            delta_i = rd_i - Rp;
            
            Omega_i = sqrt(G * Ms / pow(rd_i, 3.));
            a_dd_i = ((1. / 8.) * pow(abs(delta_i/R_H), -1.2)) + (200. * pow(abs(delta_i/R_H), -10.));
            
            tg_i = 0.35 * pow((Mp/Ms),2.) * pow(Rp,5.) * pow(Omega_p,2.) * (double)sign_val(delta_i) * rd_i * pow(1./delta_i,4.);
            t_f_nu_Omega_i = 3. * nu * Omega_i / 4.;
            H_r_term_i = pow(H_r,2.) * R_H * Rp * pow(Omega_p,2.) * rd_i * a_dd_i;
            t_f_nu_Omega_r_i = 3. * nu * Omega_i * rd_i / 2.;
            
            Sigma_log = log(Sigma[Sigma.size()-i]) + (((tg_i - t_f_nu_Omega_i) / (H_r_term_i + t_f_nu_Omega_r_i)) * (r_d[r_d.size()-(i+1)] - rd_i));
            Sigma[Sigma.size()-(i+1)] = exp(Sigma_log);
        }
    }
    
    for (int i = int(num_sig/5)-1; i < int(2*num_sig/5)-1; i++)
    {
        Mdot_maxs.push_back(Sigma[i] * sqrt(G * Ms / r_d[i]) * (r_d[i+1] - r_d[i]));
    }
    for (int i = int(3*num_sig/5)-2; i < int(4*num_sig/5)-2; i++)
    {
        Mdot_maxs.push_back(Sigma[i] * sqrt(G * Ms / r_d[i]) * (r_d[i+1] - r_d[i]));
    }
    
    double Mdot_max_sum = 0;
    for (int i = 0; i < Mdot_maxs.size(); i++)
    {
        Mdot_max_sum += Mdot_maxs[i];
    }
    std::vector<double> Mdot_max = {Mdot_max_sum};
    std::vector<double> extras = {Mp, Rp, (double)num_sig};
    
    std::vector< std::vector<double> > Max_Acc_Rate_Vals;
    
    if (Const_Sigma)
    {
        Max_Acc_Rate_Vals = {Mdot_max, Mdot_max_ana, extras};
    }
    else if (Crida_Comparison)
    {
        Max_Acc_Rate_Vals = {Sigma, r_d, extras};
    }
    else
    {
        Max_Acc_Rate_Vals = {Mdot_max, Sigma, r_d, extras};
    }
    
    std::pair< std::vector< std::vector<double> >, std::vector<bool> > Max_Acc_Rate_Result = {Max_Acc_Rate_Vals, {Const_Sigma, Crida_Comparison}};
    
    return Max_Acc_Rate_Result;
}
