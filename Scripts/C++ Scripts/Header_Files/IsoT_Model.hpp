#ifndef IsoT_Model_hpp
#define IsoT_Model_hpp

#include <vector>
#include <array>

class Iso_Temperature_Giant_Planet
{
public:
    Iso_Temperature_Giant_Planet(double mc, double me, double rp, double pd, double Td, int num_p, bool grid_given, std::vector<double> m_grid, std::array<double, 2> p_ana_c);
    ~Iso_Temperature_Giant_Planet();
    std::vector< std::vector<double> > getVariables(void);
    void printVariable(std::vector<double> variable);
    std::vector< std::vector<double> > Planet_Formation(void);
    std::vector<double> L2norm(void);
    std::array<double, 2> DensityFit(void);

private:
    // Inputs:
    const double mc, me, rp, pd, Td;                // Core mass, Envelope mass, Planetary radius, Disc density, Disc temperature
    const int num_p;                                // Number of points along mass grid m
    const bool grid_given;                          // Is a mass grid given?
    std::vector<double> m_grid;                     // Custom mass grid
    std::array<double, 2> p_ana_c;                  // Analytic density fit parameters
    
    // Constants:
    const int num_g = 1;                        // Number of ghost cells on each side
    const double AU = 1.495978707e13;           // Astronomical unit [cm]
    const double G = 6.67430e-8;                // Gravitational constant
    const double kB = 1.380649e-16;             // Boltzmann constant [erg/K]
    const double mu = 2.3;                      // Relative mass
    const double R = 8.314e7/mu;                // Gas constant [erg/mol/K]
    const double gamma = 7./5.;                 // Adiabatic constant
    const double mp = 1.67262192e-24;           // Proton mass [g]
    const double R_J = 6.9911e9;                // Jupiter radius [cm]
    const double Ma_E = 5.972e27;               // Earth mass [g]
    const double Ma_S = 1.989e33;               // Sun mass [g]
    const double pi = 3.1415926535;             // Pi
    const int lo = num_g;                       // Lower boundary for array iteration
    const int hi = num_g + num_p - 1;           // Upper boundary for array iteration
    const int double_bits = std::numeric_limits<double>::digits;
    
    // Declarations:
    double core_mass_max;       // Planet's core mass [g]
    double env_mass_max;        // Planet's envelope mass [g]
    double env_mass_min;        // [g]
    double Mp;                  // Total planet mass [g]
    double R_H;                 // Hill Radius [cm]
    double cs2;                 // Sound speed squared [cm^2/s^2]
    double rad_hi;              // Planet radius [cm]
    double den_hi;              // Disc density [g/cm^3]
    double tem_hi;              // Disc temperature [K]
    double pre_hi;              // Disc pressure [g/cm/s^2]
    std::vector<double> m, r, p, T, P, s, p_ana, dP_dr_num, dP_dr_ana;
    std::vector< std::vector<double> > variables;
    std::vector<double> L2norm_vars;
};

#endif
