#include "./Header_Files/Linspace.hpp"
#include "./Header_Files/Max_Acc_Rate_Calc.hpp"
#include "./Header_Files/IsoT_Model.hpp"
#include "./Header_Files/IsoD_Model.hpp"
#include "./Header_Files/IsoO_Model.hpp"
#include "./Header_Files/FreO_Model.hpp"
#include "./Header_Files/ComO_Model.hpp"
#include "./Header_Files/Mass_Grid_Solver.hpp"
#include "./Header_Files/Convergence_Tester.hpp"
#include "./Header_Files/Lum_Finder.hpp"
#include "./Header_Files/Evolution.hpp"
#include "./Header_Files/Data_Output.hpp"

#include <iostream>
#include <vector>
#include <string>





int main(void)
{
    std::string Run_Type = "IsoT";
    
//--------------------------------------------------------------------------------------------------------------------------------------------------
    
    if (Run_Type == "IsoT")
    {
        // Initialising test values to show the code works (not necessarily realisitic)
        //                                mc   me    rp       pd        Td    num_p acc_m
        std::array<double, 7> val_IsoT = {10., 0.01, 41.4089, 1.47e-11, 118., 1.e4, 1.e-4};
        auto [mc_IsoT, me_IsoT, rp_IsoT, pd_IsoT, Td_IsoT, n_IsoT, acc_m] = val_IsoT;
        std::array<double, 2> p_ana_c_ini_IsoT = {3.93e-11, 1.48e1};    // initial density fit parameters
        
//        // Create IsoT1 object without a mass grid
//        std::vector<double> m_grid_IsoT;
//        Iso_Temperature_Giant_Planet IsoT1(mc_IsoT, me_IsoT, rp_IsoT, pd_IsoT, Td_IsoT, n_IsoT, false, m_grid_IsoT, p_ana_c_ini_IsoT);
        
        // Calculating mass grid where density fit parameters don't change more than acc_m
        auto returned_values_IsoT = Mass_Grid_Solver_IsoT(mc_IsoT, me_IsoT, rp_IsoT, pd_IsoT, Td_IsoT, n_IsoT, acc_m, p_ana_c_ini_IsoT);
        std::vector<double> m_grid_IsoT = returned_values_IsoT.first;           // mass grid
        std::array<double, 2> p_ana_c_IsoT = returned_values_IsoT.second;       // updated density fit parameters
        
        // Create IsoT1 object with a mass grid
        Iso_Temperature_Giant_Planet IsoT1(mc_IsoT, me_IsoT, rp_IsoT, pd_IsoT, Td_IsoT, n_IsoT, true, m_grid_IsoT, p_ana_c_IsoT);
        
        // Forming the envelope from outside in
        std::vector< std::vector<double> > IsoT1_Sim_Data = IsoT1.Planet_Formation();
        // Checking we are in a region which is independent of the boundary condition
        std::cout << "Max p/p_Disc = " << *std::max_element(IsoT1_Sim_Data[2].begin(),IsoT1_Sim_Data[2].end())/pd_IsoT << std::endl;
        // Reading out all variables of envelope into .txt files to be plotted
        Read_Out(IsoT1_Sim_Data, "IsoT");
        
        // Concatenating to create a large array of of num_p values
        acc_m = 4.e-4;      // Need a larger acc_m for small num_p values
        std::vector<double> n_array_IsoT = linspace(1.e2, 1.e3, 30);
        std::vector<double> n_array_IsoT2 = linspace(1.e3, 1.e4, 50);
        std::vector<double> n_array_IsoT3 = linspace(1.e4, 1.e5, 30);
        n_array_IsoT.insert(n_array_IsoT.end(), n_array_IsoT2.begin()+1, n_array_IsoT2.end());
        n_array_IsoT.insert(n_array_IsoT.end(), n_array_IsoT3.begin()+1, n_array_IsoT3.end());
        // Performing convergence test and reading out the L2norm data to be plotted
        std::vector< std::vector<double> > IsoT1_L2norm_Data = ConvergenceIsoT(mc_IsoT, me_IsoT, rp_IsoT, pd_IsoT, Td_IsoT, n_array_IsoT, acc_m, p_ana_c_ini_IsoT);
        Read_Out(IsoT1_L2norm_Data, "ConvergenceIsoT");
    }
    
//----------------------------------------------------------------------------------------------------------------------------------------------------
    
    else if (Run_Type == "IsoD")
    {
        // Initialising test values to show the code works (not necessarily realisitic)
        //                                mc   me      rp       pd      Td    kd      lp     num_p
        std::array<double, 8> val_IsoD = {10., 0.0008, 36.3776, 7.e-11, 118., 1.5e-5, 10000, 1.e4};
        auto [mc_IsoD, me_IsoD, rp_IsoD, pd_IsoD, Td_IsoD, kd_IsoD, lp_IsoD, n_IsoD] = val_IsoD;
        
        // Create IsoD1 object without a mass grid
        Iso_Density_Giant_Planet IsoD1(mc_IsoD, me_IsoD, rp_IsoD, pd_IsoD, Td_IsoD, kd_IsoD, lp_IsoD, n_IsoD);
        
        // Forming the envelope from outside in
        std::vector< std::vector<double> > IsoD1_Sim_Data = IsoD1.Planet_Formation();
        // Checking we are in a region which is independent of the boundary condition
        std::cout << "Max T/T_Disc = " << *std::max_element(IsoD1_Sim_Data[3].begin(),IsoD1_Sim_Data[3].end())/pd_IsoD << std::endl;
        // Reading out all variables of envelope into .txt files to be plotted
        Read_Out(IsoD1_Sim_Data, "IsoD");
        
        // Concatenating to create a large array of of num_p values
        std::vector<double> n_array_IsoD = linspace(1.e2, 1.e3, 30);
        std::vector<double> n_array_IsoD2 = linspace(1.e3, 1.e4, 50);
        std::vector<double> n_array_IsoD3 = linspace(1.e4, 1.e5, 30);
        n_array_IsoD.insert(n_array_IsoD.end(), n_array_IsoD2.begin()+1, n_array_IsoD2.end());
        n_array_IsoD.insert(n_array_IsoD.end(), n_array_IsoD3.begin()+1, n_array_IsoD3.end());
        // Performing convergence test and reading out the L2norm data to be plotted
        std::vector< std::vector<double> > IsoD1_L2norm_Data = ConvergenceIsoD(mc_IsoD, me_IsoD, rp_IsoD, pd_IsoD, Td_IsoD, kd_IsoD, lp_IsoD, n_array_IsoD);
        Read_Out(IsoD1_L2norm_Data, "ConvergenceIsoD");
    }
    
//----------------------------------------------------------------------------------------------------------------------------------------------------
    
    else if (Run_Type == "IsoO")
    {
        // Initialising test values to show the code works (not necessarily realisitic)
        //                                mc   me  rp       pd        Td    kd     lp                    num_p acc_m
        std::array<double, 9> val_IsoO = {10., 1., 43.5589, 1.47e-11, 118., 6.e-5, 0.715809160364261e-4, 1.e4, 5.e-3};
        auto [mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, Td_IsoO, kd_IsoO, lp_IsoO, n_IsoO, acc_m] = val_IsoO;
        std::array<double, 2> p_ana_c_ini_IsoO = {6.14e-11, 9.20};    // initial density fit parameters
        
//        // Create IsoO1 object without a mass grid
//        std::vector<double> m_grid_IsoO;
//        Iso_Opacity_Giant_Planet IsoO1(mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, Td_IsoO, kd_IsoO, lp_IsoO, n_IsoO, false, m_grid_IsoO, p_ana_c_ini_IsoO);
        
        // Calculating mass grid where density fit parameters don't change more than acc_m
        auto returned_values_IsoO = Mass_Grid_Solver_IsoO(mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, Td_IsoO, kd_IsoO, lp_IsoO, n_IsoO, acc_m, p_ana_c_ini_IsoO);
        std::vector<double> m_grid_IsoO = returned_values_IsoO.first;         // mass grid
        std::array<double, 2> p_ana_c_IsoO = returned_values_IsoO.second;     // updated density fit parameters
        
        // Create IsoO1 object with a mass grid
        Iso_Opacity_Giant_Planet IsoO1(mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, Td_IsoO, kd_IsoO, lp_IsoO, n_IsoO, true, m_grid_IsoO, p_ana_c_IsoO);
        
        // Forming the envelope from outside in
        std::vector< std::vector<double> > IsoO1_Sim_Data = IsoO1.Planet_Formation("Both");
        // Reading out all variables of envelope into .txt files to be plotted
        Read_Out(IsoO1_Sim_Data, "IsoO");
        Read_Out(IsoO1.getInfo().second, "IsoO");
        
        // Concatenating to create a large array of of num_p values
        acc_m = 6.e-3;      // Need a larger acc_m for small num_p values
        std::vector<double> n_array_IsoO = geomspace(1.e2, 1.e5, 4, true);
        // Performing convergence test and reading out the L2norm data to be plotted
        std::vector< std::vector<double> > IsoO1_L2norm_Data = ConvergenceIsoO(mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, Td_IsoO, kd_IsoO, lp_IsoO, n_array_IsoO, acc_m, p_ana_c_ini_IsoO);
        Read_Out(IsoO1_L2norm_Data, "ConvergenceIsoO");
    }
    
//----------------------------------------------------------------------------------------------------------------------------------------------------
    
    else if (Run_Type == "FreO")
    {
        // Initialising test values to show the code works (not necessarily realisitic)
        //                                 mc   me  rp       pd        Td    kd  lp      met acc_k   num_p acc_m
        std::array<double, 11> val_FreO = {10., 1., 43.5589, 1.47e-11, 118., 1., 1.e-10, 0., 1.e-10, 1.e3, 5.e-4};
        auto [mc_FreO, me_FreO, rp_FreO, pd_FreO, Td_FreO, kd_FreO, lp_FreO, met, acc_k, n_FreO, acc_m] = val_FreO;
        std::array<double, 2> p_ana_c_ini_FreO = {5.94e-11, 9.21};      // initial density fit parameters
        
//        // Create FreO1 object without a mass grid
//        std::vector<double> m_grid_FreO;
//        Fre_Opacity_Giant_Planet FreO1(mc_FreO, me_FreO, rp_FreO, pd_FreO, Td_FreO, kd_FreO, lp_FreO, met, acc_k, n_FreO, false, m_grid_FreO, p_ana_c_ini_FreO);
        
        // Calculating mass grid where density fit parameters don't change more than acc_m
        auto returned_values_FreO = Mass_Grid_Solver_FreO(mc_FreO, me_FreO, rp_FreO, pd_FreO, Td_FreO, kd_FreO, lp_FreO, met, acc_k, n_FreO, acc_m, p_ana_c_ini_FreO);
        std::vector<double> m_grid_FreO = returned_values_FreO.first;         // mass grid
        std::array<double, 2> p_ana_c_FreO = returned_values_FreO.second;     // updated density fit parameters
        
        // Create FreO1 object with a mass grid
        Fre_Opacity_Giant_Planet FreO1(mc_FreO, me_FreO, rp_FreO, pd_FreO, Td_FreO, kd_FreO, lp_FreO, met, acc_k, n_FreO, true, m_grid_FreO, p_ana_c_FreO);
        
        // Forming the envelope from outside in
        std::vector< std::vector<double> > FreO1_Sim_Data = FreO1.Planet_Formation();
        // Reading out all variables of envelope into .txt files to be plotted
        Read_Out(FreO1_Sim_Data, "FreO");
        Read_Out(FreO1.getInfo().first.second, "FreO");
        
        // Concatenating to create a large array of of num_p values
//        acc_m = 1.e-6;      // Need a larger acc_m for small num_p values
        std::vector<double> n_array_FreO = geomspace(1.e3, 1.e4, 10, true);
        // Performing convergence test and reading out the L2norm data to be plotted
        std::vector< std::vector<double> > FreO1_L2norm_Data = ConvergenceFreO(mc_FreO, me_FreO, rp_FreO, pd_FreO, Td_FreO, kd_FreO, lp_FreO, met, acc_k, n_array_FreO, acc_m, p_ana_c_ini_FreO);
        Read_Out(FreO1_L2norm_Data, "ConvergenceFreO");
    }
    
//----------------------------------------------------------------------------------------------------------------------------------------------------
    
    else if (Run_Type == "ComO")
    {
        // Initialising test values to show the code works (not necessarily realisitic)
        //                                 mc   me  rp       pd        Td    kd  lp      ad    d_to_g acc_k   num_p acc_m
        std::array<double, 12> val_ComO = {10., 1., 43.5589, 1.47e-11, 118., 1., 1.e-10, 1.e-2, 0.01, 1.e-10, 1.e3, 1.e-3};
        auto [mc_ComO, me_ComO, rp_ComO, pd_ComO, Td_ComO, kd_ComO, lp_ComO, ad, d_to_g, acc_k, n_ComO, acc_m] = val_ComO;
        std::array<double, 2> p_ana_c_ini_ComO = {1.00e-11, 1.00};      // initial density fit parameters   [1.52e-14, 2.50e1]
        
        // Create ComO1 object without a mass grid
        std::vector<double> m_grid_ComO;
        Com_Opacity_Giant_Planet ComO1(mc_ComO, me_ComO, rp_ComO, pd_ComO, Td_ComO, kd_ComO, lp_ComO, ad, d_to_g, acc_k, n_ComO, false, m_grid_ComO, p_ana_c_ini_ComO);
        
//        // Calculating mass grid where density fit parameters don't change more than acc_m
//        auto returned_values_ComO = Mass_Grid_Solver_ComO(mc_ComO, me_ComO, rp_ComO, pd_ComO, Td_ComO, kd_ComO, lp_ComO, ad, d_to_g, acc_k, n_ComO, acc_m, p_ana_c_ini_ComO);
//        std::vector<double> m_grid_ComO = returned_values_ComO.first;         // mass grid
//        std::array<double, 2> p_ana_c_ComO = returned_values_ComO.second;     // updated density fit parameters
//        
//        // Create ComO1 object with a mass grid
//        Com_Opacity_Giant_Planet ComO1(mc_ComO, me_ComO, rp_ComO, pd_ComO, Td_ComO, kd_ComO, lp_ComO, ad, d_to_g, acc_k, n_ComO, true, m_grid_ComO, p_ana_c_ComO);
        
        // Forming the envelope from outside in
        std::vector< std::vector<double> > ComO1_Sim_Data = ComO1.Planet_Formation();
        // Reading out all variables of envelope into .txt files to be plotted
        Read_Out(ComO1_Sim_Data, "ComO");
        Read_Out(ComO1.getInfo().first.second, "ComO");
        std::cout << ComO1.DensityFit()[0] << ", " << ComO1.DensityFit()[1] << std::endl;
        
//        // Concatenating to create a large array of of num_p values
//        acc_m = 1.e-6;      // Need a larger acc_m for small num_p values
//        std::vector<double> n_array_ComO = geomspace(1.e3, 1.e4, 10, true);
//        // Performing convergence test and reading out the L2norm data to be plotted
//        std::vector< std::vector<double> > ComO1_L2norm_Data = ConvergenceComO(mc_ComO, me_ComO, rp_ComO, pd_ComO, Td_ComO, kd_ComO, lp_ComO, ad, d_to_g, acc_k, n_array_ComO, acc_m, p_ana_c_ini_ComO);
//        Read_Out(ComO1_L2norm_Data, "ConvergenceComO");
    }
    
//----------------------------------------------------------------------------------------------------------------------------------------------------
    
    else if (Run_Type == "Max_Acc")
    {
        bool Zero_mass = false;
        bool Crida_Sigma_14 = false;
        bool Diff_Alpha = false;
        
        // Initialising test values to show the code works (not necessarily realisitic)
        //                                  Sigma_outer       Ms  Mp            Tp    Rp r_in r_out alpha num_sig
        std::array<double, 9> val_Sigma = {143.365/sqrt(3.), 1., 11.*5.972e27, 118., 5.2, 2., 3., 1.e-5, 1e4};
        auto [Val_outer_Sigma, Ms_Sigma, Mp_Sigma, Tp_Sigma, Rp_Sigma, r_in_Sigma, r_out_Sigma, alpha_Sigma, num_sig_Sigma] = val_Sigma;
        std::pair< std::vector< std::vector<double> >, std::vector<bool> > Max_Acc_Rate_Result;
        
        if (Zero_mass)
        {
            // Testing 0 mass planet
            // Sigma_outer, Ms, Mp, Tp, Rp, r_in, r_out, alpha, num_sig, Const_Sigma, Crida_Comparison, L09
            Max_Acc_Rate_Result = Max_Acc_Rate(Val_outer_Sigma, Ms_Sigma, 0., Tp_Sigma, Rp_Sigma, r_in_Sigma, r_out_Sigma, alpha_Sigma, num_sig_Sigma, false, false, false);
            Read_Out(Max_Acc_Rate_Result, "Max_Acc_Rate");
            // No gap forms in disc as expected and -> perfect 1/sqrt(r) dependence in Sigma
        }
        
        else if (Crida_Sigma_14)
        {
            // Testing Crida Sigma eqn 14 Sigma solution - comparing to Figure 11
            // Crida Comparison -> Hard sets: Sigma_outer = 1/np.sqrt(3), nu = 3.e14, H/r = 0.05, Mp = 1.e-3 * Ms
            // Sigma_outer, Ms, Mp, Tp, Rp, r_in, r_out, alpha, num_sig, Const_Sigma, Crida_Comparison, L09
            Max_Acc_Rate_Result = Max_Acc_Rate(1., 1., 1., 1., 5.2, 2., 3., 1., 1e5, false, true, false);
            Read_Out(Max_Acc_Rate_Result, "Max_Acc_Rate");
        }
        
        else if (Diff_Alpha)
        {
            std::array<double, 2> alphas = {1.e-2, 1.e-5};
            std::pair< std::vector< std::vector<double> >, std::vector<bool> > Max_Acc_Rate_Result, Max_Acc_Rate_Result0, Max_Acc_Rate_Result1, Max_Acc_Rate_Result2, Max_Acc_Rate_Result3;
            std::vector< std::pair< std::vector< std::vector<double> >, std::vector<bool> > > Max_Acc_Rate_Results;

            // Sigma_outer, Ms, Mp, Tp, Rp, r_in, r_out, alpha, num_sig, Const_Sigma, Crida_Comparison, L09
            Max_Acc_Rate_Result0 = Max_Acc_Rate(Val_outer_Sigma, Ms_Sigma, 11.*5.972e27, Tp_Sigma, Rp_Sigma, r_in_Sigma, r_out_Sigma, alphas[0], num_sig_Sigma, false, false, false);
            Max_Acc_Rate_Result1 = Max_Acc_Rate(Val_outer_Sigma, Ms_Sigma, 11.*5.972e27, Tp_Sigma, Rp_Sigma, r_in_Sigma, r_out_Sigma, alphas[1], num_sig_Sigma, false, false, false);
            Max_Acc_Rate_Result2 = Max_Acc_Rate(Val_outer_Sigma, Ms_Sigma, 300.*5.972e27, Tp_Sigma, Rp_Sigma, r_in_Sigma, r_out_Sigma, alphas[0], num_sig_Sigma, false, false, false);
            Max_Acc_Rate_Result3 = Max_Acc_Rate(Val_outer_Sigma, Ms_Sigma, 300.*5.972e27, Tp_Sigma, Rp_Sigma, 1., 1., alphas[1], num_sig_Sigma, false, false, true);
            
            Max_Acc_Rate_Results = {Max_Acc_Rate_Result0, Max_Acc_Rate_Result1, Max_Acc_Rate_Result2, Max_Acc_Rate_Result3};
            Read_Out(Max_Acc_Rate_Results, "Max_Acc_Rate_Diff_Alpha");
        }
        
        // Concatenating to create a large array of of num_p values
        std::vector<double> n_array_Sigma = geomspace(1.e2, 1.e5, 4, true);
        // Performing convergence test and reading out the L2norm data to be plotted
        std::vector< std::vector<double> > Sigma_L2norm_Data = ConvergenceSigma(Val_outer_Sigma, Ms_Sigma, Mp_Sigma, Tp_Sigma, Rp_Sigma, r_in_Sigma, r_out_Sigma, alpha_Sigma, n_array_Sigma);
        Read_Out(Sigma_L2norm_Data, "Convergence_Sigma");
        
        // Performing convergence test and reading out the L2norm data to be plotted
        std::vector< std::vector<double> > MdotMax_L2norm_Data = ConvergenceMdotMax(Val_outer_Sigma, Ms_Sigma, Mp_Sigma, Tp_Sigma, Rp_Sigma, r_in_Sigma, r_out_Sigma, alpha_Sigma, n_array_Sigma);
        Read_Out(MdotMax_L2norm_Data, "Convergence_MdotMax");
    }
    
//----------------------------------------------------------------------------------------------------------------------------------------------------
    
    else if (Run_Type == "IsoO_LF")
    {
        // Initialising test values to show the code works (not necessarily realisitic)
        //                                 mc   me  rp       pd        Td    kd     lp     acc_l   num_p acc_m
        std::array<double, 10> val_IsoO = {10., 1., 43.5589, 1.47e-11, 118., 6.e-5, 1.e-4, 1.e-10, 1.e3, 1.e-3};
        auto [mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, Td_IsoO, kd_IsoO, lp_IsoO, acc_l, n_IsoO, acc_m] = val_IsoO;
        std::array<double, 2> p_ana_c_ini_IsoO = {6.14e-11, 9.20};    // initial density fit parameters
        
        // Calculating mass grid where density fit parameters don't change more than acc_m
        auto returned_values_IsoO = Mass_Grid_Solver_IsoO(mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, Td_IsoO, kd_IsoO, lp_IsoO, n_IsoO, acc_m, p_ana_c_ini_IsoO);
        std::vector<double> m_grid_IsoO = returned_values_IsoO.first;         // mass grid
        
        std::pair< std::vector< std::vector<double> >, std::vector< std::vector<double> > > IsoO_LF_Data = L_Finder_IsoO(mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, Td_IsoO, kd_IsoO, lp_IsoO, acc_l, n_IsoO, m_grid_IsoO);
        Read_Out(IsoO_LF_Data, "IsoO_LF");
        
        // Performing convergence test and reading out the L2norm data to be plotted
        std::vector< std::vector<double> > IsoO_LF_L2norm_Data = ConvergenceIsoO_LF(mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, Td_IsoO, kd_IsoO, lp_IsoO, acc_l, n_IsoO, m_grid_IsoO);
        Read_Out(IsoO_LF_L2norm_Data, "Convergence_IsoO_LF");
    }
    
//----------------------------------------------------------------------------------------------------------------------------------------------------
    
    else if (Run_Type == "FreO_LF")
    {
        // Initialising test values to show the code works (not necessarily realisitic)
        //                                 mc   me  rp       pd        Td    kd  lp      met acc_k   acc_l  num_p acc_m
        std::array<double, 12> val_FreO = {10., 1., 43.5589, 1.47e-11, 118., 1., 1.e-10, 0., 1.e-10, 5.e-4, 5.e3, 2.e-4};
        auto [mc_FreO, me_FreO, rp_FreO, pd_FreO, Td_FreO, kd_FreO, lp_FreO, met, acc_k, acc_l, n_FreO, acc_m] = val_FreO;
        std::array<double, 2> p_ana_c_ini_FreO = {5.94e-11, 9.21};    // initial density fit parameters
        // Calculating mass grid where density fit parameters don't change more than acc_m
        auto returned_values_FreO = Mass_Grid_Solver_FreO(mc_FreO, me_FreO, rp_FreO, pd_FreO, Td_FreO, kd_FreO, lp_FreO, met, acc_k, n_FreO, acc_m, p_ana_c_ini_FreO);
        std::vector<double> m_grid_FreO = returned_values_FreO.first;         // mass grid
        
        std::pair< std::vector< std::vector<double> >, std::vector< std::vector<double> > > FreO_LF_Data = L_Finder_FreO(mc_FreO, me_FreO, rp_FreO, pd_FreO, Td_FreO, kd_FreO, lp_FreO, met, acc_k, acc_l, n_FreO, m_grid_FreO);
        Read_Out(FreO_LF_Data, "FreO_LF");
        
        // Performing convergence test and reading out the L2norm data to be plotted
        std::vector< std::vector<double> > FreO_LF_L2norm_Data = ConvergenceFreO_LF(mc_FreO, me_FreO, rp_FreO, pd_FreO, Td_FreO, kd_FreO, lp_FreO, met, acc_k, acc_l, n_FreO, m_grid_FreO);
        Read_Out(FreO_LF_L2norm_Data, "Convergence_FreO_LF");
    }
    
//----------------------------------------------------------------------------------------------------------------------------------------------------
    
    else if (Run_Type == "ComO_LF")
    {
        // Initialising test values to show the code works (not necessarily realisitic)
        //                                 mc   me  rp       pd        Td    kd  lp      ad    d_to_g acc_k   acc_l   num_p acc_m
        std::array<double, 13> val_ComO = {10., 1., 43.5589, 1.47e-11, 118., 1., 1.e-10, 1.e-2, 0.01, 1.e-10, 1.e-10, 1.e3, 1.e-3};
        auto [mc_ComO, me_ComO, rp_ComO, pd_ComO, Td_ComO, kd_ComO, lp_ComO, ad, d_to_g, acc_k, acc_l, n_ComO, acc_m] = val_ComO;
        std::array<double, 2> p_ana_c_ini_ComO = {1.00e-11, 1.00};      // initial density fit parameters
        
        // Calculating mass grid where density fit parameters don't change more than acc_m
        auto returned_values_ComO = Mass_Grid_Solver_ComO(mc_ComO, me_ComO, rp_ComO, pd_ComO, Td_ComO, kd_ComO, lp_ComO, ad, d_to_g, acc_k, n_ComO, acc_m, p_ana_c_ini_ComO);
        std::vector<double> m_grid_ComO = returned_values_ComO.first;         // mass grid
        
        std::pair< std::vector< std::vector<double> >, std::vector< std::vector<double> > > ComO_LF_Data = L_Finder_ComO(mc_ComO, me_ComO, rp_ComO, pd_ComO, Td_ComO, kd_ComO, lp_ComO, ad, d_to_g, acc_k, acc_l, n_ComO, m_grid_ComO);
        Read_Out(ComO_LF_Data, "ComO_LF");
        
        // Performing convergence test and reading out the L2norm data to be plotted
        std::vector< std::vector<double> > ComO_LF_L2norm_Data = ConvergenceComO_LF(mc_ComO, me_ComO, rp_ComO, pd_ComO, Td_ComO, kd_ComO, lp_ComO, ad, d_to_g, acc_k, acc_l, n_ComO, m_grid_ComO);
        Read_Out(ComO_LF_L2norm_Data, "Convergence_ComO_LF");
    }
    
//----------------------------------------------------------------------------------------------------------------------------------------------------
    
    else if (Run_Type == "IsoO_E")
    {
        bool Evo1 = false;
        bool Evo3_diff_opacity = false;
        bool Evo3_diff_alpha = false;
        bool Evo3_Max_Acc = false;
        
        //                                 mc   me  rp       pd        Td    kd     lp     acc_l   num_p acc_m  t_tot n_t
        std::array<double, 12> val_IsoO = {10., 1., 43.5589, 1.47e-11, 118., 6.e-5, 1.e-4, 1.e-10, 1.e4, 1.e-3, 2.e6, 1.e3};
        auto [mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, Td_IsoO, kd_IsoO, lp_IsoO, acc_l, n_i_IsoO, acc_m, t_tot_IsoO, n_t_IsoO] = val_IsoO;
        std::array<double, 2> p_ana_c_ini_IsoO = {6.14e-11, 9.20};    // initial density fit parameters
        //                                  Sigma_outer      Ms  Rp r_in r_out alpha num_sig
        std::array<double, 7> val_Sigma = {143.365/sqrt(3.), 1., 5.2, 2., 3., 1.e-5, 1e4};
        auto [Val_outer_IsoO, Ms_IsoO, Rp_IsoO, r_in_IsoO, r_out_IsoO, alpha_IsoO, num_sig_IsoO] = val_Sigma;
        
        if (Evo1)
        {
            // Calculating mass grid where density fit parameters don't change more than acc_m
            auto returned_values_IsoO = Mass_Grid_Solver_IsoO(mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, Td_IsoO, kd_IsoO, lp_IsoO, n_i_IsoO, acc_m, p_ana_c_ini_IsoO);
            std::vector<double> m_grid_IsoO = returned_values_IsoO.first;         // mass grid
            
            std::pair< std::vector< std::vector<double> >, std::vector< std::vector< std::vector<double> > > > EIO1, EIO2;
            
            EIO1 = EvolverIsoO(mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, Td_IsoO, kd_IsoO, lp_IsoO, acc_l, n_i_IsoO, acc_m, p_ana_c_ini_IsoO, m_grid_IsoO, t_tot_IsoO, n_t_IsoO, Val_outer_IsoO, Ms_IsoO, Rp_IsoO, r_in_IsoO, r_out_IsoO, alpha_IsoO, num_sig_IsoO, true);
//            EIO2 = EvolverIsoO(mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, kd_IsoO, Td_IsoO, lp_IsoO, acc_l, n_i_IsoO, acc_m, p_ana_c_ini_IsoO, m_grid_IsoO, t_tot_IsoO, n_t_IsoO, Val_outer_IsoO, Ms_IsoO, Rp_IsoO, r_in_IsoO, r_out_IsoO, alpha_IsoO, num_sig_IsoO, true);
            Read_Out(EIO1, "IsoO_E");
        }
        
        if (Evo3_diff_opacity)
        {
            double kd_IsoOA = 6e-5;
            double kd_IsoOB = 1e-5;
            double kd_IsoOC = 0.5e-5;
            double lp_IsoOA = 0.70e-4;
            double lp_IsoOB = 4.14e-4;
            double lp_IsoOC = 8.29e-4;
            
            auto returned_values_IsoOA = Mass_Grid_Solver_IsoO(mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, Td_IsoO, kd_IsoOA, lp_IsoOA, n_i_IsoO, acc_m, p_ana_c_ini_IsoO);
            std::vector<double> m_grid_IsoOA= returned_values_IsoOA.first;          // mass grid A
            auto returned_values_IsoOB = Mass_Grid_Solver_IsoO(mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, Td_IsoO, kd_IsoOB, lp_IsoOB, n_i_IsoO, acc_m, p_ana_c_ini_IsoO);
            std::vector<double> m_grid_IsoOB= returned_values_IsoOB.first;          // mass grid B
            auto returned_values_IsoOC = Mass_Grid_Solver_IsoO(mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, Td_IsoO, kd_IsoOC, lp_IsoOC, 3*n_i_IsoO, acc_m, p_ana_c_ini_IsoO);
            std::vector<double> m_grid_IsoOC= returned_values_IsoOC.first;          // mass grid C
            
            std::pair< std::vector< std::vector<double> >, std::vector< std::vector< std::vector<double> > > > EIOA, EIOB, EIOC;
            EIOA = EvolverIsoO(mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, Td_IsoO, kd_IsoOA, lp_IsoOA, acc_l, n_i_IsoO, acc_m, p_ana_c_ini_IsoO, m_grid_IsoOA, t_tot_IsoO, n_t_IsoO, Val_outer_IsoO, Ms_IsoO, Rp_IsoO, r_in_IsoO, r_out_IsoO, alpha_IsoO, num_sig_IsoO, false);
            EIOB = EvolverIsoO(mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, Td_IsoO, kd_IsoOB, lp_IsoOB, acc_l, n_i_IsoO, acc_m, p_ana_c_ini_IsoO, m_grid_IsoOB, t_tot_IsoO, n_t_IsoO, Val_outer_IsoO, Ms_IsoO, Rp_IsoO, r_in_IsoO, r_out_IsoO, alpha_IsoO, num_sig_IsoO, false);
            EIOC = EvolverIsoO(mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, Td_IsoO, kd_IsoOC, lp_IsoOC, acc_l, 3*n_i_IsoO, acc_m, p_ana_c_ini_IsoO, m_grid_IsoOC, t_tot_IsoO, n_t_IsoO, Val_outer_IsoO, Ms_IsoO, Rp_IsoO, r_in_IsoO, r_out_IsoO, alpha_IsoO, num_sig_IsoO, false);
            
            Read_Out({EIOA, EIOB, EIOC}, "IsoO_E_ABC");
        }
        
        if (Evo3_diff_alpha)
        {
            double alpha_IsoOD = 1e-5;
            double alpha_IsoOE = 1e-2;
            
            auto returned_values_IsoODEF = Mass_Grid_Solver_IsoO(mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, Td_IsoO, kd_IsoO, lp_IsoO, n_i_IsoO, acc_m, p_ana_c_ini_IsoO);
            std::vector<double> m_grid_IsoODEF= returned_values_IsoODEF.first;          // mass grid DEF
            
            std::pair< std::vector< std::vector<double> >, std::vector< std::vector< std::vector<double> > > > EIOD, EIOE, EIOF;
            EIOD = EvolverIsoO(mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, Td_IsoO, kd_IsoO, lp_IsoO, acc_l, n_i_IsoO, acc_m, p_ana_c_ini_IsoO, m_grid_IsoODEF, t_tot_IsoO, n_t_IsoO, Val_outer_IsoO, Ms_IsoO, Rp_IsoO, r_in_IsoO, r_out_IsoO, alpha_IsoOD, num_sig_IsoO, false);
            EIOE = EvolverIsoO(mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, Td_IsoO, kd_IsoO, lp_IsoO, acc_l, n_i_IsoO, acc_m, p_ana_c_ini_IsoO, m_grid_IsoODEF, t_tot_IsoO, n_t_IsoO, Val_outer_IsoO, Ms_IsoO, Rp_IsoO, r_in_IsoO, r_out_IsoO, alpha_IsoOE, num_sig_IsoO, false);
            EIOF = EvolverIsoO(mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, Td_IsoO, kd_IsoO, lp_IsoO, acc_l, n_i_IsoO, acc_m, p_ana_c_ini_IsoO, m_grid_IsoODEF, t_tot_IsoO, n_t_IsoO, Val_outer_IsoO, Ms_IsoO, Rp_IsoO, r_in_IsoO, r_out_IsoO, alpha_IsoOD, num_sig_IsoO, true);
            
            Read_Out({EIOD, EIOE, EIOF}, "IsoO_E_DEF");
        }
        
        if (Evo3_Max_Acc)
        {
            double kd_IsoOA = 6e-5;
            double kd_IsoOB = 1e-5;
            double kd_IsoOC = 0.5e-5;
            double lp_IsoOA = 0.70e-4;
            double lp_IsoOB = 4.14e-4;
            double lp_IsoOC = 8.29e-4;
            double alpha_IsoOD = 1e-5;
            
            auto returned_values_IsoOA = Mass_Grid_Solver_IsoO(mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, Td_IsoO, kd_IsoOA, lp_IsoOA, n_i_IsoO, acc_m, p_ana_c_ini_IsoO);
            std::vector<double> m_grid_IsoOA= returned_values_IsoOA.first;          // mass grid A
            auto returned_values_IsoOB = Mass_Grid_Solver_IsoO(mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, Td_IsoO, kd_IsoOB, lp_IsoOB, n_i_IsoO, acc_m, p_ana_c_ini_IsoO);
            std::vector<double> m_grid_IsoOB= returned_values_IsoOB.first;          // mass grid B
            auto returned_values_IsoOC = Mass_Grid_Solver_IsoO(mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, Td_IsoO, kd_IsoOC, lp_IsoOC, 3*n_i_IsoO, acc_m, p_ana_c_ini_IsoO);
            std::vector<double> m_grid_IsoOC= returned_values_IsoOC.first;          // mass grid C
            
            std::pair< std::vector< std::vector<double> >, std::vector< std::vector< std::vector<double> > > > EIOG, EIOH, EIOI;
            EIOG = EvolverIsoO(mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, Td_IsoO, kd_IsoOA, lp_IsoOA, acc_l, n_i_IsoO, acc_m, p_ana_c_ini_IsoO, m_grid_IsoOA, t_tot_IsoO, n_t_IsoO, Val_outer_IsoO, Ms_IsoO, Rp_IsoO, r_in_IsoO, r_out_IsoO, alpha_IsoOD, num_sig_IsoO, true);
            EIOH = EvolverIsoO(mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, Td_IsoO, kd_IsoOB, lp_IsoOB, acc_l, n_i_IsoO, acc_m, p_ana_c_ini_IsoO, m_grid_IsoOB, t_tot_IsoO, n_t_IsoO, Val_outer_IsoO, Ms_IsoO, Rp_IsoO, r_in_IsoO, r_out_IsoO, alpha_IsoOD, num_sig_IsoO, true);
            EIOI = EvolverIsoO(mc_IsoO, me_IsoO, rp_IsoO, pd_IsoO, Td_IsoO, kd_IsoOC, lp_IsoOC, acc_l, 3*n_i_IsoO, acc_m, p_ana_c_ini_IsoO, m_grid_IsoOC, t_tot_IsoO, n_t_IsoO, Val_outer_IsoO, Ms_IsoO, Rp_IsoO, r_in_IsoO, r_out_IsoO, alpha_IsoOD, num_sig_IsoO, true);
            
            Read_Out({EIOG, EIOH, EIOI}, "IsoO_E_GHI");
        }
        
    }
    
//----------------------------------------------------------------------------------------------------------------------------------------------------
    
    else if (Run_Type == "FreO_E")
    {
        //
    }
    
//----------------------------------------------------------------------------------------------------------------------------------------------------
    
    else if (Run_Type == "ComO_E")
    {
        //
    }
    
    return 0;
}
