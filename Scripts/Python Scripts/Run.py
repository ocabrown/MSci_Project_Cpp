import numpy as np
import Data_Input as DI
import Class_Maker as CM
import Plotter as P
import Convergence_Plotter as CP
import Max_Acc_Rate_Plotter as MARP
import L_Finder_Plotters as LFP





#%% Loading Comparison Data

Movshovitz_Figure_3 = DI.Data_Preparer(["Figure 3 - Density.csv", "Figure 3 - Opacity.csv", "Figure 3 - Temperature.csv"])
Movshovitz_Figure_3_Opacity_Log = np.log(Movshovitz_Figure_3[1])
Movshovitz_Figure_4 = DI.Data_Preparer(["Figure 4 - Density.csv", "Figure 4 - Opacity.csv", "Figure 4 - Temperature.csv"])
Movshovitz_Figure_4_Opacity_Log = np.log(Movshovitz_Figure_4[1])

Crida_Sigma_Data = DI.Data_Extractor_Sigma("Sigma (0.05, -5.5).csv")


#%%


#                           Constant Temperature



#%% Load Data

IsoT1_Data = DI.Read_In_Var("IsoT")

#%% Main Plots

P.Plot(IsoT1_Data, "p_comp", "r", data_type = "IsoT")
P.Plot(IsoT1_Data, "dP/dr", "r", data_type = "IsoT")

#%% Extra Plots

P.Plot(IsoT1_Data, "m", "r", data_type = "IsoT")
P.Plot(IsoT1_Data, "P", "r", data_type = "IsoT")

#%% Viva plot

P.Plot(IsoT1_Data, "p+dP_dr", "r", data_type = "IsoT")

#%% Convergence Test

IsoT1_L2norm_Data = DI.Read_In_L2norm("IsoT")
CP.ConvergenceIsoT(IsoT1_L2norm_Data)

#%%


#                           Constant Density



#%% Load Data

IsoD1_Data = DI.Read_In_Var("IsoD")

#%% Main Plots

P.Plot(IsoD1_Data, "T_comp", "r", data_type = "IsoD")
P.Plot(IsoD1_Data, "dP/dr", "r", data_type = "IsoD")

#%% Extra Plots

P.Plot(IsoD1_Data, "m", "r", data_type="IsoD")
P.Plot(IsoD1_Data, "P", "r", data_type="IsoD")

#%% Viva Plot

P.Plot(IsoD1_Data, "T+dP_dr", "r", data_type = "IsoD")

#%% Convergence Test

IsoD1_L2norm_Data = DI.Read_In_L2norm("IsoD")
CP.ConvergenceIsoD(IsoD1_L2norm_Data)

#%%


#                           Constant Opacity



#%% Load Data

IsoO1_Data = DI.Read_In_Var("IsoO")
IsoO1 = CM.IsoO_Class(IsoO1_Data)

#%% Main Plots

params = [7, 2, 20, 20, 15, 20]

P.Plot([IsoO1], "m", "r", data_type="IsoO", Params=params)
P.Plot([IsoO1], "p", "r", data_type="IsoO", Params=params)
P.Plot([IsoO1], "k", "r", data_type="IsoO", Params=params) 
P.Plot([IsoO1], "T", "r", data_type="IsoO", Params=params)
P.Plot([IsoO1], "P", "r", data_type="IsoO", Params=params)
P.Plot([IsoO1], "s", "r", data_type="IsoO", Params=params)

#%% Extra Plots

params = [7, 2, 20, 20, 15, 20]

P.Plot([IsoO1], "p", "r", data_type="IsoO", Comparison_Data=Movshovitz_Figure_3, Params=params)
P.Plot([IsoO1], "k", "r", data_type="IsoO", Comparison_Data=Movshovitz_Figure_3, Params=params)
P.Plot([IsoO1], "T", "r", data_type="IsoO", Comparison_Data=Movshovitz_Figure_3, Params=params)

#%% Viva Plot

params = [7, 2, 20, 20, 15, 20]

P.Plot([IsoO1], "p+T", "r", data_type="IsoO", Params=params)

#%% Convergence Test

IsoO1_L2norm_Data = DI.Read_In_L2norm("IsoO")
CP.ConvergenceIsoO(IsoO1_L2norm_Data)

#%%


#                           Freedman Opacity



#%% Load Data

FreO1_Data = DI.Read_In_Var("FreO")
FreO1 = CM.FreO_Class(FreO1_Data)

#%% Main Plots

params = [7, 2, 20, 20, 15, 20]

P.Plot([FreO1], "m", "r", data_type="FreO", Params=params)
P.Plot([FreO1], "p", "r", data_type="FreO", Params=params)
P.Plot([FreO1], "k", "r", data_type="FreO", Params=params) 
P.Plot([FreO1], "T", "r", data_type="FreO", Params=params)
P.Plot([FreO1], "P", "r", data_type="FreO", Params=params)
P.Plot([FreO1], "s", "r", data_type="FreO", Params=params)

#%% Extra Plots

params = [7, 2, 20, 20, 15, 20]

P.Plot([FreO1], "p", "r", data_type="FreO", Comparison_Data=Movshovitz_Figure_3, Params=params)
P.Plot([FreO1], "k", "r", data_type="FreO", Comparison_Data=Movshovitz_Figure_3, Params=params)
P.Plot([FreO1], "T", "r", data_type="FreO", Comparison_Data=Movshovitz_Figure_3, Params=params)

#%% Viva Plot

params = [7, 2, 20, 20, 15, 20]

P.Plot([FreO1], "k+p+T", "r", data_type="FreO", Params=params)

#%% Convergence Test

FreO1_L2norm_Data = DI.Read_In_L2norm("FreO")
CP.ConvergenceFreO(FreO1_L2norm_Data)

#%%


#                           Combined Opacity



#%% Load Data

ComO1_Data = DI.Read_In_Var("ComO")
ComO1 = CM.ComO_Class(ComO1_Data)

#%% Main Plots

params = [7, 2, 20, 20, 15, 20]

P.Plot([ComO1], "m", "r", data_type="ComO", Params=params)
P.Plot([ComO1], "p", "r", data_type="ComO", Params=params)
P.Plot([ComO1], "k", "r", data_type="ComO", Params=params) 
P.Plot([ComO1], "T", "r", data_type="ComO", Params=params)
P.Plot([ComO1], "P", "r", data_type="ComO", Params=params)
P.Plot([ComO1], "s", "r", data_type="ComO", Params=params)

#%% Extra Plots

params = [7, 2, 20, 20, 15, 20]

P.Plot([ComO1], "p", "r", data_type="ComO", Comparison_Data=Movshovitz_Figure_3, Params=params)
P.Plot([ComO1], "k", "r", data_type="ComO", Comparison_Data=Movshovitz_Figure_3, Params=params)
P.Plot([ComO1], "T", "r", data_type="ComO", Comparison_Data=Movshovitz_Figure_3, Params=params)

#%% Viva Plot

params = [7, 2, 20, 20, 15, 20]

P.Plot([ComO1], "k+p+T", "r", data_type="ComO", Params=params)

#%% Convergence Test

ComO1_L2norm_Data = DI.Read_In_L2norm("ComO")
CP.ConvergenceComO(ComO1_L2norm_Data)

#%%


#                           Maximum Mass Accretion



#%% Plotting 0 mass planet

Sigma_Data_0_Mass = DI.Read_In_Max_Acc("Max_Acc_Rate")
MARP.Max_Acc_Rate_Plot(Sigma_Data_0_Mass)

#%% Plotting Crida Sigma eqn 14 Sigma solution - comparing to Figure 11

Sigma_Data_Crida = DI.Read_In_Max_Acc("Max_Acc_Rate")
MARP.Max_Acc_Rate_Plot(Sigma_Data_Crida, Crida_Comparison=True)

#%% Plotting different alpha values

num_sig_Sigma = 1e4
Sigma_Data_Diff_Alpha = DI.Read_In_Max_Acc("Max_Acc_Rate_Diff_Alpha")
MARP.Max_Acc_Rate_Diff_Alpha_Plot(Sigma_Data_Diff_Alpha, num_sig_Sigma)

#%% Convergence of Sigma solution

Sigma_L2norm_Data = DI.Read_In_L2norm("Max_Acc_Rate_Sigma")
CP.ConvergenceSigma(Sigma_L2norm_Data)

#%% Convergence of Mdot solution with analytic solution

MdotMax_L2norm_Data = DI.Read_In_L2norm("Max_Acc_Rate_MdotMax")
CP.ConvergenceMdotMax(MdotMax_L2norm_Data)

# Can be confident that the Mdot_max equation is suitably implemented

#%%


#                           Luminosity Finder IsoO



#%% Plotting L_Finder

IsoO_LF_Data = DI.Read_In_LF("IsoO_LF")
LFP.L_Finder_Plot(IsoO_LF_Data)

#%% Plotting L_Finder Convergence

IsoO_LF_conv_Data = DI.Read_In_L2norm("IsoO_LF")
CP.ConvergenceLum(IsoO_LF_conv_Data)

#%%


#                           Luminosity Finder FreO



#%% Plotting L_Finder

FreO_LF_Data = DI.Read_In_LF("FreO_LF")
LFP.L_Finder_Plot(FreO_LF_Data)

#%% Plotting L_Finder Convergence

FreO_LF_conv_Data = DI.Read_In_L2norm("FreO_LF")
CP.ConvergenceLum(FreO_LF_conv_Data)

#%%


#                           Luminosity Finder ComO



#%% Plotting L_Finder

ComO_LF_Data = DI.Read_In_LF("ComO_LF")
LFP.L_Finder_Plot(ComO_LF_Data)

#%% Plotting L_Finder Convergence

ComO_LF_conv_Data = DI.Read_In_L2norm("ComO_LF")
CP.ConvergenceLum(ComO_LF_conv_Data)

#%%


#                               Evolver IsoO



#%% Single Evolution

EIO1 = DI.Read_In_E("IsoO_E")[0]

params = [7, 2, 20, 20, 15, 20]
Calc1 = P.Plot([EIO1], None, "calc", Params=params)

#%% Plots for Single Evolutioin
params = [7, 2, 20, 20, 15, 20]

#P.Plot([EIO1], "M", "t", Params=params, calc=Calc1)
#P.Plot([EIO1], "L", "t", Params=params, calc=Calc1)
#P.Plot([EIO1], "L_t", "t", Params=params, calc=Calc1)
#P.Plot([EIO1], "T_t,r_t", "t", Params=params, calc=Calc1)
P.Plot([EIO1], "T_t4,r_t2", "t", Params=params, calc=Calc1)
#P.Plot([EIO1], "T_t,p_t", "r_t", Params=params, Evo=True)
#P.Plot([EIO1], "p_t comp", "r_t", Params=params, Evo=True)
P.Plot([EIO1], "m_t comp", "r_t", Params=params, Evo=True)
#P.Plot([EIO1], "T_t,p_t diff", "r_t", Params=params, Evo=True)
#P.Plot([EIO1], "T_t", "r_t", Params=params, Evo=True)
P.Plot([EIO1], "L+M", "t", Params=params, calc=Calc1)
#P.Plot([EIO1], "Mdot", "t", Params=params, calc=Calc1)

#%% Plotting different constant opacity evolutions

EIOA, EIOB, EIOC = DI.Read_In_E("IsoO_E_ABC")

params = [7, 2, 20, 20, 15, 20]

P.Plot([EIOA,EIOB,EIOC], "M", "t_pres", Params=params, Evo=True)
P.Plot([EIOA,EIOB,EIOC], "L", "t_pres", Params=params, Evo=True)

#%% Plotting different alpha values evolutions

EIOD, EIOE, EIOF = DI.Read_In_E("IsoO_E_DEF")

params = [7, 2, 20, 20, 12, 20]

P.Plot([EIOD,EIOE,EIOF], "Mdot", "t_pres", Params=params, Evo=True, acc_reg=True)

#%% Plotting different constant opacity evolutions

EIOG, EIOH, EIOI = DI.Read_In_E("IsoO_E_GHI")

params = [7, 2, 20, 20, 15, 20]

P.Plot([EIOG,EIOH,EIOI], "Mdot", "t_pres", Params=params, Evo=True)



#%% Plotting mdot and m evolution

EIOI = DI.Read_In_E("IsoO_E_GHI")[2]

params = [7, 2, 20, 20, 15, 20]

P.Plot([EIOI], "Mdot+M", "t", Params=params)

#%%


#                               Evolver FreO



#%% Single Evolution

#...

#%%


#                               Evolver ComO



#%% Single Evolution

#...


