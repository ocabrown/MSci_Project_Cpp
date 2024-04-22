import numpy as np
import scipy.optimize as so
import matplotlib.pyplot as plt
import Data_Input as DI



def fit_func(r,a,b):
    return a/np.sqrt(r) + b

def straight_func(r,b):
    return -1/2*r + b



def Sigma_Plot(data, Crida_Comparison=False):
    
    AU = 1.495978707e13
    Sig_Val = 800
    
    r, Sigma = data
    
    r /= (5.2 * AU)
    Sigma /= Sig_Val
    
    Figsize, lw, Ticksize, Labelsize, Legendfontsize = [7, 2, 20, 20, 15]
    
    plt.rc("figure", figsize=(Figsize,Figsize))
    plt.rc("xtick", labelsize=Ticksize)
    plt.rc("ytick", labelsize=Ticksize)
    plt.rc("axes", labelsize=Labelsize)
    plt.rc("legend", fontsize=Legendfontsize)
    
    plt.figure()
    plt.plot(r, Sigma, linewidth=lw, label="Analytical Model")
    
    if Crida_Comparison == True:
        Crida_Sigma_Data = DI.Data_Extractor_Sigma("Crida Data/Sigma (0.05, -5.5).csv")
        plt.plot(Crida_Sigma_Data[0], Crida_Sigma_Data[1], linewidth=lw, label="Crida Data")
    
    plt.xlabel("Radius / $R_p$")
    plt.ylabel("$\Sigma$")
    plt.legend()
    plt.show()
    plt.legend()
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()




def Max_Acc_Rate_Plot(data, Crida_Comparison=False):
    
    AU = 1.495978707e13
    M_E = 5.972e27
    R_J = 6.9911e9
    Sigma, r_d, Mp, Rp, num_sig = data
    
    darkness = 0.5
    colour = "blue"
    Figsize, lw, Ticksize, Labelsize, Legendfontsize = [7, 2, 20, 20, 15]
    
    plt.rc("figure", figsize=(Figsize,Figsize))
    plt.rc("xtick", labelsize=Ticksize)
    plt.rc("ytick", labelsize=Ticksize)
    plt.rc("axes", labelsize=Labelsize)
    plt.rc("legend", fontsize=Legendfontsize)
    
    plt.figure()
    
    if Crida_Comparison == True:
        plt.plot(r_d/AU, Sigma, color="blue", linewidth=lw, label=f"0 Myr, ~{int(np.floor(Mp/M_E))} $M_{'E'}$")
        Crida_Sigma_Data = DI.Data_Extractor_Sigma("Sigma (0.05, -5.5).csv")
        plt.plot(Crida_Sigma_Data[0]*Rp/AU, Crida_Sigma_Data[1], color="orange", linewidth=lw, linestyle="--", label="Crida Data")
        
        #plt.axvspan(r_d[int(num_sig/5)]/AU, r_d[int(2*num_sig/5)-1]/AU, alpha=darkness, facecolor=colour)
        #plt.axvspan(r_d[int(3*num_sig/5)-1]/AU, r_d[int(4*num_sig/5)-1]/AU, alpha=darkness, facecolor=colour)
        #plt.axvline(r_d[int(num_sig/5)]/AU, color="blue", linewidth=lw, linestyle="--")
        #plt.axvline(r_d[int(2*num_sig/5)-1]/AU, color="blue", linewidth=lw, linestyle="--")
        #plt.axvline(r_d[int(3*num_sig/5)-1]/AU, color="blue", linewidth=lw, linestyle="--")
        #plt.axvline(r_d[int(4*num_sig/5)-1]/AU, color="blue", linewidth=lw, linestyle="--")
        plt.xlabel("Radius / AU")
        plt.ylabel("$\Sigma$ / $\Sigma_0$")
    
    else:
        if Mp == 0.:
            plt.loglog(r_d/AU, Sigma, color="blue", linewidth=lw, label=f"{Mp/M_E} $M_\oplus$")
            coeff = so.curve_fit(straight_func,np.log(r_d),np.log(Sigma))[0]
            Sigma_ana = straight_func(np.log(r_d),coeff[0])
            Sigma_ana = np.e**(Sigma_ana)
            plt.loglog(r_d/AU, Sigma_ana, color="orange", linewidth=lw, linestyle="--", label="$1/\sqrt{r}$ Fit")
        else:
            plt.plot(r_d/AU, Sigma, color="blue", linewidth=lw, label=f"0 Myr, {Mp/M_E} $M_\oplus$")
            plt.axvspan(r_d[int(num_sig/5)]/AU, r_d[int(2*num_sig/5)-1]/AU, alpha=darkness, facecolor=colour)
            plt.axvspan(r_d[int(3*num_sig/5)-1]/AU, r_d[int(4*num_sig/5)-1]/AU, alpha=darkness, facecolor=colour)
        plt.xlabel("Radius / AU")
        plt.ylabel("$\Sigma$ / g$\cdot$cm$^{-2}$")
    
    plt.legend()
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()


def Max_Acc_Rate_Diff_Alpha_Plot(data, num_sig_Sigma):
    
    AU = 1.495978707e13
    
    Sigma0, Sigma1, Sigma2, Sigma3, r0, r1, r2, r3 = data
    
    darkness = 0.5
    
    Figsize, lw, Ticksize, Labelsize, Legendfontsize = [7, 2, 20, 20, 15]
    
    plt.rc("figure", figsize=(Figsize,Figsize))
    plt.rc("xtick", labelsize=Ticksize)
    plt.rc("ytick", labelsize=Ticksize)
    plt.rc("axes", labelsize=Labelsize)
    plt.rc("legend", fontsize=Legendfontsize)
    
    plt.figure()
    plt.plot(r0/AU, Sigma0, linestyle="-", color="blue", label=r"$\alpha$=10$^{-2}$, 0.1 $M_{\oplus}$")
    plt.plot(r2/AU, Sigma2, linestyle="--", color="blue", label=r"$\alpha$=10$^{-2}$, 300 $M_{\oplus}$")
    plt.plot(r1/AU, Sigma1, linestyle="-", color="orange", label=r"$\alpha$=10$^{-5}$, 0.1 $M_{\oplus}$")
    plt.xlabel("Radius / AU")
    plt.ylabel("$\Sigma$ / g$\cdot$cm$^{-2}$")
    plt.legend()
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    
    plt.figure()
    plt.plot(r2/AU, Sigma2, linestyle="-", color="blue", label=r"$\alpha$=10$^{-2}$, 300$M_{\oplus}$")
    plt.axvspan(r2[int(num_sig_Sigma/5)]/AU, r2[int(2*num_sig_Sigma/5)-1]/AU, alpha=darkness, facecolor="red", label="T12")
    plt.axvspan(r2[int(3*num_sig_Sigma/5)-1]/AU, r2[int(4*num_sig_Sigma/5)-1]/AU, alpha=darkness, facecolor="red")
    plt.axvspan(r3[int(num_sig_Sigma/5)]/AU, r3[int(2*num_sig_Sigma/5)-1]/AU, alpha=darkness, facecolor="orange", label="L09")
    plt.axvspan(r3[int(3*num_sig_Sigma/5)-1]/AU, r3[int(4*num_sig_Sigma/5)-1]/AU, alpha=darkness, facecolor="orange")
    plt.xlabel("Radius / AU")
    plt.ylabel("$\Sigma$ / g$\cdot$cm$^{-2}$")
    plt.legend()
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()






