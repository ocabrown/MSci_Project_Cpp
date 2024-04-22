import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit



def N1(x,a):
    return -1*x + a

def N2(x,a):
    return -2*x + a

def straight_line(x, m, c):
    return m * x + c



"""
from decimal import Decimal
F = 1e16
kappa = 1e-18
alpha = 1e-13
ntot_init = 1e3

lhs = Decimal(F)*Decimal(kappa)
a = Decimal(alpha)/lhs
b = Decimal(1)
c = -Decimal(ntot_init)
nplus = (-b + np.sqrt(b*b-Decimal(4)*a*c))/Decimal(2.)/a
fplus = nplus / Decimal(ntot_init)
fzero = (Decimal(1)-fplus)


F = 1e16
kappa = 1e-18
alpha = 1e-13
ntot_init = 1e3

lhs = F*kappa
a = alpha/lhs
b = 1
c = -ntot_init
nplus = (-b + np.sqrt(b*b-4*a*c))/2./a
fplus = nplus / ntot_init
fzero = (1-fplus)
"""



def ConvergenceIsoT(data):
    
    L2norms_p, L2norms_dP_dr, min_dP_drs, n = data
    
    vals1, val_errs1 = curve_fit(straight_line, np.log10(n), np.log10(L2norms_p))
    order1, order_std1 = vals1[0], np.sqrt(val_errs1[0][0])
    print("Order =", order1, "±", order_std1)
    
    vals2, val_errs2 = curve_fit(straight_line, np.log10(n), np.log10(L2norms_dP_dr))
    order2, order_std2 = vals2[0], np.sqrt(val_errs2[0][0])
    print("Order =", order2, "±", order_std2)
    
    p_popt1, p_pcov1 = curve_fit(N1,np.log10(n), np.log10(L2norms_p))
    P_popt1, P_pcov1 = curve_fit(N1,np.log10(n), np.log10(L2norms_dP_dr))
    p_popt2, p_pcov2 = curve_fit(N2,np.log10(n), np.log10(L2norms_p))
    P_popt2, P_pcov2 = curve_fit(N2,np.log10(n), np.log10(L2norms_dP_dr))
    
    if n[-1] > 1e3:
        ind = np.where(n==1e3)
        p_popt1, p_pcov1 = curve_fit(N1,np.log10(n[:int(ind[0][0]-5)]), np.log10(L2norms_p[:int(ind[0][0]-5)]))
        p_popt2, p_pcov2 = curve_fit(N2,np.log10(n[:int(ind[0][0]-5)]), np.log10(L2norms_p[:int(ind[0][0]-5)]))
    
    
    Figsize, lw, Ticksize, Labelsize, Legendfontsize = [7, 2, 20, 20, 15]
    
    plt.rc("figure", figsize=(Figsize,Figsize))
    plt.rc("xtick", labelsize=Ticksize)
    plt.rc("ytick", labelsize=Ticksize)
    plt.rc("axes", labelsize=Labelsize)
    plt.rc("legend", fontsize=Legendfontsize)
    
    plt.figure()
    plt.loglog(n,L2norms_p,color="blue",linewidth=lw)
    plt.loglog(n,10**N1(np.log10(n),p_popt1),color="orange",linestyle="--",linewidth=lw)
    plt.loglog(n,10**N2(np.log10(n),p_popt2),color="red",linestyle="--",linewidth=lw)
    plt.xlabel("Number of Cells")
    plt.ylabel("L$^2$-norm")
    #plt.title("Density L$^2$ Norm vs. Number of Cells")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    
    x = 1.6
    
    plt.figure()
    plt.loglog(n,L2norms_dP_dr*x,color="blue",linewidth=lw)
    plt.loglog(n,10**N1(np.log10(n),P_popt1)*x,color="orange",linestyle="--",linewidth=lw)
    plt.loglog(n,10**N2(np.log10(n),P_popt2)*x,color="red",linestyle="--",linewidth=lw)
    plt.xlabel("Number of Cells")
    plt.ylabel("L$^2$-norm")
    #plt.title("Pressure Gradient L$^2$ Norm vs. Number of Cells")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    
    plt.figure()
    plt.loglog(n, abs(min_dP_drs),linewidth=lw)
    plt.xlabel("Number of Cells")
    plt.ylabel("Minimum Pressure Gradient")
    #plt.title("Minimum Pressure Gradient vs. Number of Cells")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    
    plt.figure()
    plt.loglog(n,L2norms_p,linewidth=lw,color="red",label="Density L$^2$-norm")
    plt.loglog(n,L2norms_dP_dr,linewidth=lw,color="blue",label="dP/dr L$^2$-norm")
    plt.loglog(n,10**N1(np.log10(n),P_popt1),color="orange",linestyle="--",linewidth=lw)
    plt.xlabel("Number of Cells")
    plt.ylabel("L$^2$-norm")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.legend()
    plt.show()
    
    """
    fig, ax1 = plt.subplots()
    color1 = 'tab:red'
    ax1.set_xlabel('Number of Cells')
    ax1.set_ylabel("Density L$^2$ Norm", color=color1)
    l1 = ax1.loglog(n, L2norms_p, color=color1, linewidth=lw, label="Density L$^2$ Norm")
    ax1.tick_params(axis='y', labelcolor=color1)
    
    ax2 = ax1.twinx()
    
    color1 = 'tab:orange'
    ax2.set_ylabel('Pressure Gradient L$^2$ Norm', color=color1)
    l2 = ax2.loglog(n, L2norms_dP_dr, color=color1, linewidth=lw, label="Pressure Gradient L$^2$ Norm")
    ax2.tick_params(axis='y', labelcolor=color1)
    
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    
    leg = l1 + l2
    labs = [l.get_label() for l in leg]
    plt.legend(leg, labs, loc=9)
    """



# Constant Density Convergence Function

def ConvergenceIsoD(data):
    
    L2norms_T, L2norms_dP_dr, min_dP_drs, L2norms_dP_dr_lower, L2norms_dP_dr_upper, n = data
    
    vals, val_errs = curve_fit(straight_line, np.log10(n), np.log10(L2norms_T))
    order, order_std = vals[0], np.sqrt(val_errs[0][0])
    print(order, "±", order_std)
    
    T_popt1, T_pcov1 = curve_fit(N1,np.log10(n), np.log10(L2norms_T))
    P_popt1, P_pcov1 = curve_fit(N1,np.log10(n), np.log10(L2norms_dP_dr_upper))
    T_popt2, T_pcov2 = curve_fit(N2,np.log10(n), np.log10(L2norms_T))
    P_popt2, P_pcov2 = curve_fit(N2,np.log10(n), np.log10(L2norms_dP_dr_upper))
    
    Figsize, lw, Ticksize, Labelsize, Legendfontsize = [7, 2, 20, 20, 15]
    
    plt.rc("figure", figsize=(Figsize,Figsize))
    plt.rc("xtick", labelsize=Ticksize)
    plt.rc("ytick", labelsize=Ticksize)
    plt.rc("axes", labelsize=Labelsize)
    plt.rc("legend", fontsize=Legendfontsize)
    
    plt.figure()
    plt.loglog(n,L2norms_T, linewidth=lw)
    plt.loglog(n,10**N1(np.log10(n),T_popt1),color="orange",linestyle="--",linewidth=lw)
    plt.loglog(n,10**N2(np.log10(n),T_popt2),color="red",linestyle="--",linewidth=lw)
    plt.xlabel("Number of Cells")
    plt.ylabel("Temperature L$^2$ Norm")
    #plt.title("Temperature L$^2$ Norm vs. Number of Cells")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    
    plt.figure()
    plt.plot(n,L2norms_dP_dr, linewidth=lw)
    plt.xlabel("Number of Cells")
    plt.ylabel("Pressure Gradient L$^2$ Norm")
    #plt.title("Pressure Gradient L$^2$ Norm vs. Number of Cells")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    
    plt.figure()
    plt.plot(n,L2norms_dP_dr_lower, linewidth=lw)
    plt.xlabel("Number of Cells")
    plt.ylabel("Pressure Gradient L$^2$ Norm")
    #plt.title("Lower Pressure Gradient L$^2$ Norm vs. Number of Cells")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    
    plt.figure()
    plt.loglog(n,L2norms_dP_dr_upper, linewidth=lw)
    plt.loglog(n,10**N1(np.log10(n),P_popt1),color="orange",linestyle="--",linewidth=lw)
    plt.loglog(n,10**N2(np.log10(n),P_popt2),color="red",linestyle="--",linewidth=lw)
    plt.xlabel("Number of Cells")
    plt.ylabel("Pressure Gradient L$^2$ Norm")
    #plt.title("Upper Pressure Gradient L$^2$ Norm vs. Number of Cells")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    
    plt.figure()
    plt.plot(n, abs(min_dP_drs), linewidth=lw)
    plt.xlabel("Number of Cells")
    plt.ylabel("Minimum Pressure Gradient")
    #plt.title("Minimum Pressure Gradient vs. Number of Cells")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    
    plt.figure()
    plt.loglog(n,L2norms_T,linewidth=lw,color="red",label="Temperature L$^2$-norm")
    plt.loglog(n,10**N1(np.log10(n),T_popt1),color="orange",linestyle="--",linewidth=lw)
    plt.loglog(n,L2norms_dP_dr,linewidth=lw,color="blue",label="dP/dr L$^2$-norm")
    plt.loglog(n,10**N1(np.log10(n),P_popt1),color="orange",linestyle="--",linewidth=lw)
    plt.xlabel("Number of Cells")
    plt.ylabel("L$^2$-norm")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.legend()
    plt.show()
    
    fig, ax1 = plt.subplots()
    color1 = 'tab:red'
    ax1.set_xlabel('Number of Cells')
    ax1.set_ylabel("Temperature L$^2$ Norm", color=color1)
    l1 = ax1.loglog(n, L2norms_T, color=color1, linewidth=lw, label="Temperature L$^2$ Norm")
    ax1.tick_params(axis='y', labelcolor=color1)
    
    ax2 = ax1.twinx()
    
    color1 = 'tab:orange'
    ax2.set_ylabel('Pressure Gradient L$^2$ Norm', color=color1)
    l2 = ax2.loglog(n, L2norms_dP_dr_upper, color=color1, linewidth=lw, label="Pressure Gradient L$^2$ Norm")
    ax2.tick_params(axis='y', labelcolor=color1)
    
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    
    leg = l1 + l2
    labs = [l.get_label() for l in leg]
    plt.legend(leg, labs, loc=9)



# Constant Opacity Convergence Function

def ConvergenceIsoO(data):
    
    r0s, L2norms_dP_dr, min_dP_drs, n = data
    
    r0_high = r0s[-1]
    r0s = r0s[:-1]
    
    vals2, val_errs2 = curve_fit(straight_line, np.log10(n), np.log10(L2norms_dP_dr))
    order2, order_std2 = vals2[0], np.sqrt(val_errs2[0][0])
    print("Order =", order2, "±", order_std2)
    
    P_popt1, P_pcov1 = curve_fit(N1,np.log10(n), np.log10(L2norms_dP_dr))
    P_popt2, P_pcov2 = curve_fit(N2,np.log10(n), np.log10(L2norms_dP_dr))
    
    r_popt1, r_pcov1 = curve_fit(N1,np.log10(n), np.log10(abs((r0_high-r0s)/r0_high)))
    
    Figsize, lw, Ticksize, Labelsize = [7, 2, 20, 20]
    
    plt.rc("figure", figsize=(Figsize,Figsize))
    plt.rc("xtick", labelsize=Ticksize)
    plt.rc("ytick", labelsize=Ticksize)
    plt.rc("axes", labelsize=Labelsize)
    
    # Convergence Plots
    plt.figure()
    plt.loglog(n,L2norms_dP_dr,linewidth=lw)
    plt.loglog(n,10**N1(np.log10(n),P_popt1),color="orange",linestyle=":",linewidth=lw)
    plt.loglog(n,10**N2(np.log10(n),P_popt2),color="red",linestyle=":",linewidth=lw)
    plt.xlabel("Number of Cells")
    plt.ylabel("Pressure Gradient L$^2$ Norm")
    #plt.title("Pressure Gradient L$^2$ Norm vs. Number of Cells")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    
    plt.figure()
    plt.loglog(n, abs(min_dP_drs),linewidth=lw)
    plt.xlabel("Number of Cells")
    plt.ylabel("Minimum Pressure Gradient")
    #plt.title("Minimum Pressure Gradient vs. Number of Cells")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    
    plt.figure()
    plt.loglog(n, abs((r0_high-r0s)/r0_high),color="blue", linewidth=lw)
    plt.loglog(n,10**N1(np.log10(n),r_popt1),color="orange",linestyle="--",linewidth=lw)
    plt.xlabel("Number of Cells")
    plt.ylabel("$r_{inner}$ Relative Difference")
    #plt.title("$r_{inner}$ Relative Difference vs. Number of Cells")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    
    return abs((r0s[-1]-r0s[:-1])/r0s[-1])



# Freedman Opacity Convergence Function

def ConvergenceFreO(data):
    
    r0s, n = data
    r0_high = r0s[-1]
    r0s = r0s[:-1]
    
    Figsize, lw, Ticksize, Labelsize = [7, 2, 20, 20]
    
    plt.rc("figure", figsize=(Figsize,Figsize))
    plt.rc("xtick", labelsize=Ticksize)
    plt.rc("ytick", labelsize=Ticksize)
    plt.rc("axes", labelsize=Labelsize)
    
    # Convergence Plots
    plt.figure()
    plt.loglog(n, abs((r0_high-r0s)/r0_high), linewidth=lw)
    plt.xlabel("Number of Cells")
    plt.ylabel("$r_{inner}$ Relative Difference")
    #plt.title("$r_{inner}$ Relative Difference vs. Number of Cells")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()




# Combined Opacity Convergence Function

def ConvergenceComO(data):
    
    r0s, n = data
    r0_high = r0s[-1]
    r0s = r0s[:-1]
    
    Figsize, lw, Ticksize, Labelsize = [7, 2, 20, 20]
    
    plt.rc("figure", figsize=(Figsize,Figsize))
    plt.rc("xtick", labelsize=Ticksize)
    plt.rc("ytick", labelsize=Ticksize)
    plt.rc("axes", labelsize=Labelsize)
    
    # Convergence Plots
    plt.figure()
    plt.loglog(n, abs((r0_high-r0s)/r0_high), linewidth=lw)
    #plt.loglog(n[:-1], abs((r0s[-1]-r0s[:-1])/r0s[-1]), linewidth=lw)
    plt.xlabel("Number of Cells")
    plt.ylabel("$r_{inner}$ Relative Difference")
    #plt.title("$r_{inner}$ Relative Difference vs. Number of Cells")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()




# Sigma Convergence Function

def ConvergenceSigma(data):
    
    S0s, n = data
    S0_high = S0s[-1]
    S0s = S0s[:-1]
    
    S_popt1, S_pcov1 = curve_fit(N1, np.log10(n), np.log10(abs((S0_high-S0s)/S0_high)))
    
    Figsize, lw, Ticksize, Labelsize = [7, 2, 20, 20]
    
    plt.rc("figure", figsize=(Figsize,Figsize))
    plt.rc("xtick", labelsize=Ticksize)
    plt.rc("ytick", labelsize=Ticksize)
    plt.rc("axes", labelsize=Labelsize)
    
    # Convergence Plots
    plt.figure()
    plt.loglog(n, abs((S0_high-S0s)/S0_high), color="blue", linewidth=lw)
    plt.loglog(n,10**N1(np.log10(n),S_popt1),color="orange",linestyle="--",linewidth=lw)
    #plt.loglog(n[:-1], abs((r0s[-1]-r0s[:-1])/r0s[-1]), linewidth=lw)
    plt.xlabel("Number of Cells")
    plt.ylabel("$\Sigma_{inner}$ Relative Difference")
    #plt.title("$\Sigma_{inner}$ Relative Difference vs. Number of Cells")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()



# Sigma Convergence Function

def ConvergenceMdotMax(data):
    
    Mdot_rel_diffs, n = data
    Mdot_rel_diffs = Mdot_rel_diffs
    
    M_popt1, M_pcov1 = curve_fit(N1, np.log10(n), np.log10(Mdot_rel_diffs))

    Figsize, lw, Ticksize, Labelsize = [7, 2, 20, 20]
    
    plt.rc("figure", figsize=(Figsize,Figsize))
    plt.rc("xtick", labelsize=Ticksize)
    plt.rc("ytick", labelsize=Ticksize)
    plt.rc("axes", labelsize=Labelsize)
    
    # Convergence Plots
    plt.figure()
    plt.loglog(n, Mdot_rel_diffs, color="blue", linewidth=lw)
    plt.loglog(n,10**N1(np.log10(n),M_popt1),color="orange",linestyle="--",linewidth=lw)
    plt.xlabel("Number of Cells")
    plt.ylabel("$\dot{M}_{lim}$ Relative Difference")
    #plt.title("$\dot{M}_{max}$ Relative Difference vs. Number of Cells")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()



# Luminosity Finder IsoO Convergence Function

def ConvergenceLum(data):
    
    L_rel_diff, r_rel_diff, n_is = data
    
    Figsize, lw, Ticksize, Labelsize = [7, 2, 20, 20]
    
    plt.rc("figure", figsize=(Figsize,Figsize))
    plt.rc("xtick", labelsize=Ticksize)
    plt.rc("ytick", labelsize=Ticksize)
    plt.rc("axes", labelsize=Labelsize)
    
    plt.figure()
    plt.plot(n_is, L_rel_diff, linewidth=lw, color="blue")
    plt.yscale("log")
    plt.xlabel("Iteration Number")
    plt.ylabel("L Relative Difference")
    #plt.title("L Relative Difference Iteration")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    
    plt.figure()
    plt.plot(n_is, r_rel_diff, linewidth=lw, color="blue")
    plt.yscale("log")
    plt.xlabel("Iteration Number")
    plt.ylabel("$\epsilon_{r}$")
    #plt.title("$r_{core}$ Relative Difference Iteration")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()