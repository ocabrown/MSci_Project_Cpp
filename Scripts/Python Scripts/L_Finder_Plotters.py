import matplotlib.pyplot as plt



def L_Finder_Plot(data):
    
    R_J = 6.9911e9                  # Jupiter radius [cm]
    
    L, Ls, rc, rs, n_is = data
    
    print(L)
    
    Figsize, lw, Ticksize, Labelsize, Legendfontsize = [7, 2, 20, 20, 15]
    
    plt.rc("figure", figsize=(Figsize,Figsize))
    plt.rc("xtick", labelsize=Ticksize)
    plt.rc("ytick", labelsize=Ticksize)
    plt.rc("axes", labelsize=Labelsize)
    
    plt.figure()
    plt.plot(n_is, Ls/1e-6, linewidth=lw, color="blue")
    plt.xlabel("Iteration Number")
    plt.ylabel("Luminosity / 10$^{-6}$ L$_{\odot}$")
    #plt.title("L Iteration")
    ax = plt.gca()
    ax.set_box_aspect(1)
    plt.tight_layout()
    plt.show()
    
    plt.figure()
    plt.plot(n_is, rs/R_J, linewidth=lw, color="blue")
    plt.axhline(rc/R_J, linestyle="--", color="black")
    plt.xlabel("Iteration Number")
    plt.ylabel("$r_{inner}$ / $R_{J}$")
    #plt.title("r Iteration")
    ax = plt.gca()
    ax.set_box_aspect(1)
    ax = plt.gca()
    ax.set_yticks([0,0.1,0.2,rc/R_J,0.3,0.4])
    ax.set_yticklabels([0,0.1,0.2,"$r_{core}$",0.3,0.4])
    plt.tight_layout()
    plt.show()
    
    fig, ax1 = plt.subplots()
    color1 = 'red'
    color2 = 'orange'
    
    ax1.set_xlabel("Iteration Number")
    ax1.set_ylabel("Luminosity / 10$^{-6}$ $L_{\odot}$", color=color1)
    ax1.plot(n_is, Ls/1e-6, color=color1, linewidth=lw)
    ax1.tick_params(axis='y')
    
    ax2 = ax1.twinx()
    ax2.set_ylabel("$r_{inner}$ / $R_{J}$", color=color2)
    ax2.plot(n_is, rs/R_J, color=color2, linewidth=lw)
    ax2.axhline(rc/R_J, linestyle="--", color="black", linewidth=lw)
    ax2.tick_params(axis='y')
    fig.tight_layout()