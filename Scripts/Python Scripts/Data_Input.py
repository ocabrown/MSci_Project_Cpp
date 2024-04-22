import numpy as np



def Data_Extractor(Data, Type):
    
    X_Data = 10**Data[:,0]
    Y_Data = Data[:,1]  
    
    if Type == "Density":
        Y_Data = 10**Y_Data[:]
    elif Type == "Opacity":
        Y_Data = 10**Y_Data[:]
    elif Type == "Temperature":
        Y_Data = Y_Data
    
    return X_Data, Y_Data



def Data_Preparer(Filenames):
    
    preamble = "/Data/Movshovitz_Data/"
    DensityData = np.loadtxt(preamble+str(Filenames[0]), delimiter=",")
    OpacityData = np.loadtxt(preamble+str(Filenames[1]), delimiter=",")
    TemperatureData = np.loadtxt(preamble+str(Filenames[2]), delimiter=",")
    RD, D = Data_Extractor(DensityData, "Density")
    RO, O = Data_Extractor(OpacityData, "Opacity")
    RT, T = Data_Extractor(TemperatureData, "Temperature")
    
    return [[RD,D], [RO,O], [RT,T]]



def Data_Extractor_Sigma(Filename):
    preamble = "/Data/Crida_Data/"
    SigmaData = np.loadtxt(preamble+str(Filename), delimiter=",")
    X_Data = SigmaData[:,0]
    Y_Data = SigmaData[:,1]
    
    return X_Data, Y_Data



def Read_In_Var(Folder_Name):
    
    Folder_Name += "/"
    
    preamble = "/Desktop/Data_Outputs/"
    
    m = np.loadtxt(preamble+Folder_Name+"mass.txt", delimiter=",")
    r = np.loadtxt(preamble+Folder_Name+"radi.txt", delimiter=",")
    p = np.loadtxt(preamble+Folder_Name+"dens.txt", delimiter=",")
    T = np.loadtxt(preamble+Folder_Name+"temp.txt", delimiter=",")
    P = np.loadtxt(preamble+Folder_Name+"pres.txt", delimiter=",")
    s = np.loadtxt(preamble+Folder_Name+"entr.txt", delimiter=",")
    
    if Folder_Name == "IsoT/":
        p_ana = np.loadtxt(preamble+Folder_Name+"dens_ana.txt", delimiter=",")
        dP_dr_num = np.loadtxt(preamble+Folder_Name+"dP_dr_num.txt", delimiter=",")
        dP_dr_ana = np.loadtxt(preamble+Folder_Name+"dP_dr_ana.txt", delimiter=",")
        data = [m, r, p, T, P, s, p_ana, dP_dr_num, dP_dr_ana]
        
    elif Folder_Name == "IsoD/":
        T_ana = np.loadtxt(preamble+Folder_Name+"temp_ana.txt", delimiter=",")
        dP_dr_num = np.loadtxt(preamble+Folder_Name+"dP_dr_num.txt", delimiter=",")
        dP_dr_ana = np.loadtxt(preamble+Folder_Name+"dP_dr_ana.txt", delimiter=",")
        data = [m, r, p, T, P, s, T_ana, dP_dr_num, dP_dr_ana]
    
    elif Folder_Name == "IsoO/":
        k = np.loadtxt(preamble+Folder_Name+"opac.txt", delimiter=",")
        dP_dr_num = np.loadtxt(preamble+Folder_Name+"dP_dr_num.txt", delimiter=",")
        dP_dr_ana = np.loadtxt(preamble+Folder_Name+"dP_dr_ana.txt", delimiter=",")
        boundary_ind = np.loadtxt(preamble+Folder_Name+"boundary_i.txt", delimiter=",", dtype="int")
        boundary_rad_or_conv = np.loadtxt(preamble+Folder_Name+"boundary_rad_or_conv.txt", delimiter=",", dtype="int")
        boundary_i = []
        #for i in range(len(boundary_ind)):
        #    boundary_i.append([boundary_ind[i],boundary_rad_or_conv[i]])
        data = [m, r, p, T, k, P, s, dP_dr_num, dP_dr_ana, boundary_i]
    
    elif Folder_Name == "FreO/":
        k = np.loadtxt(preamble+Folder_Name+"opac.txt", delimiter=",")
        boundary_ind = np.loadtxt(preamble+Folder_Name+"boundary_i.txt", delimiter=",", dtype="int")
        boundary_rad_or_conv = np.loadtxt(preamble+Folder_Name+"boundary_rad_or_conv.txt", delimiter=",", dtype="int")
        boundary_i = []
        for i in range(len(boundary_ind)):
            boundary_i.append([boundary_ind[i],boundary_rad_or_conv[i]])
        boundary_i = boundary_i[1:]
        data = [m, r, p, T, k, P, s, boundary_i]
    
    elif Folder_Name == "ComO/":
        k = np.loadtxt(preamble+Folder_Name+"opac.txt", delimiter=",")
        boundary_ind = np.loadtxt(preamble+Folder_Name+"boundary_i.txt", delimiter=",", dtype="int")
        boundary_rad_or_conv = np.loadtxt(preamble+Folder_Name+"boundary_rad_or_conv.txt", delimiter=",", dtype="int")
        boundary_i = []
        for i in range(len(boundary_ind)):
            boundary_i.append([boundary_ind[i],boundary_rad_or_conv[i]])
        boundary_i = boundary_i[1:]
        data = [m, r, p, T, k, P, s, boundary_i]
    
    return data



def Read_In_L2norm(Folder_Name):
    
    Folder_Name += "/"
    
    preamble = "/Desktop/Data_Outputs/"
    
    if Folder_Name == "IsoT/":
        L2norms_p = np.loadtxt(preamble+Folder_Name+"L2norms_p.txt", delimiter=",")
        L2norms_dP_dr = np.loadtxt(preamble+Folder_Name+"L2norms_dP_dr.txt", delimiter=",")
        min_dP_drs = np.loadtxt(preamble+Folder_Name+"min_dP_drs.txt", delimiter=",")
        nums = np.loadtxt(preamble+Folder_Name+"nums.txt", delimiter=",")
        data = [L2norms_p, L2norms_dP_dr, min_dP_drs, nums]
    
    elif Folder_Name == "IsoD/":
        L2norms_T = np.loadtxt(preamble+Folder_Name+"L2norms_T.txt", delimiter=",")
        L2norms_dP_dr = np.loadtxt(preamble+Folder_Name+"L2norms_dP_dr.txt", delimiter=",")
        min_dP_drs = np.loadtxt(preamble+Folder_Name+"min_dP_drs.txt", delimiter=",")
        L2norms_dP_dr_lower = np.loadtxt(preamble+Folder_Name+"L2norms_dP_dr_lower.txt", delimiter=",")
        L2norms_dP_dr_upper = np.loadtxt(preamble+Folder_Name+"L2norms_dP_dr_upper.txt", delimiter=",")
        nums = np.loadtxt(preamble+Folder_Name+"nums.txt", delimiter=",")
        data = [L2norms_T, L2norms_dP_dr, min_dP_drs, L2norms_dP_dr_lower, L2norms_dP_dr_upper, nums]
    
    elif Folder_Name == "IsoO/":
        r0s = np.loadtxt(preamble+Folder_Name+"r0s.txt", delimiter=",")
        L2norms_dP_dr = np.loadtxt(preamble+Folder_Name+"L2norms_dP_dr.txt", delimiter=",")
        min_dP_drs = np.loadtxt(preamble+Folder_Name+"min_dP_drs.txt", delimiter=",")
        nums = np.loadtxt(preamble+Folder_Name+"nums.txt", delimiter=",")
        data = [r0s, L2norms_dP_dr, min_dP_drs, nums]
    
    elif Folder_Name == "FreO/":
        r0s = np.loadtxt(preamble+Folder_Name+"r0s.txt", delimiter=",")
        nums = np.loadtxt(preamble+Folder_Name+"nums.txt", delimiter=",")
        data = [r0s, nums]
    
    elif Folder_Name == "ComO/":
        r0s = np.loadtxt(preamble+Folder_Name+"r0s.txt", delimiter=",")
        nums = np.loadtxt(preamble+Folder_Name+"nums.txt", delimiter=",")
        data = [r0s, nums]
        
    elif Folder_Name == "Max_Acc_Rate_Sigma/":
        Folder_Name = "Max_Acc_Rate/"
        S0s = np.loadtxt(preamble+Folder_Name+"S0s.txt", delimiter=",")
        nums = np.loadtxt(preamble+Folder_Name+"nums_Sigma.txt", delimiter=",")
        data = [S0s, nums]
    
    elif Folder_Name == "Max_Acc_Rate_MdotMax/":
        Folder_Name = "Max_Acc_Rate/"
        Mdot_rel_diffs = np.loadtxt(preamble+Folder_Name+"Mdot_rel_diffs.txt", delimiter=",")
        nums = np.loadtxt(preamble+Folder_Name+"nums_MdotMax.txt", delimiter=",")
        data = [Mdot_rel_diffs, nums]
    
    elif Folder_Name == "IsoO_LF/":
        Folder_Name = "IsoO/IsoO_LF/"
        L_rel_diff = np.loadtxt(preamble+Folder_Name+"L_rel_diff.txt", delimiter=",")
        r_rel_diff = np.loadtxt(preamble+Folder_Name+"r_rel_diff.txt", delimiter=",")
        n_is = np.loadtxt(preamble+Folder_Name+"n_is.txt", delimiter=",")
        data = [L_rel_diff, r_rel_diff, n_is]
    
    elif Folder_Name == "FreO_LF/":
        Folder_Name = "FreO/FreO_LF/"
        L_rel_diff = np.loadtxt(preamble+Folder_Name+"L_rel_diff.txt", delimiter=",")
        r_rel_diff = np.loadtxt(preamble+Folder_Name+"r_rel_diff.txt", delimiter=",")
        n_is = np.loadtxt(preamble+Folder_Name+"n_is.txt", delimiter=",")
        data = [L_rel_diff, r_rel_diff, n_is]
    
    elif Folder_Name == "ComO_LF/":
        Folder_Name = "ComO/ComO_LF/"
        L_rel_diff = np.loadtxt(preamble+Folder_Name+"L_rel_diff.txt", delimiter=",")
        r_rel_diff = np.loadtxt(preamble+Folder_Name+"r_rel_diff.txt", delimiter=",")
        n_is = np.loadtxt(preamble+Folder_Name+"n_is.txt", delimiter=",")
        data = [L_rel_diff, r_rel_diff, n_is]
    
    
    return data



def Read_In_Max_Acc(Folder_Name):
    
    Folder_Name += "/"
    
    if Folder_Name == "Max_Acc_Rate/":
        preamble = "/Desktop/Data_Outputs/"
        Sigma = np.loadtxt(preamble+Folder_Name+"Sigma.txt", delimiter=",")
        r_d = np.loadtxt(preamble+Folder_Name+"r_d.txt", delimiter=",")
        extras = np.loadtxt(preamble+Folder_Name+"extras.txt", delimiter=",")
        Mp, Rp, num_sig = extras
        data = [Sigma, r_d, Mp, Rp, num_sig]
    
    elif Folder_Name == "Max_Acc_Rate_Diff_Alpha/":
        preamble = "/Desktop/Data_Outputs/Max_Acc_Rate/"
        Sigma0 = np.loadtxt(preamble+Folder_Name+"Sigma0.txt", delimiter=",")
        Sigma1 = np.loadtxt(preamble+Folder_Name+"Sigma1.txt", delimiter=",")
        Sigma2 = np.loadtxt(preamble+Folder_Name+"Sigma2.txt", delimiter=",")
        Sigma3 = np.loadtxt(preamble+Folder_Name+"Sigma3.txt", delimiter=",")
        r_d0 = np.loadtxt(preamble+Folder_Name+"r_d0.txt", delimiter=",")
        r_d1 = np.loadtxt(preamble+Folder_Name+"r_d1.txt", delimiter=",")
        r_d2 = np.loadtxt(preamble+Folder_Name+"r_d2.txt", delimiter=",")
        r_d3 = np.loadtxt(preamble+Folder_Name+"r_d3.txt", delimiter=",")
        data = [Sigma0, Sigma1, Sigma2, Sigma3, r_d0, r_d1, r_d2, r_d3]
    
    return data



def Read_In_LF(Folder_Name):
    
    preamble = "/Desktop/Data_Outputs/"
    
    if Folder_Name == "IsoO_LF":
        Folder_Name = "IsoO/IsoO_LF/"
    elif Folder_Name == "FreO_LF":
        Folder_Name = "FreO/FreO_LF/"
    elif Folder_Name == "ComO_LF":
        Folder_Name = "ComO/ComO_LF/"
    
    L = np.loadtxt(preamble+Folder_Name+"L.txt", delimiter=",")
    Ls = np.loadtxt(preamble+Folder_Name+"Ls.txt", delimiter=",")
    rc = np.loadtxt(preamble+Folder_Name+"rc.txt", delimiter=",")
    rs = np.loadtxt(preamble+Folder_Name+"rs.txt", delimiter=",")
    n_is = np.loadtxt(preamble+Folder_Name+"n_is.txt", delimiter=",")
    data = [L, Ls, rc, rs, n_is]
    
    return data



def Read_In_E(Folder_Name):
    
    preamble = "/Desktop/Data_Outputs/"
    
    if Folder_Name == "IsoO_E":
        Folder_Name = "IsoO/IsoO_E/"
        L = np.loadtxt(preamble+Folder_Name+"L.txt", delimiter=",")
        M = np.loadtxt(preamble+Folder_Name+"M.txt", delimiter=",")
        t = np.loadtxt(preamble+Folder_Name+"t.txt", delimiter=",")
        M_acc = np.loadtxt(preamble+Folder_Name+"M_acc.txt", delimiter=",")
        M_acc_max = np.loadtxt(preamble+Folder_Name+"M_acc_max.txt", delimiter=",")
        M_acc_true = np.loadtxt(preamble+Folder_Name+"M_acc_true.txt", delimiter=",")
        Mdots = np.array([M_acc, M_acc_max, M_acc_true])
        
        Folder_Name = "IsoO/IsoO_E/Varis/"
        num_i = np.loadtxt(preamble+Folder_Name+"nums.txt", delimiter=",")[0]
        num_t = np.loadtxt(preamble+Folder_Name+"nums.txt", delimiter=",")[1]
        ms = np.loadtxt(preamble+Folder_Name+"m.txt", delimiter=",")
        rs = np.loadtxt(preamble+Folder_Name+"r.txt", delimiter=",")
        ps = np.loadtxt(preamble+Folder_Name+"p.txt", delimiter=",")
        ks = np.loadtxt(preamble+Folder_Name+"k.txt", delimiter=",")
        Ts = np.loadtxt(preamble+Folder_Name+"T.txt", delimiter=",")
        Ps = np.loadtxt(preamble+Folder_Name+"P.txt", delimiter=",")
        ss = np.loadtxt(preamble+Folder_Name+"s.txt", delimiter=",")
        dP_dr_nums = np.loadtxt(preamble+Folder_Name+"dP_dr_num.txt", delimiter=",")
        dP_dr_anas = np.loadtxt(preamble+Folder_Name+"dP_dr_ana.txt", delimiter=",")
        varis = []
        for i in range(int(num_t)):
            m = ms[int((i * num_i) - i) : int(((i+1) * num_i) - i)]
            r = rs[int((i * num_i) - i) : int(((i+1) * num_i) - i)]
            p = ps[int((i * num_i) - i) : int(((i+1) * num_i) - i)]
            k = ks[int((i * num_i) - i) : int(((i+1) * num_i) - i)]
            T = Ts[int((i * num_i) - i) : int(((i+1) * num_i) - i)]
            P = Ps[int((i * num_i) - i) : int(((i+1) * num_i) - i)]
            s = ss[int((i * num_i) - i) : int(((i+1) * num_i) - i)]
            dP_dr_num = dP_dr_nums[int((i * num_i) - i) : int(((i+1) * num_i) - i)]
            dP_dr_ana = dP_dr_anas[int((i * num_i) - i) : int(((i+1) * num_i) - i)]
            vari = np.array([m, r, p, k, T, P, s, dP_dr_num, dP_dr_ana])
            varis.append(vari)
        varis = np.array(varis)
        EIO = [L, M, t, varis, Mdots]
        data = [EIO]
    
    elif Folder_Name == "IsoO_E_ABC":
        Folder_Name = "IsoO/IsoO_E/IsoO_E_ABC/"
        LA = np.loadtxt(preamble+Folder_Name+"LA.txt", delimiter=",")
        MA = np.loadtxt(preamble+Folder_Name+"MA.txt", delimiter=",")
        tA = np.loadtxt(preamble+Folder_Name+"tA.txt", delimiter=",")
        M_accA = np.loadtxt(preamble+Folder_Name+"M_accA.txt", delimiter=",")
        M_acc_maxA = np.loadtxt(preamble+Folder_Name+"M_acc_maxA.txt", delimiter=",")
        M_acc_trueA = np.loadtxt(preamble+Folder_Name+"M_acc_trueA.txt", delimiter=",")
        MdotsA = np.array([M_accA, M_acc_maxA, M_acc_trueA])
        varisA = np.array([])
        EIOA = [LA, MA, tA, varisA, MdotsA]
        LB = np.loadtxt(preamble+Folder_Name+"LB.txt", delimiter=",")
        MB = np.loadtxt(preamble+Folder_Name+"MB.txt", delimiter=",")
        tB = np.loadtxt(preamble+Folder_Name+"tB.txt", delimiter=",")
        M_accB = np.loadtxt(preamble+Folder_Name+"M_accB.txt", delimiter=",")
        M_acc_maxB = np.loadtxt(preamble+Folder_Name+"M_acc_maxB.txt", delimiter=",")
        M_acc_trueB = np.loadtxt(preamble+Folder_Name+"M_acc_trueB.txt", delimiter=",")
        MdotsB = np.array([M_accB, M_acc_maxB, M_acc_trueB])
        varisB = np.array([])
        EIOB = [LB, MB, tB, varisB, MdotsB]
        LC = np.loadtxt(preamble+Folder_Name+"LC.txt", delimiter=",")
        MC = np.loadtxt(preamble+Folder_Name+"MC.txt", delimiter=",")
        tC = np.loadtxt(preamble+Folder_Name+"tC.txt", delimiter=",")
        M_accC = np.loadtxt(preamble+Folder_Name+"M_accC.txt", delimiter=",")
        M_acc_maxC = np.loadtxt(preamble+Folder_Name+"M_acc_maxC.txt", delimiter=",")
        M_acc_trueC = np.loadtxt(preamble+Folder_Name+"M_acc_trueC.txt", delimiter=",")
        MdotsC = np.array([M_accC, M_acc_maxC, M_acc_trueC])
        varisC = np.array([])
        EIOC = [LC, MC, tC, varisC, MdotsC]
        data = [EIOA, EIOB, EIOC]
    
    elif Folder_Name == "IsoO_E_DEF":
        Folder_Name = "IsoO/IsoO_E/IsoO_E_DEF/"
        LD = np.loadtxt(preamble+Folder_Name+"LD.txt", delimiter=",")
        MD = np.loadtxt(preamble+Folder_Name+"MD.txt", delimiter=",")
        tD = np.loadtxt(preamble+Folder_Name+"tD.txt", delimiter=",")
        M_accD = np.loadtxt(preamble+Folder_Name+"M_accD.txt", delimiter=",")
        M_acc_maxD = np.loadtxt(preamble+Folder_Name+"M_acc_maxD.txt", delimiter=",")
        M_acc_trueD = np.loadtxt(preamble+Folder_Name+"M_acc_trueD.txt", delimiter=",")
        MdotsD = np.array([M_accD, M_acc_maxD, M_acc_trueD])
        varisD = np.array([])
        EIOD = [LD, MD, tD, varisD, MdotsD]
        LE = np.loadtxt(preamble+Folder_Name+"LE.txt", delimiter=",")
        ME = np.loadtxt(preamble+Folder_Name+"ME.txt", delimiter=",")
        tE = np.loadtxt(preamble+Folder_Name+"tE.txt", delimiter=",")
        M_accE = np.loadtxt(preamble+Folder_Name+"M_accE.txt", delimiter=",")
        M_acc_maxE = np.loadtxt(preamble+Folder_Name+"M_acc_maxE.txt", delimiter=",")
        M_acc_trueE = np.loadtxt(preamble+Folder_Name+"M_acc_trueE.txt", delimiter=",")
        MdotsE = np.array([M_accE, M_acc_maxE, M_acc_trueE])
        varisE = np.array([])
        EIOE = [LE, ME, tE, varisE, MdotsE]
        LF = np.loadtxt(preamble+Folder_Name+"LF.txt", delimiter=",")
        MF = np.loadtxt(preamble+Folder_Name+"MF.txt", delimiter=",")
        tF = np.loadtxt(preamble+Folder_Name+"tF.txt", delimiter=",")
        M_accF = np.loadtxt(preamble+Folder_Name+"M_accF.txt", delimiter=",")
        M_acc_maxF = np.loadtxt(preamble+Folder_Name+"M_acc_maxF.txt", delimiter=",")
        M_acc_trueF = np.loadtxt(preamble+Folder_Name+"M_acc_truef.txt", delimiter=",")
        MdotsF = np.array([M_accF, M_acc_maxF, M_acc_trueF])
        varisF = np.array([])
        EIOF = [LF, MF, tF, varisF, MdotsF]
        data = [EIOD, EIOE, EIOF]
    
    elif Folder_Name == "IsoO_E_GHI":
        Folder_Name = "IsoO/IsoO_E/IsoO_E_GHI/"
        LG = np.loadtxt(preamble+Folder_Name+"LG.txt", delimiter=",")
        MG = np.loadtxt(preamble+Folder_Name+"MG.txt", delimiter=",")
        tG = np.loadtxt(preamble+Folder_Name+"tG.txt", delimiter=",")
        M_accG = np.loadtxt(preamble+Folder_Name+"M_accG.txt", delimiter=",")
        M_acc_maxG = np.loadtxt(preamble+Folder_Name+"M_acc_maxG.txt", delimiter=",")
        M_acc_trueG = np.loadtxt(preamble+Folder_Name+"M_acc_trueG.txt", delimiter=",")
        MdotsG = np.array([M_accG, M_acc_maxG, M_acc_trueG])
        varisG = np.array([])
        EIOG = [LG, MG, tG, varisG, MdotsG]
        LH = np.loadtxt(preamble+Folder_Name+"LH.txt", delimiter=",")
        MH = np.loadtxt(preamble+Folder_Name+"MH.txt", delimiter=",")
        tH = np.loadtxt(preamble+Folder_Name+"tH.txt", delimiter=",")
        M_accH = np.loadtxt(preamble+Folder_Name+"M_accH.txt", delimiter=",")
        M_acc_maxH = np.loadtxt(preamble+Folder_Name+"M_acc_maxH.txt", delimiter=",")
        M_acc_trueH = np.loadtxt(preamble+Folder_Name+"M_acc_trueH.txt", delimiter=",")
        MdotsH = np.array([M_accH, M_acc_maxH, M_acc_trueH])
        varisH = np.array([])
        EIOH = [LH, MH, tH, varisH, MdotsH]
        LI = np.loadtxt(preamble+Folder_Name+"LI.txt", delimiter=",")
        MI = np.loadtxt(preamble+Folder_Name+"MI.txt", delimiter=",")
        tI = np.loadtxt(preamble+Folder_Name+"tI.txt", delimiter=",")
        M_accI = np.loadtxt(preamble+Folder_Name+"M_accI.txt", delimiter=",")
        M_acc_maxI = np.loadtxt(preamble+Folder_Name+"M_acc_maxI.txt", delimiter=",")
        M_acc_trueI = np.loadtxt(preamble+Folder_Name+"M_acc_trueI.txt", delimiter=",")
        MdotsI = np.array([M_accI, M_acc_maxI, M_acc_trueI])
        varisI = np.array([])
        EIOI = [LI, MI, tI, varisI, MdotsI]
        data = [EIOG, EIOH, EIOI]
    
    return data

