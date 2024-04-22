#include "./../Header_Files/Data_Output.hpp"

#include <vector>
#include <string>
#include <fstream>





void Read_Out(std::vector< std::vector<double> > data, std::string data_type)
{
    
    std::vector<std::string> filenames;
    
    if (data_type.substr(0,11) == "Convergence")
    {
        data_type = data_type.substr(11);
        
        if (data_type == "IsoT")
        {
            filenames = {"L2norms_p.txt", "L2norms_dP_dr.txt", "min_dP_drs.txt", "nums.txt"};
        }
        else if (data_type == "IsoD")
        {
            filenames = {"L2norms_T.txt", "L2norms_dP_dr.txt", "min_dP_drs.txt", "L2norms_dP_dr_lower.txt", "L2norms_dP_dr_upper.txt", "nums.txt"};
        }
        else if (data_type == "IsoO")
        {
            filenames = {"r0s.txt", "L2norms_dP_dr.txt", "min_dP_drs.txt", "nums.txt"};
        }
        else if (data_type == "FreO")
        {
            filenames = {"r0s.txt", "nums.txt"};
        }
        else if (data_type == "ComO")
        {
            filenames = {"r0s.txt", "nums.txt"};
        }
        else if (data_type == "_Sigma")
        {
            data_type = "Max_Acc_Rate";
            filenames = {"S0s.txt", "nums_Sigma.txt"};
        }
        else if (data_type == "_MdotMax")
        {
            data_type = "Max_Acc_Rate";
            filenames = {"Mdot_rel_diffs.txt", "nums_MdotMax.txt"};
        }
        else if (data_type == "_IsoO_LF")
        {
            data_type = "IsoO/IsoO_LF";
            filenames = {"L_rel_diff.txt", "r_rel_diff.txt", "nums_LF.txt"};
        }
        else if (data_type == "_FreO_LF")
        {
            data_type = "FreO/FreO_LF";
            filenames = {"L_rel_diff.txt", "r_rel_diff.txt", "nums_LF.txt"};
        }
        else if (data_type == "_ComO_LF")
        {
            data_type = "ComO/ComO_LF";
            filenames = {"L_rel_diff.txt", "r_rel_diff.txt", "nums_LF.txt"};
        }
    }
    
    else
    {
        if (data_type == "IsoT")
        {
            filenames = {"mass.txt", "radi.txt", "dens.txt", "temp.txt", "pres.txt", "entr.txt", "dens_ana.txt", "dP_dr_num.txt", "dP_dr_ana.txt"};
        }
        else if (data_type == "IsoD")
        {
            filenames = {"mass.txt", "radi.txt", "dens.txt", "temp.txt", "pres.txt", "entr.txt", "temp_ana.txt", "dP_dr_num.txt", "dP_dr_ana.txt"};
        }
        else if (data_type == "IsoO")
        {
            filenames = {"mass.txt", "radi.txt", "dens.txt", "temp.txt", "opac.txt", "pres.txt", "entr.txt", "dP_dr_num.txt", "dP_dr_ana.txt"};
        }
        else if (data_type == "FreO")
        {
            filenames = {"mass.txt", "radi.txt", "dens.txt", "temp.txt", "opac.txt", "pres.txt", "entr.txt"};
        }
        else if (data_type == "ComO")
        {
            filenames = {"mass.txt", "radi.txt", "dens.txt", "temp.txt", "opac.txt", "pres.txt", "entr.txt"};
        }
    }
    
    data_type += "/";
    
    for (int i = 0; i < data.size(); i++)
    {
        std::ofstream output_file("/Desktop/Data_Outputs/" + data_type + filenames[i]);
        for (int j = 0; j < data[i].size(); j++)
        {
            if (j == data[i].size() - 1)
            {
                output_file << std::setprecision(10) << data[i][j];
            }
            else
            {
                output_file << std::setprecision(10) << data[i][j] << ",";
            }
        }
        output_file.close();
    }

}



void Read_Out(std::vector< std::pair<int,int> > data, std::string data_type)
{
    
    std::vector<std::string> filenames;
    filenames = {"boundary_i.txt", "boundary_rad_or_conv.txt"};
    
    data_type += "/";
    
    std::ofstream output_file1("/Desktop/Data_Outputs/" + data_type + filenames[0]);
    for (int i = 0; i < data.size(); i++)
    {
        if (i == data.size() - 1)
        {
            output_file1 << data[i].first;
        }
        else
        {
            output_file1 << data[i].first << ",";
        }
    }
    output_file1.close();
    
    std::ofstream output_file2("/Desktop/Data_Outputs/" + data_type + filenames[1]);
    for (int i = 0; i < data.size(); i++)
    {
        if (i == data.size() - 1)
        {
            output_file2 << data[i].second;
        }
        else
        {
            output_file2 << data[i].second << ",";
        }
    }
    output_file2.close();
}



void Read_Out(std::pair< std::vector< std::vector<double> >, std::vector<bool> > data, std::string data_type)
{
    std::vector<std::string> filenames;
    
    if (data.second[0])
    {
        filenames = {"Mdot_max.txt", "Mdot_max_ana.txt", "extras.txt"};
    }
    else if (data.second[1])
    {
        filenames = {"Sigma.txt", "r_d.txt", "extras.txt"};
    }
    else
    {
        filenames = {"Mdot_max.txt", "Sigma.txt", "r_d.txt", "extras.txt"};
    }
    
    data_type += "/";
    
    for (int i = 0; i < data.first.size(); i++)
    {
        std::ofstream output_file("/Desktop/Data_Outputs/" + data_type + filenames[i]);
        for (int j = 0; j < data.first[i].size(); j++)
        {
            if (j == data.first[i].size() - 1)
            {
                output_file << std::setprecision(10) << data.first[i][j];
            }
            else
            {
                output_file << std::setprecision(10) << data.first[i][j] << ",";
            }
        }
        output_file.close();
    }
}



void Read_Out(std::vector< std::pair< std::vector< std::vector<double> >, std::vector<bool> > > data, std::string data_type)
{
    std::vector<std::string> filenames = {"Sigma0.txt", "Sigma1.txt", "Sigma2.txt", "Sigma3.txt", "r_d0.txt", "r_d1.txt", "r_d2.txt", "r_d3.txt"};
    
    data_type += "/";
    
    for (int i = 0; i < data.size(); i++)
    {
        std::ofstream output_file1("/Desktop/Data_Outputs/Max_Acc_Rate/" + data_type + filenames[i]);
        for (int j = 0; j < data[i].first[1].size(); j++)
        {
            if (j == data[i].first[1].size() - 1)
            {
                output_file1 << std::setprecision(10) << data[i].first[1][j];
            }
            else
            {
                output_file1 << std::setprecision(10) << data[i].first[1][j] << ",";
            }
        }
        output_file1.close();
        
        std::ofstream output_file2("/Desktop/Data_Outputs/Max_Acc_Rate/" + data_type + filenames[i+4]);
        for (int j = 0; j < data[i].first[2].size(); j++)
        {
            if (j == data[i].first[2].size() - 1)
            {
                output_file2 << std::setprecision(10) << data[i].first[2][j];
            }
            else
            {
                output_file2 << std::setprecision(10) << data[i].first[2][j] << ",";
            }
        }
        output_file2.close();
        
    }
}



void Read_Out(std::pair< std::vector< std::vector<double> >, std::vector< std::vector<double> > > data, std::string data_type)
{
    data_type = data_type.substr(0,4) + "/" + data_type + "/";
    
    std::vector<std::string> filenames = {"L.txt", "Ls.txt", "rc.txt", "rs.txt", "n_is.txt"};
    
    for (int i = 0; i < data.first.size(); i++)
    {
        std::ofstream output_file("/Desktop/Data_Outputs/" + data_type + filenames[i]);
        for (int j = 0; j < data.first[i].size(); j++)
        {
            if (j == data.first[i].size() - 1)
            {
                output_file << std::setprecision(10) << data.first[i][j];
            }
            else
            {
                output_file << std::setprecision(10) << data.first[i][j] << ",";
            }
        }
        output_file.close();
    }
}



void Read_Out(std::pair< std::vector< std::vector<double> >, std::vector< std::vector< std::vector<double> > > > data, std::string data_type)
{
    data_type = data_type.substr(0,4) + "/" + data_type + "/";
    std::vector<std::string> filenames = {"L.txt", "M.txt", "t.txt", "M_acc.txt", "M_acc_max.txt", "M_acc_true.txt"};
    
    for (int i = 0; i < data.first.size(); i++)
    {
        std::ofstream output_file("/Desktop/Data_Outputs/" + data_type + filenames[i]);
        for (int j = 0; j < data.first[i].size(); j++)
        {
            if (j == data.first[i].size() - 1)
            {
                output_file << std::setprecision(10) << data.first[i][j];
            }
            else
            {
                output_file << std::setprecision(10) << data.first[i][j] << ",";
            }
        }
        output_file.close();
    }
    
    data_type += "Varis/";
    std::vector<std::string> filenamesVaris = {"m.txt", "r.txt", "p.txt", "T.txt", "k.txt", "P.txt", "s.txt", "dP_dr_num.txt", "dP_dr_ana.txt"};
    
    for (int n = 0; n < data.second[0].size(); n++)
    {
        std::ofstream output_file("/Desktop/Data_Outputs/" + data_type + filenamesVaris[n]);
        for (int i = 0; i < data.second.size(); i++)
        {
            for (int j = 0; j < data.second[i][n].size(); j++)
            {
                if (j == data.second[i][n].size() - 1)
                {
                    output_file << std::setprecision(10) << data.second[i][n][j];
                }
                else
                {
                    output_file << std::setprecision(10) << data.second[i][n][j] << ",";
                }
            }
        }
        output_file.close();
    }
    
    std::ofstream output_file("/Desktop/Data_Outputs/" + data_type + "nums.txt");
    output_file << data.second[0][0].size() << "," << data.second.size();
    output_file.close();
}



void Read_Out(std::vector< std::pair< std::vector< std::vector<double> >, std::vector< std::vector< std::vector<double> > > > > data, std::string data_type)
{
    std::vector< std::vector<std::string> > filenames;
    
    if (data_type.substr(7) == "ABC")
    {
        filenames = {{"LA.txt", "MA.txt", "tA.txt", "M_accA.txt", "M_acc_maxA.txt", "M_acc_trueA.txt"}, {"LB.txt", "MB.txt", "tB.txt", "M_accB.txt", "M_acc_maxB.txt", "M_acc_trueB.txt"}, {"LC.txt", "MC.txt", "tC.txt", "M_accC.txt", "M_acc_maxC.txt", "M_acc_trueC.txt"}};
    }
    else if (data_type.substr(7) == "DEF")
    {
        filenames = {{"LD.txt", "MD.txt", "tD.txt", "M_accD.txt", "M_acc_maxD.txt", "M_acc_trueD.txt"}, {"LE.txt", "ME.txt", "tE.txt", "M_accE.txt", "M_acc_maxE.txt", "M_acc_trueE.txt"}, {"LF.txt", "MF.txt", "tF.txt", "M_accF.txt", "M_acc_maxF.txt", "M_acc_trueF.txt"}};
    }
    else if (data_type.substr(7) == "GHI")
    {
        filenames = {{"LG.txt", "MG.txt", "tG.txt", "M_accG.txt", "M_acc_maxG.txt", "M_acc_trueG.txt"}, {"LH.txt", "MH.txt", "tH.txt", "M_accH.txt", "M_acc_maxH.txt", "M_acc_trueH.txt"}, {"LI.txt", "MI.txt", "tI.txt", "M_accI.txt", "M_acc_maxI.txt", "M_acc_trueI.txt"}};
    }
    
    data_type = data_type.substr(0,4) + "/" + data_type.substr(0,6) + "/" + data_type + "/";
    
    for (int i = 0; i < data.size(); i++)
    {
        for (int j = 0; j < data[i].first.size(); j++)
        {
            std::ofstream output_file("/Desktop/Data_Outputs/" + data_type + filenames[i][j]);
            for (int k = 0; k < data[i].first[j].size(); k++)
            {
                if (k == data[i].first[j].size() - 1)
                {
                    output_file << std::setprecision(10) << data[i].first[j][k];
                }
                else
                {
                    output_file << std::setprecision(10) << data[i].first[j][k] << ",";
                }
            }
            output_file.close();
        }
    }
    
    
    
}
