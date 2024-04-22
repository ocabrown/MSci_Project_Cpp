#ifndef Data_Output_hpp
#define Data_Output_hpp

#include <vector>
#include <string>

void Read_Out(std::vector< std::vector<double> > data, std::string data_type);
void Read_Out(std::vector< std::pair<int,int> > data, std::string data_type);
void Read_Out(std::pair< std::vector< std::vector<double> >, std::vector<bool> > data, std::string data_type);
void Read_Out(std::vector< std::pair< std::vector< std::vector<double> >, std::vector<bool> > > data, std::string data_type);
void Read_Out(std::pair< std::vector< std::vector<double> >, std::vector< std::vector<double> > > data, std::string data_type);
void Read_Out(std::pair< std::vector< std::vector<double> >, std::vector< std::vector< std::vector<double> > > > data, std::string data_type);
void Read_Out(std::vector< std::pair< std::vector< std::vector<double> >, std::vector< std::vector< std::vector<double> > > > > data, std::string data_type);

#endif
