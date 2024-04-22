#include "./../Header_Files/Linspace.hpp"

#include <vector>





std::vector<double> linspace(double start_in, double end_in, int num_in)
{
    const double dx = (end_in - start_in) / (num_in - 1);
    std::vector<double> x(num_in);
    int iter = 0;
    std::generate(x.begin(), x.end(), [&] { return start_in + (iter++) * dx; });
    return x;
}


std::vector<double> linspace(double start_in, double end_in, int num_in, int max_stencil)
{
    const double dx = (end_in - start_in) / (num_in - 1.);
    // initialize x coords
    std::vector<double> x(num_in + 2 * max_stencil);
    int iter = -max_stencil;
    std::generate(x.begin(), x.end(), [&] { return start_in + (iter++) * dx; });
    return x;
}


std::vector<double> geomspace(double start, double stop, int num, bool endpoint)
{
    std::vector<double> result;
    
    double ratio = pow(stop / start, 1.0 / (num - 1));
    
    for (int i = 0; i < num; ++i)
    {
        double value = start * pow(ratio, i);
        result.push_back(value);
    }
    
    if (!endpoint && num > 1)
    {
        // Exclude the endpoint
        result.pop_back();
    }
    
    return result;
}
