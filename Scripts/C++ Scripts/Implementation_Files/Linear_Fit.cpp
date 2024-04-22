#include "./../Header_Files/Linear_Fit.hpp"

#include <array>
#include <vector>





std::array<double, 2> lin_fit(std::vector<double> x, std::vector<double> y)   // x-values, y-values
{
    int n_int = (int)x.size();      // array size
    double n = (double)n_int;
    
    double xsum=0., x2sum=0., ysum=0., xysum=0.;                //variables for sums/sigma of xi,yi,xi^2,xiyi etc
    for (int i = 0; i < n_int; i++)
    {
        xsum += x[i];                           //calculate sigma(xi)
        ysum += y[i];                           //calculate sigma(yi)
        x2sum += pow(x[i], 2.);                 //calculate sigma(x^2i)
        xysum += x[i]*y[i];                     //calculate sigma(xi*yi)
    }
    
    double m = ((n * xysum) - (xsum * ysum)) / ((n * x2sum) - (xsum * xsum));       //calculate slope
    double c = ((x2sum * ysum) - (xsum * xysum)) / ((x2sum * n) - (xsum * xsum));   //calculate intercept
    
    std::array<double, 2> popt = {m, c};
    
    return popt;
}
