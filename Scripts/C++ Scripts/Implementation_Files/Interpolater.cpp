#include "./../Header_Files/Interpolater.hpp"

#include <vector>
#include <cfloat>





int findNearestNeighbourIndex(double value, std::vector< double > &x)
{
    double dist = FLT_MAX;
    int idx = -1;
    for (int i = 0; i < x.size(); ++i) {
        double newDist = value - x[i];
        if (newDist > 0 && newDist < dist) {
            dist = newDist;
            idx = i;
        }
    }
    return idx;
}



std::vector< double > interp1(std::vector< double > &x, std::vector< double > &y, std::vector< double > &x_new)
{
    std::vector< double > y_new;
    y_new.reserve(x_new.size());

    std::vector< double > dx, dy, slope, intercept;
    dx.reserve(x.size());
    dy.reserve(x.size());
    slope.reserve(x.size());
    intercept.reserve(x.size());
    for(int i = 0; i < x.size(); ++i){
        if(i < x.size()-1)
        {
            dx.push_back(x[i+1] - x[i]);
            dy.push_back(y[i+1] - y[i]);
            slope.push_back(dy[i] / dx[i]);
            intercept.push_back(y[i] - x[i] * slope[i]);
        }
        else
        {
            dx.push_back(dx[i-1]);
            dy.push_back(dy[i-1]);
            slope.push_back(slope[i-1]);
            intercept.push_back(intercept[i-1]);
        }
    }

    for (int i = 0; i < x_new.size(); ++i)
    {
        int idx = findNearestNeighbourIndex(x_new[i], x);
        y_new.push_back(slope[idx] * x_new[i] + intercept[idx]);

    }
    return y_new;
}
