#include "./../Header_Files/Sign.hpp"





int sign_val(double x)
{
    int result;
    if (x < 0.)
    {
        result = -1;
    }
    else if (x > 0.)
    {
        result = 1;
    }
    else
    {
        result = 0;
    }
    return result;
}
