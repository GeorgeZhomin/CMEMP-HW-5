#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include "linfuncs.h"
#include "spline.h"

int main()
{
    const double a = 0, b = M_PI;

    std::ofstream file("output/output.txt");
    for (int n = 5; n < 50000; n*=2)
    {
        //make spline
        Grid *grid = make_chebyshev_grid(n, a, b);
        Spline spline(grid, map(grid->knots, grid->len, [](double x){return std::sin(x);}));

        //estimate error norm
        double norm = 0, diff, x;
        for (int i = 0; i < grid->len-1; i++)
        {
            x = (grid->knots[i] + grid->knots[i+1])/2;
            diff = std::abs(spline.f(x) - std::sin(x));
            if (diff > norm)
            {
                norm = diff;
            }
        }

        //write in file
        file << n << '\t' << grid->norm << '\t' << norm << '\n';
    }
    file.close();

    return 0;
}