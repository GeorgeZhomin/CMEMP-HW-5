#include "spline.h"
#include "linfuncs.h"
#define _USE_MATH_DEFINES
#include <cmath>

Grid::Grid(int len, double *knots, double norm, double *steps)
{
    this->len = len;
    this->knots = knots;
    this->norm = norm;
    this->steps = steps;
}

Grid::~Grid()
{
    delete[] knots;
    delete[] steps;
}

Grid* make_chebyshev_grid(int n, double a, double b)
{
    double *knots = new double[n];
    double mid = (a+b)/2, delta = (b-a)/2;
    for (int i = 0; i < n; i++)
    {
        knots[n-1-i] = mid + delta*std::cos((M_PI*(2*i+1))/(2*n));
    }
    double *steps = new double[n-1];
    double coeff = (b-a)*std::sin(M_PI/(2*n));
    for (int i = 1; i < n; i++)
    {
        steps[n-1-i] = coeff*std::sin((M_PI*i)/n);
    }
    double norm = steps[(n-1)/2];
    return new Grid(n, knots, norm, steps);
}

double* Spline::calc_coeffs()
{
    double *m = new double[grid->len-2];
    double *l = new double[grid->len-3];
    double *b = new double[grid->len-2];
    for (int i = 0; i < grid->len-3; i++)
    {
        l[i] = grid->steps[i+1]/6;
    }
    for (int i = 0; i < grid->len-2; i++)
    {
        b[i] = (u[i+2] - u[i+1])/grid->steps[i+1] - (u[i+1]-u[i])/grid->steps[i];
        m[i] = (grid->steps[i]+grid->steps[i+1])/3;
    }
    double *result = solve_tma(l, m, l, b, grid->len-2);
    delete[] m;
    delete[] l;
    delete[] b;
    return result;
}
    
int Spline::get_segment_index(double x)
{
    int l = 0, h = grid->len-1, m;
    while(l+1 < h)
    {
        m = (l+h)/2;
        if (grid->knots[m] < x)
        {
            l = m;
        }
        else
        {
            h = m;
        }
    }
    return l;
}

double Spline::get_value_on_segment(double x, int i)
{
    double Dx0 = (x-grid->knots[i]), Dx1 = (grid->knots[i+1]-x), h = grid->steps[i];
    return (Dx1*u[i] + Dx0*u[i+1] + coeffs[i]*(Dx1 - h)*(Dx1 + h)*Dx1/6 + coeffs[i+1]*(Dx0 - h)*(Dx0+h)*Dx0/6)/h;
}

Spline::Spline(Grid *grid, double *u)
{
    this->grid = grid;
    this->u = u;

    double *result = calc_coeffs();
    coeffs = new double[grid->len];
    coeffs[0] = 0;
    coeffs[grid->len-1] = 0;
    for (int i = 1; i < grid->len-1; i++)
    {
        coeffs[i] = result[i-1];
    }
    delete[] result;
}

Spline::~Spline()
{
    delete[] coeffs;
    delete[] u;
    delete grid;
}