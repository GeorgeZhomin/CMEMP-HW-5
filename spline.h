#pragma once

class Grid
{
public:
    int len; //number of knots in grid
    double *knots;
    double norm;
    double *steps;

    Grid(int len, double *knots, double norm, double *steps);
    ~Grid();
};

//Чебышевская сетка с равномерной нормой
Grid* make_chebyshev_grid(int n, double a, double b);

class Spline
{
private:
    double* calc_coeffs();
    int get_segment_index(double x);
    double get_value_on_segment(double x, int i);

public:
    double *coeffs;
    double *u;
    Grid *grid;
    
    Spline(Grid *grid, double *u);
    ~Spline();

    double f(double x) {return get_value_on_segment(x, get_segment_index(x));}
};