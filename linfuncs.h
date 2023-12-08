#pragma once

//нижняя диагональ, главная диагональ, верхняя диагональ, вектор значений, число переменных
double* solve_tma(double *l, double *m, double *u, double *b, int n_vars);
double* map(double* x, int len, double (*f)(double));