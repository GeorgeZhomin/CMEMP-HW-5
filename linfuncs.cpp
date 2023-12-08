double* solve_tma(double *l, double *m, double *u, double *b, int n_vars)
{
    double *alpha = new double[n_vars-1];
    double *beta = new double[n_vars];

    alpha[0] = -u[0]/m[0];
    beta[0] = b[0]/m[0];
    double denom;
    for (int i = 1; i < n_vars-1; i++)
    {
        denom = m[i]+alpha[i-1]*l[i-1];
        alpha[i] = -u[i]/denom;
        beta[i] = (b[i]-l[i-1]*beta[i-1])/denom;
    }
    beta[n_vars-1] = (b[n_vars-1]-l[n_vars-2]*beta[n_vars-2])/(m[n_vars-1]+alpha[n_vars-2]*l[n_vars-2]);

    double* result = new double[n_vars];
    result[n_vars-1] = beta[n_vars-1];
    for (int i = n_vars-2; i >= 0; i--)
    {
        result[i] = alpha[i]*result[i+1] + beta[i];
    }
    delete[] alpha;
    delete[] beta;
    return result;
}

double* map(double* x, int len, double (*f)(double))
{
    double* result = new double[len];
    for (int i = 0; i < len; i++)
    {
        result[i] = f(x[i]);
    }
    return result;
}