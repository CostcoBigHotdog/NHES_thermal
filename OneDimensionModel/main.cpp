#include <iostream>
#include <fstream>
#include <cmath>

// Time & space steps
const int nx = 100;  // Number of spatial points
const int nt = 40000; // Number of time points
const double dt = 1e-3; // Time step size
const double dx = 400.0 / (nx - 1); // Spatial step size

// Kinetics parameters
const double D = 1.43;
const double v = 88.46e6;
const double sum_a = 0.004835;
const double sum_f = 0.001782;
const double rho_0 = 1.7141;
double rho = rho_0;
const double lifetime = 4.48274e-7;
const double fisnum = 2.92;
const double ln = 2.95630e-7;

const double lamb1 = 0.0127095;
const double lamb2 = 0.0301051;
const double lamb3 = 0.112269;
const double lamb4 = 0.327525;
const double lamb5 = 1.23099;
const double lamb6 = 8.1784;

const double beta1 = 6.75e-5;
const double beta2 = 6.66e-4;
const double beta3 = 5.11e-4;
const double beta4 = 1.43e-3;
const double beta5 = 6.43e-4;
const double beta6 = 1.7e-4;
const double beta = beta1 + beta2 + beta3 + beta4 + beta5 + beta6;

//
// Core geometry
//
const double L = 400;  // Core height
const int n_pin = 11382;

//
// Function to solve tridiagonal matrix using Thomas Algorithm
//
void thomasAlgorithm(double A[], double B[], double C[], double D[], double X[]) {
    double* c_star = new double[nx];
    double* d_star = new double[nx];

    c_star[0] = C[0] / B[0];
    d_star[0] = D[0] / B[0];

    for (int i = 1; i < nx; ++i) {
        double m = 1.0 / (B[i] - A[i] * c_star[i - 1]);
        c_star[i] = C[i] * m;
        d_star[i] = (D[i] - A[i] * d_star[i - 1]) * m;
    }

    X[nx - 1] = d_star[nx - 1];

    for (int i = nx - 2; i >= 0; --i) {
        X[i] = d_star[i] - c_star[i] * X[i + 1];
    }

    delete[] c_star;
    delete[] d_star;
}

int main() {
    double p_0[nx] = { 0 }; // initial power
    double p_real[nx] = { 0 };  // real power
    double p_fractional = 0;    // total relatice power
    double current_power = 0;

    double a[nx], b[nx], c[nx], d[nx];
    double h_clad = 2 * 3.14 * 0.2 * 0.34; // Example value
    
    //
    // Initialization
    //
    // Relative power
    double p[nx] = { 0 };  
    for (int i = 0; i < nx; ++i) {
        p[i] = 1;
    }
    // C 
    double c1[nx];
    for (int i = 0; i < nx; ++i) {
        c1[i] = beta1 * fisnum * sum_f / lamb1 * 1;
    }
    double c2[nx];
    for (int i = 0; i < nx; ++i) {
        c2[i] = beta2 * fisnum * sum_f / lamb2 * 1;
    }
    double c3[nx];
    for (int i = 0; i < nx; ++i) {
        c3[i] = beta3 * fisnum * sum_f / lamb3 * 1;
    }
    double c4[nx];
    for (int i = 0; i < nx; ++i) {
        c4[i] = beta4 * fisnum * sum_f / lamb4 * 1;
    }
    double c5[nx];
    for (int i = 0; i < nx; ++i) {
        c5[i] = beta5 * fisnum * sum_f / lamb5 * 1;
    }
    double c6[nx];
    for (int i = 0; i < nx; ++i) {
        c6[i] = beta6 * fisnum * sum_f / lamb6 * 1;
    }

    // Time iteration
    for (int k = 0; k < nt; ++k) {
        // Core kinetics update
        for (int i = 0; i < nx; ++i) {
            b[i] = (1 - ((rho - beta) / lifetime - (1 / ln) - sum_a * v) * dt + 2 * D * v * dt / (dx * dx));
            a[i] = -D * v * dt / (dx * dx);
            c[i] = -D * v * dt / (dx * dx);
            d[i] = p[i] + dt * v * (lamb1 * c1[i] + lamb2 * c2[i] + lamb3 * c3[i] + lamb4 * c4[i] + lamb5 * c5[i] + lamb6 * c6[i]);
        }
        thomasAlgorithm(a, b, c, d, p);

        for (int i = 0; i < nx; ++i) {
            p_real[i] = p[i] * p_0[i] * n_pin;
        }

        for (int i = 0; i < nx; ++i) {
            current_power += p[i];
        }
        p_fractional = current_power / nx;

        //// Fuel temperature update
        //for (int i = 0; i < nx; ++i) {
        //    T_f[i][k] = T_f[i][k - 1] + dt * (p[i][k - 1] * p_0[i] + h_clad * (T_cl[i][k - 1] - T_f[i][k - 1])) / (16.2 * 0.34 * A_f);
        //}

        //// Reactivity change
        //if (k <= 50000) {
        //    rho += a_ext * (current_power - nx);
        //}
        //else if (k > 50000 && k <= 60000) {
        //    // Update reactivity based on conditions
        //}
        //else if (k > 60000) {
        //    // Update reactivity based on conditions
        //}
    }

    return 0;
}
