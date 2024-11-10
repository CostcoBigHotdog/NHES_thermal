#include <iostream>
#include <fstream>
#include <cmath>

const double pi = 3.141592653589793;

//
// Function to solve tridiagonal matrix using Thomas Algorithm
//
void thomasAlgorithm(double A[], double B[], double C[], double D[], double X[], int nx) {
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


//
// Time steps
//
const int nt = 40000;                   // Number of time points
const double dt = 1e-3;                 // Time step size


//
// Kinetics parameters
//
const double D = 1.43;
const double v = 88.46e6;
const double sum_a = 0.004835;
const double sum_f = 0.001782;
const double rho_0 = 1.748567697835726;
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
const double L = 47 * 2.54;  // Core height
const int nx = 100;                     // Number of spatial points
const double dx = L / (nx - 1);     // Spatial step size

//
// Fuel pin & cladding parameters
//
// Geometry
const double r = 94.2 / 2 * 2.54;
const double A_core = pi * r * r;
const double r_f = 0.272;
const double r_cl = 0.328;
const double k_cl = 0.11;
const double A_f = pi * r_f * r_f;
const double A_cl = pi * r_cl * r_cl - pi * r_f * r_f;
const int n_pin = 11382;
const double W_na = 5.4e6;
const double v_na = 600;
const double A_na = W_na / v_na / n_pin;
// PDE
const double rhoCA_fuel = 16.2 * 0.34 * A_f * (840 / 425);
const double rhoCA_clad = 6.551 * 0.33 * A_cl * (840 / 425);
const double h_clad = 2 * pi * r_f * 0.4 * (840 / 425);


//
// Main function
//
int main() {
    //
    // Initialization
    //
    // Reactor core
    double p_0[nx] = { 0 };          // initial power
    double p_real[nx] = { 0 };       // Real power, initial power * n_pin * relative power
    double current_power = 0;        // Middle varible to calculate p_totalFractional
    double p_totalFractional = 0;    // Total relative power
    double a[nx], b[nx], c[nx], d[nx];
    double p[nx] = {
    0.0774436968891907, 0.124536481969205, 0.171511789696009, 0.218325307473560, 0.264932875324990,
    0.311290527549453, 0.357354534195706, 0.403081442313309, 0.448428116942519, 0.493351781804225,
    0.537810059651525, 0.581761012244901, 0.625163179913258, 0.667975620663528, 0.710157948801944,
    0.751670373030537, 0.792473733982939, 0.832529541164067, 0.871800009258855, 0.910248093775777,
    0.947837525991535, 0.984532847163953, 1.02029944198081, 1.05510357121304, 1.08891240354150,
    1.12169404652735, 1.15341757669666, 1.18405306871109, 1.21357162359694, 1.24194539600608,
    1.26914762048284, 1.29515263671243, 1.31993591372667, 1.34347407304446, 1.36574491072515,
    1.38672741831384, 1.40640180265902, 1.42474950458381, 1.44175321639309, 1.45739689820018,
    1.47166579305754, 1.48454644087719, 1.49602669112792, 1.50609571429702, 1.51474401210597,
    1.52196342647032, 1.52774714719537, 1.53208971840029, 1.53498704366481, 1.53643638989338,
    1.53643638989338, 1.53498704366481, 1.53208971840029, 1.52774714719537, 1.52196342647032,
    1.51474401210597, 1.50609571429702, 1.49602669112792, 1.48454644087719, 1.47166579305753,
    1.45739689820018, 1.44175321639309, 1.42474950458381, 1.40640180265902, 1.38672741831384,
    1.36574491072515, 1.34347407304446, 1.31993591372666, 1.29515263671243, 1.26914762048284,
    1.24194539600607, 1.21357162359694, 1.18405306871108, 1.15341757669666, 1.12169404652735,
    1.08891240354150, 1.05510357121304, 1.02029944198081, 0.984532847163951, 0.947837525991532,
    0.910248093775775, 0.871800009258853, 0.832529541164064, 0.792473733982937, 0.751670373030535,
    0.710157948801942, 0.667975620663526, 0.625163179913255, 0.581761012244899, 0.537810059651523,
    0.493351781804223, 0.448428116942518, 0.403081442313308, 0.357354534195705, 0.311290527549452,
    0.264932875324990, 0.218325307473559, 0.171511789696008, 0.124536481969205, 0.0774436968891905
    };


    // neutron diffusion: C_i
    double c1[nx];
    for (int i = 0; i < nx; ++i) {
        c1[i] = beta1 * fisnum * sum_f / lamb1 * p[i];
    }
    double c2[nx];
    for (int i = 0; i < nx; ++i) {
        c2[i] = beta2 * fisnum * sum_f / lamb2 * p[i];
    }
    double c3[nx];
    for (int i = 0; i < nx; ++i) {
        c3[i] = beta3 * fisnum * sum_f / lamb3 * p[i];
    }
    double c4[nx];
    for (int i = 0; i < nx; ++i) {
        c4[i] = beta4 * fisnum * sum_f / lamb4 * p[i];
    }
    double c5[nx];
    for (int i = 0; i < nx; ++i) {
        c5[i] = beta5 * fisnum * sum_f / lamb5 * p[i];
    }
    double c6[nx];
    for (int i = 0; i < nx; ++i) {
        c6[i] = beta6 * fisnum * sum_f / lamb6 * p[i];
    }

    // Fuel pin & cladding temperature
    double T_f[100] = {
    338.867144056664, 350.500355435206, 362.160537154101, 373.857213207554, 385.579351240069,
    397.315893562218, 409.055768896691, 420.787902822078, 432.501228219585, 444.184695712805,
    455.827284090803, 467.418010704590, 478.945941827258, 490.400202967936, 501.769989129882,
    513.044575003128, 524.213325081410, 535.265703695476, 546.191284951116, 556.979762564194,
    567.620959582742, 578.104837987049, 588.421508158697, 598.561238209591, 608.514463162229,
    618.271793972491, 627.824026386493, 637.162149623110, 646.277354873997, 655.161043613080,
    663.804835707702, 672.200577323735, 680.340348617227, 688.216471205329, 695.821515409440,
    703.148307263745, 710.189935282545, 716.939756979966, 723.391405135936, 729.538793802490,
    735.376124044743, 740.897889411124, 746.098881127707, 750.974193011732, 755.519226099684,
    759.729692985588, 763.601621865376, 767.131360283566, 770.315578578675, 773.151273024155,
    775.635768661849, 777.766721825331, 779.542122350727, 780.960295472925, 782.019903405428,
    782.719946602284, 783.059764701005, 783.039037145474, 782.657783488340, 781.916363372579,
    780.815476192226, 779.356160432632, 777.539792690840, 775.368086377023, 772.843090098183,
    769.967185725684, 766.743086148381, 763.173832713503, 759.262792357714, 755.013654431004,
    750.430427216482, 745.517434149292, 740.279309738234, 734.720995193962, 728.847733767860,
    722.665065805981, 716.178823522765, 709.395125499398, 702.320370912059, 694.961233495465,
    687.324655247434, 679.417839880369, 671.248246025885, 662.823580198958, 654.151789528225,
    645.241054259346, 636.099780038426, 626.736589982829, 617.160316546855, 607.379993189927,
    597.404845855186, 587.244284266501, 576.907893052116, 566.405422703311, 555.746780376605,
    544.942020548153, 534.001335529219, 522.935045851361, 511.753587851104, 500.431569467341
    }; // Copied from .mat file
    double h_clad = 2 * 3.14 * 0.2 * 0.34;
  
    //
    // Time iteration
    //
    for (int k = 0; k < nt; ++k) {
        // Core kinetics update
        for (int i = 0; i < nx; ++i) {
            double f_c1 = beta1 * fisnum * sum_f * p[i] - lamb1 * c1[i];
            c1[i] = c1[i] + dt * f_c1;

            double f_c2 = beta2 * fisnum * sum_f * p[i] - lamb2 * c2[i];
            c2[i] = c2[i] + dt * f_c2;

            double f_c3 = beta3 * fisnum * sum_f * p[i] - lamb3 * c3[i];
            c3[i] = c3[i] + dt * f_c3;

            double f_c4 = beta4 * fisnum * sum_f * p[i] - lamb4 * c4[i];
            c4[i] = c4[i] + dt * f_c4;

            double f_c5 = beta5 * fisnum * sum_f * p[i] - lamb5 * c5[i];
            c5[i] = c5[i] + dt * f_c5;

            double f_c6 = beta6 * fisnum * sum_f * p[i] - lamb6 * c6[i];
            c6[i] = c6[i] + dt * f_c6;
        }

        for (int i = 0; i < nx; ++i) {
            b[i] = (1 - ((rho - beta) / lifetime - (1 / ln) - sum_a * v) * dt + 2 * D * v * dt / (dx * dx));
            a[i] = -D * v * dt / (dx * dx);
            c[i] = -D * v * dt / (dx * dx);
            d[i] = p[i] + dt * v * (lamb1 * c1[i] + lamb2 * c2[i] + lamb3 * c3[i] + lamb4 * c4[i] + lamb5 * c5[i] + lamb6 * c6[i]);
        }
        thomasAlgorithm(a, b, c, d, p, nx);

        for (int i = 0; i < nx; ++i) {
            p_real[i] = p[i] * p_0[i] * n_pin;
        }
        for (int i = 0; i < nx; ++i) {
            current_power += p[i];
        }
        p_totalFractional = current_power / nx;

        // Fuel temperature update


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
