//This code computes the stability time for a given set of parameters and initial conditions
// The stability time is the number of iterations before the distance between the points exceeds a threshold
//An initial condition is considered stable if its distance from the origin is less than a certain control radius (r_c) when N=N_iterations
//Since the computation is time-consuming, dynamical indicators are computed and it is possible to compare all the results


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Function prototypes
void henon(double *x, double *p_x, double *y, double *p_y, double nu_x, double nu_y, int i, double eps);
double stabilityTime(double x0, double y0, double nu_x, double nu_y, double eps, double r_c, int N_iterations);

int main(int argc, char *argv[])
{
    // Check if the correct number of arguments is provided
    if (argc != 7)
    {
        printf("Usage: %s nu_x nu_y eps N_steps lattice_length N_iterations\n", argv[0]);
        return (1);
    }

    // Parse command line arguments
    double nu_x = atof(argv[1]);
    double nu_y = atof(argv[2]);
    double eps = atof(argv[3]);
    int N_steps = atoi(argv[4]);
    double lattice_length = atof(argv[5]);
    int N_iterations = atoi(argv[6]);
    double x0, px_0, y0, py_0;
    double r_c = 10000.;  // Threshold distance for stability check
    int i, j;

    // Create a filename based on the number of iterations
    char filename[50];
    sprintf(filename, "stabtime%s.txt", argv[6]);
    FILE *stabtime = fopen(filename, "w");
    if (stabtime == NULL) {
        perror("Error opening file");
        return 1;
    }

    // Generate a grid of points and calculate stability time for each point
    for (int i = 0; i < N_steps; i++) // Loop over x
    {
        for (int j = 0; j < N_steps; j++) // Loop over y
        {
            double x0 = (double)i * (lattice_length / (double)N_steps);
            double y0 = (double)j * (lattice_length / (double)N_steps);
            double px_0 = 0.;
            double py_0 = 0.;
            
            // Calculate and write stability time for the current point
            fprintf(stabtime, "%lf\t%lf\t%lf\n", x0, y0, stabilityTime(x0, y0, nu_x, nu_y, eps, r_c, N_iterations));
            printf("Processing lattice positions: x0=%d\t y0=%d\n", i, j);
        }
    }
	printf("Stability time calculation completed\n");
    fclose(stabtime);
}

// Henon map function
void henon(double *x, double *p_x, double *y, double *p_y, double nu_x, double nu_y, int i, double eps)
{
    double omega_x;
    double omega_y;
    double omega1 = 2. * M_PI / 868.12;  // Fundamental frequency
    double tk[7] = {omega1, 2. * omega1, 3. * omega1, 6. * omega1, 7. * omega1, 10. * omega1, 12. * omega1}; // Harmonic frequencies
    double epsilon[7] = {1.E-4, 0.218E-4, 0.708E-4, 0.254E-4, 0.100E-4, 0.078E-4, 0.218E-4};  // Epsilon parameters
    double t0_x = 2. * M_PI * nu_x;
    double t0_y = 2. * M_PI * nu_y;

    omega_x = t0_x * (1. + eps * (epsilon[0] * cos(tk[0] * i) + epsilon[1] * cos(tk[1] * i) + epsilon[2] * cos(tk[2] * i) + epsilon[3] * cos(tk[3] * i) + epsilon[4] * cos(tk[4] * i) + epsilon[5] * cos(tk[5] * i) + epsilon[6] * cos(tk[6] * i)));
    omega_y = t0_y * (1. + eps * (epsilon[0] * cos(tk[0] * i) + epsilon[1] * cos(tk[1] * i) + epsilon[2] * cos(tk[2] * i) + epsilon[3] * cos(tk[3] * i) + epsilon[4] * cos(tk[4] * i) + epsilon[5] * cos(tk[5] * i) + epsilon[6] * cos(tk[6] * i)));

    double x_old;
    double y_old;
    double px_old;
    double py_old;

    // Store the current values
    x_old = *x;
    px_old = *p_x;
    y_old = *y;
    py_old = *p_y;

    // Update values according to the Henon map equations
    *x = x_old * cos(omega_x) + (px_old + (x_old * x_old - y_old * y_old)) * sin(omega_x);
    *p_x = -sin(omega_x) * x_old + cos(omega_x) * (px_old + (x_old * x_old - y_old * y_old));
    *y = y_old * cos(omega_y) + (py_old - (2. * x_old * y_old)) * sin(omega_y);
    *p_y = -sin(omega_y) * y_old + cos(omega_y) * (py_old - (2. * x_old * y_old));
}

// Function to calculate the stability time
double stabilityTime(double x0, double y0, double nu_x, double nu_y, double eps, double r_c, int N_iterations)
{
    double x, p_x, y, p_y;
    double dist;
    int k;

    x = x0;
    y = y0;
    p_x = 0.;
    p_y = 0.;

    // Iterate over the number of iterations
    for (k = 1; k <= N_iterations; k++)
    {
        henon(&x, &p_x, &y, &p_y, nu_x, nu_y, k, eps);
        dist = (x * x + y * y + p_x * p_x + p_y * p_y);

        // Check if the distance exceeds the threshold
        if (dist >= r_c)
        {
            return log10((double)k);  // Return the log of the iteration count
        }
    }

    return log10(N_iterations);  // If not exceeded, return the log of the max iterations
}
