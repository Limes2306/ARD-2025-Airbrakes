#include <stdio.h>
#include <math.h>
// #include "rk4_solver/rk4.h"

// Predict_Apogee()

double dvdt(float v, float t){

    double rho = 1.01065; // Avg air density across 1-3km altitude [kg/m^3]
    double CD = 0.52; // Avg drag coefficient of rocket
    double A = 0.0201062; // Cross-sectional area of rocket body [m^2]
    double m = 16.318; // Mass of rocket after burnout [kg]
    double g = 9.81; // Acceleration due to gravity [m/s^2]

    double res = (-1 * rho * CD * A)*v*v/(2*m) - 9.81;

    // printf("%.f\n", res);

    return res;

    // return (-0.000647*v*v - 9.81);
}

double dhdt(float v){
    return v; //I know this is dumb but its just here to make the process clearer
}

// takes 1 Runge-Kutta step
double rk4(double t0, double v0, double dt){

    double k1 = dvdt(v0, t0);
    
    double t1 = t0 + dt/2.0;
    double v1 = v0 + k1*dt/2.0;
    double k2 = dvdt(v1, t1);

    double t2 = t0 + dt/2.0;
    double v2 = v0 + k2*dt/2.0;
    double k3 = dvdt(v2, t2);
    
    double t3 = t0 + dt;
    double v3 = v0 + dt*k3;
    double k4 = dvdt(v3, t3);

    double v_step = v0 + dt*(k1 + 2*k2 + 2*k3 + k4)/6.0;

    return v_step;

}

// // takes 1 Runge-Kutta step
// double rk4_dhdt(double t0, double v0, double dt){

//     double k1 = dhdt(v0, t0);
    
//     double t1 = t0 + dt/2.0;
//     double h1 = v0 + k1*dt/2.0;
//     double k2 = dhdt(v1, t1);

//     double t2 = t0 + dt/2.0;
//     double h2 = v0 + k2*dt/2.0;
//     double k3 = dhdt(v2, t2);
    
//     double t3 = t0 + dt;
//     double h3 = v0 + dt*k3;
//     double k4 = dhdt(v3, t3);

//     double v_step = v0 + dt*(k1 + 2*k2 + 2*k3 + k4)/6.0;

//     return v_step;

// }

// double 

// double rk4_step()

int main(void){
    double t0 = 4.729;
    double v0 = 294.06;
    double dt = 0.1;
    // double h0 = 
    double t_max = 25;

    double res = dvdt(v0, t0);

    double v1;
    double t1;

    while(1){

        if(t0 >= t_max){
            break;
        } else if (v0 <= 0){
            break;
        }
        
        v1 = rk4(t0, v0, dt);
        t1 = t0 + dt;

        v0 = v1;
        t0 = t1;

        printf ( "  %g  %g\n", t0, v0);
    }

    // printf("%.2f\n", );

    // rk4()
}