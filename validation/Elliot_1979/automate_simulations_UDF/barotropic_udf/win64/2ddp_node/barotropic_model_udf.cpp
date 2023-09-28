// Global coefficients for the barotropic_model
double pressure_scaling = 2.0e6;
double rhomass_coefficients[] = {4.756162e-01, 1.562659e+03, -2.284117e+03, 2.861341e+03, -2.736432e+03, 1.820571e+03, -7.756865e+02, 1.882134e+02, -1.966771e+01};
double viscosity_coefficients[] = {1.848772e-05, 1.483233e-03, -2.244149e-03, 2.906750e-03, -2.861206e-03, 1.946811e-03, -8.436345e-04, 2.073393e-04, -2.188098e-05};
double speed_sound_coefficients[] = {4.207661e+01, 6.656170e+01, -1.770000e+00, 5.997488e+00, -8.561204e+00, 7.301848e+00, -3.667316e+00, 9.975987e-01, -1.132715e-01};


/**
 * Evaluate a polynomial using Horner's rule.
 *
 * @param x The variable value.
 * @param coefficients The polynomial coefficients.
 * @param degree The polynomial degree.
 * @return The evaluated polynomial value for x.
 */
double horner_eval(double x, const double coefficients[], int degree) {
    double result = coefficients[degree];
    for (int i = degree-1; i >= 0; --i) {
        result = result * x + coefficients[i];
    }
    return result;
}


/**
 * Helper function to evaluate a property based on a given pressure and coefficients.
 * This function scales the pressure, computes the degree of the polynomial from the coefficients,
 * and then evaluates the polynomial using Horner's method.
 *
 * @param pressure The given pressure value.
 * @param coefficients The polynomial coefficients.
 * @return The evaluated property value for the given pressure.
 */
double evaluate_property(double pressure, const double coefficients[], int degree) {
    double x = pressure / pressure_scaling; // Scale the pressure
    return horner_eval(x, coefficients, degree); // Evaluate the polynomial for the scaled pressure
}


/**
 * Compute the density based on a given pressure.
 *
 * @param pressure The given pressure value.
 * @return The density value for the given pressure.
 */
double compute_density(double pressure) {
    int degree = sizeof(rhomass_coefficients) / sizeof(rhomass_coefficients[0]) - 1;
    return evaluate_property(pressure, rhomass_coefficients, degree);
}


/**
 * Compute the viscosity based on a given pressure.
 *
 * @param pressure The given pressure value.
 * @return The viscosity value for the given pressure.
 */
double compute_viscosity(double pressure) {
    int degree = sizeof(viscosity_coefficients) / sizeof(viscosity_coefficients[0]) - 1;
    return evaluate_property(pressure, viscosity_coefficients, degree);
}


/**
 * Compute the speed of sound based on a given pressure.
 *
 * @param pressure The given pressure value.
 * @return The speed of sound value for the given pressure.
 */
double compute_speed_sound(double pressure) {
    int degree = sizeof(speed_sound_coefficients) / sizeof(speed_sound_coefficients[0]) - 1;
    return evaluate_property(pressure, speed_sound_coefficients, degree);
}


/* 
 * ----------------------------
 * ANSYS Fluent UDF Definitions
 * ----------------------------
 * 
 * This section defines the User Defined Functions (UDFs) for ANSYS Fluent.
 * These functions allow Fluent to access custom property calculations based on 
 * barotropic relations. Each DEFINE_PROPERTY function retrieves the pressure
 * from the Fluent solver, computes the respective property (density, viscosity, 
 * or speed of sound) using the provided polynomial relations, and then returns
 * the computed value to the Fluent solver.
 * 
 */

#include "udf.h"

DEFINE_PROPERTY(density, cell, thread) {
    real pressure = C_P(cell, thread);
    return (real) compute_density(pressure);
}

DEFINE_PROPERTY(viscosity, cell, thread) {
    real pressure = C_P(cell, thread);
    return (real) compute_viscosity(pressure);
}

DEFINE_PROPERTY(speed_sound, cell, thread) {
    real pressure = C_P(cell, thread);
    return (real) compute_speed_sound(pressure);
}




/*
 * ==========================
 * Exporting User-Defined Scalars (UDS) for Postprocessing
 * --------------------------
 * This section includes definitions and functions responsible for exporting 
 * user-defined scalars to Fluent. These scalars, originating from the barotropic 
 * model calculations, can be used for visualization within the Fluent.
 * ==========================
 */

/* 
 * Enumerated constants representing user-defined scalars for the barotropic model.
 * Each scalar corresponds to a specific property of the model.
 * `N_REQUIRED_UDS` is used to count the total number of user-defined scalars 
 * defined before it in the enumeration.
 */
enum
{
  PRESSURE_UDS,
  DENSITY_UDS,
  VISCOSITY_UDS,
  SPEED_SOUND_UDS,
  N_REQUIRED_UDS
};


DEFINE_EXECUTE_AT_END(eval_barotropic)
{
    Domain *domain;
    Thread *t;
    cell_t c;
    face_t f;

    // Obtain the domain pointer
    domain = Get_Domain(1); 
    
    // Check the number of UDS
    if (n_uds < N_REQUIRED_UDS)
    {
        Message("Error: Number of allocated user-defined scalars (%d) is less than required (%d).\n", n_uds, N_REQUIRED_UDS);
        Message("Please allocate the required number of user-defined scalars under 'User-Defined -> Scalars' in the GUI.\n");
        // Internal_Error("User-defined scalars allocation mismatch.");
    }


    /* Loop through all cells in the domain. */
    thread_loop_c(t,domain)
    {
        /* Check if the UDS storage for pressure is available for the current thread. */
        if (NULL != THREAD_STORAGE(t,SV_UDS_I(PRESSURE_UDS)))
        {
            /* Begin iterating over cells in the current thread. */
            begin_c_loop (c,t)
            {
                /* Obtain the pressure value of the current cell. */
                real P_cell = C_P(c,t);
                
                /* Populate the User Defined Scalars (UDS) for the cell based on the pressure. */
                C_UDSI(c,t,PRESSURE_UDS) = P_cell;
                C_UDSI(c,t,DENSITY_UDS) = compute_density(P_cell);
                C_UDSI(c,t,VISCOSITY_UDS) = compute_viscosity(P_cell);
                C_UDSI(c,t,SPEED_SOUND_UDS) = compute_speed_sound(P_cell);
            }
            /* End iteration over cells. */
            end_c_loop (c,t)
        }
    }

    /* Loop through all faces in the domain. */
    thread_loop_f(t, domain) {
        /* Check if the UDS storage for pressure is available for the current face thread. */
        if (NULL != THREAD_STORAGE(t, SV_UDS_I(PRESSURE_UDS))) {
            /* Begin iterating over faces in the current thread. */
            begin_f_loop(f, t) {
                /* Initialize face pressure to 0. */
                real P_face = 0.;
                
                /* If the face thread has pressure storage, obtain the pressure value. 
                   Otherwise, obtain it from the adjacent cell. */
                if (NULL != THREAD_STORAGE(t, SV_P)) {
                    P_face = F_P(f, t);
                } else if (NULL != THREAD_STORAGE(t->t0, SV_P)) {
                    P_face = C_P(F_C0(f, t), t->t0);
                }

                /* Populate the User Defined Scalars (UDS) for the face based on the pressure. */
                F_UDSI(f, t, PRESSURE_UDS) = P_face;
                F_UDSI(f, t, DENSITY_UDS) = compute_density(P_face);
                F_UDSI(f, t, VISCOSITY_UDS) = compute_viscosity(P_face);
                F_UDSI(f, t, SPEED_SOUND_UDS) = compute_speed_sound(P_face);
            }
            /* End iteration over faces. */
            end_f_loop(f, t)
        }
    }


    // Print a message in the Fluent console indicating the postprocessing execution
    // Message("Barotropic model postprocessing executed successfully. (Process %d)\n", myid);


}





// #include <iostream>
// #include <iomanip>
// int main() {
//     // Test pressures
//     double pressures[] = {0.5e5, 1e5, 5e5, 10e5, 20e5, 40e5};
//     std::cout << "-----------------------------------------------------------------------------------" << std::endl;
//     std::cout << "Testing the barotropic model properties at different pressures:" << std::endl;
//     std::cout << "-----------------------------------------------------------------------------------" << std::endl;
//     std::cout << std::setw(15) << "Pressure (Pa)"
//               << std::setw(18) << "Density (kg/m^3)"
//               << std::setw(20) << "Viscosity (Pa.s)"
//               << std::setw(22) << "Speed of Sound (m/s)" << std::endl;
//     std::cout << "-----------------------------------------------------------------------------------" << std::endl;

//     // std::cout << std::fixed;
//     std::cout << std::scientific; // Set scientific notat

//     for (double pressure : pressures) {
//         double rho = compute_density(pressure);
//         double mu = compute_viscosity(pressure);
//         double speed = compute_speed_sound(pressure);

//         std::cout << std::setw(15) << std::setprecision(4) << pressure
//                   << std::setw(18) << std::setprecision(4) << rho
//                   << std::setw(20) << std::setprecision(4) << mu
//                   << std::setw(22) << std::setprecision(4) << speed << std::endl;
//     }
//     std::cout << "-----------------------------------------------------------------------------------" << std::endl;

//     return 0;
// }

