/*
UDFs for modeling centrifugal reverse osmosis
*/

// Don't worry about this these imports, Fluent does all the compiling
#include "udf.h"
#include "sg.h"
#include "mem.h"
#include <stdlib.h>
#include <math.h>

DEFINE_ON_DEMAND(species_list)
{
    /*
    NOTE: NOT NEEDED FOR SIMULATIONS, SIMPLY A TOOL THAT IS HELPFUL TO USE
    Function to list the different species in a mixture template defined in Fluent and their corresponding index number
    */
    int i = -1;
    char *spe_name;
    Domain *d;
    Material *mix_mat;
    d=Get_Domain(1);
    mix_mat= mixture_material(d);
    mixture_species_loop_i(mix_mat,i)
        {
        spe_name=MIXTURE_SPECIE_NAME(mix_mat,i);
        Message("species index=%d,species name=%s\n",i,spe_name);
        }
}

DEFINE_ON_DEMAND(face_conc)
{
    /* 
    NOTE: NOT NEEDED FOR SIMULATIONS, SIMPLY A TOOL THAT IS HELPFUL TO USE
    Useful UDF that you can call whenever you want in Fluent to print out data you are interested in 
    */
    face_t f; cell_t c;
    Domain *d;
    d = Get_Domain(1);
    int zone_mem = 13; // Change to zone of interest
    Thread *t = Lookup_Thread(d, zone_mem);
    begin_c_loop(c,t)
    {
        Message("C_YI: %16.13f\t C_V: %16.13f\n", C_YI(c,t,0), C_V(c,t));
    }
    end_c_loop(c,t)
}

DEFINE_INIT(init_conc, d)
{
    /*
    Initialize the concentration of salt in the solution
    Inputs:
        - init_conc: UDF name
        - d: pointer to the domain which the initialization function is to be applied
    Return:
        - void
    */
    cell_t c;
    Thread *t;
    real C0, rho;
    // Loop over all cell threads
    thread_loop_c(t,d)
    {
        // Loop over all cells in this cell thread
        begin_c_loop_all(c,t)
        {
            if (THREAD_ID(t) != 16) // Change this to whatever cell zone id your porous zone is
            {
                // Non-zero concentration for the feed side
                C0 = 6; // kg/m^3
                // Density function of concentration comes from reference paper (need to derive new formula for ocean water)
                rho = 0.909*C0 + 997.1; // kg/m^3
                C_YI(c,t,0) = C0/rho; // mass fraction of salt species
            }
            else 
            {
                // Initialize the concentration in the membrane cells to 0
                C_YI(c,t,0) = 0;
            }
        }
        end_c_loop_all(c,t)
    }
}

DEFINE_PROFILE(outlet_flux, t, yi)
{
    /*
        Approximating zero diffusive flux at the feed outlet by assuming the concentration on outlet face is equal to the adjacent cell concentration
        Inputs: 
            - outlet_flux: UDF name
            - t: pointer to thread on which boundary condition is to be applied
            - yi: int index that identifies the variable that is to be defined
        Return:
            - void
    */

    real A[ND_ND], A_tmp[ND_ND];
    real x[ND_ND], x_tmp[ND_ND];
    cell_t cell_m;
    face_t f;
    Thread *t_cm, *t_f;
    int n, i;

    // Loop through faces in face thread
    begin_f_loop(f, t)
    {
        // Get area vector and face centroid
        F_AREA(A, f, t);
        F_CENTROID(x, f, t);

        // Find adjacent cell connectivity
        cell_m = F_C0(f, t);
        t_cm = THREAD_T0(t);

        // Assign the face value to that of the cell center
        F_PROFILE(f,t,yi) = C_YI(cell_m,t_cm,0);
    }
    end_f_loop(f, t)
}

DEFINE_PROPERTY(density, c, t)
{
    /* 
        Update the cell density (kg/m^3)
        Inputs:
            - density: UDF name
            - c: int cell id
            - t: pointer to thread where property is to be applied
        Return:
            - density as a function of concentration
    */
    return 997.1 + 0.909*C_R(c, t)*C_YI(c, t, 0);
}

DEFINE_DIFFUSIVITY(diffusivity, c, t, i)
{
    /* 
        Update the diffusivity of the solution (m^2/s)
        Inputs: 
            - diffusivity: UDF name
            - c: int cell id
            - t: pointer to thread where property is to be applied
            - i: index that identifies the species (not using) 
        Return:
            - diffusivity for species mixture as function of concentration
    */
    return 1.16e-9 - 3.9e-12*C_R(c, t)*C_YI(c, t, 0);
}

DEFINE_PROPERTY(viscosity, c, t)
{
    /* 
        Update the viscosity of the solution (Pa.s) 
        Inputs: 
            - viscosity: UDF name
            - c: int cell id
            - t: pointer to thread where property is to be applied
        Return:
            - viscosity as a function of concentration
    */
    return 8.9e-4 + 3.133e-6*C_R(c, t)*C_YI(c, t, 0);
}

DEFINE_PROFILE(porous_zone, t, i)
{
    /* 
        Make the permeability of the membrane a function of concentration 
        Inputs: 
            - porous_zone: UDF name
            - t: pointer to cell thread on which boundary condition will be applied
            - i: int index that indicates the variable that is to be defined
        Return:
            - void
    */
    cell_t c;
    real K = 1.1467e-11, dz = 0.04e-3; // Hard coded will have to update
    
    // Loop through cells in cell thread
    begin_c_loop(c, t)
    {
        C_PROFILE(c,t,i) = 1/(K*C_MU_EFF(c,t)*dz);
    }
    end_c_loop(c, t)
}

DEFINE_PROFILE(membrane_concentration, t, y_i)
{
    /*
        Calculate the mass fraction of salt on the membrane surface as a Dirichlet BC
        Inputs:
            - membrane_concentration: UDF name
            - t: pointer to thread on which boundary condition is to be applied
            - y_i: int index specifying the variable to be changed
        Return:
            - void
    */
    real dz = 0.04e-3, dp, p_osm, p_atm = 101325, K = 1.1467e-11; // dz is hard coded, you can easily change to get that number with code
    real d2, d2_tmp;
    cell_t c, cell_a, cell_b, a_id, b_id;
    real C_m, C_a, C_b, D, J, rho_avg, yi_avg;
    real x[ND_ND], x_tmp[ND_ND];
    // Get neighboring cell zones
    Domain *d = Get_Domain(1);
    int zone_a = 14, zone_b = 15;
    Thread *t_a = Lookup_Thread(d, zone_a);
    Thread *t_b = Lookup_Thread(d, zone_b);

    // Loop over cells in the membrane adjacent cell zone
    begin_c_loop(c,t)
    {
        // Find cell centroid
        C_CENTROID(x,c,t);

        // Loop over cells in next zone
        d2 = -1;
        begin_c_loop(cell_a,t_a)
        {
            // Find cell centroid
            C_CENTROID(x_tmp,cell_a,t_a);

            // Check if we have already defined our starting distance (measure of distance d2 = dx^2 + dy^2)
            // If not, initialize it
            if (d2 != -1) 
            {
                // Get tmp distance
                d2_tmp = pow(x_tmp[0] - x[0], 2.) + pow(x_tmp[1] - x[1], 2.);
                // Check if this distance is smaller
                if (d2_tmp < d2)
                {
                    // Update d2 and cell id
                    d2 = d2_tmp;
                    a_id = cell_a;
                }
            }
            else 
            {
                d2 = pow(x_tmp[0] - x[0], 2.) + pow(x_tmp[1] - x[1], 2.);
                a_id = cell_a;      
            }
        }
        end_c_loop(cell_a,t_a)

        // Found the closest cell in cell zone A
        // Access C_a value
        C_a = C_YI(a_id, t_a, 0)*C_R(a_id, t_a);

        d2 = -1;
        begin_c_loop(cell_b,t_b)
        {
            // Find cell centroid
            C_CENTROID(x_tmp,cell_b,t_b);

            // Check if we have already defined our starting distance
            // If not, initialize it
            if (d2 != -1) 
            {
                // Get tmp distance
                d2_tmp = pow(x_tmp[0] - x[0], 2.) + pow(x_tmp[1] - x[1], 2.);
                // Check if this distance is smaller
                if (d2_tmp < d2)
                {
                    // Update d2 and cell id
                    d2 = d2_tmp;
                    b_id = cell_b;
                }
            }
            else 
            {
                d2 = pow(x_tmp[0] - x[0], 2.) + pow(x_tmp[1] - x[1], 2.);
                b_id = cell_b;      
            }
        }
        end_c_loop(cell_b,t_b)

        // Found the closest cell in cell zone B
        // Access C_b value
        C_b = C_YI(b_id, t_b, 0)*C_R(b_id, t_b);

        // Calculate mass fraction of salt on membrane
        dp = C_P(c,t) - p_atm; // units: Pa
        p_osm = 0.523*C_R(c, t)*C_YI(c, t, 0)*1e5; // units: Pa
        // Calculate the volume flux / velocity through the membrane
        J = K*(dp - p_osm); // v = -J   

        // Technically you should be able to replace the last three lines with the one below but I haven't tested to see if you get the same result
        // This can be tested by simply printing out or writing the two different values to a file
        // J = -C_V(c,t); // v = -J only negative because of the orientation of the system

        // Should be way to access the diffusivity but haven't looked for it, instead just calculating
        D = 1.16e-9 - 3.9e-12*C_R(c, t)*C_YI(c, t, 0);

        C_m  = D*(4*C_a - C_b)/(3*D - 2*J*dz);
        
        C_PROFILE(c,t,y_i) = C_m/C_R(c,t);
    }
    end_c_loop(c,t)
}


DEFINE_PROFILE(permeate_pressure, t, p)
{
    /*
        UDF that assigns the pressure profile at the backside of the membrane (permeate side) such that the flux through the membrane is J = K(dp - p_osm)
        Inputs:
            - permeate_pressure: UDF name
            - t: pointer to thread on which boundary condition is to be applied
            - p: int index specifying the variable to be changed
        Return:
            - void
    */
    real p_osm;
    real x[ND_ND], x_tmp[ND_ND];
    face_t f;
    cell_t c, c_id;
    real d2, d2_tmp;
    Domain *d = Get_Domain(1);
    int zone = 13;
    Thread *t_m = Lookup_Thread(d, zone);

    begin_f_loop(f,t)
    {
        // Find face centroid
        F_CENTROID(x,f,t);

        // Loop over cells in next zone
        d2 = -1;
        begin_c_loop(c,t_m)
        {
            // Find cell centroid
            C_CENTROID(x_tmp,c,t_m);

            // Check if we have already defined our starting distance
            // If not, initialize it
            if (d2 != -1) 
            {
                // Get tmp distance
                d2_tmp = pow(x_tmp[0] - x[0], 2.) + pow(x_tmp[1] - x[1], 2.);
                // Check if this distance is smaller
                if (d2_tmp < d2)
                {
                    // Update d2 and cell id
                    d2 = d2_tmp;
                    c_id = c;
                }
            }
            else 
            {
                d2 = pow(x_tmp[0] - x[0], 2.) + pow(x_tmp[1] - x[1], 2.);
                c_id = c;      
            }
        }
        end_c_loop(c,t_m)

        p_osm = 0.523*C_YI(c_id, t_m, 0)*C_R(c_id, t_m)*1e5; // units: Pa

        F_PROFILE(f,t,p) = p_osm;
    }
    end_f_loop(f,t)
}

DEFINE_EXECUTE_AT_END(vol)
{
    /*
        UDF executed at the end of each time step to write the volume of fresh water produced
        Inputs:
            - vol: UDF name
        Return:
            - void
    */
    int mem_id = 24;
    Domain *d = Get_Domain(1);
    Thread *t = Lookup_Thread(d, mem_id);
    real A[ND_ND], J;
    face_t f;
    real dt = CURRENT_TIMESTEP, vol = 0; 

    // Loop through faces on membrane
    begin_f_loop(f,t)
    {
        J = fabs(F_V(f,t));
        F_AREA(A,f,t);
        // Calculate volume produced
        vol += J*NV_MAG(A)*dt;
    }
    end_f_loop(f,t)
    
    FILE *fp;
    fp = fopen("C:/Users/scant/fall23/CRO/flux.txt", "a"); // Change file to whatever you want
    // Remember to delete file when starting new simulation
    fprintf(fp, "%e, %d\n", vol, N_TIME); // N_TIME is the time step
    fclose(fp);
}