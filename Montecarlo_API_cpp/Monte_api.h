/*
 * Monte_api.h
 *
 *  Created on: Sep 22, 2014
 *      Author: ganeshfx
 */

#ifndef MONTE_API_H_
#define MONTE_API_H_


/*  function prototypes */

#ifdef __cplusplus
extern "C" {
#endif

    // prepare Montecarlo simulation setup from MD
    void init_montecarlo(int* md_cpu_dim, int md_tot_types,int md_real_types, int* md_restriction,double* md_simbox,
    		double md_temperature,double rsweep,double wall_thick,int sample_seed,MPI_Comm comm_name);

    // import data required for Montecarlo
    void pack_config_to_montecarlo(long md_mc_tatoms_cpu,long *md_mc_atomnumber,int *md_mc_atomtypes,
		double  *md_mc_atommass,double  *md_mc_positions , double *md_mc_epot);

    // perform montecarlo sampling and export updated data
    void do_montecarlo(int* md_pid,long *md_tatoms_cpu,long **md_atomnumber,int **md_atomtypes,
    		double **md_atommass,double **md_positions);

    // export current velocities from MOntecarlo
    void get_velocities(double **md_velocities);

    // export Potential energy to MD
    void get_pot_energy(double **md_epot);

    // clear all data structures created in MC scope
    void clean_montecarlo();

    // clear all MC simulation setup
    void shut_down_montecarlo();

#ifdef __cplusplus
}
#endif


#endif /* MONTE_API_H_ */
