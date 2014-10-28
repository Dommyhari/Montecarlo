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


    void pack_config_to_montecarlo(int md_mc_pid,long md_mc_tatoms_cpu,long *md_mc_atomnumber,int *md_mc_atomtypes,
		double  *md_mc_atommass,double  *md_mc_positions);

    void do_montecarlo(int* md_pid,long *md_tatoms_cpu,long **md_atomnumber,int **md_atomtypes,double **md_atommass,double **md_positions);


    void clean_montecarlo();

#ifdef __cplusplus
}
#endif

// utility methods prototypes

// to be verified
//void make_particles();



#endif /* MONTE_API_H_ */
