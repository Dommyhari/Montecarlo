/*
 * Monte_globals.h
 *
 *  Created on: Oct 2, 2014
 *      Author: ganeshfx
 */



#ifndef MONTE_GLOBALS_H_
#define MONTE_GLOBALS_H_

#include<random>
//#include "Monte_classes.h"


// some type definitions

/* 3d Vector real */
typedef struct {double x; double y; double z; } vector3d;

typedef vector3d vec3d;

/* 3d Vector integer */
typedef struct {int x; int y; int z; } ivector3d;

typedef ivector3d ivec3d;


// some flags
extern int mc_parallel_flag; // to be initialized during init_mc method call
extern int mc_ncpus;       // =4 hard coded for testing

// variables for cells
// M: ivec3d get_cpu_gcoord(int myrank)
extern ivec3d cpu_gcoord;    // cpu global position vector

// M:  make_mc_tbox()
extern vec3d mc_tbox_x,mc_tbox_y,mc_tbox_z;                                  // mc_transformation box                                  MC
// M:  calc_cell_dim(rcut,rsample)
extern vec3d mc_cell_dim;                                                    // cell dimensions for MC                                 MC

// M:make_mc_tbox() -- implicit definition
extern double mc_volume;                                                    // some internal attributes for mapping
extern vec3d mc_height;                                                      // used for global cell coordinates


//**********************************************************************************
//     NOTE: values to be initialized using new methods
//**********************************************************************************
// have to be assigned from MD call using Init method callable from C

extern ivec3d mc_cpu_dim;                                                    // get and assign cpu dim array                           MD-init
extern int mc_ncpu;                                                          // no of cpus get & assign from                           MD-init
extern vec3d mc_simbox_x,mc_simbox_y,mc_simbox_z;                            // get & assign box physical dimension from               MD-init

extern int mc_tot_types;                                                     // get total types from MD
extern double mc_temp;                                                      // get simulation temperature (T=0: static and T>0: dynamic)
//extern ivec3d  mc_restriction[mc_tot_types];                                 // get restriction vectors from MD all element types
extern ivec3d  mc_restriction[];                                 // get restriction vectors from MD all element types

// ---------------------------------------------------------------

// M: calc_mc_cpu_cell_array()
extern ivec3d mc_cpu_cell_dim;                                               // cpu cell array dim                                     MC

// M: calc_mc_global_cell_array()
extern ivec3d mc_global_cell_dim;                                            // global cell array dim computation                      MC

// M: calc_mc_cpu_box()
extern vec3d  mc_cpu_box_x,mc_cpu_box_y,mc_cpu_box_z;                        // cpu_box physical dimension method                      MC


// same for all cpu--  domain decomposition -- in preliminary implementation (later changes for load balancing)

//to be included suitable declaration to assign the following quantities

extern int stack_total;  // total no of stacks
extern int ncells_cpu;   // no of cells per cpu
extern int ncells_stack; // no of cells per stack in a cpu


// M: calc_block_dim() from rank 0

extern int rows_block;     // no of rows in global domain
extern int cols_block;     // no of cols in global domain
extern int stacks_block;   // no of stacks in global domain

// some global variables in c++ scope using containers-Vector

extern long mc_tatoms_cpu;             //  total particles per cpu
extern int  mc_pid;                    //  process id
extern int  mc_prank;                  //  process rank
extern int  mc_nbcells ;             // =6  NB cells count
extern double mc_rsweep ;              //  r_cut + r_sample


#endif /* MONTE_GLOBALS_H_ */
