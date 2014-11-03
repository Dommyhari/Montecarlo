/*
 * Monte_globals.h
 *
 *  Created on: Oct 2, 2014
 *      Author: ganeshfx
 */



#ifndef MONTE_GLOBALS_H_
#define MONTE_GLOBALS_H_

#include<random>
#include<vector>
#include<mpi.h>

using namespace std;


// STL definitions
vector<int> mc_atomtypes;        //  vector container for atom types
vector<long> mc_atomnumber;      //  vector container for atom number
vector<double> mc_atommass;      //  vector container for atom mass
vector<double> mc_positions;     //  vector container for atom positions
vector<double> mc_velocities;    //  vector container for atom velocities
vector<double> mc_epot;          //  vector container for atom potential energy


// some user type definitions

/* 3d Vector real */
typedef struct {double x; double y; double z; } vector3d;

typedef vector3d vec3d;

/* 3d Vector integer */
typedef struct {int x; int y; int z; } ivector3d;

typedef ivector3d ivec3d;


// some flags
int mc_parallel_flag;       // to be initialized during init_mc method call

int mc_get_velocity = 0;    // 0: default, 1: MD need velocity data

char* file_name = "mc_sphere.chkpt";

string md_binary = "imd_eam_glok_homdef_nbl_virtual_atoms";

string md_param = "local_md.param";

int accept_flag= 0;

//**********************************************************************************
//                           Initialization from MD part
//**********************************************************************************

ivec3d mc_cpu_dim;                              // get and assign cpu dim array                           MD-init
int mc_ncpus;                                   // no of cpus get & assign from                           MD-init
vec3d mc_simbox_x,mc_simbox_y,mc_simbox_z;      // get & assign box physical dimension from               MD-init
int mc_tot_types;                               // get total types from MD
int mc_real_types;                              // get real types from MD
ivec3d  mc_restriction[mc_tot_types];           // get restriction vectors from MD all element types

// to be initialized in Monte_init methods
long   mc_tatoms_cpu;                           //  total particles per cpu
int    mc_prank;                                //  process id

double mc_temp;                                 //  get simulation temperature (T=0: static and T>0: dynamic)
double mc_rsweep ;                              //  r_cut + r_sample
double mc_sphere_wall ;                         //  sphere boundary thickness
int    mc_sample_seed;                          //  seed for zone sampler


// Montecarlo variables
int    mc_nbcells = 26 ;                        // 26  NB cells count
long   mc_real_cpu;                             // total real particles in this system
long   mc_phold_cpu;                            // total placeholders in cpu
long   mc_carb_cpu;                             // total carbons in cpu

// ---------------------------------------------------------------

// M: calc_mc_cpu_cell_array()
ivec3d mc_cpu_cell_dim;                                  // cpu cell array dim                                     MC

// M: calc_mc_global_cell_array()
ivec3d mc_global_cell_dim;                               // global cell array dim computation                      MC

//: calc_mc_cpu_box()
vec3d  mc_cpu_box_x,mc_cpu_box_y,mc_cpu_box_z;           // cpu_box physical dimension method                      MC

// same for all cpu--  domain decomposition -- in preliminary implementation (later changes for load balancing)

// variables for cells

// M: ivec3d get_cpu_gcoord(int myrank)
ivec3d cpu_gcoord;    // cpu global position vector

//:  make_mc_tbox()
vec3d mc_tbox_x,mc_tbox_y,mc_tbox_z;                                  // mc_transformation box                                  MC

//:  calc_cell_dim(rsample)
vec3d mc_cell_dim;                                                    // cell dimensions for MC                                 MC

//:make_mc_tbox() -- implicit definition
double mc_volume;                                                    // some internal attributes for mapping
vec3d mc_height;                                                      // used for global cell coordinates



//to be included suitable declaration to assign the following quantities

int stack_total;  // total no of stacks
int ncells_cpu;   // no of cells per cpu
int ncells_stack; // no of cells per stack in a cpu


// M: calc_block_dim() from rank 0

int rows_block;     // no of rows in global computation domain
int cols_block;     // no of cols in global computation domain
int stacks_block;   // no of stacks in global computation domain



#endif /* MONTE_GLOBALS_H_ */
