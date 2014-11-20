/*
 * Monte_globals.h
 *
 *  Created on: Oct 2, 2014
 *      Author: ganeshfx
 */



#ifndef MONTE_GLOBALS_H_
#define MONTE_GLOBALS_H_


#include<iomanip>
#include<fstream>
#include<string>
#include<iostream>
#include<ctime>
#include<cmath>
#include<cstdlib>
#include<random>
#include<vector>
#include "Monte_classes.h"
#include "mpi.h"

using namespace std;


// STL definitions
vector<int> mc_atomtypes;        //  vector container for atom types
vector<long> mc_atomnumber;      //  vector container for atom number
vector<double> mc_atommass;      //  vector container for atom mass
vector<double> mc_positions;     //  vector container for atom positions
vector<double> mc_velocities;    //  vector container for atom velocities
vector<double> mc_epot;          //  vector container for atom potential energy


// some flags
int mc_parallel_flag;       // to be initialized during init_mc method call
int mc_get_velocity = 0;    // 0: default, 1: MD need velocity data


// should be passed as a parameter file // param file reader to be developed
char* file_name = "mc_sphere.chkpt";
string md_binary = "imd_eam_glok_homdef_nbl_virtual_atoms";
string md_param = "local_md.param";

int accept_flag= 0;

//**********************************************************************************
//                           Initialization from MD part
//**********************************************************************************

// MPI handles
MPI_Comm comm_name; // MPI Cartesian communicator
MPI_Status status;


ivec3d mc_cpu_dim;                              // get and assign cpu dim array                           MD-init
int mc_ncpus;                                   // no of cpus get & assign from                           MD-init
vec3d mc_simbox_x,mc_simbox_y,mc_simbox_z;      // get & assign box physical dimension from               MD-init
int mc_tot_types;                               // get total types from MD
int mc_real_types;                              // get real types from MD
ivec3d  mc_restriction[mc_tot_types];           // get restriction vectors from MD all element types

// to be initialized in Monte_init methods
long   mc_tatoms_cpu;                           //  total particles per cpu
int    mc_prank;                                //  process rank

double mc_temp;                                 //  get simulation temperature (T=0: static and T>0: dynamic)
double mc_rsweep ;                              //  r_cut + r_sample
double mc_sphere_wall ;                         //  sphere boundary thickness
int    mc_sample_seed;                          //  seed for zone sampler


// Montecarlo variables
int    mc_nbcells = 27 ;                        // 27  NB cells count
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


// windows sampler

// all possible window neighbors

int windows_count = 8;

// possible neighbor index as per window location

int window_x[8] = {-1,-1,+1,+1,-1,-1,+1,+1};
int window_y[8] = {-1,-1,-1,-1,+1,+1,+1,+1};
int window_z[8] = {-1,+1,-1,+1,-1,+1,-1,+1};


// sample cell ranges (HARDCODED FOR 8 CENTRAL CELLS)
ivec6d zone_limit_0 = {0,1,0,1,0,1}; // window 0
ivec6d zone_limit_1 = {0,1,0,1,1,2}; // window 1
ivec6d zone_limit_2 = {1,2,0,1,0,1}; // window 2
ivec6d zone_limit_3 = {1,2,0,1,1,2}; // window 3
ivec6d zone_limit_4 = {0,1,1,2,0,1}; // window 4
ivec6d zone_limit_5 = {0,1,1,2,1,2}; // window 5
ivec6d zone_limit_6 = {1,2,1,2,0,1}; // window 6
ivec6d zone_limit_7 = {1,2,1,2,1,2}; // window 7

ivec6d zone_limit[8]={zone_limit_0,zone_limit_1,zone_limit_2,zone_limit_3,zone_limit_4,zone_limit_5,zone_limit_6,zone_limit_7};

// windows position

ivec3d win_pos_0 = {0,0,0};
ivec3d win_pos_1 = {0,0,1};
ivec3d win_pos_2 = {1,0,0};
ivec3d win_pos_3 = {1,0,1};
ivec3d win_pos_4 = {0,1,0};
ivec3d win_pos_5 = {0,1,1};
ivec3d win_pos_6 = {1,1,0};
ivec3d win_pos_7 = {1,1,1};

ivec3d window_position[8]={win_pos_0,win_pos_1,win_pos_2,win_pos_3,win_pos_4,win_pos_5,win_pos_6,win_pos_7};

// some MPI information (NEED ALLOCATION)
long buffer_size; // allocated as double x 10 x buffer_size

//to be included suitable declaration to assign the following quantities

int stack_total;  // total no of stacks
int ncells_cpu;   // no of cells per cpu
int ncells_stack; // no of cells per stack in a cpu


// M: calc_block_dim() from rank 0

int rows_block;     // no of rows in global computation domain
int cols_block;     // no of cols in global computation domain
int stacks_block;   // no of stacks in global computation domain



#endif /* MONTE_GLOBALS_H_ */
