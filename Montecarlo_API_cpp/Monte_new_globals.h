/*
 * Monte_new_globals.h
 *
 *  Created on: Mar 6, 2015
 *      Author: ganeshfx
 */

#ifndef MONTE_NEW_GLOBALS_H_
#define MONTE_NEW_GLOBALS_H_




#include "Monte_new_classes.h" // class declarations


//using namespace std;


// Standard Template Library definitions
vector<int> mc_atomtypes;        //  vector container for atom types
vector<long> mc_atomnumber;      //  vector container for atom number
vector<double> mc_atommass;      //  vector container for atom mass
vector<double> mc_positions;     //  vector container for atom positions
vector<double> mc_velocities;    //  vector container for atom velocities
vector<double> mc_epot;          //  vector container for atom potential energy

// ------------------------------------------------------------------------
// list of particle id exported for sphere construction for each neighbor
//-------------------------------------------------------------------------
vector<long> list_nb1,list_nb2,list_nb3,list_nb4,list_nb5,list_nb6,list_nb7;

vector<long> list_ids[7] = {list_nb1,list_nb2,list_nb3,list_nb4,list_nb5,list_nb6,list_nb7};

// ***************************************************************************
//                           Check if required
// some flags
int mc_parallel_flag;       // to be initialized during init_mc method call
int mc_get_velocity = 0;    // 0: default, 1: MD need velocity data back from MC
// ***************************************************************************

// IMPORTANT:  NEED DEVELOPMENT: param file reader to be developed
// should be passed as a parameter file
// Dev Note:  following keys have to be passed for writing param file

string str_fname         = "mc_sphere.chkpt";
string str_fname_relax   = "mc_sphere_relax";
const char* file_name   =  str_fname.c_str();
const char* file_relax  =  str_fname_relax.c_str();

string md_binary  = "imd_eam_fire_nbl_virtual_atoms_spline";
string md_param   = "local_md";

//**********************************************************************************
//                           Initialization from MD part
//**********************************************************************************

// MPI handles
MPI_Comm comm_name; // MPI Cartesian communicator
MPI_Status status;

// MONTECARLO SIMULATION PARAMETERS
ivec3d mc_cpu_dim;                              // get and assign cpu dim array                           MD-init
int mc_ncpus;                                   // no of cpus get & assign from                           MD-init
vec3d mc_simbox_x,mc_simbox_y,mc_simbox_z;      // get & assign box physical dimension from               MD-init
int mc_tot_types;                               // get total types from MD
int mc_real_types;                              // get real types from MD
// for compilation : mc_tot_types = 9
// Dnote: Check whether it conflict with sphere wall type??
ivec3d  mc_restriction[9];           // get restriction vectors from MD all element types

// to be initialized in Monte_init methods
long   mc_tatoms_cpu;                           //  total particles per cpu
int    mc_prank;                                //  process rank

// start from here
double mc_temp;                                 //  get simulation temperature (T=0: static and T>0: dynamic)
double mc_rsweep ;                              //  r_cut + r_sample
double mc_sphere_wall ;                         //  sphere boundary thickness
int    mc_sample_seed;                          //  seed for zone sampler


// Montecarlo variables
int    mc_nbcells = 27 ;                        // 27  NB cells count
long   mc_real_cpu=0 ;                            // total real particles per CPU
long   mc_phold_cpu=0;                            // total placeholders in CPU
long   mc_carb_cpu=0 ;                            // total carbons in CPU

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

// Energy computation for acceptance condition

double energy_ref;
double energy_new;


//// possible neighbor index as per window location
//
//int window_x[8] = {-1,-1,+1,+1,-1,-1,+1,+1};
//int window_y[8] = {-1,-1,-1,-1,+1,+1,+1,+1};
//int window_z[8] = {-1,+1,-1,+1,-1,+1,-1,+1};
//
//// windows position
//
//ivec3d win_pos_0 = {0,0,0};
//ivec3d win_pos_1 = {0,0,1};
//ivec3d win_pos_2 = {1,0,0};
//ivec3d win_pos_3 = {1,0,1};
//ivec3d win_pos_4 = {0,1,0};
//ivec3d win_pos_5 = {0,1,1};
//ivec3d win_pos_6 = {1,1,0};
//ivec3d win_pos_7 = {1,1,1};
//
//ivec3d window_position[8]={win_pos_0,win_pos_1,win_pos_2,win_pos_3,win_pos_4,win_pos_5,win_pos_6,win_pos_7};


// some MPI information (NEED ALLOCATION)
long buffer_size; // total no of particles + threshold value


int ncells_cpu;   // no of cells per cpu
int ncells_stack; // no of cells per stack in a cpu


// M: calc_block_dim() from rank 0

// if possible replace methods (with MPI Standard functions)
int rows_block;     // no of rows in global computation domain
int cols_block;     // no of cols in global computation domain
int stacks_block;   // no of stacks in global computation domain

celltype old_sphere_config;
celltype new_sphere_config;



#endif /* MONTE_NEW_GLOBALS_H_ */
