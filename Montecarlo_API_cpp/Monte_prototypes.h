/*
 * Monte_prototypes.h
 *
 *  Created on: Oct 23, 2014
 *      Author: ganeshfx
 */

#ifndef MONTE_PROTOTYPES_H_
#define MONTE_PROTOTYPES_H_

//#include "Monte_globals.h"  // global variable definitions

#include "Monte_classes.h"
#include<string>

// utility methods definition
int get_cpu_rank(int inde0,int inde1,int inde2);

void setup_config(void);

void update_particle(cellblock bobj,int  cell_id, particle ref_atom );

// not yet defined
//void remove_particles(vector<long> list); // to be defined

void do_local_mdrun(string bin_name,string param_name);

particle sample_zone(cellblock bobj,int win_id);

// argument definition: (type of acceptance check, ref_sphere config before trial, new_sphere config after trial)
int acceptance_check(int type,celltype old_sphere,celltype new_sphere);

// split this method into read and update
//void read_update_config (int accep_tag,int win_id,char* fname, particle pobj,cellblock bobj);

void read_update_config (int win_id,const char* fname, particle pobj,cellblock bobj); // updated once the acceptance condition is met.

//celltype read_config(int win_id,char* fname, particle pobj,cellblock bobj); // read new config and return new sphere cell

//void update_config(int accep_tag);

//void construct_sphere (particle pobj, cellblock bobj, int win_id,char* filename); // calling method - sampler would have signatures

celltype construct_sphere (particle pobj, cellblock bobj, int win_id,const char* filename); // calling method - sampler would have signatures

void make_mc_nblist(celltype cobj);

void create_maxwell_velocities(cellblock loc_obj,double mc_temp, ivec3d*  mc_restriction);

ivec3d get_cpu_gcoord(int myrank);

ivec3d get_cell_loc_coord(ivec3d glob_coord, ivec3d cpu_glob_pos);

void make_particles(cellblock loc_obj);

void fill_mc_container(cellblock loc_obj);

void make_cells(cellblock loc_obj);

void count_particles(cellblock loc_obj); // verify if required

void  calc_mc_cpu_box(void);

void calc_mc_global_cell_array(void);

void calc_mc_cpu_cell_array(void);

void make_mc_tbox(void);

ivec3d cell_coordinate(double x, double y, double z);

//void calc_cpu_block_limits();

void calc_cell_dim(double rsample);

int calc_ncells_cpu(void);

int calc_ncells_stack(void);

int calc_cpu_col_total(void);

int calc_cpu_row_total(void);

int calc_cpu_stack_total(void);

void calc_block_dim(void);

ivec3d calc_cpu_coord(ivec3d cell_coord);

double get_gaussian(double sigma);

double scalar_prod(vec3d u, vec3d v);

int iscalar_prod(ivec3d u, ivec3d v);

vec3d vector_prod(vec3d u, vec3d v);

ivec3d ivector_prod(ivec3d u, ivec3d v);

double distance_vect(vec3d v1,vec3d v2);

#endif /* MONTE_PROTOTYPES_H_ */
