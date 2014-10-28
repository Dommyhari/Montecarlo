/*
 * Monte_prototypes.h
 *
 *  Created on: Oct 23, 2014
 *      Author: ganeshfx
 */

#ifndef MONTE_PROTOTYPES_H_
#define MONTE_PROTOTYPES_H_

#include "Monte_globals.h"
#include "Monte_classes.h"

// utility methods definition

particle get_target_particle(long particle_id, int  cell_id);

void read_update_config (char* fname);

void construct_sphere (particle pobj, celltype cobj, char *fname); // calling method - sampler would have signatures

void make_mc_nblist(celltype cobj);

void create_maxwell_velocities(double mc_temp, ivec3d*  mc_restriction);

ivec3d get_cpu_gcoord(int myrank);

ivec3d get_cell_loc_coord(ivec3d glob_coord, ivec3d cpu_glob_pos);

void make_particles(void);

void make_cells(void);

void  calc_mc_cpu_box(void);

void calc_mc_global_cell_array(void);

void calc_mc_cpu_cell_array(void);

void make_mc_tbox(void);

ivec3d cell_coordinate(double x, double y, double z);

void calc_cpu_block_limits();

vec3d calc_cell_dim(double rcut,double rsample);

int calc_ncells_cpu(void);

int calc_ncells_stack(void);

int calc_cpu_col_total(void);

int calc_cpu_row_total(void);

int calc_cpu_stack_total(void);

void calc_block_dim(void);

int get_cpu_rank(ivec3d cell_coord);

double get_gaussian(double sigma);

double scalar_prod(vec3d u, vec3d v);

int iscalar_prod(ivec3d u, ivec3d v);

vec3d vector_prod(vec3d u, vec3d v);

ivec3d ivector_prod(ivec3d u, ivec3d v);

double distance_vect(vec3d v1,vec3d v2);



#endif /* MONTE_PROTOTYPES_H_ */
