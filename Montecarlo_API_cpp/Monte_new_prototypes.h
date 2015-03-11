/*
 * Monte_new_prototypes.h
 *
 *  Created on: Mar 5, 2015
 *      Author: ganeshfx
 */

#ifndef MONTE_NEW_PROTOTYPES_H_
#define MONTE_NEW_PROTOTYPES_H_


#include "Monte_new_classes.h"

vec3d calc_cell_dim(double rsample, vec3d mc_simbox_dim );

vec3d calc_mc_cpu_box_vector(vec3d simbox_vec, int dim);

vec3d* make_mc_tbox_vector(vec3d* simbox,int prank);

ivec3d calc_mc_global_cell_array(vec3d simbox_diag,vec3d cell_dim);

ivec3d calc_mc_cpu_cell_array(vec3d cpubox_diag,vec3d cell_dim, ivec3d cpu_dim,ivec3d global_cell_dim);

ivec3d get_cpu_gcoord(int myrank,MPI_Comm c_name);

ivec3d calc_block_dim(ivec3d global_cell_dim,ivec3d cpu_cell_dim);

int calc_ncells_cpu(vec3d cpu_box_diag,vec3d cell_dim);

ivec3d calc_cpu_cell_division(vec3d cpu_cd_diag, vec3d cell_dim );

ivec3d get_cell_loc_coord(ivec3d glob_coord, ivec3d cpu_glob_pos, ivec3d cpu_cell_dim );

ivec3d cell_coordinate(double* pos,double* tbox,ivec3d global_cell_dim);

cellblock make_cells(cellblock loc_obj,vec3d cpu_box_diag, vec3d cell_dim,ivec3d gcoord,ivec3d cpu_cell_dim, int prank);

cellblock make_particles(cellblock loc_obj, long tatoms_cpu, double* tbox_dim, ivec3d global_cell_dim, ivec3d loc_cpu_gcoord,ivec3d cpu_cell_dim,
		vector<long> atomnumber,vector<int> atomtypes,vector<double> atommass,vector<double> positions,vector<double> epot,int prank,int debug);

cellblock create_maxwell_velocities(cellblock loc_obj,double temp, ivec3d*  restriction, int prank,int debug);

double get_gaussian(double sigma);

double scalar_prod(vec3d u, vec3d v);

int iscalar_prod(ivec3d u, ivec3d v);

vec3d vector_prod(vec3d u, vec3d v);

ivec3d ivector_prod(ivec3d u, ivec3d v);

double distance_vect(vec3d v1,vec3d v2);

#endif /* MONTE_NEW_PROTOTYPES_H_ */