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

ivec3d calc_mc_global_cell_array(vec3d simbox_diag,vec3d cell_dim);

ivec3d calc_mc_cpu_cell_array(vec3d cpubox_diag,vec3d cell_dim, ivec3d cpu_dim,ivec3d global_cell_dim);

ivec3d get_cpu_gcoord(int myrank,MPI_Comm c_name);

ivec3d calc_block_dim(ivec3d global_cell_dim,ivec3d cpu_cell_dim);

int calc_ncells_cpu(vec3d cpu_box_diag,vec3d cell_dim);

ivec3d calc_cpu_cell_division(vec3d cpu_cd_diag, vec3d cell_dim );

ivec3d get_cell_loc_coord(ivec3d glob_coord, ivec3d cpu_glob_pos, ivec3d cpu_cell_dim );

ivec3d get_particle_glob_coordinate(double* pos,vec3d mc_cell_phy_dim);

void make_cells(cellblock** loc_obj,vec3d cpu_box_diag, vec3d cell_dim,ivec3d gcoord,ivec3d cpu_cell_dim, int prank);

void make_particles(cellblock** loc_obj, long tatoms_cpu, ivec3d global_cell_dim, ivec3d loc_cpu_gcoord,ivec3d cpu_cell_dim,
		vector<long> atomnumber,vector<int> atomtypes,vector<double> atommass,vector<double> positions,vector<double> epot,int prank,int debug,vec3d cell_dim);

void create_maxwell_velocities(cellblock** loc_obj,double temp, ivec3d*  restriction, int prank,int debug);

double get_gaussian(double sigma);

int get_cpu_rank(int inde0,int inde1,int inde2, MPI_Comm c_name);

particle sample_zone(cellblock** bobj,int win_id,ivec3d cpu_cell_dim,int prank);

void construct_sphere(particle pobj, cellblock** bobj, int win_id,int prank,MPI_Comm comm_name, MPI_Status stat
		,int test_rank,ivec3d cpu_dim,double* data_list,ivec3d cpu_cell_dim,celltype** catalog);

void do_local_mdrun(string bin_name,string param_name,int prank);

int acceptance_check(int type,celltype old_sphere,celltype new_sphere);

void read_update_config (int win_id,particle pobj,cellblock** bobj,double* data_list,int prank,int test_rank,
		vec3d cell_dim,ivec3d cpu_cell_dim,MPI_Status stat,MPI_Comm comm_name,celltype** catalog);

double scalar_prod(vec3d u, vec3d v);

int iscalar_prod(ivec3d u, ivec3d v);

vec3d vector_prod(vec3d u, vec3d v);

ivec3d ivector_prod(ivec3d u, ivec3d v);

double distance_vect(vec3d v1,vec3d v2);

#endif /* MONTE_NEW_PROTOTYPES_H_ */
