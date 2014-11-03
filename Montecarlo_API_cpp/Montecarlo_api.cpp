/*****************************************************************************************************
 *
 *                                 Generic Montecarlo API
 *  specifications to be included
 *  Created on: Sep 4, 2014
 *      Author: ganeshfx
 ****************************************************************************************************/

/*  MC routine */

#include<iostream>
#include<ctime>
#include<cmath>
#include<cstdlib>
#include<vector>
#include "Monte_classes.h" // class definitions
#include "Monte_globals.h" // method definitions with extern
#include "Monte_prototypes.h"

using namespace std;


extern "C" void init_montecarlo(int* md_cpu_dim, int* md_types, int* md_restriction,double *md_simbox,double* md_data,MPI_Comm comm_name){

          // think of arguments required for setting up Montecarlo
	      // declare here globally or in Monte_globals.h
          // from xcode
    int ierror;
    
    // getting cpu dimension
    mc_cpu_dim.x = *(md_cpu_dim);
    mc_cpu_dim.y = *(md_cpu_dim++);
    mc_cpu_dim.z = *(md_cpu_dim++);
    
    // total types
    mc_tot_types = *(md_types);

    // real types
    mc_real_types = *(md_types++);

    // mc simbox dimensions
    mc_simbox_x.x = *(md_simbox); mc_simbox_x.y = 0.0; mc_simbox_x.z = 0.0;
    mc_simbox_y.x = 0.0; mc_simbox_y.y = *(md_simbox++); mc_simbox_y.z = 0.0;
    mc_simbox_z.x = 0.0; mc_simbox_z.y = 0.0; mc_simbox_z.z = *(md_simbox++);

    // restriction vector for all element types
    int count=0;
    for(int i=0; i<mc_tot_types; i++){

    	mc_restriction[i].x = *(md_restriction+count); count++;
    	mc_restriction[i].y = *(md_restriction+count); count++;
    	mc_restriction[i].z = *(md_restriction+count); count++;

    }

    // MC simulation and sampling data
    mc_temp         = *(md_data);           // temperature T
    mc_rsweep       = *(md_data++);         // sphere radius
    mc_sphere_wall  = *(md_data++);         // sphere boundary thickness
    mc_sample_seed  = (int) *(md_data++);   // seed for sampling random generator

    // no of process
    ierror = MPI_Comm_size(comm_name, &mc_ncpus);
    
    // process rank
    ierror = MPI_Comm_rank(comm_name, &mc_prank);
    
}

extern "C" void pack_config_to_montecarlo(long md_mc_tatoms_cpu,long *md_mc_atomnumber,int *md_mc_atomtypes,
		double  *md_mc_atommass,double  *md_mc_positions, double *md_mc_epot){

	    long rand_no,m=0; // random number variable

		mc_tatoms_cpu = md_mc_tatoms_cpu; // get total particles from MD

		for(long i=0; i<mc_tatoms_cpu; i++){

			// filling data container with corresponding values
			mc_atomnumber.push_back(md_mc_atomnumber[i]);
			mc_atomtypes.push_back(md_mc_atomtypes[i]);
			mc_atommass.push_back(md_mc_atommass[i]);

			mc_positions.push_back(md_mc_positions[m++]); // holds x
			mc_positions.push_back(md_mc_positions[m++]); // holds y
			mc_positions.push_back(md_mc_positions[m++]); // holds z

			// similarly add epot
			mc_epot.push_back(md_mc_epot[i]);

		}

		// check for unallocated vector
		if(mc_atomnumber.empty()||mc_atomtypes.empty()||mc_atommass.empty()||mc_positions.empty()){
			std::cerr<<"Montecarlo: problem allocating datas from MD "<<endl;
		}
		else if((mc_atomnumber.size()==mc_tatoms_cpu)&&(mc_atomtypes.size()==mc_tatoms_cpu)&&(mc_atommass.size()==mc_tatoms_cpu)
				&&((mc_positions.size()/3)==mc_tatoms_cpu)){

			std::cout<<" ------------------------------------------- "<<endl;
			std::cout<<" Process id : "<<mc_prank<<endl;
			std::cout<<" Total particles received : " <<mc_tatoms_cpu<<endl;
			std::cout<<"Montecarlo: datas are allocated correctly " <<endl;
			std::cout<<" ------------------------------------------- "<<endl;
	    }
		else {
			std::cerr<<"Montecarlo: unknown error with allocation" <<endl;
		}


}

extern "C" void do_montecarlo(int* md_pid,long *md_tatoms_cpu,long **md_atomnumber,int **md_atomtypes,
		double **md_atommass,double **md_positions,double **md_velocities ,double **md_epot)
{



    // compute Montecarlo cell dimension
	calc_cell_dim(mc_rsweep);

    // compute cpu box physical dimensions
	calc_mc_cpu_box();

	// construct transformation box
	make_mc_tbox();

	// compute global cell array
	calc_mc_global_cell_array();

	// compute cpu cell array
	calc_mc_cpu_cell_array();

	// all cpu global position
	get_cpu_gcoord(mc_prank);

	// compute domain block dimension
	calc_block_dim();

	// creating cellblock object
	cellblock c_obj;


	// construct cells
    make_cells(c_obj);

    // create particles
    make_particles(c_obj);

    // create velocities (T != 0K)
    if(!mc_temp) create_maxwell_velocities(c_obj,mc_temp,mc_restriction);

    // select random cell ex: cell-0 // for moving window cell after cell
    celltype cell_obj = c_obj.get_cell(0);

    // construct neighbor list for chosen cell
    make_mc_nblist(cell_obj);

    // particle instance
    particle sam_particle = sample_zone(c_obj.get_cell(0));


    //run local MD
    do_local_mdrun(md_binary,md_param);

    // if accept update
    update_particle(cell_obj.get_cell_id(),sam_particle);

	// Data mirroring procedure - Assigning pointers

	md_pid=&mc_prank;
	md_tatoms_cpu = &mc_tatoms_cpu;

	// preliminary pointer assignment

	*md_atomnumber = mc_atomnumber.data(); //error part
    *md_atomtypes  = mc_atomtypes.data();
	*md_atommass   = mc_atommass.data();  //they do have initialized pointers and memory
    *md_positions  = mc_positions.data();
    *md_epot       = mc_epot.data();

    // based on velocity request flag
    if(mc_get_velocity) *md_velocities = mc_velocities.data();


    // clear STL container
    mc_atomtypes.clear();


//    std::cout<<" **********************************************"<<endl;
    //    std::cout<<" Hello from do_montecarlo with process : "<<*md_pid<<endl;
    //std::cout<<" **********************************************"<<endl;
    //std::cout<<" some data check "<<endl;
    //std::cout<<" created random no : "<<rand_no<<endl;
    // std::cout<<" md_pid  :" <<*md_pid<<" mc_pid : "<<mc_prank<<endl;
    //std::cout<<" md_tatoms_cpu..  :" << *md_tatoms_cpu <<" mc_tatoms_cpu : "<<mc_tatoms_cpu<<endl;
    //std::cout<<" mc_atomnumber[rand_no]  :"<<mc_atomnumber[rand_no]<<endl;
    //std::cout<<" mc_atomtypes[rand_no]  :"<<mc_atomtypes[rand_no]<<endl;
    //std::cout<<" mc_atommass[rand_no]  :" <<mc_atommass[rand_no]<<endl;
    //std::cout<<" pos_x : mc_positions[rand_no]  :"<<mc_positions[rand_no*3]<<endl;
    //std::cout<<" pos_y : mc_positions[rand_no]  :"<<mc_positions[rand_no*3+1]<<endl;
    //std::cout<<" pos_z : mc_positions[rand_no]  :"<<mc_positions[rand_no*3+2]<<endl;
    //std::cout<<" ------------------------------------------- "<<endl;


}

extern "C" void clean_montecarlo(){

	// clear or delete all data structures created in MC scope

	std::cout<<" **********************************************"<<endl;
	std::cout<<"  Montecarlo: clearing data structures "<<endl;
	std::cout<<" **********************************************"<<endl;
	mc_atomnumber.clear();
	mc_atomtypes.clear();
	mc_atommass.clear();
	mc_positions.clear();

	std::cout<< "   Some Post check       " <<endl;
	std::cout<<  " size : mc_atomnumber  "<<mc_atomnumber.size()<<endl;
	std::cout<<  " size : mc_atomtypes  "<<mc_atomtypes.size()<<endl;
	std::cout<<  " size : mc_atommass  "<<mc_atommass.size()<<endl;
	std::cout<<  " size : mc_positions  "<<mc_positions.size()<<endl;
	std::cout<<" ------------------------------------------- "<<endl;
}
