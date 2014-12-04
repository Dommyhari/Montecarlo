/*****************************************************************************************************
 *
 *                                 Generic Montecarlo API
 *  specifications to be included
 *  Created on: Sep 4, 2014
 *      Author: ganeshfx
 ****************************************************************************************************/

/*  MC routine */

#include "Monte_classes.h"     // class declarations
#include "Monte_globals.h"     // global variable definitions
#include "Monte_prototypes.h"  // Internal method declarations
#include "api.h"               // Interface method declarations

using namespace std;

// seems OK
extern "C" void init_montecarlo(int* md_cpu_dim, int md_tot_types,int md_real_types, int* md_restriction,double* md_simbox,
		double md_temperature,double rsweep,double wall_thick,int sample_seed,MPI_Comm comm){

    // getting cpu dimension
    mc_cpu_dim.x = *(md_cpu_dim);
    mc_cpu_dim.y = *(md_cpu_dim++);
    mc_cpu_dim.z = *(md_cpu_dim++);
    
    // total types
    mc_tot_types = md_tot_types;

    // real types
    mc_real_types = md_real_types;

    // mc simbox dimensions
    mc_simbox_x.x = *(md_simbox);   mc_simbox_x.y = 0.0;              mc_simbox_x.z = 0.0;
    mc_simbox_y.x = 0.0;            mc_simbox_y.y = *(md_simbox++);   mc_simbox_y.z = 0.0;
    mc_simbox_z.x = 0.0;            mc_simbox_z.y = 0.0;              mc_simbox_z.z = *(md_simbox++);

    // restriction vector for all element types
    int count=0;
    for(int i=0; i<mc_tot_types; i++){
    	mc_restriction[i].x = *(md_restriction+count); count++;
    	mc_restriction[i].y = *(md_restriction+count); count++;
    	mc_restriction[i].z = *(md_restriction+count); count++;
    }

    // MC simulation and sampling data
    mc_temp         =  md_temperature;      // temperature T
    mc_rsweep       =  rsweep;              // sphere radius
    mc_sphere_wall  =  wall_thick;          // sphere boundary thickness
    mc_sample_seed  =  sample_seed;         // seed for sampling random generator

    // assign communicator
    comm_name = comm;

    // no of process
    MPI_Comm_size(comm_name, &mc_ncpus);
    
    // process rank
    MPI_Comm_rank(comm_name, &mc_prank);
    
    // prepare Montecarlo configuration
    setup_config();
}

// seems OK
extern "C" void pack_config_to_montecarlo(long md_mc_tatoms_cpu,long *md_mc_atomnumber,int *md_mc_atomtypes,
		double  *md_mc_atommass,double  *md_mc_positions, double *md_mc_epot){

	    long m=0; // index

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
		if(mc_atomnumber.empty()||mc_atomtypes.empty()||mc_atommass.empty()||mc_positions.empty()||mc_epot.empty()){
			std::cerr<<"Montecarlo: problem allocating datas from MD "<<endl;
		}
		else if((mc_atomnumber.size()==mc_tatoms_cpu)&&(mc_atomtypes.size()==mc_tatoms_cpu)&&(mc_atommass.size()==mc_tatoms_cpu)
				&&((mc_positions.size()/3)==mc_tatoms_cpu) && (mc_epot.size()==mc_tatoms_cpu)){
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

// NEED TESTING ( acceptance condition part)

extern "C" void do_montecarlo(int* md_pid,long *md_tatoms_cpu,long **md_atomnumber,int **md_atomtypes,
    		double **md_atommass,double **md_positions){


	int accept_flag = 0; // acceptance flag

	// creating cellblock object
	cellblock c_obj;

	// construct cells
    make_cells(c_obj);

    // create particles
    make_particles(c_obj);

    // create velocities (T != 0K)
    if(!mc_temp) create_maxwell_velocities(c_obj,mc_temp,mc_restriction);


    //**************************************************************
    //!!!!!!!!!!!!!!!!!!!!!!!!!! check if required
    //**************************************************************

    // select random cell ex: cell-0 // for moving window cell after cell
    celltype cell_obj = c_obj.get_cell(0);

    // construct neighbor list for chosen cell
    make_mc_nblist(cell_obj);

    // sample window id
    int sample_id= 0;
    // particle instance
    particle sam_particle = sample_zone(c_obj,sample_id);


    //run local MD
    do_local_mdrun(md_binary,md_param);

    // if accept update
    //update_particle(cell_obj.get_cell_id(),sam_particle);

    //*****************************************
    //          some energy computation as per ensemble definitions (to initiate acceptance flag!!)
    //*****************************************


    // reading updated configuration after simulation
    read_update_config (accept_flag,sample_id,file_name,sam_particle,c_obj);

    // fill mc container
    fill_mc_container(c_obj);

    // clear created objects

    for (long count=0; count < c_obj.get_ncells(); count++){
        // clear all particles in cell
    	c_obj.get_cell(count).clear_all_particles();
    }
    // clear all cells in cellblock
    c_obj.clear_all_cells();

	// Data mirroring procedure - Assigning pointers

	md_pid         = &mc_prank;              // process rank
	md_tatoms_cpu  = &mc_tatoms_cpu;         // current total particles

	// preliminary pointer assignment
	*md_atomnumber = mc_atomnumber.data();   // atom id
    *md_atomtypes  = mc_atomtypes.data();    // atom types
	*md_atommass   = mc_atommass.data();     // atom mass
    *md_positions  = mc_positions.data();    // atom positions

}

// seems OK
extern "C" void get_velocities(double **md_velocities){
    // based on velocity request flag

	//if(mc_get_velocity) *md_velocities = mc_velocities.data();

	*md_velocities = mc_velocities.data();
}

// seems OK
extern "C" void get_pot_energy(double **md_epot){
    // get potential energy
    *md_epot       = mc_epot.data();
}

// seems OK
extern "C" void clean_montecarlo(){
	// clear or delete all data structures created in MC scope

	std::cout<<" **********************************************"<<endl;
	std::cout<<"  Montecarlo: clearing data structures "<<endl;
	std::cout<<" **********************************************"<<endl;

	mc_atomnumber.clear();
	mc_atomtypes.clear();
	mc_atommass.clear();
	mc_positions.clear();
	mc_velocities.clear();
	mc_epot.clear();

	cout<<  " Some Post check       " <<endl;
	cout<<  " size : mc_atomnumber  " << mc_atomnumber.size() << endl;
	cout<<  " size : mc_atomtypes   " << mc_atomtypes.size()  << endl;
	cout<<  " size : mc_atommass    " << mc_atommass.size()   << endl;
	cout<<  " size : mc_positions   " << mc_positions.size()  << endl;
	cout<<  " size : mc_positions   " << mc_velocities.size() << endl;
	cout<<  " size : mc_epot        " << mc_epot.size()       << endl;
	cout<<" ------------------------------------------- "<<endl;
}

// seems OK
extern "C" void shut_down_montecarlo(){

    // reinitialize cpu dimension
    mc_cpu_dim.x = 0; mc_cpu_dim.y = 0; mc_cpu_dim.z = 0;

    // total types           // real types
    mc_tot_types = 0;       mc_real_types = 0;

    // mc simbox dimensions
    mc_simbox_x.x = 0.0;   mc_simbox_x.y = 0.0; mc_simbox_x.z = 0.0;
    mc_simbox_y.x = 0.0;   mc_simbox_y.y = 0.0; mc_simbox_y.z = 0.0;
    mc_simbox_z.x = 0.0;   mc_simbox_z.y = 0.0; mc_simbox_z.z = 0.0;

    // restriction vector for all element types
    int count=0;

    for(int i=0; i<mc_tot_types; i++){

    	mc_restriction[i].x = 0; count++;
    	mc_restriction[i].y = 0; count++;
    	mc_restriction[i].z = 0; count++;

    }

    // MC simulation and sampling data
    mc_temp         =  0.0;        // temperature T
    mc_rsweep       =  0.0;        // sphere radius
    mc_sphere_wall  =  0.0;        // sphere boundary thickness
    mc_sample_seed  =  0;          // seed for sampling random generator

    // NOTE: look for further possibilities

}
