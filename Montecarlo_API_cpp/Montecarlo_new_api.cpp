/*
 * Montecarlo_new_api.cpp
 *
 *  Created on: Mar 6, 2015
 *      Author: ganeshfx
 */

#include "Monte_new_prototypes.h"  // Internal method declarations
#include "Monte_new_globals.h"  // global variable definitions
#include "Monte_api.h"         // Interface method declarations

using namespace std;

// Tested
extern "C" void init_montecarlo(int* md_cpu_dim, int md_tot_types,int md_real_types, int* md_restriction,double* md_simbox,
		double md_temperature,double rsweep,double wall_thick,int sample_seed,MPI_Comm comm){

	MPI_Barrier(comm); // for ordered printing

    // getting cpu dimension
    mc_cpu_dim.x = *(md_cpu_dim++);
    mc_cpu_dim.y = *(md_cpu_dim++);
    mc_cpu_dim.z = *(md_cpu_dim++);

    // total types
    mc_tot_types = md_tot_types;

    // real types
    mc_real_types = md_real_types;

    // mc simbox dimensions
    mc_simbox_x.x = *(md_simbox++); mc_simbox_x.y = *(md_simbox++);              mc_simbox_x.z = *(md_simbox++);
    mc_simbox_y.x = *(md_simbox++); mc_simbox_y.y = *(md_simbox++);              mc_simbox_y.z = *(md_simbox++);
    mc_simbox_z.x = *(md_simbox++); mc_simbox_z.y = *(md_simbox++);              mc_simbox_z.z = *(md_simbox++);

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

    // compute Montecarlo cell dimension

    vec3d loc_cell_dim;
    loc_cell_dim.x = mc_simbox_x.x;
    loc_cell_dim.y = mc_simbox_y.y;
    loc_cell_dim.z = mc_simbox_z.z;


    // MonteCarlo cell dimension
    mc_cell_dim.x = mc_rsweep + mc_sphere_wall ; mc_cell_dim.y = mc_rsweep  + mc_sphere_wall; mc_cell_dim.z = mc_rsweep  + mc_sphere_wall;

    //mc_cell_dim.x = mc_rsweep; mc_cell_dim.y = mc_rsweep; mc_cell_dim.z = mc_rsweep;

    // compute cpu box physical dimensions

    mc_cpu_box_x = calc_mc_cpu_box_vector(mc_simbox_x, mc_cpu_dim.x); // CPU box x vector
    mc_cpu_box_y = calc_mc_cpu_box_vector(mc_simbox_y, mc_cpu_dim.y); // CPU box y vector
    mc_cpu_box_z = calc_mc_cpu_box_vector(mc_simbox_z, mc_cpu_dim.z); // CPU box z vector

	// compute global cell array
    vec3d simbox_diag;
    simbox_diag.x = mc_simbox_x.x; simbox_diag.y = mc_simbox_y.y; simbox_diag.z = mc_simbox_z.z;

    mc_global_cell_dim = calc_mc_global_cell_array(simbox_diag,mc_cell_dim);

    // CPU cell array computation
    vec3d cpubox_diag;
    cpubox_diag.x = mc_cpu_box_x.x; cpubox_diag.y = mc_cpu_box_y.y; cpubox_diag.z = mc_cpu_box_z.z;

    mc_cpu_cell_dim = calc_mc_cpu_cell_array(cpubox_diag,mc_cell_dim, mc_cpu_dim,mc_global_cell_dim);

	// change it accordingly to compute from Virtual topology
	cpu_gcoord = get_cpu_gcoord(mc_prank,comm_name);

    // computing block dimension -- could be replaced by MPI Cart option
    ivec3d bloc_dim;

    bloc_dim = calc_block_dim(mc_global_cell_dim,mc_cpu_cell_dim);

    cols_block = bloc_dim.x; rows_block = bloc_dim.y; stacks_block = bloc_dim.z;


    if (mc_prank == 0) {

       cout << "*****************************************************" << endl;
       cout << "             Values  pass  check   "  << endl;
       cout << " mc_tot_types :  " << mc_tot_types   << endl;
       cout << " mc_real_types : " << mc_real_types  << endl;
       cout << "*****************************************************" << endl;

       cout << "***********************************"<< endl;
       cout << " Parameter passing check -- init_montecarlo (Montecarlo_new_api.cpp) " << endl;
       cout << "***********************************"<< endl;
       cout << "mc_CPU dimension [ :"<< mc_cpu_dim.x << mc_cpu_dim.y << mc_cpu_dim.z << "]" <<endl;
       cout << " simbox dimension x : [" << mc_simbox_x.x << " "<<mc_simbox_x.y <<" "<< mc_simbox_x.z << " " << " ]" <<endl;
       cout << " simbox dimension y : [" << mc_simbox_y.x << " "<<mc_simbox_y.y <<" "<< mc_simbox_y.z << " " << " ]" <<endl;
       cout << " simbox dimension z : [" << mc_simbox_z.x << " "<<mc_simbox_z.y <<" "<< mc_simbox_z.z << " " << " ]" <<endl;

       for(int i=0; i<mc_tot_types; i++){
          	cout<<"Restriction [ "<<i<<" ] : "<<mc_restriction[i].x<<" "<<mc_restriction[i].y<<" "<<mc_restriction[i].z<<endl;
       }
       cout << " mc_temperature : " << mc_temp << endl;
       cout << " mc_sweep_radii : " << mc_rsweep << endl;
       cout << " wall thickness : " << mc_sphere_wall << endl;
       cout << " random no : " << mc_sample_seed << endl;

       cout << " ==================================="<<endl;
       cout << " Montecarlo_new_api.cpp (calc_cell_dim -- after method call)    " << endl;
       cout << " mc_cell_dim check " << endl;
       cout << " mc_cell_dim.x    :  " <<  mc_cell_dim.x << endl;
       cout << " mc_cell_dim.y    :  " <<  mc_cell_dim.y << endl;
       cout << " mc_cell_dim.z    :  " <<  mc_cell_dim.z << endl;
       cout << " ==================================="<<endl;

       cout << " ==================================="<<endl;
       cout << " Montecarlo_new_api.cpp (calc_mc_cpu_box_vector -- after method call)    " << endl;
       cout << " mc_cpu_box check " << endl;
       cout << " mc_cpu_box_x    :  [ " <<  mc_cpu_box_x.x <<" "<<mc_cpu_box_x.y<<" "<<mc_cpu_box_x.z<<" "<<"]"<< endl;
       cout << " mc_cpu_box_y    :  [ " <<  mc_cpu_box_y.x <<" "<<mc_cpu_box_y.y<<" "<<mc_cpu_box_y.z<<" "<<"]"<< endl;
       cout << " mc_cpu_box_z    :  [ " <<  mc_cpu_box_z.x <<" "<<mc_cpu_box_z.y<<" "<<mc_cpu_box_z.z<<" "<<"]"<< endl;
       cout << " ==================================="<<endl;
       cout << " ==================================="<<endl;
       cout << " Montecarlo_new_api.cpp (calc_mc_global_cell_array -- after method call)    " << endl;
       cout << " mc_global_cell_dim check " << endl;
       cout << " mc_global_cell_dim.x    :  " <<  mc_global_cell_dim.x << endl;
       cout << " mc_global_cell_dim.y    :  " <<  mc_global_cell_dim.y << endl;
       cout << " mc_global_cell_dim.z    :  " <<  mc_global_cell_dim.z << endl;
       cout << " ==================================="<<endl;
       cout << " Montecarlo_new_api.cpp (calc_mc_cpu_cell_array -- after method call)    " << endl;
       cout << " mc_cpu_cell_dim check " << endl;
       cout << " mc_cpu_cell_dim.x    :  " <<  mc_cpu_cell_dim.x << endl;
       cout << " mc_cpu_cell_dim.y    :  " <<  mc_cpu_cell_dim.y << endl;
       cout << " mc_cpu_cell_dim.z    :  " <<  mc_cpu_cell_dim.z << endl;
       cout << " ==================================="<<endl;
       cout << " Montecarlo_new_api.cpp (get_cpu_gcoord -- after method call)    " << endl;
       cout << " cpu_gcoord dim check " << endl;
       cout << " cpu_gcoord.x    :  " <<  cpu_gcoord.x << endl;
       cout << " cpu_gcoord.y    :  " <<  cpu_gcoord.y << endl;
       cout << " cpu_gcoord.z    :  " <<  cpu_gcoord.z << endl;
       cout << " ==================================="<<endl;
       cout << " Montecarlo_new_api.cpp (calc_block_dim -- after method call)    " << endl;
       cout << " bloc_dim  check " << endl;
       cout << " cols_block    :  " <<  cols_block << endl;
       cout << " rows_block    :  " <<  rows_block << endl;
       cout << " stacks_block  :  " <<  stacks_block << endl;
       cout << " ==================================="<<endl;
    }


    MPI_Barrier(comm); // for ordered printing

}

// Tested
extern "C" void pack_config_to_montecarlo(long md_mc_tatoms_cpu,long *md_mc_atomnumber,int *md_mc_atomtypes,
		double  *md_mc_atommass,double  *md_mc_positions, double *md_mc_epot){

	    MPI_Barrier(comm_name); // for ordered printing

	    long m=0; // index

		mc_tatoms_cpu = md_mc_tatoms_cpu; // get total particles from MD

	    // allocate buffer size
	    buffer_size = mc_tatoms_cpu; // total no of particles + threshold value

		for(long i=0; i<mc_tatoms_cpu; i++){

			// filling data container with corresponding values
			mc_atomnumber.push_back(md_mc_atomnumber[i]);
			mc_atomtypes.push_back(md_mc_atomtypes[i]);
			mc_atommass.push_back(md_mc_atommass[i]);
            /*
			if ((mc_prank == 0)  && (i==0)){
				long n =0;
				n=m;
				cout << "***********************************************************" << endl;
				cout << "   check while packing " << endl;
				cout << " check for value i :  " << i << endl;
				cout << " md_mc_atomnumber[i] : " << md_mc_atomnumber[i] << endl;
				cout << " md_mc_atomtypes[i]  : " << md_mc_atomtypes[i]  << endl;
				cout << " md_mc_atommass[i]   : " << md_mc_atommass[i]   << endl;
				cout << " md_mc_positions[n++] : " << md_mc_positions[n++]  << endl;
				cout << " md_mc_positions[n++] : " << md_mc_positions[n++]  << endl;
				cout << " md_mc_positions[n++] : " << md_mc_positions[n++]  << endl;
				cout << " md_mc_epot[i]        : " << md_mc_epot[i] << endl;
				cout << "***********************************************************" << endl;
			}
            */

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
			std::cout<<" buffer_size : "<< buffer_size << endl;
			std::cout<<" Total particles received : " <<mc_tatoms_cpu<<endl;
			std::cout<<"Montecarlo: datas are allocated correctly " <<endl;
			std::cout<<" ------------------------------------------- "<<endl;
	    }
		else {
			std::cerr<<"Montecarlo: unknown error with allocation" <<endl;
		}
		MPI_Barrier(comm_name); // for ordered printing
}

// NEED TESTING ( acceptance condition part)

extern "C" void do_montecarlo(int md_pid,long *md_tatoms_cpu,long **md_atomnumber,int **md_atomtypes,
    		double **md_atommass,double **md_positions){

	MPI_Barrier(comm_name); // for ordered printing

	int accept_flag = 0; // acceptance flag
	int win_id = 0;  // Hardcoded for testing
	int test_cpu = 0;
	int dbug_flag = 0;   // debug flag -- print check statements


	// creating cellblock object
	cellblock c_obj, r_cell,r_partic, r_vel;

	// CPU box dimension
	vec3d cpu_box_diag;

	cpu_box_diag.x = mc_cpu_box_x.x; cpu_box_diag.y = mc_cpu_box_y.y; cpu_box_diag.z = mc_cpu_box_z.z;

	//ncells_cpu = calc_ncells_cpu(cpu_box_diag,mc_cell_dim);
	c_obj.set_mycpu(mc_prank);

	// construct cells
	r_cell = make_cells(c_obj,cpu_box_diag, mc_cell_dim,cpu_gcoord,mc_cpu_cell_dim,mc_prank);

	// create particles
    r_partic = make_particles(r_cell,mc_tatoms_cpu, mc_global_cell_dim,cpu_gcoord,mc_cpu_cell_dim,
        		mc_atomnumber,mc_atomtypes,mc_atommass,mc_positions,mc_epot,mc_prank,dbug_flag,mc_cell_dim);

    // clear STL container (once particle objects are created)

    // (NOT REQUIRED DURING TESTING PHASE -- INCLUDE LATER -- in sync with fill_mc_container)
	mc_atomnumber.clear();
	mc_atomtypes.clear();
    mc_atommass.clear();
    mc_positions.clear();
    mc_epot.clear();


    // create velocities (T != 0K)
	if(mc_temp>0){
		//cellblock  r_vel = new cellblock;

		r_vel = create_maxwell_velocities(r_partic,mc_temp,mc_restriction,mc_prank,dbug_flag);
		r_partic = r_vel;

		//delete r_vel;
	}


    // extend the fragment with random selector
    int sample_id= win_id;

    // particle instance -- sample zone method
    particle rand_particle = sample_zone(r_partic,sample_id,mc_cpu_cell_dim,mc_prank);

    // data_list for sphere construction routine
    double* dat_list;
    dat_list[0] = (mc_cpu_cell_dim.x/2 * mc_cpu_cell_dim.y/2 * mc_cpu_cell_dim.z/2 ); // no of cells
    dat_list[1] = mc_rsweep;        // sphere radius
    dat_list[2] = mc_sphere_wall;   // sphere wall thickness
    dat_list[3] = mc_simbox_x.x;
    dat_list[4] = mc_simbox_y.y;
    dat_list[5] = mc_simbox_z.z;
    dat_list[6] = (double) mc_real_types;


    // some test case examples

    if((mc_prank==0)) {
    	rand_particle.set_myposition(2.8617,2.8651,1.4384); // 80x80x80 -win 0
//    	//rand_particle.set_myposition(11.4372, 2.8651, 2.8677); // 300x300x300 -win 0
    }
//
    if( (mc_prank==1)  ){
    	rand_particle.set_myposition(2.8617,2.8651,41.4575);
    }
//
    if((mc_prank==2)){
    	rand_particle.set_myposition(4.2909,41.4549,1.4384);
    }
//
    if((mc_prank==3)){
    	rand_particle.set_myposition(2.8617,42.8842,41.4575);
    }
//
    if((mc_prank==4)){
    	rand_particle.set_myposition(44.3100,2.8651,2.8677);
    }
//
    if((mc_prank==5)){
    	rand_particle.set_myposition(42.8808,2.8651,41.4575);
//    	rand_particle.set_myposition(60.0318,21.4454,61.4671); // critical check
    }
//
    if((mc_prank==6)){
    	rand_particle.set_myposition(42.8808,42.8842,1.4384);
    }
//
    if((mc_prank==7)){
    	rand_particle.set_myposition(42.8808,42.8842,41.4575);
    }

	// debug/control checking part

    if(mc_prank == test_cpu){


       long tot_particles=0; // check variable
       cout << "****************************************************" << endl;

       for(int i=0; i<r_partic.get_cell_list_size(); i++){
    	     //cout << " No of particles in cell id  [ " << i << " ] : "<< r_partic.get_cell(i).get_nparticles()  << endl;
    	     tot_particles += r_partic.get_cell(i).get_nparticles();

       }

       cout << " Total particles allocated after MAKE CELLS in CPU 0 " << tot_particles << endl;
       cout << "****************************************************" << endl;

   	   tot_particles=0; // check variable
       cout << "****************************************************" << endl;

       for(int i=0; i<r_vel.get_cell_list_size(); i++){
   	     //cout << " No of particles in cell id  [ " << i << " ] : "<< r_partic.get_cell(i).get_nparticles()  << endl;
   	     tot_particles += r_vel.get_cell(i).get_nparticles();

       }
       cout << " Total particles allocated after VELOCITY COMPUTATION in CPU 0 " << tot_particles << endl;
       cout << "****************************************************" << endl;

       cout << " After sample zone test " << endl;
	   cout << "  Randomly selected placeholder scope check  " << endl;
	   cout << "  Random particle id  : " << rand_particle.get_mynumber() << endl;
	   cout << "  Random particle type: " << rand_particle.get_mytype() << endl;
	   cout << "  Random particle mass  " << rand_particle.get_mymass() << endl;
	   vec3d loc_pos,loc_vel;
	   loc_pos = rand_particle.get_myposition();
	   loc_vel = rand_particle.get_myvelocity();
	   cout << "  Random particle position :  [ " << loc_pos.x <<" "<<loc_pos.y<<" "<<loc_pos.z<<" "<< " ]"<<endl;
       cout << "  Random particle velocity :  [ " << loc_vel.x <<" "<<loc_vel.y<<" "<<loc_vel.z<<" "<< " ]"<<endl;
       cout << "  Random particle Epot       : " << rand_particle.get_myepot() << endl;

	}


    celltype sphr, nbr_1,nbr_2,nbr_3,nbr_4,nbr_5,nbr_6,nbr_7,my_list;

    celltype* ptr_list[9]={&sphr,&nbr_1,&nbr_2,&nbr_3,&nbr_4,&nbr_5,&nbr_6,&nbr_7,&my_list};

    // construct sphere

    cellblock aft_sphere = construct_sphere(rand_particle, r_partic, win_id,file_name,mc_prank,comm_name,status,test_cpu,mc_cpu_dim,dat_list,
   		mc_cpu_cell_dim,ptr_list);

    // run Local MD
    do_local_mdrun(md_binary,md_param,mc_prank);

    // read or update config

    cellblock loc_obj = read_update_config(win_id,rand_particle,aft_sphere,dat_list,mc_prank,test_cpu,mc_cell_dim,mc_cpu_cell_dim,status,comm_name,ptr_list);

    // fill mc container

    // edited // fill_mc_container(c_obj); // (this method could be defined locally here!!)

    vector<int> mod_mc_atomtypes;        //  vector container for atom types
    vector<long> mod_mc_atomnumber;      //  vector container for atom number
    vector<double> mod_mc_atommass;      //  vector container for atom mass
    vector<double> mod_mc_positions;     //  vector container for atom positions
    vector<double> mod_mc_velocities;    //  vector container for atom velocities
    vector<double> mod_mc_epot;          //  vector container for atom potential energy


//	// loop over cells
//
//    celltype main_cell;
//    particle main_part;
//
//	for (long i=0; i<loc_obj.get_cell_list_size();i++){
//
//		par_count = loc_obj.get_cell(i).get_nparticles();
//
//		for(long j=0; j<par_count;j++ ){
//
//			main_part = loc_obj.get_cell(i).get_particle(j);
//
//            // add particle attributes
//			mod_mc_atomnumber.push_back(main_part.get_mynumber());    // atom id
//			mod_mc_atomtypes.push_back(main_part.get_mytype());       // atom type
//            mod_mc_atommass.push_back(main_part.get_mymass());        // atom mass
//			mod_mc_positions.push_back(main_part.get_myposition().x);
//			mod_mc_positions.push_back(main_part.get_myposition().y); // atom positions
//			mod_mc_positions.push_back(main_part.get_myposition().z);
//
//			mod_mc_velocities.push_back(main_part.get_myvelocity().x);
//			mod_mc_velocities.push_back(main_part.get_myvelocity().y); // atom velocities
//			mod_mc_velocities.push_back(main_part.get_myvelocity().z);
//
//			mod_mc_epot.push_back(main_part.get_myepot());            // atom epot
//
//		}// loop over particle
//
//		// not required here -- declared in api using clear_all
//		// delete all particles in cell
//		loc_obj.get_cell(i).clear_all_particles();
//
//	}// loop over cells
////
//	// delete all cells
//	loc_obj.clear_all_cells();

    // clear created objects

    // edited //for (long count=0; count < c_obj.get_ncells(); count++){
        // clear all particles in cell
    // edited //	c_obj.get_cell(count).clear_all_particles();
    // edited // }

    // clear all cells in cellblock

    // edited // c_obj.clear_all_cells();

	// Data mirroring procedure - Assigning pointers

	md_pid         = mc_prank;              // process rank
	md_tatoms_cpu  = &mc_tatoms_cpu;         // current total particles

	// preliminary pointer assignment
	*md_atomnumber = mc_atomnumber.data();   // atom id
    *md_atomtypes  = mc_atomtypes.data();    // atom types
	*md_atommass   = mc_atommass.data();     // atom mass
    *md_positions  = mc_positions.data();    // atom positions

    MPI_Barrier(comm_name); // for ordered printing

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

	//MPI_Barrier(comm_name);
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

	//MPI_Barrier(comm_name);
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
