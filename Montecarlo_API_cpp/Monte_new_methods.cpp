/*
 * Monte_new_methods.cpp
 *
 *  Created on: Mar 5, 2015
 *      Author: ganeshfx
 */

#include<iostream>
#include "Monte_new_prototypes.h"

using namespace std;

vec3d calc_cell_dim(double rsample, vec3d mc_simbox_dim ){

	//    compute MonteCarlo cell dimension

	vec3d mc_cell_loc;

	mc_cell_loc.x = (mc_simbox_dim.x) / (rsample);
	mc_cell_loc.y = (mc_simbox_dim.y) / (rsample);
	mc_cell_loc.z = (mc_simbox_dim.z) / (rsample);

	return mc_cell_loc;
}


vec3d calc_mc_cpu_box_vector(vec3d simbox_vec, int dim){

	vec3d cpu_box_vec;
	// cpu box  vector
	cpu_box_vec.x = simbox_vec.x / dim;
	cpu_box_vec.y = simbox_vec.y / dim;
	cpu_box_vec.z = simbox_vec.z / dim;

	return cpu_box_vec;

}

ivec3d calc_mc_global_cell_array(vec3d simbox_diag,vec3d cell_dim){

	//           MonteCarlo global cell array computation
	// get simbox physical dimension (simbox_diag = mc_simbox_x.x,mc_simbox_y.y,mc_simbox_z.z)

	ivec3d global_cell_dim;

	global_cell_dim.x = (int) floor((simbox_diag.x / cell_dim.x));
	global_cell_dim.y = (int) floor((simbox_diag.y / cell_dim.y));
	global_cell_dim.z = (int) floor((simbox_diag.z / cell_dim.z));

	return global_cell_dim;
}

ivec3d calc_mc_cpu_cell_array(vec3d cpubox_diag,vec3d cell_dim, ivec3d cpu_dim,ivec3d global_cell_dim ){

      //          MonteCarlo cpu cell array computation

	ivec3d cpu_cell_dim;

	cpu_cell_dim.x = (int) floor(cpubox_diag.x/ cell_dim.x);
	cpu_cell_dim.y = (int) floor(cpubox_diag.y/ cell_dim.y);
	cpu_cell_dim.z = (int) floor(cpubox_diag.z/ cell_dim.z);

	if( cpu_cell_dim.x != (global_cell_dim.x/cpu_dim.x)) std::cerr<<" Incompatible global and cell array dimension "<<endl;
	if( cpu_cell_dim.y != (global_cell_dim.y/cpu_dim.y)) std::cerr<<" Incompatible global and cell array dimension "<<endl;
	if( cpu_cell_dim.z != (global_cell_dim.z/cpu_dim.z)) std::cerr<<" Incompatible global and cell array dimension "<<endl;

	return cpu_cell_dim;
}


ivec3d get_cpu_gcoord(int myrank,MPI_Comm c_name){


	// compute CPU global coordinate based on process rank

	ivec3d my_coord;
	int grid_coord[3];

	MPI_Cart_coords(c_name,myrank,3,grid_coord);
	my_coord.x = grid_coord[0];
	my_coord.y = grid_coord[1];
	my_coord.z = grid_coord[2];

	return my_coord;
}

ivec3d calc_block_dim(ivec3d global_cell_dim,ivec3d cpu_cell_dim){

	ivec3d block_dim;
	// check if really required??
	//  method for computing domain block dimensions

	block_dim.x   = global_cell_dim.x / cpu_cell_dim.x;  // cols_block
	block_dim.y   = global_cell_dim.y / cpu_cell_dim.y;  // rows_block
	block_dim.z   = global_cell_dim.z / cpu_cell_dim.z;  // stacks_block

    return block_dim;
}

// do_montecarlo step

int calc_ncells_cpu(vec3d cpu_box_diag,vec3d cell_dim){

	 //   compute no of cells in CPU

	int ret = (int) floor((cpu_box_diag.x * cpu_box_diag.y * cpu_box_diag.z )/(cell_dim.x * cell_dim.y * cell_dim.z));
    return ret;
}

ivec3d calc_cpu_cell_division(vec3d cpu_cd_diag, vec3d cell_dim ){

	ivec3d cell_div;

	// compute no of cpu cell columns
	cell_div.x = (int) floor(cpu_cd_diag.x/ cell_dim.x);

	// compute no of cpu cell rows
	cell_div.y = (int) floor(cpu_cd_diag.y/ cell_dim.y);

	// compute no of cpu cell stacks
	cell_div.z = (int) floor(cpu_cd_diag.z/ cell_dim.z);

	return cell_div;
}

ivec3d get_cell_loc_coord(ivec3d glob_coord, ivec3d cpu_glob_pos, ivec3d cpu_cell_dim ){

	// compute cell local coordinates from cell global coordinate and cpu glob position

	ivec3d cell_loc_coord;

	cell_loc_coord.x  = glob_coord.x - (cpu_glob_pos.x * cpu_cell_dim.x);
	cell_loc_coord.y  = glob_coord.y - (cpu_glob_pos.y * cpu_cell_dim.y);
	cell_loc_coord.z  = glob_coord.z - (cpu_glob_pos.z * cpu_cell_dim.z);

	return cell_loc_coord;
}


ivec3d get_particle_glob_coordinate(double* pos,vec3d mc_cell_phy_dim){

	//  cell_coord computes the (global) cell coordinates of a position
	ivec3d coord;
	double x,y,z;                   // holds position
	int x_quo,y_quo,z_quo;       // holds quotient
	double x_rem,y_rem,z_rem;       // holds remainder

	// assignment
	x = pos[0]; y = pos[1]; z = pos[2];

	// quotient computation
	x_quo = (int) x/mc_cell_phy_dim.x ; y_quo = (int) y/mc_cell_phy_dim.y ; z_quo = (int) z/mc_cell_phy_dim.z ;

	// remainder computation
	//x_rem = x%mc_cell_phy_dim.x ; y_rem = y%mc_cell_phy_dim.y ; z_rem = z%mc_cell_phy_dim.z ;

	x_rem = fmod (x,mc_cell_phy_dim.x) ; y_rem = fmod(y,mc_cell_phy_dim.y) ; z_rem = fmod(z,mc_cell_phy_dim.z) ;

	if(x_rem > 0.0) { x_quo++; }
	if(y_rem > 0.0) { y_quo++; }
	if(z_rem > 0.0) { z_quo++; }


	coord.x = x_quo; coord.y = y_quo; coord.z = z_quo;

	return coord;

}

// Tested
cellblock make_cells(cellblock loc_obj,vec3d cpu_box_diag, vec3d cell_dim,ivec3d gcoord,ivec3d cpu_cell_dim, int prank){

    //    method for creating cell objects and fill in cellblock container
	//    should assign cell boundary and ids correctly from the field

	int cell_no     = 0;
	ivec3d cell_total;

	cellblock ret_obj;

	int line_counter=0;

	//  (could be hardcoded for Moving window approach as nxn )

	//int ncells_tot    = loc_obj.get_ncells();

	cell_total = calc_cpu_cell_division(cpu_box_diag,cell_dim);


	int col_total     = cell_total.x;
	int row_total     = cell_total.y;
	int stack_total   = cell_total.z;

    ivec3d cell_gcoord = {0,0,0}; // initialization of cell global coordinate
    ivec3d cell_lcoord = {0,0,0}; // initialization of cell local coordinate

    int cells_stack=0;
    int cell_counter = 0;

    //*****************************************************************
    // create cell objects and fill into cell block container eventually
    // NOTE: hierarchy constraint: before make_particles() make_cells() should be called
    //*****************************************************************

    // loop over stack


	for (int col_count=0; col_count<col_total; col_count++){

		for(int row_count=0; row_count<row_total; row_count++){

			for(int stack_no=0; stack_no<stack_total; stack_no++){

				  celltype cell_obj;

                  // cell id computation
        		  //cell_no = cell_counter + ( ncells_tot * prank ); // cell id computation

				  cell_no = cell_counter; // cell id computation

        		  // numbering sequence: x--> y --> z

        		  // cell global coordinates
        		  cell_gcoord.x = col_count +  (col_total * gcoord.x);
        		  cell_gcoord.y = row_count +  (row_total * gcoord.y);
        		  cell_gcoord.z = stack_no  +  (stack_total * gcoord.z);

        		  // cell local coordinate
        		  cell_lcoord = get_cell_loc_coord(cell_gcoord,gcoord,cpu_cell_dim);

        		  // assign values
       		   	  cell_obj.set_cell_id(cell_no);

       		      cell_obj.set_glob_coord(cell_gcoord); // setting global cell coordinate

       		      cell_obj.set_loc_coord(cell_lcoord);  // setting local cell coordinate

       		      loc_obj.add_cell(cell_obj); // adding cell into cel_block

       		      cell_counter++; //  incrementing cell counter

			}//column_count

		}//row_count

	}//stack loop

//    if(prank != 10) {
//      cout << " ========================================================================" << endl;
//      cout << " my process id :  " << prank << endl;
//      cout << " CPU coordinate  : " << gcoord.x <<" " << gcoord.y <<" "<<gcoord.z <<" "<< endl;
//      //cout << " Cell allocation check :  no of cells  :  " <<  loc_obj.get_cell_list_size() << endl;
//      //cout << " cell ID :" << loc_obj.get_cell(cell_no).get_cell_id();
//      //cout << "   cell global coordinate        "<< loc_obj.get_cell(cell_no).get_cell_glob_coord().x<<" "<<loc_obj.get_cell(cell_no).get_cell_glob_coord().y
//      //		  <<" "<<loc_obj.get_cell(cell_no).get_cell_glob_coord().z<<endl;
//      //cout << "   cell local coordinate        "<< loc_obj.get_cell(cell_no).get_cell_loc_coord().x<<" "<<loc_obj.get_cell(cell_no).get_cell_loc_coord().y
//      //		  <<" "<<loc_obj.get_cell(cell_no).get_cell_loc_coord().z<<endl;
//      cout << " ========================================================================" << endl;
//    }



	ret_obj = loc_obj;

	return ret_obj;

} // make_cells


// Tested
cellblock make_particles(cellblock loc_obj, long tatoms_cpu, ivec3d global_cell_dim, ivec3d loc_cpu_gcoord,ivec3d cpu_cell_dim,
		vector<long> atomnumber,vector<int> atomtypes,vector<double> atommass,vector<double> positions,vector<double> epot,int prank,int debug,vec3d cell_dim){


     // method for creating particle objects and fill in cell container

	 //*************************************************
	 // NOTE:before make_particles():--  cell_block()--> make_cells()--> with cell id and glob coord
	 //*************************************************

	 long m=0; long part_counter;
	 ivec3d cell_glob_coord,nb_cpu_gcoord;
	 int loc_rank;
     double inst_pos[3];
     cellblock ret_obj;

     int n_counter = 0;
     int l_counter = 0;

     int cell_index = 0;

     long total_cells = loc_obj.get_cell_list_size();

	 for(long i=0;i<tatoms_cpu;i++){

		 particle atom;

		 // assigning particle attributes
		 atom.set_mynumber(atomnumber.at(i));
 		 atom.set_mytype(atomtypes.at(i));
         atom.set_mymass(atommass.at(i));
         // better split as position x,y & z
         atom.set_myposition(positions.at((i*3)+0),positions.at((i*3)+1),positions.at((i*3)+2));
 		 atom.set_myepot(epot.at(i));

 		 // get global cell coordinate from particle position
 		 inst_pos[0] = positions.at((i*3)+0); inst_pos[1] = positions.at((i*3)+1); inst_pos[2] = positions.at((i*3)+2);

 		 // counting particles

 		 // include later (March) (with in scope of Montecarlo_new_api.cpp)

// 		 if((mc_atomtypes.at(i)%mc_real_types) == 0) mc_real_cpu++;   // add real Fe particle
// 		 if((mc_atomtypes.at(i)%mc_real_types) == 1) mc_carb_cpu++;   // add real C particle
// 		 if((mc_atomtypes.at(i)%mc_real_types) == 2) mc_phold_cpu++;  // add placeholder particle


         /*
 		 if (i==0){
 		      cout << "===================================================================================" << endl;
 	          cout << "         Inside make_particles  check while reading and assigning              " << endl;
 		      cout << " particle 0 -- id   :" << atom.get_mynumber() << endl;
 	          cout << " particle 0 -- type :" << atom.get_mytype() << endl;
 	          cout << " particle 0 -- mass :" << atom.get_mymass() << endl;
 	          cout << " particle 0 -- pos_x :" << atom.get_myposition().x << endl;
 	          cout << " particle 0 -- pos_y :" << atom.get_myposition().y << endl;
 	          cout << " particle 0 -- pos_z :" << atom.get_myposition().z << endl;
 	          //cout << " particle 0 -- epot  :" << atom.get_myepot() << endl;
 	          cout << "===================================================================================" << endl;
 		 }
 		 */


 		 cell_glob_coord = get_particle_glob_coordinate(inst_pos,cell_dim); // New try with my method

 		 // decrement operations
 		 cell_glob_coord.x = cell_glob_coord.x-1; cell_glob_coord.y = cell_glob_coord.y-1; cell_glob_coord.z = cell_glob_coord.z-1;

 		 //cell_glob_coord = cell_coordinate(inst_pos,tbox_dim,global_cell_dim); // Earlier implementation

 		 ivec3d cpu_fact = loc_cpu_gcoord;

 		 // get cell local coordinate
 		 ivec3d cell_loc_coord = get_cell_loc_coord(cell_glob_coord,cpu_fact,cpu_cell_dim);

 		 // to get memory index of cell in cell list
 		 //cell_index = cell_loc_coord.x + (cell_loc_coord.y * cpu_cell_dim.x) + (cell_loc_coord.z * cpu_cell_dim.y * cpu_cell_dim.x);


 		 cell_index = loc_obj.cell_with_lcoord(cell_loc_coord,total_cells).get_cell_id(); // expensive!!!!

      	 // just for testing add all particles to cell 1

 		 celltype cell0 = loc_obj.get_cell(cell_index);

 		 cell0.add_particle(atom);

         loc_obj.set_cell(cell0,cell_index);  // update cell object with included particle

//         if (debug == 1){
//
//        	   if((prank == 6) &&  (cell_index==63) && (l_counter == 0)){
//
//        		         cout << "  Super check: Inside make particles " << endl;
//                         cout << "  cell index  : " << cell_index << endl;
//                         cout << "  cell size :  " <<  loc_obj.get_cell(cell_index).get_nparticles() << endl;
//                         cout << "  my no :" << loc_obj.get_cell(cell_index).get_particle(0).get_mynumber() << endl;
//    	                 cout << "  my type :" << loc_obj.get_cell(cell_index).get_particle(0).get_mytype() << endl;
//        	             cout << "  my mass :" << loc_obj.get_cell(cell_index).get_particle(0).get_mymass() << endl;
//        	             cout << "  my pos:x " << loc_obj.get_cell(cell_index).get_particle(0).get_myposition().x << endl;
//        	             cout << "  my pos:y " << loc_obj.get_cell(cell_index).get_particle(0).get_myposition().y << endl;
//        	             cout << "  my pos:z " << loc_obj.get_cell(cell_index).get_particle(0).get_myposition().z << endl;
//    	                 cout << "  epot      "<< loc_obj.get_cell(cell_index).get_particle(0).get_myepot() << endl;
//
//    	                 l_counter++;
//           }
//         }

 		 }// end of for loop

	     // returning updated cellblock
	     ret_obj = loc_obj;

	     return ret_obj;

} //make particles

// Tested
cellblock create_maxwell_velocities(cellblock loc_obj,double temp, ivec3d*  restriction, int prank,int debug){

	 //     this code fragment/logic is inspired and inherited from IMD- imd_maxwell.c
	 //     the impulse computation is re-implemented in terms of velocity .(velocity = impulse/mass)

	  cellblock ret_obj;

	  //   create and fill initial velocities for particles

      double vx,vy,vz;                      // velocity holders
      double imp_x,imp_y,imp_z;             // impulse holders
      double sum_x, sum_y, sum_z;           // sum holders
      double rtemp;                         // reduced temp variable

      int   dof_x,dof_y,dof_z;               // degrees of freedom holders
      long  tot_dof_x,tot_dof_y,tot_dof_z ;

      long cell_count = loc_obj.get_cell_list_size(); // no of cells in cellblock
      long particles_count = 0;
      int mytype=0;
      int cell_index=0, l_counter=0;

      for(long i=0; i<cell_count; i++){ // loop over all cells

    	  celltype cobj;
    	  cobj = loc_obj.get_cell(i);

    	  particles_count = cobj.get_nparticles();

    	  //particles_count = loc_obj.get_cell(i).get_nparticles();

    	  for(long j=0; j<particles_count; j++){ // loop over all particles

        	  particle atom;

        	  atom = cobj.get_particle(j);
        	  mytype = atom.get_mytype();

        	  dof_x  = restriction[mytype].x;
        	  dof_y  = restriction[mytype].y;
        	  dof_z  = restriction[mytype].z;

        	  // reduced temp
        	  rtemp = sqrt(atom.get_mymass() * temp);

        	  imp_x = get_gaussian(rtemp) * dof_x ;
        	  imp_y = get_gaussian(rtemp) * dof_y ;
        	  imp_z = get_gaussian(rtemp) * dof_z ;

        	  // initially storing impulse
        	  atom.set_myvelocity(imp_x,imp_y,imp_z);

        	  //update atom particle object
              cobj.set_particle(atom,j);

        	  // summing total dof
        	  tot_dof_x += dof_x;
        	  tot_dof_y += dof_y;
        	  tot_dof_z += dof_z;

        	  // summing total impulse
        	  sum_x += imp_x;
        	  sum_y += imp_y;
        	  sum_z += imp_z;

    	  } // loop particles -1

    	  //update cell object
    	  loc_obj.set_cell(cobj,i);

      } // loop cells-1

      // this snippet has to be verified
//	  sum_x = tot_dof_x == 0 ? 0.0 : sum_x / tot_dof_x;
//	  sum_y = tot_dof_y == 0 ? 0.0 : sum_y / tot_dof_y;
//	  sum_z = tot_dof_z == 0 ? 0.0 : sum_z / tot_dof_z;


	  sum_x =  sum_x / tot_dof_x;
	  sum_y =  sum_y / tot_dof_y;
	  sum_z =  sum_z / tot_dof_z;

      double new_vx,new_vy,new_vz;

	  for(long i=0; i<cell_count; i++){ // loop over all cells

    	  celltype cobj;
    	  cobj = loc_obj.get_cell(i);
    	  particles_count = cobj.get_nparticles();

    	  //particles_count = loc_obj.get_cell(cell_count).get_nparticles();


		  for(long j=0; j<particles_count; j++){ // loop over all particles

        	  particle atom;
        	  int mytype;
        	  atom = cobj.get_particle(j);
        	  mytype = atom.get_mytype();

        	  dof_x  = restriction[mytype].x;
        	  dof_y  = restriction[mytype].y;
        	  dof_z  = restriction[mytype].z;

        	  // new velocity computation

        	  new_vx = atom.get_myvelocity().x - ( (sum_x /atom.get_mymass()) * dof_x ) ;
        	  new_vy = atom.get_myvelocity().y - ( (sum_y /atom.get_mymass()) * dof_y ) ;
        	  new_vz = atom.get_myvelocity().z - ( (sum_z /atom.get_mymass()) * dof_z ) ;


//        	  new_vx = ((atom.get_myvelocity().x - sum_x )/atom.get_mymass()) * dof_x;
//        	  new_vy = ((atom.get_myvelocity().y - sum_y )/atom.get_mymass()) * dof_y;
//        	  new_vz = ((atom.get_myvelocity().z - sum_z )/atom.get_mymass()) * dof_z;

        	  //updating velocity
        	  atom.set_myvelocity(new_vx,new_vy,new_vz);

        	  //update atom particle
              cobj.set_particle(atom,j);

		  } // loop particles -2

    	  //update cell object
    	  loc_obj.set_cell(cobj,i);

          cell_index = i;

          if((prank == 0) && (l_counter == 0) && (debug==1)){

         	  cout << "  Super check: create Maxwell velocities " << endl;
              cout << "  cell index  : " << cell_index << endl;
              cout << "  cell size :  " <<  loc_obj.get_cell(cell_index).get_nparticles() << endl;
              cout << "  my no :" << loc_obj.get_cell(cell_index).get_particle(0).get_mynumber() << endl;
     	      cout << "  my type :" << loc_obj.get_cell(cell_index).get_particle(0).get_mytype() << endl;
         	  cout << "  my mass :" << loc_obj.get_cell(cell_index).get_particle(0).get_mymass() << endl;
         	  cout << "  my pos:x " << loc_obj.get_cell(cell_index).get_particle(0).get_myposition().x << endl;
         	  cout << "  my pos:y " << loc_obj.get_cell(cell_index).get_particle(0).get_myposition().y << endl;
         	  cout << "  my pos:z " << loc_obj.get_cell(cell_index).get_particle(0).get_myposition().z << endl;
              cout << "  my vel:x " << loc_obj.get_cell(cell_index).get_particle(0).get_myvelocity().x << endl;
              cout << "  my vel:y " << loc_obj.get_cell(cell_index).get_particle(0).get_myvelocity().y << endl;
              cout << "  my vel:z " << loc_obj.get_cell(cell_index).get_particle(0).get_myvelocity().z << endl;
         	  cout << "  epot      "<< loc_obj.get_cell(cell_index).get_particle(0).get_myepot() << endl;

     	      l_counter++;
          }

	  } // loop cells -2


	  //updating block object
	  ret_obj = loc_obj;

	  // check msg
	  cout<<" initial velocities are computed for cpu id:"<< prank <<endl;

	  return ret_obj;

}
// Tested
double get_gaussian(double sigma){

	      //     this code fragment/logic is inspired and inherited from IMD- imd_maxwell.c
	      //     Polar (Box-Mueller) method; See Knuth v2, 3rd ed, p122

		  double x, y, r2;

		  do{
		      /* choose x,y in uniform square (-1,-1) to (+1,+1) */
			  // NOTE: c++ random generator to be updated check later

		      x = -1 + 2 * drand48();
		      y = -1 + 2 * drand48();

		      /* see if it is in the unit circle */
		      r2 = x * x + y * y;

		  }while (r2 > 1.0 || r2 == 0);


		  /* Box- Muller transform */
		  return (double) (sigma * y * sqrt (-2.0 * log (r2) / r2));

}

// Tested
particle sample_zone(cellblock bobj,int win_id,ivec3d cpu_cell_dim,int prank){

	celltype cobj;             // cell object

	int sam_cell_counter = (cpu_cell_dim.x/2 * cpu_cell_dim.y/2 * cpu_cell_dim.z/2 );   // no of sample cells belonging to window
	celltype sample_cells[sam_cell_counter];  // sample cells for given window type
	long rand_no,n_particles;
	particle atom;

	// local cell coordinates zone
    ivec3d test;
    int count=0,rand_cell;
    long ncells = bobj.get_cell_list_size();


    // defining zone limits
    int x_start=0,  x_mid = cpu_cell_dim.x/2,    x_end = cpu_cell_dim.x;   // zone limit parameters in x direction
    int y_start=0,  y_mid = cpu_cell_dim.y/2,    y_end = cpu_cell_dim.y;   // zone limit parameters in y direction
    int z_start=0,  z_mid = cpu_cell_dim.z/2,    z_end = cpu_cell_dim.z;   // zone limit parameters in z direction


    // sample cell ranges
    //ivec6d zone_limit_0 = {0,4,0,4,0,1}; // window 0

    ivec6d zone_limit_0 = {x_start,x_mid-1,y_start,y_mid-1,z_start,z_mid-1}; // window 0
    ivec6d zone_limit_1 = {x_start,x_mid-1,y_start,y_mid-1,z_mid,z_end-1};   // window 1
    ivec6d zone_limit_2 = {x_mid,x_end-1,y_start,y_mid-1,z_start,z_mid-1};   // window 2
    ivec6d zone_limit_3 = {x_mid,x_end-1,y_start,y_mid-1,z_mid,z_end-1};     // window 3

    ivec6d zone_limit_4 = {x_start,x_mid-1,y_mid,y_end-1,z_start,z_mid-1}; // window 4
    ivec6d zone_limit_5 = {x_start,x_mid-1,y_mid,y_end-1,z_mid,z_end-1};   // window 5
    ivec6d zone_limit_6 = {x_mid,x_end-1,y_mid,y_end-1,z_start,z_mid-1};   // window 6
    ivec6d zone_limit_7 = {x_mid,x_end-1,y_mid,y_end-1,z_mid,z_end-1};     // window 7


    ivec6d zone_limit[8]={zone_limit_0,zone_limit_1,zone_limit_2,zone_limit_3,zone_limit_4,zone_limit_5,zone_limit_6,zone_limit_7};

    // zone_limit: define sample window positions (x,y,z)

    // prepare cell list as per sample window id

    for(int i=zone_limit[win_id].xmin;i<=zone_limit[win_id].xmax;i++){
    	for(int j=zone_limit[win_id].ymin;j<=zone_limit[win_id].ymax;j++){
    		for(int k=zone_limit[win_id].zmin;k<=zone_limit[win_id].zmax;k++){

    			test.x = i; test.y = j; test.z = k;
    			sample_cells[count] = bobj.cell_with_lcoord(test,ncells);

                count++;
    		}
    	}
    }


//    if(prank == 6){
//    	cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
//    	cout << "---                Inside sample zone method                      ---------" << endl;
//    	cout << " My process CPU  : " << prank << endl;
//    	cout << " Chosen window position : " << win_id << endl;
//    	cout << " computed no of sample cell -- Sample_cell_counter : " << sam_cell_counter << endl;
//    	cout << " Included no of cells from loop counter :   " << count << endl;
//    	// print cell local coordinate and check??
//    	ivec3d loc_coord;
//    	ivec3d glo_coord;
//    	vec3d pos_ph;
//    	cout << "=======================================================" << endl;
//    	cout << "                       sample cells check     "  << endl;
//    	cout << "=======================================================" << endl;
//    	for( int i=0; i<sam_cell_counter;i++){
//    		loc_coord = sample_cells[i].get_cell_loc_coord();
//    		glo_coord = sample_cells[i].get_cell_glob_coord();
//    		cout << " sample get cell id : " << sample_cells[i].get_cell_id()  << endl;
//    		cout << " sample cell local coordinate :  [ " <<loc_coord.x<<" "<<loc_coord.y<<" "<<loc_coord.z<<"  ]"<<endl;
//    		cout << " sample cell global coordinate :  [ " <<glo_coord.x<<" "<<glo_coord.y<<" "<<glo_coord.z<<"  ]"<<endl;
//    		cout << " no of particles : " << sample_cells[i].get_nparticles() << endl;
//
//
//    		cout << " Inside sample get cell id : " << sample_cells[i].get_cell_id()  << endl;
//
//    		for(long j=0;j<sample_cells[i].get_nparticles();j++){
//				pos_ph = sample_cells[i].get_particle(j).get_myposition();
//				cout << " particle id : " << j << " "<< " position : "<< pos_ph.x <<" "<<pos_ph.y<<" "<<pos_ph.z<<" "<<endl;
//			}
//
//    	}
//    	cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
//
//    	cout << " CPU Cell local coordinates assignments check  " << endl;
//    	cout << " =====================================================" << endl;
//
//
//    }


    // choose one cell in random (belong to sample window)

    srand(prank);
	rand_cell = rand()%sam_cell_counter;
	cobj = sample_cells[rand_cell];



	do{
		//no of particles
		n_particles= cobj.get_nparticles();
		// choosing random placeholder
		rand_no=(long) rand()%n_particles;
	}while(cobj.get_particle(rand_no).get_mytype() !=2);


	atom = cobj.get_particle(rand_no);

//	cout << "************************************************************" << endl;
//	cout << " Random generator check    "  << endl;
//	cout << " my process rank : " << prank << endl;
//	cout << " Randomly chosen cell no:  "  << rand_cell<< endl;
//	cout << " Random cell local coordinate :  " <<cobj.get_cell_loc_coord().x<<" "<<cobj.get_cell_loc_coord().y<<" "<<cobj.get_cell_loc_coord().z<<endl;
//	cout << "  Random cell global coordinate :" <<cobj.get_cell_glob_coord().x<<" "<<cobj.get_cell_glob_coord().y<<" "<<cobj.get_cell_glob_coord().z<<endl;
//	cout << " Randomly chosen particle index : "<< rand_no << endl;
//    cout << " Randomly chosen particle id:"  << atom.get_mynumber() << endl;
//	cout << " Random particle position x :"  << atom.get_myposition().x << endl;
//	cout << " Random particle position y :"  << atom.get_myposition().y << endl;
//	cout << " Random particle position z :"  << atom.get_myposition().z << endl;
//
//    cout << "************************************************************" << endl;

	return atom;
}

celltype construct_sphere(particle pobj, cellblock bobj, int win_id,const char* filename,int prank,MPI_Comm comm_name, MPI_Status stat,
		int test_rank,ivec3d cpu_dim,double* data_list){

	//########################################################################################################################################################
	//
    //                                      PHASE 1 - Random particle selection and communication with neighbors
	//
	//                   Random particle is chosen from sample cell list defined by appropriate window position
	//########################################################################################################################################################

	celltype sphere_old;

	// possible neighbor index as per window location

	int window_x[8] = {-1,-1,+1,+1,-1,-1,+1,+1};
	int window_y[8] = {-1,-1,-1,-1,+1,+1,+1,+1};
	int window_z[8] = {-1,+1,-1,+1,-1,+1,-1,+1};

	// TOO BE VERIFIED

	// 8 windows position in CPU (fixed positions -- window 0 to window 7)

	// Old assignment
	ivec3d win_pos_0 = {0,0,0}; ivec3d win_pos_1 = {0,0,1}; ivec3d win_pos_2 = {1,0,0}; ivec3d win_pos_3 = {1,0,1};
	ivec3d win_pos_4 = {0,1,0}; ivec3d win_pos_5 = {0,1,1}; ivec3d win_pos_6 = {1,1,0}; ivec3d win_pos_7 = {1,1,1};

	// New assignment
//	ivec3d win_pos_0 = {0,0,0}; ivec3d win_pos_1 = {1,0,0}; ivec3d win_pos_2 = {0,0,1}; ivec3d win_pos_3 = {1,0,1};
//	ivec3d win_pos_4 = {0,1,0}; ivec3d win_pos_5 = {1,1,0}; ivec3d win_pos_6 = {0,1,1}; ivec3d win_pos_7 = {1,1,1};

	ivec3d window_position[8] = {win_pos_0,win_pos_1,win_pos_2,win_pos_3,win_pos_4,win_pos_5,win_pos_6,win_pos_7};

	// neighbor CPU global coordinate
	int nb_gcoord [7][3]; // [neighbor] [gcoord_x gcoord_y gcoord_z]

	ivec3d neig_coord;   // variable to fetch and hold neighbor coordinate


	// ==========================================================
	// Detect neighbor processes (maximum 7 if periodicity is enabled in all direction) as per window position
	// ==========================================================
    // select neighbors as per sample window position
	int xfact = window_x[win_id]; int yfact = window_y[win_id]; int zfact = window_z[win_id];

	// send neighbors (to whom I send)        ---- process communications (in z direction)
	int x_send_phase[4] = { 0,xfact,xfact, 0}; int y_send_phase[4] = { 0, 0,yfact,yfact}; int z_send_phase[4] = {zfact,zfact,zfact,zfact};

	// Receive neighbors (to whom I receive) ----- process communications (in z direction)
	int x_recv_phase[4] = { 0,-xfact,-xfact, 0}; int y_recv_phase[4] = { 0, 0,-yfact,-yfact}; int z_recv_phase[4] = {-zfact,-zfact,-zfact,-zfact};

    // send neighbors (to whom I send) -----  process communications (in x/y direction)
	int x_send_next[3] = {xfact,xfact, 0 }; int y_send_next[3] = { 0,yfact,yfact }; int z_send_next[3] = { 0, 0, 0 };

	// Receive neighbors (to whom I send) -----  process communications (in x/y direction)
	int x_recv_next[3] = {-xfact,-xfact, 0 }; int y_recv_next[3] = { 0,-yfact,-yfact }; int z_recv_next[3] = { 0, 0, 0 };

	// variables to store the position of selected particle on my CPU and neighbor CPU's
	// randomly selected particle positions (my CPU and rec_neighbor CPU)

	double my_pos[3], rec_pos_0[3],rec_pos_1[3],rec_pos_2[3],rec_pos_3[3],rec_pos_4[3],rec_pos_5[3],rec_pos_6[3];

	double* rec_pos[7]={rec_pos_0,rec_pos_1,rec_pos_2,rec_pos_3,rec_pos_4,rec_pos_5,rec_pos_6}; // get and store positions from all seven neighbors

	//  get chosen random particle position
	my_pos[0]=pobj.get_myposition().x; my_pos[1]=pobj.get_myposition().y; my_pos[2]=pobj.get_myposition().z;


	//***************************************************************************
	// communication part -- send/request neighbors with chosen particle
	//***************************************************************************

	// get my process global coordinate
	ivec3d gcoord = get_cpu_gcoord(prank,comm_name);
	int x = gcoord.x; int y = gcoord.y; int z = gcoord.z; // [ x y z] -- this process coordinate from MPI Cartesian system

	int send_id,rec_id;
	double* posit;

	// synchronize communication
	MPI_Barrier(comm_name);

	if(prank == test_rank){
		cout << "  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$   " << endl;
		cout << "   Selected Window position : " << win_id << endl;
		cout << "   MPI Communication in Z direction     " << endl;
		cout << "      "<< endl;
	}

	// z direction communication
	for (int ind=0;ind<4;ind++){

	     if (z%2==0){ // if z is even

	       send_id = get_cpu_rank(x+x_send_phase[ind],y+y_send_phase[ind],z+z_send_phase[ind],comm_name);

	       MPI_Send(my_pos,3,MPI_DOUBLE,send_id,0,comm_name);

	       if(prank == test_rank){
	    	   cout << " I sent :  [ " << my_pos[0] << " " << my_pos[1] << ""<<" " <<my_pos[2] << endl;
	    	   cout << " Process :  " << prank <<" sent data to front target process : " << send_id << endl;
	    	   cout << "      "<< endl;
	       }

	       rec_id = get_cpu_rank(x+x_recv_phase[ind],y+y_recv_phase[ind],z+z_recv_phase[ind],comm_name);

	       MPI_Recv(rec_pos[ind],3,MPI_DOUBLE,rec_id,1,comm_name,&stat);

	       neig_coord = get_cpu_gcoord(rec_id,comm_name);
	       nb_gcoord[ind][0] = neig_coord.x;  nb_gcoord[ind][1] = neig_coord.y; nb_gcoord[ind][2] = neig_coord.z;

	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" received data from back target process : " << rec_id << endl;
	    	   posit = rec_pos[ind];
	    	   cout << " I received :  [ " << posit[0] << " " << posit[1] << ""<<" " <<posit[2] << endl;
	    	   cout << "      "<< endl;
	       }

	     }
	     else{       // if z is odd

	       rec_id = get_cpu_rank(x+x_recv_phase[ind],y+y_recv_phase[ind],z+z_recv_phase[ind],comm_name);

	       MPI_Recv(rec_pos[ind],3,MPI_DOUBLE,rec_id,0,comm_name,&stat);

	       neig_coord = get_cpu_gcoord(rec_id,comm_name);
	       nb_gcoord[ind][0] = neig_coord.x;  nb_gcoord[ind][1] = neig_coord.y; nb_gcoord[ind][2] = neig_coord.z;


	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" received data from back target process : " << rec_id << endl;
	    	   posit = rec_pos[ind];
	    	   cout << " I received :  [ " << posit[0] << " " << posit[1] << ""<<" " <<posit[2] << endl;
	    	   cout << "      "<< endl;
	       }

	       send_id = get_cpu_rank(x+x_send_phase[ind],y+y_send_phase[ind],z+z_send_phase[ind],comm_name);

	       MPI_Send(my_pos,3,MPI_DOUBLE,send_id,1,comm_name);

	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" sent data to front target process : " << send_id << endl;
	    	   cout << "      "<< endl;
	       }

	     }
	}

	// Next phase even-even or odd-odd communications

	if(prank == test_rank){
		cout << "  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$   " << endl;
		cout << "   MPI Communication in X direction     " << endl;
		cout << "      "<< endl;

	}

	// x direction communication
	if (x%2 == 0){    // if x is even
	         // right & left comm

		   send_id = get_cpu_rank(x+x_send_next[0],y+y_send_next[0],z+z_send_next[0],comm_name);

		   MPI_Send(my_pos,3,MPI_DOUBLE,send_id,2,comm_name);

	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" sent data to right target process : " << send_id << endl;
	    	   cout << "      "<< endl;
	       }

	       rec_id = get_cpu_rank(x+x_recv_next[0],y+y_recv_next[0],z+z_recv_next[0],comm_name);

	       MPI_Recv(rec_pos[4],3,MPI_DOUBLE,rec_id,3,comm_name,&stat);

	       neig_coord = get_cpu_gcoord(rec_id,comm_name);
	       nb_gcoord[4][0] = neig_coord.x;  nb_gcoord[4][1] = neig_coord.y; nb_gcoord[4][2] = neig_coord.z;

	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" received data from left target process : " << rec_id << endl;
	    	   posit = rec_pos[4];
	    	   cout << " I received :  [ " << posit[0] << " " << posit[1] << ""<<" " <<posit[2] << endl;
	    	   cout << "      "<< endl;
	       }

	         // top-right & bottom-left comm
	       send_id = get_cpu_rank(x+x_send_next[1],y+y_send_next[1],z+z_send_next[1],comm_name);

	       MPI_Send(my_pos,3,MPI_DOUBLE,send_id,4,comm_name);

	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" sent data to top - right target process : " << send_id << endl;
	    	   cout << "      "<< endl;
	       }

	       rec_id = get_cpu_rank(x+x_recv_next[1],y+y_recv_next[1],z+z_recv_next[1],comm_name);

	       MPI_Recv(rec_pos[5],3,MPI_DOUBLE,rec_id,5,comm_name,&stat);

	       neig_coord = get_cpu_gcoord(rec_id,comm_name);
	       nb_gcoord[5][0] = neig_coord.x;  nb_gcoord[5][1] = neig_coord.y; nb_gcoord[5][2] = neig_coord.z;


	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" received data from bottom - left target process : " << rec_id << endl;
	    	   posit = rec_pos[5];
	    	   cout << " I received :  [ " << posit[0] << " " << posit[1] << ""<<" " <<posit[2] << endl;
	    	   cout << "      "<< endl;
	       }

	}
	else{            // if x is odd

		   rec_id = get_cpu_rank(x+x_recv_next[0],y+y_recv_next[0],z+z_recv_next[0],comm_name);
		   // right & left comm
		   MPI_Recv(rec_pos[4],3,MPI_DOUBLE,rec_id,2,comm_name,&stat);

	       neig_coord = get_cpu_gcoord(rec_id,comm_name);
	       nb_gcoord[4][0] = neig_coord.x;  nb_gcoord[4][1] = neig_coord.y; nb_gcoord[4][2] = neig_coord.z;

	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" received data from left target process : " << rec_id << endl;
	    	   posit = rec_pos[4];
	    	   cout << " I received :  [ " << posit[0] << " " << posit[1] << ""<<" " <<posit[2] << endl;
	    	   cout << "      "<< endl;
	       }

	       send_id = get_cpu_rank(x+x_send_next[0],y+y_send_next[0],z+z_send_next[0],comm_name);

	       MPI_Send(my_pos,3,MPI_DOUBLE,send_id,3,comm_name);

	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" sent data to right target process : " << send_id << endl;
	    	   cout << "      "<< endl;
	       }

	       rec_id = get_cpu_rank(x+x_recv_next[1],y+y_recv_next[1],z+z_recv_next[1],comm_name);

	       // top-right & bottom-left comm
	       MPI_Recv(rec_pos[5],3,MPI_DOUBLE,rec_id,4,comm_name,&stat);

	       neig_coord = get_cpu_gcoord(rec_id,comm_name);
	       nb_gcoord[5][0] = neig_coord.x;  nb_gcoord[5][1] = neig_coord.y; nb_gcoord[5][2] = neig_coord.z;

	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" received data from bottom - left target process : " << rec_id << endl;
	    	   posit = rec_pos[5];
	    	   cout << " I received :  [ " << posit[0] << " " << posit[1] << ""<<" " <<posit[2] << endl;
	    	   cout << "      "<< endl;
	       }

	       send_id = get_cpu_rank(x+x_send_next[1],y+y_send_next[1],z+z_send_next[1],comm_name);

	       MPI_Send(my_pos,3,MPI_DOUBLE,send_id,5,comm_name);

	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" sent data to top - right target process : " << send_id << endl;
	    	   cout << "      "<< endl;
	       }

	}

	if(prank == test_rank){
		cout << "  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$   " << endl;
		cout << "   MPI Communication in Y direction     " << endl;
		cout << "      "<< endl;

	}

	// y direction communication

	if (y%2 == 0){      // if y is even

		   send_id = get_cpu_rank(x+x_send_next[2],y+y_send_next[2],z+z_send_next[2],comm_name);
		   // top & bottom comm

		   MPI_Send(my_pos,3,MPI_DOUBLE,send_id,4,comm_name);

	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" sent data to bottom target process : " << send_id << endl;
	    	   cout << "      "<< endl;
	       }

	       rec_id = get_cpu_rank(x+x_recv_next[2],y+y_recv_next[2],z+z_recv_next[2],comm_name);

	       MPI_Recv(rec_pos[6],3,MPI_DOUBLE,rec_id,5,comm_name,&stat);

	       neig_coord = get_cpu_gcoord(rec_id,comm_name);
	       nb_gcoord[6][0] = neig_coord.x;  nb_gcoord[6][1] = neig_coord.y; nb_gcoord[6][2] = neig_coord.z;

	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" received data from top target process : " << rec_id << endl;
	    	   posit = rec_pos[6];
	    	   cout << " I received :  [ " << posit[0] << " " << posit[1] << ""<<" " <<posit[2] << endl;
	    	   cout << "      "<< endl;
	       }
	}
	else{              // if y is odd

		   rec_id = get_cpu_rank(x+x_recv_next[2],y+y_recv_next[2],z+z_recv_next[2],comm_name);
		   MPI_Recv(rec_pos[6],3,MPI_DOUBLE,rec_id,4,comm_name,&stat);

	       neig_coord = get_cpu_gcoord(rec_id,comm_name);
	       nb_gcoord[6][0] = neig_coord.x;  nb_gcoord[6][1] = neig_coord.y; nb_gcoord[6][2] = neig_coord.z;


	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" received data from top target process : " << rec_id << endl;
	    	   posit = rec_pos[6];
	    	   cout << " I received :  [ " << posit[0] << " " << posit[1] << ""<<" " <<posit[2] << endl;
	    	   cout << "      "<< endl;
	       }

	       send_id = get_cpu_rank(x+x_send_next[2],y+y_send_next[2],z+z_send_next[2],comm_name);
	       MPI_Send(my_pos,3,MPI_DOUBLE,send_id,5,comm_name);

	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" sent data to bottom target process : " << send_id << endl;
	    	   cout << "      "<< endl;
	       }
	}

	if(prank == test_rank){
		cout << "  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$   " << endl;
		// neighbor CPU global coordinate check
		cout <<"***********************************************************"<<endl;
		for(int i=0; i<7; i++){
			cout << " nb_gcoord[ "<<i<<" ]  : "<< nb_gcoord[i][0] <<" "<< nb_gcoord[i][1] <<" "<< nb_gcoord[i][2] <<" "<<endl;
		}
		cout << " data_list[0]  : " << data_list[0] << endl;
		cout << " no of cells in CPU  :   " << bobj.get_cell_list_size() << endl;
		cout << " no of particles in one cell 0  :" << bobj.get_cell(0).get_nparticles()<<endl;


		cout <<"***********************************************************"<<endl;
	}


	// synchronize communication
	MPI_Barrier(comm_name);

//########################################################################################################################################################
//
//                                      PHASE 2 - Sweep test and buffer construction
//
//
//########################################################################################################################################################

		// flags for coordinate shifting (if they are on edges!!)
	    int flag [7][3]; // [neighbor] [x y z]

	    // global cpu array dimension range
	  	ivec3d cpu_max = {cpu_dim.x-1,cpu_dim.y-1,cpu_dim.z-1};            // dimension - Max
	    ivec3d cpu_min = {0,0,0};                                          // dimension - Min

		// buffer allocation for receive neighbors (NEED REDEFINITION)

		int buf_count = int(data_list[0]) * 10; // buffer_size


		// N1 N2 N4 N6 N3 N5 N7
		double *nb_0, *nb_1, *nb_2, *nb_3, *nb_4, *nb_5, *nb_6; // to be allocated and send

		// allocate memory (to send to neighbors)
		nb_0 = new double [buf_count]; nb_1= new double [buf_count]; nb_2 = new double [buf_count];
		nb_3 = new double [buf_count]; nb_4= new double [buf_count]; nb_5 = new double [buf_count];
		nb_6 = new double [buf_count];

		double* nb_buffer[7] = {nb_0,nb_1,nb_2,nb_3,nb_4,nb_5,nb_6};

		double *my_0, *my_1, *my_2, *my_3, *my_4, *my_5, *my_6; // to be received and filled in sphere cells

		// allocate memory (to receive from neighbors)
		my_0 = new double [buf_count]; my_1= new double [buf_count]; my_2 = new double [buf_count];
		my_3 = new double [buf_count]; my_4= new double [buf_count]; my_5 = new double [buf_count];
		my_6 = new double [buf_count];

		double* my_buffer[7] = {my_0,my_1,my_2,my_3,my_4,my_5,my_6};


	    // boundary checks and flag assignments (ENSURE COORDINATES SHIFT IF REQUIRED)

	    // z check  with  N1 ( Neighbor -1)
	    if( (nb_gcoord[0][2] == cpu_min.z) && (window_position[win_id].z==0) ) { flag[0][2] = -1;}
	    if( (nb_gcoord[0][2] == cpu_max.z) && (window_position[win_id].z==1) ) { flag[0][2] = +1;}

	    // x-z check with N2 ( Neighbor -2)
	    if( (nb_gcoord[1][0] == cpu_min.x) && (window_position[win_id].x==0) ) { flag[1][0] = -1;}
	    if( (nb_gcoord[1][0] == cpu_max.x) && (window_position[win_id].x==1) ) { flag[1][0] = +1;}
	    if( (nb_gcoord[1][2] == cpu_min.z) && (window_position[win_id].z==0) ) { flag[1][2] = -1;}
	    if( (nb_gcoord[1][2] == cpu_max.z) && (window_position[win_id].z==1) ) { flag[1][2] = +1;}

	    // x-y-z check with N4 ( Neighbor -3)
	    if( (nb_gcoord[2][0] == cpu_min.x) && (window_position[win_id].x==0) ) { flag[2][0] = -1;}
	    if( (nb_gcoord[2][0] == cpu_max.x) && (window_position[win_id].x==1) ) { flag[2][0] = +1;}
	    if( (nb_gcoord[2][1] == cpu_min.y) && (window_position[win_id].y==0) ) { flag[2][1] = -1;}
	    if( (nb_gcoord[2][1] == cpu_max.y) && (window_position[win_id].y==1) ) { flag[2][1] = +1;}
	    if( (nb_gcoord[2][2] == cpu_min.z) && (window_position[win_id].z==0) ) { flag[2][2] = -1;}
	    if( (nb_gcoord[2][2] == cpu_max.z) && (window_position[win_id].z==1) ) { flag[2][2] = +1;}

	    // y-z check with N6 ( Neighbor -4)
	    if( (nb_gcoord[3][1] == cpu_min.y) && (window_position[win_id].y==0) ) { flag[3][1] = -1;}
	    if( (nb_gcoord[3][1] == cpu_max.y) && (window_position[win_id].y==1) ) { flag[3][1] = +1;}
	    if( (nb_gcoord[3][2] == cpu_min.z) && (window_position[win_id].z==0) ) { flag[3][2] = -1;}
	    if( (nb_gcoord[3][2] == cpu_max.z) && (window_position[win_id].z==1) ) { flag[3][2] = +1;}

	    // x check with N3  ( Neighbor -5)
	    if( (nb_gcoord[4][0] == cpu_min.x) && (window_position[win_id].x==0) ) { flag[4][0] = -1;}
	    if( (nb_gcoord[4][0] == cpu_max.x) && (window_position[win_id].x==1) ) { flag[4][0] = +1;}

	    // x-y check with N5 ( Neighbor -6)
	    if( (nb_gcoord[5][0] == cpu_min.x) && (window_position[win_id].x==0) ) { flag[5][0] = -1;}
	    if( (nb_gcoord[5][0] == cpu_max.x) && (window_position[win_id].x==1) ) { flag[5][0] = +1;}
	    if( (nb_gcoord[5][1] == cpu_min.y) && (window_position[win_id].y==0) ) { flag[5][1] = -1;}
	    if( (nb_gcoord[5][1] == cpu_max.y) && (window_position[win_id].y==1) ) { flag[5][1] = +1;}

	    // y check with N7 ( Neighbor -7)
	    if( (nb_gcoord[6][1] == cpu_min.y) && (window_position[win_id].y==0) ) { flag[6][1] = -1;}
	    if( (nb_gcoord[6][1] == cpu_max.y) && (window_position[win_id].y==1) ) { flag[6][1] = +1;}

		// local cell coordinates zone
	//    ivec3d test, nb_test[19];
	//
	//    // window sampling cells (my own cells)
	//    // later generalize for more than 8 cells
	//
	//    celltype sample_cells[8];
	//
	//    // ***************************************************************************
	//    //window neighbor cells (to be sent)
	//    // no of each neighbor cells ( N1 N2 N4 N6 N3 N5 N7 )
	//    // to be generalized later
	//    int neig_cells_size [7]={4,2,1,2,4,2,4};
	//
	//    // (Check this for generality)
	//    //****************************************************************************
	//
	//    // list of neighbor receive cells
	//    celltype neighb[19];
	//
	//    // total cells in cpu block
	//    long ncells= bobj.get_ncells();
	//
	//    // variables for neighbor cell ranges
	//    int nb_xmin, nb_ymin, nb_zmin, nb_xmax, nb_ymax, nb_zmax, nb_x, nb_y, nb_z;
	//
	//    // assign neighbor cell ranges
	//    nb_xmin = zone_limit[win_id].xmin; nb_xmax = zone_limit[win_id].xmax;
	//    nb_ymin = zone_limit[win_id].ymin; nb_ymax = zone_limit[win_id].ymax;
	//    nb_zmin = zone_limit[win_id].zmin; nb_zmax = zone_limit[win_id].zmax;
	//
	//    // Neighbor cells to be communicated as per sample window position
	//    if(nb_xmax ==1) {nb_x = 2;}
	//    else{ nb_x=0; }
	//
	//    if(nb_ymax ==1) {nb_y = 2;}
	//    else{ nb_y=0; }
	//
	//    if(nb_zmax ==1) {nb_z = 2;}
	//    else{ nb_z=0; }
	//
	//    // prepare NEIGHBOR cell list as per sample window id
	//    //Neighbor-1 (N1)
	//    nb_test[0].x =nb_xmin; nb_test[0].y = nb_ymin; nb_test[0].z = nb_z;    neighb[0] = bobj.cell_with_lcoord(nb_test[0],ncells);
	//    nb_test[1].x =nb_xmax; nb_test[1].y = nb_ymin; nb_test[1].z = nb_z;    neighb[1] = bobj.cell_with_lcoord(nb_test[1],ncells);
	//    nb_test[2].x =nb_xmin; nb_test[2].y = nb_ymax; nb_test[2].z = nb_z;    neighb[2] = bobj.cell_with_lcoord(nb_test[2],ncells);
	//    nb_test[3].x =nb_xmax; nb_test[3].y = nb_ymax; nb_test[3].z = nb_z;    neighb[3] = bobj.cell_with_lcoord(nb_test[3],ncells);
	//
	//
	////      conceptual idea
	////	    //updating cell
	////	    bobj.delete_cell(bobj.cell_with_lcoord(nb_test[0],ncells).get_cell_id());
	////
	////	    bobj.add_cell(neighb[0]);
	//
	//    //Neighbor-2 (N2)
	//    nb_test[4].x = nb_x; nb_test[4].y = nb_ymin; nb_test[4].z = nb_z;      neighb[4] = bobj.cell_with_lcoord(nb_test[4],ncells);
	//    nb_test[5].x = nb_x; nb_test[5].y = nb_ymax; nb_test[5].z = nb_z;      neighb[5] = bobj.cell_with_lcoord(nb_test[5],ncells);
	//
	//    //Neighbor-3 (N4)
	//    nb_test[6].x = nb_x; nb_test[6].y = nb_y;    nb_test[6].z = nb_z;      neighb[6] = bobj.cell_with_lcoord(nb_test[6],ncells);
	//
	//    //Neighbor-4 (N6)
	//    nb_test[7].x = nb_xmin; nb_test[7].y = nb_y; nb_test[7].z = nb_z;      neighb[7] = bobj.cell_with_lcoord(nb_test[7],ncells);
	//    nb_test[8].x = nb_xmax; nb_test[8].y = nb_y; nb_test[8].z = nb_z;      neighb[8] = bobj.cell_with_lcoord(nb_test[8],ncells);
	//
	//    //Neighbor-5 (N3)
	//    nb_test[9].x =nb_x;   nb_test[9].y = nb_ymin; nb_test[9].z = nb_zmin;  neighb[9] = bobj.cell_with_lcoord(nb_test[9],ncells);
	//    nb_test[10].x =nb_x; nb_test[10].y = nb_ymin; nb_test[10].z = nb_zmax; neighb[10]= bobj.cell_with_lcoord(nb_test[10],ncells);
	//    nb_test[11].x =nb_x; nb_test[11].y = nb_ymax; nb_test[11].z = nb_zmin; neighb[11]= bobj.cell_with_lcoord(nb_test[11],ncells);
	//    nb_test[12].x =nb_x; nb_test[12].y = nb_ymax; nb_test[12].z = nb_zmax; neighb[12]= bobj.cell_with_lcoord(nb_test[12],ncells);
	//
	//    //Neighbor-6 (N5)
	//    nb_test[13].x = nb_x; nb_test[13].y = nb_y; nb_test[13].z = nb_zmin;   neighb[13] = bobj.cell_with_lcoord(nb_test[13],ncells);
	//    nb_test[14].x = nb_x; nb_test[14].y = nb_y; nb_test[14].z = nb_zmax;   neighb[14] = bobj.cell_with_lcoord(nb_test[14],ncells);
	//
	//    //Neighbor-7 (N7)
	//    nb_test[15].x =nb_xmin; nb_test[15].y = nb_y; nb_test[15].z = nb_zmin; neighb[15] = bobj.cell_with_lcoord(nb_test[15],ncells);
	//    nb_test[16].x =nb_xmin; nb_test[16].y = nb_y; nb_test[16].z = nb_zmax; neighb[16] = bobj.cell_with_lcoord(nb_test[16],ncells);
	//    nb_test[17].x =nb_xmax; nb_test[17].y = nb_y; nb_test[17].z = nb_zmin; neighb[17] = bobj.cell_with_lcoord(nb_test[17],ncells);
	//    nb_test[18].x =nb_xmax; nb_test[18].y = nb_y; nb_test[18].z = nb_zmax; neighb[18] = bobj.cell_with_lcoord(nb_test[18],ncells);
	//
	//    //***************************************************************************
	//    //                  My received neighbor part - fill buffers for my neighbors
	//    //***************************************************************************
	//
	//    // sweep distance computation (with squares)
	//    double r_full   = pow((mc_rsweep + mc_sphere_wall),2.0);
	//    double r_sphere = pow(mc_rsweep,2.0);
	//
	//    vec3d ref_pos, curr_pos;                       // reference/current particle position
	//
	//    // temporary holder for particle attributes
	//	double temp_id; double temp_type;
	//	double temp_mass, temp_epot, dist_check;
	//	vec3d temp_pos, temp_vel;
	//
	//	long total_particles=0;                       // counter for total particle no in buffers
	//    int ind =0;                                   // neighbor cell counter
	//    long buf_counter,part_count;
	//
	//    for(int j=0; j<7;j++ ){ // loop over all neighbor cpu
	//    	    buf_counter=1;
	//    	    // get received neighbor chosen particle positions
	//    	    ref_pos.x = *rec_pos[j+0]; ref_pos.y = *rec_pos[j+1]; ref_pos.z = *rec_pos[j+2];
	//
	//    	    for (int k=0; k<neig_cells_size[j]; k++){ // loop over each neighbor sub total
	//    	         part_count = neighb[ind].get_nparticles(); //get total particles
	//
	//    	         for (long val=0;val<part_count;val++){ // loop over particles in each neighbor cell
	//	    		          curr_pos =neighb[ind].get_particle(val).get_myposition();
	//	    		          // coordinate shifting as per flag initialization(for cells on edges of CPU)
	//	    		          curr_pos.x = curr_pos.x + (flag[j][0] * mc_simbox_x.x);
	//	    		          curr_pos.y = curr_pos.y + (flag[j][1] * mc_simbox_y.y);
	//	    		          curr_pos.z = curr_pos.z + (flag[j][2] * mc_simbox_z.z);
	//
	//	    		          dist_check = distance_vect(ref_pos,curr_pos);
	//	    		          // ----- cutoff check
	//	    		          if(dist_check <= r_full ){
	// 		    		            // declare virtual particles on sphere boundary
	//     	    		            if( dist_check > r_sphere ){
	//     	    			             temp_type = (double) neighb[ind].get_particle(val).get_mytype();
	//
	//     	    			             // ignore placeholders on sphere wall
	//	    			                 if (temp_type != (double) 2 ) { // HC: now hardcoded for placeholders
	//
	//	    			        	           temp_id     = (double) (neighb[ind].get_particle(val).get_mynumber());
	//
	//	    			        	           // add particle id that are sent for sphere constructions
	//	    			        	           list_ids[j].push_back(temp_id);
	//
	//	    			        	           temp_type   = (double) (temp_type + mc_real_types);
	//                                           temp_mass   =  neighb[ind].get_particle(val).get_mymass();
	//
	//                                           // assign temporary position (coordinates shifted if neighbors on cpu boundary)
	//                                           temp_pos  =  curr_pos;
	//                                           temp_vel  =  neighb[ind].get_particle(val).get_myvelocity();
	//                                           temp_epot =  neighb[ind].get_particle(val).get_myepot();
	//
	//                                           // *************************************************
	//                                           // CRITICAL PART
	//                                           // *************************************************
	//
	//                                           // ADD TO BUFFER nb_buffer
	//                                           *(nb_buffer[j]+buf_counter++)  = temp_id;
	//                                           *(nb_buffer[j]+buf_counter++)  = temp_type;
	//                                           *(nb_buffer[j]+buf_counter++)  = temp_mass;
	//                                           *(nb_buffer[j]+buf_counter++)  = temp_pos.x;
	//                                           *(nb_buffer[j]+buf_counter++)  = temp_pos.y;
	//                                           *(nb_buffer[j]+buf_counter++)  = temp_pos.z;
	//                                           *(nb_buffer[j]+buf_counter++)  = temp_vel.x;
	//                                           *(nb_buffer[j]+buf_counter++)  = temp_vel.y;
	//                                           *(nb_buffer[j]+buf_counter++)  = temp_vel.z;
	//                                           *(nb_buffer[j]+buf_counter++)  = temp_epot;
	//	    			                 }
	//	    		                }
	//	    		                else{ // include all particles inside mc_rsweep
	//
	// 			        	            temp_id     = (double) (neighb[ind].get_particle(val).get_mynumber());
	// 			        	            temp_type   = (double) (neighb[ind].get_particle(val).get_mytype());
	//                                    temp_mass   =  neighb[ind].get_particle(val).get_mymass();
	//                                    // assign temporary position (coordinates shifted if neighbors on cpu boundary)
	//                                    temp_pos  =  curr_pos;
	//                                    temp_vel  =  neighb[ind].get_particle(val).get_myvelocity();
	//                                    temp_epot =  neighb[ind].get_particle(val).get_myepot();
	//
	//                                    // ADD TO BUFFER nb_nb_buffer_0
	//                                    *(nb_buffer[j]+buf_counter++) = temp_id;
	//                                    *(nb_buffer[j]+buf_counter++) = temp_type;
	//                                    *(nb_buffer[j]+buf_counter++) = temp_mass;
	//                                    *(nb_buffer[j]+buf_counter++) = temp_pos.x;
	//                                    *(nb_buffer[j]+buf_counter++) = temp_pos.y;
	//                                    *(nb_buffer[j]+buf_counter++) = temp_pos.z;
	//                                    *(nb_buffer[j]+buf_counter++) = temp_vel.x;
	//                                    *(nb_buffer[j]+buf_counter++) = temp_vel.y;
	//                                    *(nb_buffer[j]+buf_counter++) = temp_vel.z;
	//                                    *(nb_buffer[j]+buf_counter++) = temp_epot;
	//
	//	             	            }
	//
	//	                      }  // cut-off check
	//	    	     } //  end of particles loop
	//	    	     ind++; // neighbor cells counter
	//
	//    	    } // loop over each neighbor sub total
	//    	    // no of particles in each neighbor buffer
	//    	    total_particles = (buf_counter-1) * 0.1 ;     //10 attributes for each particle
	//            *(nb_buffer[j])  = (double) total_particles;  // each buffer first location has particles count
	//
	//	} // end of Neighbor cpu loop




	return sphere_old;
}

int get_cpu_rank(int inde0,int inde1,int inde2, MPI_Comm c_name){
	// gives my cpu rank with my own coordinates
	int nb_grid_coord[3],nb_rank;
    nb_grid_coord[0] = inde0;
    nb_grid_coord[1] = inde1;
    nb_grid_coord[2] = inde2;

    // get process rank from grid coord
    MPI_Cart_rank(c_name,nb_grid_coord,&nb_rank);
    return nb_rank;

}

// some utility methods

double scalar_prod(vec3d u, vec3d v){

	    // double vector - scalar product

	   double val = (u.x * v.x ) + (u.y * v.y) + (u.z * v.z);
       return val;
}

int iscalar_prod(ivec3d u, ivec3d v){

	   // integer vector - scalar product

       int val = (u.x * v.x ) + (u.y * v.y) + (u.z * v.z);
       return val;
}

vec3d vector_prod(vec3d u, vec3d v) {

	   //  float vector product

        vec3d w;
        w.x = u.y * v.z - u.z * v.y;
        w.y = u.z * v.x - u.x * v.z;
        w.z = u.x * v.y - u.y * v.x;
        return w;
}

ivec3d ivector_prod(ivec3d u, ivec3d v) {

	  //  integer vector - vector product

	   ivec3d w;
       w.x = u.y * v.z - u.z * v.y;
       w.y = u.z * v.x - u.x * v.z;
       w.z = u.x * v.y - u.y * v.x;
       return w;
}

double distance_vect(vec3d v1,vec3d v2){

   //  distance between two vectors

	double dist;
	vec3d temp_v;

	temp_v.x = v2.x - v1.x ;
	temp_v.y = v2.y - v1.y ;
	temp_v.z = v2.z - v1.z ;

	// compute distance as sum of square
	dist = (scalar_prod(temp_v,temp_v));

	return dist;
}

// end of utility methods
