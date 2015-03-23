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
		int test_rank,ivec3d cpu_dim,double* data_list,ivec3d cpu_cell_dim){

	//########################################################################################################################################################
	//
    //                                      PHASE 1 - Random particle selection and communication with neighbors
	//
	//                   Random particle is chosen from sample cell list defined by appropriate window position
	//########################################################################################################################################################



	// possible neighbor index as per window location

	int window_x[8] = {-1,-1,+1,+1,-1,-1,+1,+1};
	int window_y[8] = {-1,-1,-1,-1,+1,+1,+1,+1};
	int window_z[8] = {-1,+1,-1,+1,-1,+1,-1,+1};


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

    // defining zone limits
    int x_start=0,  x_mid = cpu_cell_dim.x/2,    x_end = cpu_cell_dim.x;   // zone limit parameters in x direction
    int y_start=0,  y_mid = cpu_cell_dim.y/2,    y_end = cpu_cell_dim.y;   // zone limit parameters in y direction
    int z_start=0,  z_mid = cpu_cell_dim.z/2,    z_end = cpu_cell_dim.z;   // zone limit parameters in z direction

    // sample cell ranges

    ivec6d zone_limit_0 = {x_start,x_mid-1,y_start,y_mid-1,z_start,z_mid-1}; // window 0
    ivec6d zone_limit_1 = {x_start,x_mid-1,y_start,y_mid-1,z_mid,z_end-1};   // window 1
    ivec6d zone_limit_2 = {x_mid,x_end-1,y_start,y_mid-1,z_start,z_mid-1};   // window 2
    ivec6d zone_limit_3 = {x_mid,x_end-1,y_start,y_mid-1,z_mid,z_end-1};     // window 3

    ivec6d zone_limit_4 = {x_start,x_mid-1,y_mid,y_end-1,z_start,z_mid-1}; // window 4
    ivec6d zone_limit_5 = {x_start,x_mid-1,y_mid,y_end-1,z_mid,z_end-1};   // window 5
    ivec6d zone_limit_6 = {x_mid,x_end-1,y_mid,y_end-1,z_start,z_mid-1};   // window 6
    ivec6d zone_limit_7 = {x_mid,x_end-1,y_mid,y_end-1,z_mid,z_end-1};     // window 7

    ivec6d zone_limit[8]={zone_limit_0,zone_limit_1,zone_limit_2,zone_limit_3,zone_limit_4,zone_limit_5,zone_limit_6,zone_limit_7};

	// flags for coordinate shifting (if they are on edges!!)
	int flag [7][3]; // [neighbor] [x y z]
	for(int i=0;i<7;i++){
		flag [i][0] = 0; flag [i][1] = 0; flag [i][2] = 0;
	}

	// global cpu array dimension range
	ivec3d cpu_max = {cpu_dim.x-1,cpu_dim.y-1,cpu_dim.z-1};            // dimension - Max
	ivec3d cpu_min = {0,0,0};                                          // dimension - Min

	// boundary checks and flag assignments (ENSURE COORDINATES SHIFT IF REQUIRED)

	//flag[0][2] = window_position[win_id].z + (nb_gcoord[0][2] cpu_max.z)
	// NEED BETTER IMPLEMENTATION
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

	// variables for neighbor cell ranges
	int nb_xmin, nb_ymin, nb_zmin, nb_xmax, nb_ymax, nb_zmax, nb_x, nb_y, nb_z;

    ivec3d test;

    // Neighbor cells to be communicated as per sample window position
    // assign neighbor cell ranges (as per window position)
	nb_xmin = zone_limit[win_id].xmin; nb_xmax = zone_limit[win_id].xmax;
	nb_ymin = zone_limit[win_id].ymin; nb_ymax = zone_limit[win_id].ymax;
	nb_zmin = zone_limit[win_id].zmin; nb_zmax = zone_limit[win_id].zmax;

	// as window chosen on one corner of CPU; then cells on other corner for particular direction should be selected.

	// Maximum of CPU cell dimension
	ivec3d glob_c_dim = {cpu_cell_dim.x-1,cpu_cell_dim.y-1,cpu_cell_dim.z-1};

	nb_x = glob_c_dim.x - (window_position[win_id].x * glob_c_dim.x );
	nb_y = glob_c_dim.y - (window_position[win_id].y * glob_c_dim.y );
	nb_z = glob_c_dim.z - (window_position[win_id].z * glob_c_dim.z );

    // List of neighbor cells to be exported
	vector<celltype> nb_cells;
	long ncells = bobj.get_cell_list_size();

	// cell counters for neighbor communications
	int count_1=0, count_2=0, count_3=0, count_4=0, count_5=0, count_6=0, count_7=0;

	// particle counters for neighbor communications
    long tot_part_1 =0, tot_part_2 =0, tot_part_3 =0, tot_part_4 =0, tot_part_5 =0, tot_part_6 =0, tot_part_7 =0;


	// prepare NEIGHBOR cell list (to export) as per sample window id

	// Neighbor-1 (N1) (z direction)
    for(int i=nb_xmin;i<=nb_xmax;i++){
    	for(int j=nb_ymin;j<=nb_ymax;j++){
    			test.x = i; test.y = j; test.z = nb_z;
    			nb_cells.push_back(bobj.cell_with_lcoord(test,ncells));
    			tot_part_1 += bobj.cell_with_lcoord(test,ncells).get_nparticles(); // add particle count
    			count_1++;
    	}
    }

    // Neighbor-2 (N2) ( x-z direction)
	for(int j=nb_ymin;j<=nb_ymax;j++){
			test.x = nb_x; test.y = j; test.z = nb_z;
			nb_cells.push_back(bobj.cell_with_lcoord(test,ncells));
			tot_part_2 += bobj.cell_with_lcoord(test,ncells).get_nparticles(); // add particle count
			count_2++;
	}

	//Neighbor-3 (N4) (x-y-z direction)
	test.x = nb_x; test.y = nb_y; test.z = nb_z;
	nb_cells.push_back(bobj.cell_with_lcoord(test,ncells));
	tot_part_3 += bobj.cell_with_lcoord(test,ncells).get_nparticles(); // add particle count
	count_3++;

	//Neighbor-4 (N6) (y-z direction)
	for(int i=nb_xmin;i<=nb_xmax;i++){
			test.x = i; test.y = nb_y; test.z = nb_z;
			nb_cells.push_back(bobj.cell_with_lcoord(test,ncells));
			tot_part_4 += bobj.cell_with_lcoord(test,ncells).get_nparticles(); // add particle count
			count_4++;
	}

	//    //Neighbor-5 (N3) (x direction)
    for(int j=nb_ymin;j<=nb_ymax;j++){
    	for(int k=nb_zmin; k<=nb_zmax; k++){
    			test.x = nb_x; test.y = j; test.z = k;
    			nb_cells.push_back(bobj.cell_with_lcoord(test,ncells));
    			tot_part_5 += bobj.cell_with_lcoord(test,ncells).get_nparticles(); // add particle count
    			count_5++;
    	}
    }

	//    //Neighbor-6 (N5) (x-y direction)
	for(int k=nb_zmin;k<=nb_zmax;k++){
			test.x = nb_x; test.y = nb_y; test.z = k;
			nb_cells.push_back(bobj.cell_with_lcoord(test,ncells));
			tot_part_6 += bobj.cell_with_lcoord(test,ncells).get_nparticles(); // add particle count
			count_6++;
	}

    //    //Neighbor-7 (N7) (y direction)
    for(int i=nb_xmin;i<=nb_xmax;i++){
    	for(int k=nb_zmin; k<=nb_zmax; k++){
    			test.x = i; test.y = nb_y; test.z = k;
    			nb_cells.push_back(bobj.cell_with_lcoord(test,ncells));
    			tot_part_7 += bobj.cell_with_lcoord(test,ncells).get_nparticles(); // add particle count
    			count_7++;
    	}
    }

    int count_list[7]     = {count_1,count_2,count_3,count_4,count_5,count_6,count_7};                      // list of neighbor export cells counter
    long* to_send_list[7] = {&tot_part_1,&tot_part_2,&tot_part_3,&tot_part_4,&tot_part_5,&tot_part_6,&tot_part_7}; // list of neighbor export particles counter


//    if (prank == test_rank){
//       ivec3d nb_g_coord;
//       cout << " **************************************************** " << endl;
//       cout << " After allocation check : Neighbor cells  " << endl;
//       cout << " nb_cells.size() " << nb_cells.size()<< endl;
//
//       cout << " neighbor cells counter check " << endl;
//
//       for(int i=0;i<7;i++){
//    	   cout << " Neighbor "<< i<<" cell export count " << count_list[i] << "  particle export count   "<<tot_part_list[i]<< endl;
//       }
//
//       cout << "  Neighbor cells coordinate check " << endl;
//
//       for(int i=0; i<nb_cells.size(); i++){
//    	    nb_g_coord = nb_cells[i].get_cell_glob_coord();
//    	    cout <<" cell count " << i <<" " << "glob coordinate "<< nb_g_coord.x <<" "<< nb_g_coord.y<<" "<< nb_g_coord.z << endl;
//
//       }
//
//       cout << " **************************************************** " << endl;
//    }

	//***************************************************************************
	//                  My received neighbor part - fill buffers for my neighbors
	//***************************************************************************


	// buffer allocation for receive neighbors (NEED REDEFINITION)

	int buf_count = int(data_list[0]) * 10; // buffer_size

	// N1 N2 N4 N6 N3 N5 N7
	double *nb_0, *nb_1, *nb_2, *nb_3, *nb_4, *nb_5, *nb_6; // to be allocated and send

	// allocate memory (to send to neighbors)

	nb_0 = new double [tot_part_1*10]; nb_1= new double [tot_part_2*10]; nb_2 = new double [tot_part_3*10];
	nb_3 = new double [tot_part_4*10]; nb_4= new double [tot_part_5*10]; nb_5 = new double [tot_part_6*10];
	nb_6 = new double [tot_part_7*10];

	double* nb_buffer[7] = {nb_0,nb_1,nb_2,nb_3,nb_4,nb_5,nb_6};


    double rsweep      = data_list[1];  // sphere radius
    double sphere_wall = data_list[2];  // wall thickness

    // sweep distance computation (with squares)
	double r_full   = (rsweep + sphere_wall) * (rsweep + sphere_wall);
	double r_sphere = (rsweep) * (rsweep);

	vec3d ref_pos, curr_pos;                       // reference/current particle position

	// temporary holder for particle attributes
	double temp_id; double temp_type;
	double temp_mass, temp_epot, dist_check;
	vec3d temp_pos, temp_vel,simbox_size;
    int real_types;

	simbox_size.x = data_list[3];         // simbox x dimension
	simbox_size.y = data_list[4];         // simbox y dimension
	simbox_size.z = data_list[5];         // simbox z dimension
	real_types    = (int) data_list[6];   // no of real types

	long total_particles;                       // counter for total particle no in buffers
	int ind =0;                                   // neighbor cell counter
	long buf_counter,part_count;
    double* pos_holder, buffer_holder;
	int nb_ptr=0;

//	if(prank == test_rank){
//	    cout << "================================================" << endl;
//		cout << " Received attributes check data_list " << endl;
//	    cout << " rsweep   :" << rsweep << endl;
//	    cout << " sphere_wall  :" <<sphere_wall << endl;
//	    cout << " simbox_size.x :" << simbox_size.x << endl;
//	    cout << " simbox_size.y :" << simbox_size.y << endl;
//	    cout << " simbox_size.z :" << simbox_size.z << endl;
//	    cout << " real_types  : "<<real_types << endl;
//
//		cout << "================================================" << endl;
//	}

	for(int j=0; j<7;j++){ // loop over all neighbor CPU

		buf_counter=0;

	    // get received neighbor chosen particle position
	    pos_holder = rec_pos[j];
	    ref_pos.x = pos_holder[0]; ref_pos.y = pos_holder[1]; ref_pos.z = pos_holder[2];

	    if(prank == test_rank){
	    	cout<<"  "<<endl;
	    	cout << " my neighbor id : " << j << " & chosen particle position : " <<ref_pos.x<<" "<<ref_pos.y<<" "<<ref_pos.z << endl;
            cout << " flag[j][0]  :" << flag[j][0] << "flag[j][1]  :" << flag[j][1] << " flag[j][2] : "<< flag[j][2] << endl;
	    }

	    //total_particles=0;

	    for(int k=nb_ptr; k<(nb_ptr+count_list[j]); k++){ // loop over each neighbor  total export cells

	    	    part_count = nb_cells[k].get_nparticles(); //get total particles

//	    	    if(prank == test_rank){
//	    	    	cout << " Neigh ID : "<< j << " Export cell ID :  "<< k << "  particles count :   "<< part_count << endl;
//	    	    	cout << " Export cell global coordinate : " << nb_cells[k].get_cell_glob_coord().x<<" "<<nb_cells[k].get_cell_glob_coord().y
//	    	    			<<" "<<nb_cells[k].get_cell_glob_coord().y<< endl;
//	    	    }

	    	    total_particles=0;

	    	    for (long val=0; val<part_count; val++){ // loop over particles in each neighbor cell

		    		          curr_pos = nb_cells[k].get_particle(val).get_myposition();

		    		          // coordinate shifting as per flag initialization(for cells on edges of CPU)
		    		          curr_pos.x = curr_pos.x + (flag[j][0] * simbox_size.x);
		    		          curr_pos.y = curr_pos.y + (flag[j][1] * simbox_size.y);
		    		          curr_pos.z = curr_pos.z + (flag[j][2] * simbox_size.z);

		    		          // compute distance
		    		          dist_check = distance_vect(ref_pos,curr_pos);

		    		          // ----- cutoff check
		    		          if(dist_check <= r_full ){

	 		    		            // declare virtual particles on sphere boundary
	     	    		            if( dist_check > r_sphere ){
	     	    			             temp_type = (double) nb_cells[k].get_particle(val).get_mytype();

	     	    			             // ignore placeholders on sphere wall
		    			                 if (temp_type != (double) 2 ) { // HC: now hardcoded for placeholders

		    			        	           temp_id     = (double) (nb_cells[k].get_particle(val).get_mynumber());

		    			        	           // add particle id that are sent for sphere constructions
		    			        	           //list_ids[j].push_back(temp_id);

		    			        	           temp_type   = (double) (temp_type + real_types);
	                                           temp_mass   =  nb_cells[k].get_particle(val).get_mymass();

	                                           // assign temporary position (coordinates shifted if neighbors on CPU boundary)
	                                           temp_pos  =  curr_pos;
	                                           temp_vel  =  nb_cells[k].get_particle(val).get_myvelocity();
	                                           temp_epot =  nb_cells[k].get_particle(val).get_myepot();


	                                           // ADD TO BUFFER nb_buffer

	                                           *(nb_buffer[j]+buf_counter++)  = temp_id;
	                                           *(nb_buffer[j]+buf_counter++)  = temp_type;
	                                           *(nb_buffer[j]+buf_counter++)  = temp_mass;
	                                           *(nb_buffer[j]+buf_counter++)  = temp_pos.x;
	                                           *(nb_buffer[j]+buf_counter++)  = temp_pos.y;
	                                           *(nb_buffer[j]+buf_counter++)  = temp_pos.z;
	                                           *(nb_buffer[j]+buf_counter++)  = temp_vel.x;
	                                           *(nb_buffer[j]+buf_counter++)  = temp_vel.y;
	                                           *(nb_buffer[j]+buf_counter++)  = temp_vel.z;
	                                           *(nb_buffer[j]+buf_counter++)  = temp_epot;
	                                           total_particles++;
		    			                 }
		    		                }
		    		                else{ // include all particles inside mc_rsweep

	 			        	            temp_id     = (double) (nb_cells[k].get_particle(val).get_mynumber());
	 			        	            temp_type   = (double) (nb_cells[k].get_particle(val).get_mytype());
	                                    temp_mass   =  nb_cells[k].get_particle(val).get_mymass();

	                                    // assign temporary position (coordinates shifted if neighbors on cpu boundary)
	                                    temp_pos  =  curr_pos;
	                                    temp_vel  =  nb_cells[k].get_particle(val).get_myvelocity();
	                                    temp_epot =  nb_cells[k].get_particle(val).get_myepot();

	                                    // ADD TO BUFFER nb_nb_buffer_0

	                                    *(nb_buffer[j]+buf_counter++) = temp_id;
	                                    *(nb_buffer[j]+buf_counter++) = temp_type;
	                                    *(nb_buffer[j]+buf_counter++) = temp_mass;
	                                    *(nb_buffer[j]+buf_counter++) = temp_pos.x;
	                                    *(nb_buffer[j]+buf_counter++) = temp_pos.y;
	                                    *(nb_buffer[j]+buf_counter++) = temp_pos.z;
	                                    *(nb_buffer[j]+buf_counter++) = temp_vel.x;
	                                    *(nb_buffer[j]+buf_counter++) = temp_vel.y;
	                                    *(nb_buffer[j]+buf_counter++) = temp_vel.z;
	                                    *(nb_buffer[j]+buf_counter++) = temp_epot;

	                                    total_particles++;
		             	            }

		                      }  // cut-off check
		    	     } //  end of particles loop

	    	         //total_particles = buf_counter * 0.1;

	    	         if(prank == test_rank){
	    	         //if(j==0){
	    	          cout<<"  "<<endl;
			    	  cout << " my rank : "<< prank << " "<< " Neigh ID : "<< j << " Added  particles to buffer count :   "<< buf_counter * 0.1 << endl;

			    	 }

	    	    } // loop over each neighbor sub total
	            nb_ptr += count_list[j]; // updating neigbor pointer

	            *to_send_list[j] =  buf_counter * 0.1; //10 attributes for each particle


	} // end of Neighbor cpu loop

	if(prank == test_rank){

	    cout << " //////////////////////////////////////////////////////////" << endl;
   	    cout << "    Particle counted for export  "  << endl;
	    cout << "==========================================================0" << endl;

	    for(int i=0;i<7;i++){
		     cout<<"  "<<endl;
			 cout << "my_rank  : "<< prank << "  Neighbor : " << i << " total particles to send   :" << *to_send_list[i] << endl;
	    }
	    cout << " //////////////////////////////////////////////////////////" << endl;
	}

    //========================================================================================================
    //                                communicating buffer_size to neighbors
	//========================================================================================================

	// allocate list of send neighbors

	// allocate list of receive neighbors
	// particle counters for neighbor communications
    long to_recv_1 =0, to_recv_2 =0, to_recv_3 =0, to_recv_4 =0, to_recv_5 =0, to_recv_6 =0, to_recv_7 =0;

    long* to_recv_list[7] ={&to_recv_1, &to_recv_2, &to_recv_3, &to_recv_4, &to_recv_5, &to_recv_6, &to_recv_7 };

	// synchronize communication
	MPI_Barrier(comm_name);


	// z direction communication
	for (int ind=0;ind<4;ind++){

	     if (z%2==0){ // if z is even

	       //send_id = get_cpu_rank(x+x_send_phase[ind],y+y_send_phase[ind],z+z_send_phase[ind],comm_name);
	       send_id = get_cpu_rank(x+x_recv_phase[ind],y+y_recv_phase[ind],z+z_recv_phase[ind],comm_name);

	       MPI_Send(to_send_list[ind],1,MPI_LONG,send_id,0,comm_name);

	       if(prank == test_rank){
	    	   cout << " I sent :   " << *to_send_list[ind]  << endl;
	    	   cout << " Process :  " << prank <<" sent data to front target process : " << send_id << endl;
	    	   cout << "      "<< endl;
	       }

	       //rec_id = get_cpu_rank(x+x_recv_phase[ind],y+y_recv_phase[ind],z+z_recv_phase[ind],comm_name);
	       rec_id = get_cpu_rank(x+x_send_phase[ind],y+y_send_phase[ind],z+z_send_phase[ind],comm_name);


	       MPI_Recv(to_recv_list[ind],1,MPI_LONG,rec_id,1,comm_name,&stat);

	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" received data from back target process : " << rec_id << endl;
	    	   cout << " I received : " << *to_recv_list[ind] << endl;
	    	   cout << "      "<< endl;
	       }

	     }
	     else{       // if z is odd

	       //rec_id = get_cpu_rank(x+x_recv_phase[ind],y+y_recv_phase[ind],z+z_recv_phase[ind],comm_name);
	       rec_id = get_cpu_rank(x+x_send_phase[ind],y+y_send_phase[ind],z+z_send_phase[ind],comm_name);

	       MPI_Recv(to_recv_list[ind],1,MPI_LONG,rec_id,0,comm_name,&stat);

	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" received data from back target process : " << rec_id << endl;
	    	   cout << " I received : " << *to_recv_list[ind] << endl;
	    	   cout << "      "<< endl;
	       }

	       //send_id = get_cpu_rank(x+x_send_phase[ind],y+y_send_phase[ind],z+z_send_phase[ind],comm_name);
	       send_id = get_cpu_rank(x+x_recv_phase[ind],y+y_recv_phase[ind],z+z_recv_phase[ind],comm_name);

	       MPI_Send(to_send_list[ind],1,MPI_LONG,send_id,1,comm_name);

	       if(prank == test_rank){
	    	   cout << " I sent :  " << *to_send_list[ind]  << endl;
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

		   //send_id = get_cpu_rank(x+x_send_next[0],y+y_send_next[0],z+z_send_next[0],comm_name);
		   send_id = get_cpu_rank(x+x_recv_next[0],y+y_recv_next[0],z+z_recv_next[0],comm_name);

		   MPI_Send(to_send_list[4],1,MPI_LONG,send_id,2,comm_name);

	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" sent data to right target process : " << send_id << endl;
	    	   cout << " I sent :   " << *to_send_list[4]  << endl;
	    	   cout << "      "<< endl;
	       }

	       //rec_id = get_cpu_rank(x+x_recv_next[0],y+y_recv_next[0],z+z_recv_next[0],comm_name);
	       rec_id = get_cpu_rank(x+x_send_next[0],y+y_send_next[0],z+z_send_next[0],comm_name);

	       MPI_Recv(to_recv_list[4],1,MPI_LONG,rec_id,3,comm_name,&stat);

	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" received data from left target process : " << rec_id << endl;
	    	   cout << " I received :   " << *to_recv_list[4]<< endl;
	    	   cout << "      "<< endl;
	       }

	         // top-right & bottom-left comm
	       //send_id = get_cpu_rank(x+x_send_next[1],y+y_send_next[1],z+z_send_next[1],comm_name);
	       send_id = get_cpu_rank(x+x_recv_next[1],y+y_recv_next[1],z+z_recv_next[1],comm_name);


	       //MPI_Send(my_pos,3,MPI_DOUBLE,send_id,4,comm_name);
		   MPI_Send(to_send_list[5],1,MPI_LONG,send_id,4,comm_name);

	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" sent data to top - right target process : " << send_id << endl;
	    	   cout << " I sent :   " << *to_send_list[5]  << endl;
	    	   cout << "      "<< endl;
	       }

	       //rec_id = get_cpu_rank(x+x_recv_next[1],y+y_recv_next[1],z+z_recv_next[1],comm_name);
	       rec_id = get_cpu_rank(x+x_send_next[1],y+y_send_next[1],z+z_send_next[1],comm_name);

	       MPI_Recv(to_recv_list[5],1,MPI_LONG,rec_id,5,comm_name,&stat);

	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" received data from bottom - left target process : " << rec_id << endl;
	    	   cout << " I received :   " <<*to_recv_list[5] << endl;
	    	   cout << "      "<< endl;
	       }

	}
	else{            // if x is odd

		   //rec_id = get_cpu_rank(x+x_recv_next[0],y+y_recv_next[0],z+z_recv_next[0],comm_name);
		   rec_id = get_cpu_rank(x+x_send_next[0],y+y_send_next[0],z+z_send_next[0],comm_name);

		   // right & left comm
		   MPI_Recv(to_recv_list[4],1,MPI_LONG,rec_id,2,comm_name,&stat);

	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" received data from left target process : " << rec_id << endl;
	    	   cout << " I received :  " << *to_recv_list[4]  << endl;
	    	   cout << "      "<< endl;
	       }

	       //send_id = get_cpu_rank(x+x_send_next[0],y+y_send_next[0],z+z_send_next[0],comm_name);
	       send_id = get_cpu_rank(x+x_recv_next[0],y+y_recv_next[0],z+z_recv_next[0],comm_name);

	       MPI_Send(to_send_list[4],1,MPI_LONG,send_id,3,comm_name);

	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" sent data to right target process : " << send_id << endl;
	    	   cout << " I sent :   " << *to_send_list[4]  << endl;
	    	   cout << "      "<< endl;
	       }

	       //rec_id = get_cpu_rank(x+x_recv_next[1],y+y_recv_next[1],z+z_recv_next[1],comm_name);
	       rec_id = get_cpu_rank(x+x_send_next[1],y+y_send_next[1],z+z_send_next[1],comm_name);

	       // top-right & bottom-left comm
	       MPI_Recv(to_recv_list[5],1,MPI_LONG,rec_id,4,comm_name,&stat);

	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" received data from bottom - left target process : " << rec_id << endl;
	    	   cout << " I received :  " << *to_recv_list[5] << endl;
	    	   cout << "      "<< endl;
	       }

	       //send_id = get_cpu_rank(x+x_send_next[1],y+y_send_next[1],z+z_send_next[1],comm_name);
	       send_id = get_cpu_rank(x+x_recv_next[1],y+y_recv_next[1],z+z_recv_next[1],comm_name);

	       MPI_Send(to_send_list[5],1,MPI_LONG,send_id,5,comm_name);

	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" sent data to top - right target process : " << send_id << endl;
	    	   cout << " I sent :   " << *to_send_list[5]  << endl;
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

		   //send_id = get_cpu_rank(x+x_send_next[2],y+y_send_next[2],z+z_send_next[2],comm_name);
		   send_id = get_cpu_rank(x+x_recv_next[2],y+y_recv_next[2],z+z_recv_next[2],comm_name);


		   // top & bottom comm
		   MPI_Send(to_send_list[6],1,MPI_LONG,send_id,6,comm_name);

	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" sent data to bottom target process : " << send_id << endl;
	    	   cout << " I sent :   " << *to_send_list[6]  << endl;
	    	   cout << "      "<< endl;
	       }

	       //rec_id = get_cpu_rank(x+x_recv_next[2],y+y_recv_next[2],z+z_recv_next[2],comm_name);
	       rec_id = get_cpu_rank(x+x_send_next[2],y+y_send_next[2],z+z_send_next[2],comm_name);


	       MPI_Recv(to_recv_list[6],1,MPI_LONG,rec_id,7,comm_name,&stat);

	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" received data from top target process : " << rec_id << endl;
	    	   cout << " I received :" << *to_recv_list[6] << endl;
	    	   cout << "      "<< endl;
	       }
	}
	else{              // if y is odd

		   //rec_id = get_cpu_rank(x+x_recv_next[2],y+y_recv_next[2],z+z_recv_next[2],comm_name);
		   rec_id = get_cpu_rank(x+x_send_next[2],y+y_send_next[2],z+z_send_next[2],comm_name);

		   MPI_Recv(to_recv_list[6],1,MPI_LONG,rec_id,6,comm_name,&stat);

	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" received data from top target process : " << rec_id << endl;
	    	   cout << " I received : " <<*to_recv_list[6] << endl;
	    	   cout << "      "<< endl;
	       }

	       //send_id = get_cpu_rank(x+x_send_next[2],y+y_send_next[2],z+z_send_next[2],comm_name);
	       send_id = get_cpu_rank(x+x_recv_next[2],y+y_recv_next[2],z+z_recv_next[2],comm_name);

	       MPI_Send(to_send_list[6],1,MPI_LONG,send_id,7,comm_name);

	       if(prank == test_rank){
	    	   cout << " Process :  " << prank <<" sent data to bottom target process : " << send_id << endl;
	    	   cout << " I sent :  " << *to_send_list[6]  << endl;
	    	   cout << "      "<< endl;
	    	   cout << "  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$   " << endl;
	       }
	}

	// synchronize communication
	MPI_Barrier(comm_name);

	//========================================================================================================

	double *my_0, *my_1, *my_2, *my_3, *my_4, *my_5, *my_6; // to be received and filled in sphere cells

	// allocate memory (to receive from neighbors)
	my_0 = new double [*to_recv_list[0]*10]; my_1= new double [*to_recv_list[1]*10]; my_2 = new double [*to_recv_list[2]*10];
	my_3 = new double [*to_recv_list[3]*10]; my_4= new double [*to_recv_list[4]*10]; my_5 = new double [*to_recv_list[5]*10];
	my_6 = new double [*to_recv_list[6]*10];

	double* my_buffer[7] = {my_0,my_1,my_2,my_3,my_4,my_5,my_6};

	//*********************************************************************************************************************************
    //                                            send/receive buffers to/from corresponding neighbors
	//*********************************************************************************************************************************

	// synchronize communication
	MPI_Barrier(comm_name);

    // z direction

	for (int i=0;i<4;i++){
        if (z%2==0){
          // SWAP ORDER: receive from neighbors to whom I requested
          //           : send to neighbors from whom I received request
          // front comm
          MPI_Send(nb_buffer[i],(*to_send_list[i]*10),MPI_DOUBLE,get_cpu_rank(x+x_recv_phase[i],y+y_recv_phase[i],z+z_recv_phase[i],comm_name),7,comm_name);
          MPI_Recv(my_buffer[i],(*to_recv_list[i]*10),MPI_DOUBLE,get_cpu_rank(x+x_send_phase[i],y+y_send_phase[i],z+z_send_phase[i],comm_name),6,comm_name,&stat);
        }
        else{

          MPI_Recv(my_buffer[i],(*to_recv_list[i]*10),MPI_DOUBLE,get_cpu_rank(x+x_send_phase[i],y+y_send_phase[i],z+z_send_phase[i],comm_name),7,comm_name,&stat);
          MPI_Send(nb_buffer[i],(*to_send_list[i]*10),MPI_DOUBLE,get_cpu_rank(x+x_recv_phase[i],y+y_recv_phase[i],z+z_recv_phase[i],comm_name),6,comm_name);
        }
    } // loop  over neighbors (N1 N2 N4 N6)

    // Next phase even-even or odd-odd communications

    // x direction communication
    if (x%2 == 0){
            // right & left comm
          MPI_Send(nb_buffer[4],(*to_send_list[4]*10),MPI_DOUBLE,get_cpu_rank(x+x_recv_next[0],y+y_recv_next[0],z+z_recv_next[0],comm_name),9,comm_name);
          MPI_Recv(my_buffer[4],(*to_recv_list[4]*10),MPI_DOUBLE,get_cpu_rank(x+x_send_next[0],y+y_send_next[0],z+z_send_next[0],comm_name),8,comm_name,&stat);

          // top-right & bottom-left comm
          MPI_Send(nb_buffer[5],(*to_send_list[5]*10),MPI_DOUBLE,get_cpu_rank(x+x_recv_next[1],y+y_recv_next[1],z+z_recv_next[1],comm_name),11,comm_name);
          MPI_Recv(my_buffer[5],(*to_recv_list[5]*10),MPI_DOUBLE,get_cpu_rank(x+x_send_next[1],y+y_send_next[1],z+z_send_next[1],comm_name),10,comm_name,&stat);
    }
    else{

          MPI_Recv(my_buffer[4],(*to_recv_list[4]*10),MPI_DOUBLE,get_cpu_rank(x+x_send_next[0],y+y_send_next[0],z+z_send_next[0],comm_name),9,comm_name,&stat);
          MPI_Send(nb_buffer[4],(*to_send_list[4]*10),MPI_DOUBLE,get_cpu_rank(x+x_recv_next[0],y+y_recv_next[0],z+z_recv_next[0],comm_name),8,comm_name);

          // top-right & bottom-left comm
          MPI_Recv(my_buffer[5],(*to_recv_list[5]*10),MPI_DOUBLE,get_cpu_rank(x+x_send_next[1],y+y_send_next[1],z+z_send_next[1],comm_name),11,comm_name,&stat);
          MPI_Send(nb_buffer[5],(*to_send_list[5]*10),MPI_DOUBLE,get_cpu_rank(x+x_recv_next[1],y+y_recv_next[1],z+z_recv_next[1],comm_name),10,comm_name);
    }

    // y direction communication
    if (y%2 == 0){
    	  // top & bottom comm
          MPI_Send(nb_buffer[6],(*to_send_list[6]*10),MPI_DOUBLE,get_cpu_rank(x+x_recv_next[2],y+y_recv_next[2],z+z_recv_next[2],comm_name),13,comm_name);
          MPI_Recv(my_buffer[6],(*to_recv_list[6]*10),MPI_DOUBLE,get_cpu_rank(x+x_send_next[2],y+y_send_next[2],z+z_send_next[2],comm_name),12,comm_name,&stat);
    }
    else{
          MPI_Recv(my_buffer[6],(*to_recv_list[6]*10),MPI_DOUBLE,get_cpu_rank(x+x_send_next[2],y+y_send_next[2],z+z_send_next[2],comm_name),13,comm_name,&stat);
          MPI_Send(nb_buffer[6],(*to_send_list[6]*10),MPI_DOUBLE,get_cpu_rank(x+x_recv_next[2],y+y_recv_next[2],z+z_recv_next[2],comm_name),12,comm_name);
    }

//	if(prank == 7){
//
//		int j=2; // neighbor order
//		int par_coun = *to_send_list[j];
//		cout<<"=====================================" <<endl;
//		cout<< " Send particles check : " << endl;
//		cout<< " My rank : " << prank << endl;
//		cout<<"=====================================" <<endl;
//
//		for(int i=0;i<par_coun;i++){
//			cout << "  " <<endl;
//			cout << "Particle id   : " <<*(nb_buffer[j]+ (10*i+0)) << endl;
//			cout << "Particle type : " <<*(nb_buffer[j]+ (10*i+1)) << endl;
//			cout << "Particle mass : " <<*(nb_buffer[j]+ (10*i+2)) << endl;
//			cout << "Particle pos_x: " <<*(nb_buffer[j]+ (10*i+3)) << endl;
//			cout << "Particle pos_y: " <<*(nb_buffer[j]+ (10*i+4)) << endl;
//			cout << "Particle pos_z: " <<*(nb_buffer[j]+ (10*i+5)) << endl;
//			cout << "Particle vel_x: " <<*(nb_buffer[j]+ (10*i+6)) << endl;
//			cout << "Particle vel_y: " <<*(nb_buffer[j]+ (10*i+7)) << endl;
//			cout << "Particle vel_z: " <<*(nb_buffer[j]+ (10*i+8)) << endl;
//			cout << "Particle E_pot: " <<*(nb_buffer[j]+ (10*i+9)) << endl;
//            cout << "  " <<endl;
//
//			cout<<"=====================================" <<endl;
//		}
//
//	}
//
//
//	if(prank == test_rank){
//
//		int j=2; // neighbor order
//		int par_coun = *to_recv_list[j];
//		cout<<"=====================================" <<endl;
//		cout<< " Received particles check : " << endl;
//		cout<< " My rank : " << prank << endl;
//		cout<<"=====================================" <<endl;
//
//		for(int i=0;i<par_coun;i++){
//			cout << "  " <<endl;
//			cout << "Particle id   : " <<*(my_buffer[j]+ (10*i+0)) << endl;
//			cout << "Particle type : " <<*(my_buffer[j]+ (10*i+1)) << endl;
//			cout << "Particle mass : " <<*(my_buffer[j]+ (10*i+2)) << endl;
//			cout << "Particle pos_x: " <<*(my_buffer[j]+ (10*i+3)) << endl;
//			cout << "Particle pos_y: " <<*(my_buffer[j]+ (10*i+4)) << endl;
//			cout << "Particle pos_z: " <<*(my_buffer[j]+ (10*i+5)) << endl;
//			cout << "Particle vel_x: " <<*(my_buffer[j]+ (10*i+6)) << endl;
//			cout << "Particle vel_y: " <<*(my_buffer[j]+ (10*i+7)) << endl;
//			cout << "Particle vel_z: " <<*(my_buffer[j]+ (10*i+8)) << endl;
//			cout << "Particle E_pot: " <<*(my_buffer[j]+ (10*i+9)) << endl;
//            cout << "  " <<endl;
//
//			cout<<"=====================================" <<endl;
//		}
//
//	}

	// synchronize communication
	MPI_Barrier(comm_name);

    // **************************************************************************

    //========================================================================================================================
	//                                                  Phase-3
    //                                          allocate received buffers if any
    //                analyze received buffers and fill into sphere cells if necessary
	//========================================================================================================================

    // **************************************************************************
    // **************************************************************************
    celltype sphere_cell, sphere_old;    // contains particles in spherical domain for LOCAL MD simulation
    particle atom,my_atom;   // particle object

	//int sam_cell_counter = (cpu_cell_dim.x/2 * cpu_cell_dim.y/2 * cpu_cell_dim.z/2 );   // no of sample cells belonging to window

    // New definition
	int sam_cell_counter = (cpu_cell_dim.x* cpu_cell_dim.y* cpu_cell_dim.z );   // no of sample cells belonging to window

	celltype sample_cells[sam_cell_counter];  // sample cells for given window type

	long ncells_3 = bobj.get_cell_list_size();

    int count=0;

    // prepare SAMPLE cell list as per sample window id for my own CPU

//    for(int i=zone_limit[win_id].xmin;i<=zone_limit[win_id].xmax;i++){
//    	for(int j=zone_limit[win_id].ymin;j<=zone_limit[win_id].ymax;j++){
//    		for(int k=zone_limit[win_id].zmin;k<=zone_limit[win_id].zmax;k++){
//    			test.x = i; test.y = j; test.z = k;
//    			sample_cells[count] = bobj.cell_with_lcoord(test,ncells_3);
//    			// add to sample factor (INCLUDE LATER!!)
//    			//bobj.cell_with_lcoord(test,ncells_3).add_sample();
//                count++;
//    		}// k loop
//    	}// j loop
//    }// i loop



    for(int i=0;i<cpu_cell_dim.x;i++){
    	for(int j=0;j<cpu_cell_dim.y;j++){
    		for(int k=0;k<cpu_cell_dim.z;k++){
    			test.x = i; test.y = j; test.z = k;
    			sample_cells[count] = bobj.cell_with_lcoord(test,ncells_3);
    			// add to sample factor (INCLUDE LATER!!)
    			//bobj.cell_with_lcoord(test,ncells_3).add_sample();
                count++;
    		}// k loop
    	}// j loop
    }// i loop




    // Reference particle position
    vec3d ref_position = pobj.get_myposition(); // my own selected particle
    vec3d my_curr; // my current position

    // loop over my sample zone and add to sphere cells with in sweep
    // get reference particle sweep boundary
    double ref_xmin = ref_position.x - (rsweep + sphere_wall);
    double ref_ymin = ref_position.y - (rsweep + sphere_wall);   // minimum boundary
    double ref_zmin = ref_position.z - (rsweep + sphere_wall);

    //*************************************************************
    //                  fill sphere cells - My neighbor part
    //*************************************************************

    // loop over my_received buffers and add particles to sphere cell
    // shift to origin

    long p_count=0, buf_ind;
    double temp_position_x,temp_position_y,temp_position_z;

    for(int j=0; j<7; j++){ //loop over neighbors

    	p_count = *to_recv_list[j]; // no of particles in received buffer

    	buf_ind = 0;                      // initialize for each buffer list

    	if(p_count !=0){                  // check for empty buffer
         	for (long i=0; i< p_count; i++){

         		  atom.set_mynumber((long) *(my_buffer[j]+buf_ind++));
                  atom.set_mytype((int) *(my_buffer[j]+buf_ind++));
                  atom.set_mymass((double) *(my_buffer[j]+buf_ind++));

                  // shifting origin of coordinate system
                  temp_position_x = (double) (*(my_buffer[j]+buf_ind++) - ref_xmin);
                  temp_position_y = (double) (*(my_buffer[j]+buf_ind++) - ref_ymin);
                  temp_position_z = (double) (*(my_buffer[j]+buf_ind++) - ref_zmin);

                  atom.set_myposition(temp_position_x,temp_position_y,temp_position_z);
                  atom.set_myvelocity((double) *(my_buffer[j]+buf_ind++),(double) *(my_buffer[j]+buf_ind++),(double) *(my_buffer[j]+buf_ind++));
                  atom.set_myepot((double) *(my_buffer[j]+buf_ind++));

                  // appending particle instance to sphere cell list
                  sphere_cell.add_particle(atom);
    	    }
    	}
    }

    int rec_counter=0;
    if(prank == test_rank){

    	for(int i=0;i<7;i++){
    		rec_counter+= *to_recv_list[i];
    	}

    	cout<< "===============================================" << endl;
    	cout<< " " << endl;
        cout<< " Total received particles from my neighbors :  " << rec_counter << endl;
        cout<< " Total particles from neighbor in sphere cell : "<< sphere_cell.get_nparticles()<<endl;
        cout<< " Total no of sample cells  : " << sam_cell_counter << endl;
        cout<< " minimum boundary  :  " << ref_xmin <<"  "<<ref_ymin <<"  "<< ref_zmin<<"  "<< endl;
        cout<< " " << endl;
    	cout<< "===============================================" << endl;
    }



    //*************************************************************
    //                   My particle part -fill sphere cells
    //*************************************************************

    // loop over sample_cells -- perform distance check and add particles to sphere_cell
    // shift to origin
    int my_temp_type;

    for(int j=0; j<sam_cell_counter; j++){ // loop over sample cells

    	p_count = sample_cells[j].get_nparticles(); // no of particles in each sample cell

    	for(long i=0; i< p_count; i++){ // loop over particles

    		atom = sample_cells[j].get_particle(i);
            my_curr = atom.get_myposition();

	        dist_check = distance_vect(ref_position,my_curr);

	        // ----- cutoff check
	        if(dist_check <= r_full ){
		            // declare virtual particles on sphere boundary
		            if( dist_check > r_sphere ){

		            	 my_temp_type = atom.get_mytype();

		                 // ignore placeholders on sphere wall
	                     if (my_temp_type !=  2 ) { // HC: now hardcoded for placeholders

	                	        my_temp_type = my_temp_type + real_types ;

	                	        // shifting origin of coordinate system
	                            temp_position_x = atom.get_myposition().x - ref_xmin;
	                            temp_position_y = atom.get_myposition().y - ref_ymin;
	                            temp_position_z = atom.get_myposition().z - ref_zmin;

	                            atom.set_myposition(temp_position_x,temp_position_y,temp_position_z); // update shifted position
	                	        atom.set_mytype(my_temp_type); // update virtual type

                                sphere_cell.add_particle(atom);
	                     }
		            }
		            else{ // include all particles inside mc_rsweep

                	     // shifting origin of coordinate system
                         temp_position_x = atom.get_myposition().x - ref_xmin;
                         temp_position_y = atom.get_myposition().y - ref_ymin;
                         temp_position_z = atom.get_myposition().z - ref_zmin;

                         atom.set_myposition(temp_position_x,temp_position_y,temp_position_z); // update shifted position

		            	 sphere_cell.add_particle(atom);
		            }
	        } // end of cut-off check

    	} // loop over particles
    } // loop over sample cells

    // adding randomly chosen particle of own cpu
    my_atom = pobj;


    // shifting origin of coordinate system
    temp_position_x = my_atom.get_myposition().x - ref_xmin;
    temp_position_y = my_atom.get_myposition().y - ref_ymin;
    temp_position_z = my_atom.get_myposition().z - ref_zmin;

    my_atom.set_myposition(temp_position_x,temp_position_y,temp_position_z);

    // flip its type to real Carbon (SWAPPING PLACEHOLDER INTO CARBON)

    //sphere_cell.get_particle(my_atom.get_mynumber()).set_mytype(1);
    my_atom.set_mytype(1);

    // adding chosen particle
    sphere_cell.add_particle(my_atom);

    // creating reference sphere configuration
    sphere_old = sphere_cell;

    // ***********************************************************************
    // INCLUDE OTHER TWO TRIAL MOVES HERE
    // ***********************************************************************

    if(prank == test_rank){

    	cout<< "===============================================" << endl;
    	cout<< " " << endl;
        cout<< " Total particles for sphere :  " << rec_counter << endl;
        cout<< " Total particles (neighbor + my part ) sphere_cell : "<< sphere_cell.get_nparticles()<<endl;
        cout<< " Total particles (neighbor + my part ) sphere_old: "<< sphere_old.get_nparticles()<<endl;
        cout<< " " << endl;
    	cout<< "===============================================" << endl;
    }

    // *********************************************************************
    //                          file writing part
    // *********************************************************************

    //if (prank == test_rank){
    long file_number;
    int file_type;
    double file_mass, file_epot;
    vec3d file_pos, file_vel;

    // opening file stream object
    //stringstream convert;
    //convert << prank;
    string filename_local = "MC_sphere_config_p_"+to_string(prank)+".chkpt";

    //ofstream fout(filename, ios_base::out);
    ofstream fout(filename_local, ios_base::out);

    cout<<" Sphere constructor writing to file : " << filename << endl;

    //********** loop over revised sphere cell list - 2 **************
    long fp_total = sphere_cell.get_nparticles();


    // print IMD header -- check

    fout <<"#F A 1 1 1 3 0 0"<<endl;
    fout<<"#C number type mass x y z vx vy vz epot"<<endl;
    fout<<"#X 40.009552 0 0 "<<endl;
    fout<<"#Y 0 40.009552 0 "<<endl;
    fout<<"#Z 0 0 40.009552 "<<endl;
    fout<<"#E "<<endl;


    for(long fp=0; fp<fp_total; fp++){

        	file_number  = sphere_cell.get_particle(fp).get_mynumber();
        	file_type    = sphere_cell.get_particle(fp).get_mytype();
        	file_mass    = sphere_cell.get_particle(fp).get_mymass();
        	file_pos     = sphere_cell.get_particle(fp).get_myposition();
        	file_vel     = sphere_cell.get_particle(fp).get_myvelocity(); // optional with #ifdef or so
        	file_epot    = sphere_cell.get_particle(fp).get_myepot();

            	  // file flush

            fout <<file_number
            << "  " << file_type
            << "  " << setprecision(6) << file_mass
            << "  " << setprecision(6) << file_pos.x
            << "  " << setprecision(6) << file_pos.y
            << "  " << setprecision(6) << file_pos.z
            << "  " << setprecision(6) << file_vel.x
            << "  " << setprecision(6) << file_vel.y
            << "  " << setprecision(6) << file_vel.z
            << "  " << setprecision(6) << file_epot
            << endl;


    }

    fout.close(); // closing outfile connection

    //}



    // DELETE BUFFER MEMORY ACCORDINGLY
    delete[] nb_0; delete[] nb_1;delete[] nb_2; delete[] nb_3;
    delete[] nb_4;delete[] nb_5;delete[] nb_6;

    delete[] my_0; delete[] my_1;delete[] my_2; delete[] my_3;
    delete[] my_4;delete[] my_5;delete[] my_6;

	return sphere_old;
}

// pass IMD_executable and param file as arguments
void do_local_mdrun(string bin_name,string param_name,int prank){

	       // method to run local MD simulation
	       // later extend format for Lammps,etc..

	// make binary and parameter file clone

	string srank    =  to_string(prank);
	string spath    =  "./" ;
	string space    =  "  " ;
    string verflag  =  " -p ";
    string nline    =  " \n ";
    string uscore   =  "_";
    string f_extent =  ".param";
    string copy_cmd =  "cp"; // bash command
    string remov_cmd=  "rm"; // bash command

    string ensem       = "glok";
    string coord_name  = "MC_sphere_config_p_"+srank+".chkpt";
    string out_file    = "./sphere_out/MC_sphere_config_p_"+srank+"_relax";

    string filename_local = "local_md_part1_"+to_string(prank)+".param";

    // writing part1 param file
    ofstream fout(filename_local, ios_base::out);

    fout << "ensemble" <<"           "<< ensem << endl;
    fout << " " << endl;
    fout << "coordname"<<"           "<< coord_name << endl;
    fout << " " << endl;
    fout << "outfiles"<<"           "<< out_file << endl;
    fout << " " << endl;
    fout.close(); // closing outfile connection


    string param_1    =   filename_local;
    string param_2    =   "local_md_part2.param";
    string param_out  =   "local_md_"+srank+".param";

    string merge_cmd = "cat" + space + param_1 + space + param_2 + space+
    		            ">>" + space + param_out;
    system( merge_cmd.c_str());

	string new_bin   = bin_name + uscore + srank ;

	string new_param = param_out;

    // copy commands
    string cmd1 = copy_cmd + space + bin_name +space + new_bin;

    //string cmd2 = copy_cmd + space + param_name+f_extent +space + new_param;

    // remove commands
    string cmd3 = remov_cmd + space + new_bin;
    string cmd4 = remov_cmd + space + new_param;

    // IMD command
    string exec_string = spath + new_bin + space + verflag + space + new_param + nline;

    // getting char* from string
    const char * cmd1_ptr = cmd1.c_str();
    const char * cmd3_ptr = cmd3.c_str();
    const char * cmd4_ptr = cmd4.c_str();
    const char * expr = exec_string.c_str();

    // file copying
    system(cmd1_ptr); // copy binary clone

	system(expr);     // run MD simulation

	//file removing
    system(cmd3_ptr); // remove binary clone
    system(cmd4_ptr); // remove param clone
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
	//dist = (scalar_prod(temp_v,temp_v));

	dist = temp_v.x*temp_v.x + temp_v.y*temp_v.y +temp_v.z*temp_v.z ;

	return dist;
}

// end of utility methods
