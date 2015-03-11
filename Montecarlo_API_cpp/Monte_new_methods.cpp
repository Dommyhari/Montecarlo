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

vec3d* make_mc_tbox_vector(vec3d* simbox,int prank){

 // transformation box concept is mandatory to implement PBC effect correctly
 //  compute box transformation matrix;
 //  NOTE: need to be extended or verified

 vec3d 	tbox_x, tbox_y,tbox_z;
 double mc_volume;                                                    // some internal attributes for mapping
 vec3d mc_height;                                                      // used for global cell coordinates
 vec3d v_list[3];

 /* first unnormalized */

 tbox_x = vector_prod( simbox[1], simbox[2] );
 tbox_y = vector_prod( simbox[2], simbox[0] );
 tbox_z = vector_prod( simbox[0], simbox[1] );

 /* volume */
 mc_volume = scalar_prod( simbox[0], tbox_x );

 if ((0==prank) && (0==mc_volume))
	  std::cerr <<"Box Edges are parallel.";

 /* normalization */
 tbox_x.x /= mc_volume;  tbox_x.y /= mc_volume;  tbox_x.z /= mc_volume;
 tbox_y.x /= mc_volume;  tbox_y.y /= mc_volume;  tbox_y.z /= mc_volume;
 tbox_z.x /= mc_volume;  tbox_z.y /= mc_volume;  tbox_z.z /= mc_volume;

 /* squares of the box heights perpendicular to the faces */
 mc_height.x = 1.0 / scalar_prod(tbox_x,tbox_x);
 mc_height.y = 1.0 / scalar_prod(tbox_y,tbox_y);
 mc_height.z = 1.0 / scalar_prod(tbox_z,tbox_z);

 v_list[0] = tbox_x; v_list[1] = tbox_y; v_list[2] = tbox_z;

 return v_list;
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

ivec3d cell_coordinate(double* pos,double* tbox,ivec3d global_cell_dim){

  //  cell_coord computes the (global) cell coordinates of a position
  //   NOTE: need to be extended or verified

  double x,y,z;
  ivec3d coord;
  vec3d tbox_x,tbox_y,tbox_z;
  int m=0;
  // assignment
  x = pos[0]; y = pos[1]; z = pos[2];


  tbox_x.x = tbox[0]; tbox_x.y = tbox[1]; tbox_x.z = tbox[2];
  tbox_y.x = tbox[3]; tbox_y.y = tbox[4]; tbox_y.z = tbox[5];
  tbox_z.x = tbox[6]; tbox_z.y = tbox[7]; tbox_z.z = tbox[8];

  /* Map positions to boxes */
  coord.x = (int)(global_cell_dim.x * (x*tbox_x.x + y*tbox_x.y + z*tbox_x.z));
  coord.y = (int)(global_cell_dim.y * (x*tbox_y.x + y*tbox_y.y + z*tbox_y.z));
  coord.z = (int)(global_cell_dim.z * (x*tbox_z.x + y*tbox_z.y + z*tbox_z.z));


  /* rounding errors may put atoms slightly outside the simulation cell */
  /* in the case of no pbc they may even be far outside */

  if      (coord.x >= global_cell_dim.x) coord.x = global_cell_dim.x - 1;
  else if (coord.x < 0)                  coord.x = 0;
  if      (coord.y >= global_cell_dim.y) coord.y = global_cell_dim.y - 1;
  else if (coord.y < 0)                  coord.y = 0;
  if      (coord.z >= global_cell_dim.z) coord.z = global_cell_dim.z - 1;
  else if (coord.z < 0)                  coord.z = 0;

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

    cout << " ========================================================================" << endl;
    cout << " my process id :  " << prank << endl;
    cout << " CPU coordinate  : " << gcoord.x <<" " << gcoord.y <<" "<<gcoord.z <<" "<< endl;
    cout << " Cell allocation check :  no of cells  :  " <<  loc_obj.get_cell_list_size() << endl;
    cout << " cell ID :" << loc_obj.get_cell(cell_no).get_cell_id();
    cout << "   cell global coordinate        "<< loc_obj.get_cell(cell_no).get_cell_glob_coord().x<<" "<<loc_obj.get_cell(cell_no).get_cell_glob_coord().y
      		  <<" "<<loc_obj.get_cell(cell_no).get_cell_glob_coord().z<<endl;
    cout << "   cell local coordinate        "<< loc_obj.get_cell(cell_no).get_cell_loc_coord().x<<" "<<loc_obj.get_cell(cell_no).get_cell_loc_coord().y
      		  <<" "<<loc_obj.get_cell(cell_no).get_cell_loc_coord().z<<endl;
    cout << " ========================================================================" << endl;




	ret_obj = loc_obj;

	return ret_obj;

} // make_cells



cellblock make_particles(cellblock loc_obj, long tatoms_cpu, double* tbox_dim, ivec3d global_cell_dim, ivec3d loc_cpu_gcoord,ivec3d cpu_cell_dim,
		vector<long> atomnumber,vector<int> atomtypes,vector<double> atommass,vector<double> positions,vector<double> epot,int prank){

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


 		 cell_glob_coord = cell_coordinate(inst_pos,tbox_dim,global_cell_dim);

 		 ivec3d cpu_fact = loc_cpu_gcoord;

 		 // get cell local coordinate
 		 ivec3d cell_loc_coord = get_cell_loc_coord(cell_glob_coord,cpu_fact,cpu_cell_dim);

 		 // to get memory index of cell in cell list

 		 cell_index = cell_loc_coord.x + (cell_loc_coord.y * cpu_cell_dim.x) + (cell_loc_coord.z * cpu_cell_dim.y * cpu_cell_dim.x);


 		 //int cell_index = loc_obj.cell_with_lcoord(cell_loc_coord,total_cells).get_cell_id(); // expensive!!!!

      	 // just for testing add all particles to cell 1

 		 celltype cell0 = loc_obj.get_cell(cell_index);

 		 cell0.add_particle(atom);

         loc_obj.set_cell(cell0,cell_index);

//       loc_obj.get_cell(cell_index).add_particle(atom,cell_index);


         if((prank == 0) &&  (cell_index==100) && (l_counter == 0)){
             cout << "  Super check: Inside make particles " << endl;
             cout << "  cell index  : " << cell_index << endl;
             cout << "  cell size :  " <<  loc_obj.get_cell(cell_index).get_nparticles() << endl;
             cout << "  my no :" << loc_obj.get_cell(cell_index).get_particle(0).get_mynumber() << endl;
    	     cout << "  my type :" << loc_obj.get_cell(cell_index).get_particle(0).get_mytype() << endl;
        	 cout << "  my mass :" << loc_obj.get_cell(cell_index).get_particle(0).get_mymass() << endl;
        	 cout << "  my pos:x " << loc_obj.get_cell(cell_index).get_particle(0).get_myposition().x << endl;
        	 cout << "  my pos:y " << loc_obj.get_cell(cell_index).get_particle(0).get_myposition().y << endl;
        	 cout << "  my pos:z " << loc_obj.get_cell(cell_index).get_particle(0).get_myposition().z << endl;
    	     cout << "  epot      "<< loc_obj.get_cell(cell_index).get_particle(0).get_myepot() << endl;

    	     l_counter++;

          }


 		 }// end of for loop

	     // returning updated cellblock
	     ret_obj = loc_obj;

	     // clear STL container (once particle objects are created)

	     // (NOT REQUIRED DURING TESTING PHASE -- INCLUDE LATER -- in sync with fill_mc_container)
//	       mc_atomnumber.clear();
//	       mc_atomtypes.clear();
//         mc_atommass.clear();
//         mc_positions.clear();
//         mc_epot.clear();

	     return ret_obj;

} //make particles



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
