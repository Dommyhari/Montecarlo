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

    if(prank == 0) {
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
    }



	ret_obj = loc_obj;

	return ret_obj;

} // make_cells


// Tested
cellblock make_particles(cellblock loc_obj, long tatoms_cpu, double* tbox_dim, ivec3d global_cell_dim, ivec3d loc_cpu_gcoord,ivec3d cpu_cell_dim,
		vector<long> atomnumber,vector<int> atomtypes,vector<double> atommass,vector<double> positions,vector<double> epot,int prank,int debug){

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

         loc_obj.set_cell(cell0,cell_index);  // update cell object with included particle

         if (debug == 1){

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
         }

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

          if((prank == 0) &&  (cell_index==100) && (l_counter == 0) && (debug==1)){

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


    if(prank == 0){

    	cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    	cout << "---                Inside sample zone method                      ---------" << endl;
    	cout << " My process CPU  : " << prank << endl;
    	cout << " Chosen window position : " << win_id << endl;
    	cout << " computed no of sample cell -- Sample_cell_counter : " << sam_cell_counter << endl;
    	cout << " Included no of cells from loop counter :   " << count << endl;

    	// print cell local coordinate and check??

    	ivec3d loc_coord;

    	cout << "=======================================================" << endl;
    	cout << "                       sample cells check     "  << endl;
    	cout << "=======================================================" << endl;
    	for( int i=0; i<sam_cell_counter;i++){

    		loc_coord = sample_cells[i].get_cell_loc_coord();
    		cout << " sample cell id : " << i  << endl;
    		cout << " sample cell local coordinate :  [ " <<loc_coord.x<<" "<<loc_coord.y<<" "<<loc_coord.z<<"  ]"<<endl;
    		cout << " no of particles : " << sample_cells[i].get_nparticles() << endl;

    	}


    	cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;

    }

    // choose one cell in random (belong to sample window)

	rand_cell = rand()%sam_cell_counter;
	cobj = sample_cells[rand_cell];

	do{
		//no of particles
		n_particles= cobj.get_nparticles();
		// choosing random placeholder
		rand_no=(long) rand()%n_particles;
	}while(cobj.get_particle(rand_no).get_mytype() !=2);

	atom = cobj.get_particle(rand_no);

	//--------------------- Unwanted replication - refer construct_sphere ---------------------//
	// make change only in particle object instance
	// swapping placeholder into carbon
	// HC: hardcoded to 1
    //atom.set_mytype(1);
    //-----------------------------------------------------------------------------------------//

    // construct sphere around selected particle and get the reference sphere before any trial move

	//old_sphere_config = construct_sphere(atom,bobj,win_id,file_name);

	return atom;
}

//celltype construct_sphere(particle pobj, cellblock bobj, int win_id,const char* filename,int prank,MPI_Comm comm_name){
//
//	celltype sphere_old;
//
//	// possible neighbor index as per window location
//
//	int window_x[8] = {-1,-1,+1,+1,-1,-1,+1,+1};
//	int window_y[8] = {-1,-1,-1,-1,+1,+1,+1,+1};
//	int window_z[8] = {-1,+1,-1,+1,-1,+1,-1,+1};
//
//	// 8 windows position in CPU (fixed positions -- window 0 to window 7)
//
//	ivec3d win_pos_0 = {0,0,0}; ivec3d win_pos_1 = {0,0,1}; ivec3d win_pos_2 = {1,0,0}; ivec3d win_pos_3 = {1,0,1};
//	ivec3d win_pos_4 = {0,1,0}; ivec3d win_pos_5 = {0,1,1}; ivec3d win_pos_6 = {1,1,0}; ivec3d win_pos_7 = {1,1,1};
//
//	ivec3d window_position[8] = {win_pos_0,win_pos_1,win_pos_2,win_pos_3,win_pos_4,win_pos_5,win_pos_6,win_pos_7};
//
//	// ==========================================================
//	// Detect neighbor processes as per window position
//	// ==========================================================
//    // select neighbors as per sample window position
//	int xfact = window_x[win_id]; int yfact = window_y[win_id]; int zfact = window_z[win_id];
//
//	// send neighbors (to whom I send)        ---- process communications (in z direction)
//	int x_send_phase[4] = { 0,xfact,xfact, 0}; int y_send_phase[4] = { 0, 0,yfact,yfact}; int z_send_phase[4] = {zfact,zfact,zfact,zfact};
//
//	// Receive neighbors (to whom I receive) ----- process communications (in z direction)
//	int x_recv_phase[4] = { 0,-xfact,-xfact, 0}; int y_recv_phase[4] = { 0, 0,-yfact,-yfact}; int z_recv_phase[4] = {-zfact,-zfact,-zfact,-zfact};
//
//    // send neighbors (to whom I send) -----  process communications (in x/y direction)
//	int x_send_next[3] = {xfact,xfact, 0 }; int y_send_next[3] = { 0,yfact,yfact }; int z_send_next[3] = { 0, 0, 0 };
//
//	// Receive neighbors (to whom I send) -----  process communications (in x/y direction)
//	int x_recv_next[3] = {-xfact,-xfact, 0 }; int y_recv_next[3] = { 0,-yfact,-yfact }; int z_recv_next[3] = { 0, 0, 0 };
//
//	// get my process global coordinate
//	ivec3d gcoord = get_cpu_gcoord(prank,comm_name);
//	int x = gcoord.x; int y = gcoord.y; int z = gcoord.z;
//
//	// variables to store the position of selected particle on my CPU and neighbor CPU's
//	// randomly selected particle positions (my CPU and rec_neighbor CPU)
//
//	double my_pos[3], rec_pos_0[3],rec_pos_1[3],rec_pos_2[3],rec_pos_3[3],rec_pos_4[3],rec_pos_5[3],rec_pos_6[3];
//
//	double* rec_pos[7]={rec_pos_0,rec_pos_1,rec_pos_2,rec_pos_3,rec_pos_4,rec_pos_5,rec_pos_6};
//
//	my_pos[0]=pobj.get_myposition().x; my_pos[1]=pobj.get_myposition().y; my_pos[2]=pobj.get_myposition().z;
//
//	return sphere_old;
//}

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
