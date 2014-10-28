/********************************************************************************
 *                        Monte_methods.cpp
 *
 *              contain all utility methods for Montecarlo
 *******************************************************************************/

#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<cstdlib>
#include<cmath>
#include "Monte_classes.h"
#include "Monte_globals.h"
#include "Monte_prototypes.h"

using namespace std;

void update_particle( int  cell_id, particle ref_atom ){


	   // method for updating particle attributes for the
	   // given cell and particle id

	   cellblock loc_obj;

	   // get cell
	   celltype  cell_obj = loc_obj.get_cell(cell_id) ;

	   particle p_obj;

	   // get particles count
	   long p_count = cell_obj.get_nparticles();

	   long loc_id = ref_atom.get_mynumber();

	   // loop over particles in cell
	   for ( long i=0; i<p_count; i++ ){
		   if ( cell_obj.get_particle(i).get_mynumber() == loc_id) { p_obj = cell_obj.get_particle(i); } // particle id matching
	   }

	   // updating atom attributes
	   p_obj.set_mytype(ref_atom.get_mytype());
	   p_obj.set_myposition(ref_atom.get_myposition().x,ref_atom.get_myposition().y,ref_atom.get_myposition().z);
	   p_obj.set_myvelocity(ref_atom.get_myvelocity().x,ref_atom.get_myvelocity().y,ref_atom.get_myvelocity().z);
	   p_obj.set_myepot(ref_atom.get_myepot());

}

void read_update_config (char* fname){

	// method for reading current file configuration and updating particle attributes
	// with in and across cpu process

	//----------------------------------------------//
	//            1 - File reading part             //
    //----------------------------------------------//

	const int LIMIT=1000;

	// file buffers to store atom data's read from last configuration
    vector<long>   fbuffer_id;
    vector<int>    fbuffer_type;
    vector<double> fbuffer_mass,fbuffer_epot;
    vector<double> fbuffer_pos_x,fbuffer_pos_y,fbuffer_pos_z;
    vector<double> fbuffer_vel_x,fbuffer_vel_y,fbuffer_vel_z;


    ifstream fin(fname,std::ios_base::in);

    // local holders
    long f_id; int f_type;
    double f_mass,f_posx,f_posy,f_posz,f_velx,f_vely,f_velz,f_epot;

    char headerline[LIMIT];

    cout<< " Reading simulation sphere file " << fname << endl;

    // crunch header part
    while (fin.get() == '#'){
          fin.getline(headerline,LIMIT,'\n');
          continue;
    }

    while( !fin.eof()){ // end of file check

    	// reading         // pushing
        fin>>f_id;        fbuffer_id.push_back(f_id);
        fin>>f_type;      fbuffer_type.push_back(f_type);
        fin>>f_mass;      fbuffer_mass.push_back(f_mass);
        fin>>f_posx;      fbuffer_pos_x.push_back(f_posx);
        fin>>f_posy;      fbuffer_pos_y.push_back(f_posy);
        fin>>f_posz;      fbuffer_pos_z.push_back(f_posz);
        fin>>f_velx;      fbuffer_vel_x.push_back(f_velx);
        fin>>f_vely;      fbuffer_vel_y.push_back(f_vely);
        fin>>f_velz;      fbuffer_vel_z.push_back(f_velz);
        fin>>f_epot;      fbuffer_epot.push_back(f_epot);
    }

    fin.clear(); // resetting bit states

    fin.close();

    //----------------------------------------------------------//
    //                     2- Updating part                     //
    //----------------------------------------------------------//

    // NOTE: later could be changed with n - particles in sphere count

    cellblock loc_obj;

    particle atom;

    cout<< "Updating particle attributes " << endl;

    for (int i=0; i< (fbuffer_id.size()-1); i++){

    	// get global cell coordinate from particle position
    	ivec3d glob_coord = cell_coordinate(fbuffer_pos_x.at(i), fbuffer_pos_y.at(i), fbuffer_pos_z.at(i));

		// get cpu rank from particle physical coordinate
		int loc_rank = get_cpu_rank(glob_coord);


		//   NOTE ********************
		// if loc_rank != myrank send particle to target cpu
		//****************************


		// get cpu glob position
		ivec3d cpu_fact = get_cpu_gcoord(loc_rank);

		// get cell local coordinate
		ivec3d cell_loc_coord = get_cell_loc_coord(glob_coord,cpu_fact);

		// to get memory index of cell in cell list
		int cell_index = cell_loc_coord.x + (cell_loc_coord.y * mc_cpu_cell_dim.x) + (cell_loc_coord.z * mc_cpu_cell_dim.x * mc_cpu_cell_dim.x);

		long particle_id = fbuffer_id.at(i);

		// defining particle attributes
		atom.set_mynumber(fbuffer_id.at(i));
		atom.set_mytype(fbuffer_type.at(i));
		atom.set_myposition(fbuffer_pos_x.at(i), fbuffer_pos_y.at(i), fbuffer_pos_z.at(i));
		atom.set_myvelocity(fbuffer_vel_x.at(i),fbuffer_vel_y.at(i),fbuffer_vel_z.at(i));
		atom.set_myepot(fbuffer_epot.at(i));

		// update particle
		update_particle(cell_index,atom);


    } // for loop

}

void construct_sphere(particle pobj, celltype cobj, char *fname){

	// Construct sphere around chosen particle and write the particles in (r_cut + r_sample) to configuration file
	//NOTE: now hard coded for mc_nbcells=6 neighbors

	ivec3d temp_glob, temp_loc, temp_block ; // temporary holder for cell global coordinates

	vec3d temp_pos, temp_vel, vec1, vec2;
	long temp_id; int temp_type;            // temporary holder for particle attributes
	double temp_mass, temp_epot, dist_check;

	char * filename = fname;                // file name as per simulation state from signature

	// object definitions

	particle temp_part;
	celltype sphere_cells[mc_nbcells]; // list of neighbor cells
	cellblock loc_obj;

	ofstream fout(filename, ios_base::out);

	// manipulators settings
	cout << fixed << right;

    //********** loop over Neighbor cells list **************

	for(int i=0; i<mc_nbcells; i++){

          // getting nbl cells, cpu and cell address

	      // get cell global coordinate
	      temp_glob = cobj.get_nbl_id(i);

	      // get cpu rank
	      int loc_rank = get_cpu_rank(temp_glob);

	      // get cpu block coordinate
	      temp_block = get_cpu_gcoord(loc_rank);

	      // cpu Foreign cell
	      if(mc_prank != loc_rank) {

	    	   // preliminary notifier
		       cout<<"Note from sphere construction pid: "<<mc_pid <<"\n"<<" Require process communication with : "<< loc_rank <<endl;

		       // unpack datas from MPI packed buffers

		       // after receiving and unpacking

		       long p_total; // should get from MPI receive
		       long p_c,v_c; // value holders

		       // loop over n_particles received

		       for( long k=0; k<p_total; k++ ){

		    	  //  p_c =  k+(size*3);
		    	  //  v_c =  k+(size*3*3);
		    	  //  temp_type   = buffer[k];
		    	  //  temp_id     = buffer[k+size];
		    	  //  temp_mass   = buffer[k+(size*2)];
		    	  //  temp_epot   = buffer[k+(size*3)]; // NOTE: need to be verified
		    	  //  temp_pos.x  = buffer[p_c++];
		    	  //  temp_pos.y  = buffer[p_c++];
		    	  //  temp_pos.z  = buffer[p_c++];
		    	  //  temp_vel.x  = buffer[v_c++];
		    	  //  temp_vel.y  = buffer[v_c++];
		    	  //  temp_vel.z  = buffer[v_c++];

		          temp_part.set_mytype(temp_type);
		          temp_part.set_mynumber(temp_id);
		          temp_part.set_mymass(temp_mass);
		          temp_part.set_myposition(temp_pos.x,temp_pos.y,temp_pos.z);
		          temp_part.set_myvelocity(temp_vel.x,temp_vel.y,temp_vel.z);
		          temp_part.set_myepot(temp_epot);

		          // ----- cutoff check

		          vec1 = pobj.get_myposition();
		          vec2 = temp_part.get_myposition();
		          dist_check = distance_vect(vec1,vec2);

		          //**************************************************************
		          //NOTE: to be updated in callee process to make dist_check
		          //**************************************************************

		          //if(dist_check <= mc_rsweep)    sphere_cells[i].add_particle(temp_part);

		       }

	      }//end of cpu foreign cells


	           // cpu native cells
	      else{

		       // get cell local coordinate
		       ivec3d cell_loc_coord = get_cell_loc_coord(temp_glob,temp_block);

		       // to get memory index of cell in cell list
		       int cell_index = cell_loc_coord.x + (cell_loc_coord.y * mc_cpu_cell_dim.x) + (cell_loc_coord.z * mc_cpu_cell_dim.x * mc_cpu_cell_dim.x);

		       long p_total = loc_obj.get_cell(cell_index).get_nparticles();

		       for( long pcount; pcount<p_total; pcount++ ){

		    	   // reading values
		    	   temp_id    = loc_obj.get_cell(cell_index).get_particle(pcount).get_mynumber();
		    	   temp_type  = loc_obj.get_cell(cell_index).get_particle(pcount).get_mytype();
		    	   temp_mass  = loc_obj.get_cell(cell_index).get_particle(pcount).get_mymass();
		    	   temp_pos   = loc_obj.get_cell(cell_index).get_particle(pcount).get_myposition();
		    	   temp_vel   = loc_obj.get_cell(cell_index).get_particle(pcount).get_myvelocity();
		    	   temp_epot  = loc_obj.get_cell(cell_index).get_particle(pcount).get_myepot();

		    	   // pre-check before assignment



			       vec1 = pobj.get_myposition();
			       vec2 = temp_pos;


			       dist_check = distance_vect(vec1,vec2);

			       // ----- cutoff check

			       if(dist_check <= mc_rsweep){

			    	   temp_part.set_mynumber(temp_id);
			    	   temp_part.set_mytype(temp_type);
			    	   temp_part.set_mymass(temp_mass);
			    	   temp_part.set_myposition(temp_pos.x,temp_pos.y,temp_pos.z);
			    	   temp_part.set_myvelocity(temp_vel.x,temp_vel.y,temp_vel.z);
			    	   temp_part.set_myepot(temp_epot);


			    	   sphere_cells[i].add_particle(temp_part);
			       }
		       }



          }// end of native cell


	      // file writing methods
	      long fp_total = sphere_cells[i].get_nparticles();

	      long my_number;
	      int my_type;
          double my_mass, my_epot;
          vec3d my_pos, my_vel;

          cout<<" Sphere constructor writing to file : " << filename<<endl;

          for(long fp; fp<fp_total; fp++){

        	  my_number  = sphere_cells[i].get_particle(fp).get_mynumber();
        	  my_type    = sphere_cells[i].get_particle(fp).get_mytype();
        	  my_mass    = sphere_cells[i].get_particle(fp).get_mymass();
        	  my_pos     = sphere_cells[i].get_particle(fp).get_myposition();
        	  my_vel     = sphere_cells[i].get_particle(fp).get_myvelocity(); // optional with #ifdef or so
              my_epot    = sphere_cells[i].get_particle(fp).get_myepot();



        	  // file flush
        	  fout << setw(8) << my_number
        		   << setw(6) << my_type
        		   << setw(6) << setprecision(8) << my_mass
        		   << setw(6) << setprecision(8) << my_pos.x
        		   << setw(6) << setprecision(8) << my_pos.y
        		   << setw(6) << setprecision(8) << my_pos.z
        		   << setw(6) << setprecision(8) << my_vel.x
        		   << setw(6) << setprecision(8) << my_vel.y
        		   << setw(6) << setprecision(8) << my_vel.z
        		   << setw(6) << setprecision(8) << my_epot
        		   << endl;

          }

          fout.close(); // closing outfile connection

	}// loop over sphere cells


}

void make_mc_nblist(celltype cobj){

	 // create neighbor list for the given cell

	 ivec3d g_max = mc_global_cell_dim;    // global cell array dimension - Max
	 ivec3d g_min = {0,0,0};               // global cell array dimension - Min

	 ivec3d ref = cobj.get_cell_glob_coord();
	 ivec3d ip_left,ip_right,ip_top,ip_bottom,ip_front,ip_back; // possible neighbors (NEED TO BE EXTENDED IF NECESSARY)


	 // general nbl assignment
	 ip_left.x  = --ref.x;       ip_left.y  = ref.y;      ip_left.z    = ref.z;
	 ip_right.x = ++ref.x;       ip_right.y = ref.y;      ip_right.z   = ref.z;
	 ip_front.x = ref.x;         ip_front.y = ++ref.y;    ip_front.z   = ref.z;
	 ip_back.x  = ref.x;         ip_back.y  = --ref.y;    ip_back.z    = ref.z;
	 ip_top.x   = ref.x;        ip_top.y    = ref.y;      ip_top.z     = ++ref.z;
	 ip_bottom.x = ref.x;       ip_bottom.y = ref.y;      ip_bottom.z  = --ref.z;

	 // cell on x boundary
	 if(ref.x == g_min.x) ip_left.x  = g_max.x;         // extreme left    - xmax
	 if(ref.x == g_max.x) ip_right.x  = g_min.x;        // extreme right   - xmin

	 // cell on y boundary
	 if(ref.y == g_min.y) ip_back.y   = g_max.y;        // extreme back    - ymax
	 if(ref.y == g_max.y) ip_front.y  = g_min.y;        // extreme front   - ymin

	 // cell on z boundary
	 if(ref.z == g_min.z) ip_bottom.z = g_max.z;        // extreme bottom  - zmax
	 if(ref.z == g_max.z) ip_top.z  = g_min.z;          // extreme top     - zmin

     // assigning nbl coordinates to cell object

	 cobj.set_my_left(ip_left);
	 cobj.set_my_right(ip_right);
	 cobj.set_my_front(ip_front);
	 cobj.set_my_back(ip_back);
	 cobj.set_my_bottom(ip_bottom);
	 cobj.set_my_top(ip_top);

}

ivec3d get_cpu_gcoord(int myrank){

	// compute cpu global coordinate based on process rank

	ivec3d cpu_array[mc_ncpus];  // CHECK: have to be initialized via corresponding method
	int count=0;

	for(int i=0;i<stacks_block;i++ ){
		for(int j=0;j<rows_block;j++){
			for(int k=0;k<cols_block;k++){

				cpu_array[count].x = k;
				cpu_array[count].y = j;
				cpu_array[count].z = i;

				count++;
			}
		}
	}
	return cpu_array[myrank];
}

ivec3d get_cell_loc_coord(ivec3d glob_coord, ivec3d cpu_glob_pos){

	// compute cell local coordinates from cell global coordinate and cpu glob position

	ivec3d cell_loc_coord;

	cell_loc_coord.x  = glob_coord.x - (cpu_glob_pos.x * mc_cpu_cell_dim.x);
	cell_loc_coord.y  = glob_coord.y - (cpu_glob_pos.y * mc_cpu_cell_dim.y);
	cell_loc_coord.z  = glob_coord.z - (cpu_glob_pos.z * mc_cpu_cell_dim.z);

	return cell_loc_coord;
}

void make_particles(){

     // method for creating particle objects and fill in cell container


	 std::cout<<" make particles for process : "<< mc_pid << endl;
	 std::cout<<" total particle objects created : "<< mc_tatoms_cpu <<endl;

	 //*************************************************
	 // NOTE:before make_particles():--  cell_block()--> make_cells()--> with cell id and glob coord

	 cellblock loc_obj;

	 //*************************************************

//	 particle* atom;

	 long m=0;
	 for(long i=0;i<mc_tatoms_cpu;i++){

		 particle atom;

		 // IGNORE THIS PART IF NECESSARY
//		 atom = new particle;
//		 atom->set_mynumber(mc_atomnumber.at(i));
//		 atom->set_mytype(mc_atomtypes.at(i));
//		 atom->set_mymass(mc_atommass.at(i));
//		 atom->set_myposition(mc_positions.at(m++),mc_positions.at(m++),mc_positions.at(m++));

		 // assigning attributes - old implementation

		 atom.mynumber = mc_atomnumber.at(i);
 		 atom.mytype   = mc_atomtypes.at(i);
 		 atom.mymass   = mc_atommass.at(i);

 		 atom.myposition.x = mc_positions.at(m++);
 		 atom.myposition.y = mc_positions.at(m++);
 		 atom.myposition.z = mc_positions.at(m++);


 		 //**************************************************
 		 // NOTE: generate velocity
 		 //**************************************************
 		 // by now should have constructed cell with boundary

 		 // also method to return cell_id from particle position


 		 // get global cell coordinate from particle position

   		 ivec3d glob_coord = cell_coordinate(atom.get_myposition().x, atom.get_myposition().y, atom.get_myposition().z);

 		 // get cpu rank from particle physical coordinate

 		 int loc_rank = get_cpu_rank(glob_coord);

 		 // get cpu glob position

 		 ivec3d cpu_fact = get_cpu_gcoord(loc_rank);

 		 // get cell local coordinate

 		 ivec3d cell_loc_coord = get_cell_loc_coord(glob_coord,cpu_fact);

 		 // to get memory index of cell in cell list
 		 int cell_index = cell_loc_coord.x + (cell_loc_coord.y * mc_cpu_cell_dim.x) + (cell_loc_coord.z * mc_cpu_cell_dim.x * mc_cpu_cell_dim.x);


         // add particle
 		 loc_obj.get_cell(cell_index).add_particle(atom);

 		 }

} //make particles

void create_maxwell_velocities(double mc_temp, ivec3d*  mc_restriction){

	  //   create and fill initial velocities for particles

      double vx,vy,vz;                      // velocity holders
      double imp_x,imp_y,imp_z;             // impulse holders
      double sum_x, sum_y, sum_z;           // sum holders
      double rtemp;                         // reduced temp variable

      int   dof_x,dof_y,dof_z;               // degrees of freedom holders
      long  tot_dof_x,tot_dof_y,tot_dof_z ;

      cellblock loc_obj;

      long cell_count = loc_obj.get_ncells();
      long m=0;

      for(long i=0; i<cell_count; i++){

    	  celltype cobj;
    	  long particles_count = loc_obj.get_cell(cell_count).get_nparticles();
    	  cobj = loc_obj.get_cell(i);

    	  for(long j=0; j<particles_count; j++){

        	  particle atom;
        	  int mytype;

        	  atom = cobj.get_particle(j);

        	  mytype = atom.get_mytype();
        	  dof_x  = mc_restriction[mytype].x;
        	  dof_y  = mc_restriction[mytype].y;
        	  dof_z  = mc_restriction[mytype].z;

        	  // reduced temp
        	  rtemp = sqrt(atom.get_mymass() * mc_temp);

        	  imp_x = get_gaussian(rtemp) * dof_x ;
        	  imp_y = get_gaussian(rtemp) * dof_y ;
        	  imp_z = get_gaussian(rtemp) * dof_z ;

        	  // initially storing impulse
        	  atom.set_myvelocity(imp_x,imp_y,imp_z);


        	  // summing total dof
        	  tot_dof_x += dof_x;
        	  tot_dof_y += dof_y;
        	  tot_dof_z += dof_z;

        	  // summing total impulse
        	  sum_x += imp_x;
        	  sum_y += imp_y;
        	  sum_z += imp_z;

    	  } // loop particles -1

//    	  sum_x = tot_dof_x == 0 ? 0.0 : sum_x / tot_dof_x;
//    	  sum_y = tot_dof_y == 0 ? 0.0 : sum_y / tot_dof_y;
//    	  sum_z = tot_dof_z == 0 ? 0.0 : sum_z / tot_dof_z;

      } // loop cells-1

	  sum_x =  sum_x / tot_dof_x;
	  sum_y =  sum_y / tot_dof_y;
	  sum_z =  sum_z / tot_dof_z;

      double new_vx,new_vy,new_vz;

	  for(long i=0; i<cell_count; i++){

    	  celltype cobj;
    	  long particles_count = loc_obj.get_cell(cell_count).get_nparticles();
    	  cobj = loc_obj.get_cell(i);

		  for(long j=0; j<particles_count; j++){

        	  particle atom;
        	  int mytype;
        	  atom = cobj.get_particle(j);
        	  mytype = atom.get_mytype();

        	  dof_x  = mc_restriction[mytype].x;
        	  dof_y  = mc_restriction[mytype].y;
        	  dof_z  = mc_restriction[mytype].z;

        	  // new velocity computation
        	  new_vx = ((atom.get_myvelocity().x - sum_x )/atom.get_mymass()) * dof_x;
        	  new_vy = ((atom.get_myvelocity().y - sum_y )/atom.get_mymass()) * dof_y;
        	  new_vz = ((atom.get_myvelocity().z - sum_z )/atom.get_mymass()) * dof_z;


        	  //updating velocity
        	  atom.set_myvelocity(new_vx,new_vy,new_vz);


		  } // loop particles -2

	  } // loop cells -2

	  // check msg
	  cout<<" initial velocities are computed for cpu id:"<< mc_pid <<endl;

}

void make_cells(){

    //    method for creating cell objects and fill in cellblock container
	//    should assign cell boundary and ids correctly from the field

	std::cout<<" make cells for process : "<< mc_pid << endl;
	std::cout<<" total cell objects created : "<< calc_ncells_cpu() <<endl;

	int cell_no     = 0;

	// I: some method calls
	stack_total   = calc_cpu_stack_total();
	ncells_stack  = calc_ncells_stack(); // NOTE: verify requirement
	ncells_cpu    = calc_ncells_cpu();
	int row_total = calc_cpu_row_total();
	int col_total = calc_cpu_col_total();

//	vec3d min_bound = {0.0,0.0,0.0};
//  vec3d max_bound = {0.0,0.0,0.0};

    ivec3d cell_gcoord = {0,0,0}; // initialization

    int cells_stack=0;

    int x_fac = cpu_gcoord.x;
    int y_fac = cpu_gcoord.y;
    int z_fac = cpu_gcoord.z;

    //*****************************************************************
    // create cell objects and fill into cell block container eventually
    // NOTE: hierarchy constraint: before make_particles() make_cells() should be called
    //*****************************************************************

    cellblock loc_obj;

    // loop over stack

	for (int stack_no=0; stack_no<stack_total; stack_no++){

		for(int row_count=0; row_count<row_total; row_count++){

			for(int col_count=0; col_count<col_total; col_count++){

				  celltype cell_obj;

        		  cell_no += (stack_no * cells_stack) + ( ncells_cpu * mc_pid ); // cell id computation

        		  cell_gcoord.x = col_count +  (col_total * x_fac);

        		  cell_gcoord.y = row_count +  (row_total * y_fac);

        		  cell_gcoord.z = stack_no  +  (stack_total * z_fac);

        		  // assign values
       		   	  cell_obj.set_cell_id(cell_no);

       		      cell_obj.set_glob_coord(cell_gcoord); // setting global cell coordinate

       		      loc_obj.add_cell(cell_obj); // adding cell into cel_block

			      cell_no++; //  updating cell id

			}//column_count

		}//row_count
		//z dim to be updated here
	}//stack loop

} // make_cells

void calc_mc_cpu_box(){

	//  computes cpu_box physical dimensions
    // NOTE: CHECK IF REQUIRED Scope verification

	// cpu box x vector
	mc_cpu_box_x.x = mc_simbox_x.x / mc_cpu_dim.x;
	mc_cpu_box_x.y = mc_simbox_x.y / mc_cpu_dim.x;
	mc_cpu_box_x.z = mc_simbox_x.z / mc_cpu_dim.x;

	// cpu box y vector
	mc_cpu_box_y.x = mc_simbox_y.x / mc_cpu_dim.y;
	mc_cpu_box_y.y = mc_simbox_y.y / mc_cpu_dim.y;
	mc_cpu_box_y.z = mc_simbox_y.z / mc_cpu_dim.y;

	// cpu box z vector
	mc_cpu_box_z.x = mc_simbox_z.x / mc_cpu_dim.z;
	mc_cpu_box_z.y = mc_simbox_z.y / mc_cpu_dim.z;
	mc_cpu_box_z.z = mc_simbox_z.z / mc_cpu_dim.z;

}

void calc_mc_global_cell_array(){

    //          MonteCarlo global cell array computation

	mc_global_cell_dim.x = (int) std::floor((mc_simbox_x.x / mc_cell_dim.x));
	mc_global_cell_dim.y = (int) std::floor((mc_simbox_y.y / mc_cell_dim.y));
	mc_global_cell_dim.z = (int) std::floor((mc_simbox_z.z / mc_cell_dim.z));
}

void calc_mc_cpu_cell_array(){

      //          MonteCarlo cpu cell array computation

	mc_cpu_cell_dim.x = (int) std::floor(mc_cpu_box_x.x/ mc_cell_dim.x);
	mc_cpu_cell_dim.y = (int) std::floor(mc_cpu_box_y.y/ mc_cell_dim.y);
	mc_cpu_cell_dim.z = (int) std::floor(mc_cpu_box_z.z/ mc_cell_dim.z);

	if( mc_cpu_cell_dim.x != (mc_global_cell_dim.x/mc_cpu_dim.x)) std::cerr<<" Incompatible global and cell array dimension "<<endl;
	if( mc_cpu_cell_dim.y != (mc_global_cell_dim.y/mc_cpu_dim.y)) std::cerr<<" Incompatible global and cell array dimension "<<endl;
	if( mc_cpu_cell_dim.z != (mc_global_cell_dim.z/mc_cpu_dim.z)) std::cerr<<" Incompatible global and cell array dimension "<<endl;
}

 void make_mc_tbox(){

  // transformation box concept is mandatory to implement PBC effect correctly
  //  compute box transformation matrix;
  //  NOTE: need to be extended or verified

  /* first unnormalized */

  mc_tbox_x = vector_prod( mc_simbox_y, mc_simbox_z );
  mc_tbox_y = vector_prod( mc_simbox_z, mc_simbox_x );
  mc_tbox_z = vector_prod( mc_simbox_x, mc_simbox_y );

  /* volume */
  mc_volume = scalar_prod( mc_simbox_x, mc_tbox_x );
  if ((0==mc_pid) && (0==mc_volume)) std::cerr("Box Edges are parallel.");

  /* normalization */
  mc_tbox_x.x /= mc_volume;  mc_tbox_x.y /= mc_volume;  mc_tbox_x.z /= mc_volume;
  mc_tbox_y.x /= mc_volume;  mc_tbox_y.y /= mc_volume;  mc_tbox_y.z /= mc_volume;
  mc_tbox_z.x /= mc_volume;  mc_tbox_z.y /= mc_volume;  mc_tbox_z.z /= mc_volume;

  /* squares of the box heights perpendicular to the faces */
  mc_height.x = 1.0 / scalar_prod(mc_tbox_x,mc_tbox_x);
  mc_height.y = 1.0 / scalar_prod(mc_tbox_y,mc_tbox_y);
  mc_height.z = 1.0 / scalar_prod(mc_tbox_z,mc_tbox_z);
}

ivec3d cell_coordinate(double x, double y, double z){

  //  cell_coord computes the (global) cell coordinates of a position
  //   NOTE: need to be extended or verified

  ivec3d coord;

  /* Map positions to boxes */
  coord.x = (int)(mc_global_cell_dim.x * (x*mc_tbox_x.x + y*mc_tbox_x.y + z*mc_tbox_x.z));
  coord.y = (int)(mc_global_cell_dim.y * (x*mc_tbox_y.x + y*mc_tbox_y.y + z*mc_tbox_y.z));
  coord.z = (int)(mc_global_cell_dim.z * (x*mc_tbox_z.x + y*mc_tbox_z.y + z*mc_tbox_z.z));

  /* rounding errors may put atoms slightly outside the simulation cell */
  /* in the case of no pbc they may even be far outside */
  if      (coord.x >= mc_global_cell_dim.x) coord.x = mc_global_cell_dim.x - 1;
  else if (coord.x < 0)                  coord.x = 0;
  if      (coord.y >= mc_global_cell_dim.y) coord.y = mc_global_cell_dim.y - 1;
  else if (coord.y < 0)                  coord.y = 0;
  if      (coord.z >= mc_global_cell_dim.z) coord.z = mc_global_cell_dim.z - 1;
  else if (coord.z < 0)                  coord.z = 0;

  return coord;
}

void calc_cpu_block_limits(){

	// NOTE : COULD BE REMOVED
    //*********************************************************
    // seemingly not important
	// compute CPU block limits in global dimension - called from rank 0
	//*********************************************************

	int x_min_val=0,y_min_val=0,z_min_val=0;
	int x_max_val=mc_cpu_cell_dim.x,y_max_val=mc_cpu_cell_dim.y,z_max_val=mc_cpu_cell_dim.z;

	for(auto i=0; i<stacks_block;i++){

		for(auto j=0; j<rows_block; j++){

			for(auto k=0; k<cols_block; k++){

                  cellblock block_obj;

                  block_obj.set_x_min(x_min_val);
                  block_obj.set_x_max(x_max_val-1);

                  block_obj.set_y_min(y_min_val);
                  block_obj.set_y_max(y_max_val-1);

                  block_obj.set_z_min(z_min_val);
                  block_obj.set_z_max(z_max_val-1);

                  x_min_val = x_max_val*(k+1);
                  x_max_val += x_max_val ;

			}// Column loop

			x_min_val = 0;
			x_max_val = mc_cpu_cell_dim.x;

			y_min_val =  y_max_val * (j+1);
			y_max_val += y_max_val;

		}// Row loop

		x_min_val = 0;
		x_max_val = mc_cpu_cell_dim.x;

		y_min_val =  0;
		y_max_val =  mc_cpu_cell_dim.y;

		z_min_val  =  z_max_val * (i+1);
		z_max_val +=  z_max_val;


	}// Stack loop



} //calc_cpu_block_limits

vec3d calc_cell_dim(double rcut,double rsample){

	//    compute MonteCarlo cell dimension

	mc_cell_dim.x = (mc_simbox_x.x) / (rcut + rsample);
	mc_cell_dim.y = (mc_simbox_y.y) / (rcut + rsample);
	mc_cell_dim.z = (mc_simbox_z.z) / (rcut + rsample);

	return mc_cell_dim;
}

int calc_ncells_cpu(){

	 //   compute no of cells in cpu

	int ret = (int) std::floor((mc_cpu_box_x.x * mc_cpu_box_y.y * mc_cpu_box_z.z )/(mc_cell_dim.x * mc_cell_dim.y * mc_cell_dim.z));
    return ret;
}

int calc_ncells_stack(){
	// get no of cells per stack
	// default stack pan (X * Y)
	// CHECK IF NECESSARY
	int ret = (int) std::floor((mc_cpu_box_x.x * mc_cpu_box_y.y)/(mc_cell_dim.x * mc_cell_dim.y));
	return ret;
}
//: D NOTE try to merge row, column and stack into single method
int calc_cpu_col_total(void){
    // compute no of cpu cell columns
	int ret = (int) std::floor(mc_cpu_box_x.x/ mc_cell_dim.x);
	return ret;
}


int calc_cpu_row_total(void){
    // compute no of cpu cell rows
	int ret = (int) std::floor(mc_cpu_box_y.y/ mc_cell_dim.y);
	return ret;
}

int calc_cpu_stack_total(){
	// compute no of cpu cell stacks
	// default cell stack direction z
	int	ret = (int) std::floor((mc_cpu_box_z.z/mc_cell_dim.z)); //to be checked and included
	return ret;
}

void calc_block_dim(){

	//  method for computing domain block dimensions

	cols_block   = mc_global_cell_dim.x / mc_cpu_cell_dim.x;
	rows_block   = mc_global_cell_dim.y / mc_cpu_cell_dim.y;
	stacks_block = mc_global_cell_dim.z / mc_cpu_cell_dim.z;
}

int get_cpu_rank(ivec3d cell_coord){

	// method to compute cpu_rank from global cell coordinates

	ivec3d red_vec;
	int proc_rank;

	red_vec.x = cell_coord.x / mc_cpu_cell_dim.x;
	red_vec.y = cell_coord.y / mc_cpu_cell_dim.y;
	red_vec.z = cell_coord.z / mc_cpu_cell_dim.z;

	proc_rank = (cols_block * rows_block * red_vec.z ) + (cols_block * red_vec.y) + (red_vec.x);
	return proc_rank;
}

// essential methods

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

	   vec3d w;
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

	dist = sqrt(scalar_prod(temp_v,temp_v));

	return dist;
}

// end of utility methods
