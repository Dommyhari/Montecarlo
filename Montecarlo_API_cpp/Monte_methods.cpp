/********************************************************************************
 *                        Monte_methods.cpp
 *
 *              contain all utility methods for Montecarlo
 *******************************************************************************/

#include "Monte_classes.h"
#include "Monte_globals.h"
#include "Monte_prototypes.h"

using namespace std;

int get_cpu_neighbor(int inde0,int inde1,int inde2){
	// gives my cpu rank with my own coordinates
	int nb_grid_coord[3],nb_rank;
    nb_grid_coord[0] = inde0;
    nb_grid_coord[1] = inde1;
    nb_grid_coord[2] = inde2;

    // get process rank from grid coord
    MPI_Cart_rank(comm_name,nb_grid_coord,&nb_rank);
    return nb_rank;

}

void setup_config(){

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

	// change it accordingly to compute from Vitrual topology
	get_cpu_gcoord(mc_prank);

	// compute domain block dimension
	calc_block_dim();

}

void update_particle( cellblock bobj,int  cell_id, particle ref_atom ){


	   // method for updating particle attributes for the
	   // given cell block,cell index and particle id
	   // get cell

	   // C++:  make use of copy constructor after checking

	   // get particles count
	   long p_count = bobj.get_cell(cell_id).get_nparticles();

	   long loc_id = ref_atom.get_mynumber();

	   // loop over particles in cell
	   for ( long i=0; i<p_count; i++ ){
		   if ( bobj.get_cell(cell_id).get_particle(i).get_mynumber() == loc_id) { // particle id matching
			   // updating atom attributes
		       bobj.get_cell(cell_id).get_particle(i).set_mytype(ref_atom.get_mytype());
			   bobj.get_cell(cell_id).get_particle(i).set_myposition(ref_atom.get_myposition().x,ref_atom.get_myposition().y,ref_atom.get_myposition().z);
			   bobj.get_cell(cell_id).get_particle(i).set_myvelocity(ref_atom.get_myvelocity().x,ref_atom.get_myvelocity().y,ref_atom.get_myvelocity().z);
			   bobj.get_cell(cell_id).get_particle(i).set_myepot(ref_atom.get_myepot());
		   }
	   }

}

void do_local_mdrun(string bin_name,string param_name){

	       // method to run local MD simulation
	       // later extend format for Lammps,etc..

	// make binary and parameter file clone

	string srank = to_string(mc_prank);
	string spath   =  "./" ;
	string space   =  "  " ;
    string verflag =  " -p ";
    string nline   =  " \n ";
    string uscore  =  "_";
    string copy_cmd  = "cp"; // bash command
    string remov_cmd = "rm"; // bash command

	string new_bin   = bin_name + uscore + srank ;
    string new_param = param_name + uscore + srank;

    // copy commands
    string cmd1 = copy_cmd + space + bin_name +space + new_bin;
    string cmd2 = copy_cmd + space + param_name +space + new_param;

    // remove commands
    string cmd3 = remov_cmd + space + new_bin;
    string cmd4 = remov_cmd + space + new_param;

    // IMD command
    string exec_string = spath + new_bin + space + verflag + space + new_param + nline;

    // getting char* from string
    const char * cmd1_ptr = cmd1.c_str();
    const char * cmd2_ptr = cmd2.c_str();
    const char * cmd3_ptr = cmd3.c_str();
    const char * cmd4_ptr = cmd4.c_str();
    const char * expr = exec_string.c_str();

    // file copying
    system(cmd1_ptr); // copy binary clone
    system(cmd2_ptr); // copy param clone

	system(expr);     // run MD simulation

	//file removing
    system(cmd3_ptr); // remove binary clone
    system(cmd4_ptr); // remove param clone

}

// sample_seed parameter have to be used at callee function/method
particle sample_zone(cellblock bobj,int win_id){

	// NEED CHANGES

	celltype cobj;
	celltype sample_cells[8];

	// local cell coordinates zone
    ivec3d test;
    int count=0,rand_cell;

    ivec6d zone_limit[8];

    // NOTE: HARDCODED FOR 8 CELLS IN A WINDOW
    zone_limit[0] = {0,1,0,1,0,1};
    zone_limit[1] = {0,1,0,1,1,2};
    zone_limit[2] = {1,2,0,1,0,1};
    zone_limit[3] = {1,2,0,1,1,2};
    zone_limit[4] = {0,1,1,2,0,1};
    zone_limit[5] = {0,1,1,2,1,2};
    zone_limit[6] = {1,2,1,2,0,1};
    zone_limit[7] = {1,2,1,2,1,2};


    // prepare cell list as per sample window id

    for(int i=zone_limit[win_id].xmin;i<=zone_limit[win_id].xmax;i++){
    	for(int j=zone_limit[win_id].ymin;j<=zone_limit[win_id].ymax;j++){
    		for(int k=zone_limit[win_id].zmin;k<=zone_limit[win_id].zmax;k++){

    			test.x = i; test.y = j; test.z = k;
    			sample_cells[count] = bobj.cell_with_lcoord(test);

                count++;
    		}
    	}
    }

	long rand_no,n_particles;
	particle atom;

	//srand(time(NULL));
	//srand(sample_seed);

	// NOTE: HARDCODED FOR 8 CELLS IN A WINDOW
	rand_cell = rand()%8;

	cobj = sample_cells[rand_cell];

	do{
		//no of particles
		n_particles= cobj.get_nparticles();
		// choosing random placeholder
		rand_no=(long) rand()%n_particles;
	}while(cobj.get_particle(rand_no).get_mytype() !=2);

	atom = cobj.get_particle(rand_no);



	// make change only in particle object
	// swapping placeholder into carbon
	// HC: hardcoded to 1
    atom.set_mytype(1);


    // construct sphere around selected particle
    construct_sphere(atom,bobj,win_id,file_name);

	return atom;
}

void read_update_config (int accep_tag,int win_id,char* fname,particle pobj,cellblock bobj){

	// NEED CHANGES

	// method for reading current file configuration and updating particle attributes
	// with in and across cpu process

	//----------------------------------------------//
	//            1 - File reading part             //
    //----------------------------------------------//

	const int LIMIT=1000;

	// file buffers to store atom data's read from last configuration
    vector<long>    fbuffer_id;
    vector<int>     fbuffer_type;
    vector<double> fbuffer_mass,fbuffer_epot;
    vector<double> fbuffer_pos_x,fbuffer_pos_y,fbuffer_pos_z;
    vector<double> fbuffer_vel_x,fbuffer_vel_y,fbuffer_vel_z;


    //******************************************************************************
    //                           Flags
    //******************************************************************************

    // select neighbors as per sample window position
	int xfact = window_x[win_id]; int yfact = window_y[win_id]; int zfact = window_z[win_id];

	// send neighbors ----- odd -even process communications (in z direction)
	int x_send_phase[4] = { 0,xfact,xfact, 0}; int y_send_phase[4] = { 0, 0,yfact,yfact}; int z_send_phase[4] = {zfact,zfact,zfact,zfact};

	// Receive neighbors ----- odd -even process communications (in z direction)
	int x_recv_phase[4] = { 0,-xfact,-xfact, 0}; int y_recv_phase[4] = { 0, 0,-yfact,-yfact}; int z_recv_phase[4] = {-zfact,-zfact,-zfact,-zfact};

    // send neighbors ----- odd-odd / even-even process communications (in x/y direction)
	int x_send_next[3] = {xfact,xfact, 0 }; int y_send_next[3] = { 0,yfact,yfact }; int z_send_next[3] = { 0, 0, 0 };

	// Receive neighbors ----- odd-odd / even-even process communications (in x/y direction)
	int x_recv_next[3] = {-xfact,-xfact, 0 }; int y_recv_next[3] = { 0,-yfact,-yfact }; int z_recv_next[3] = { 0, 0, 0 };

	// rank list as per window id
    int rank_list[7];

    long buff_size[7]; // buffer size corresponding to neighbor CPU rank

    // neighbors to whom should I send updated configuration (from those I requested sphere portions (if any))
    for (int i=0; i<4; i++){
    	rank_list[i] = get_cpu_neighbor(x_send_phase[i],y_send_phase[i],z_send_phase[i]);
    }
    for (int i=4; i<7; i++){
    	rank_list[i] = get_cpu_neighbor(x_send_next[i],y_send_next[i],z_send_next[i]);
    }

    //******************************************************************************
    //                           Buffer allocation
    //******************************************************************************

	// buffer allocation for receive neighbors
	int buf_count = buffer_size * 10; // buffer_size (to be allocated via suitable param file)

	// N1 N2 N4 N6 N3 N5 N7
	double* nb_buffer[7],nb_0,nb_1,nb_2,nb_3,nb_4,nb_5,nb_6; // to be allocated and send

	// allocate memory (to send to neighbors)
	nb_0 = new double [buf_count]; nb_1= new double [buf_count]; nb_2 = new double [buf_count];
	nb_3 = new double [buf_count]; nb_4= new double [buf_count]; nb_5 = new double [buf_count];
	nb_6 = new double [buf_count];

	nb_buffer[7]={nb_0,nb_1,nb_2,nb_3,nb_4,nb_5,nb_6};

	double* my_buffer [7],my_0,my_1,my_2,my_3,my_4,my_5,my_6; // to be received and filled in sphere cells

	// allocate memory (to receive for neighbors)
	my_0 = new double [buf_count]; my_1= new double [buf_count]; my_2 = new double [buf_count];
	my_3 = new double [buf_count]; my_4= new double [buf_count]; my_5 = new double [buf_count];
	my_6 = new double [buf_count];

	my_buffer={my_0,my_1,my_2,my_3,my_4,my_5,my_6};

    //********************************************************************************

	// open and read file only if configuration get accepted

	if(accep_tag){

	        // input file stream instance
            ifstream fin(fname,std::ios_base::in);

            // local value holders
            long f_id; int f_type;
            double f_mass,f_posx,f_posy,f_posz,f_velx,f_vely,f_velz,f_epot;

	        // initialize reference particle
	        vec3d vec1 = pobj.get_myposition();

	        // get reference particle sweep boundary
	        double ref_xmin = vec1.x - (mc_rsweep + mc_sphere_wall);
	        double ref_ymin = vec1.y - (mc_rsweep + mc_sphere_wall);   // minimum boundary
	        double ref_zmin = vec1.z - (mc_rsweep + mc_sphere_wall);

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

                      // set back virtual particles to real particles on sphere boundary
                      fin>>f_type;
                      if(f_type > mc_real_types){ f_type-=mc_real_types; }

                      fbuffer_type.push_back(f_type);
                      fin>>f_mass;      fbuffer_mass.push_back(f_mass);

                      // shift back to native coordinate origin system
                      fin>>f_posx;      fbuffer_pos_x.push_back(f_posx + ref_xmin);
                      fin>>f_posy;      fbuffer_pos_y.push_back(f_posy + ref_ymin);
                      fin>>f_posz;      fbuffer_pos_z.push_back(f_posz + ref_zmin);

                      fin>>f_velx;      fbuffer_vel_x.push_back(f_velx);
                      fin>>f_vely;      fbuffer_vel_y.push_back(f_vely);
                      fin>>f_velz;      fbuffer_vel_z.push_back(f_velz);
                      fin>>f_epot;      fbuffer_epot.push_back(f_epot);
            }
            fin.clear(); // resetting bit states
            fin.close();

            //----------------------------------------------------------------//
            //           2- Updating and shifting part                        //
            //----------------------------------------------------------------//

            // NOTE: later could be changed with n - particles in sphere count

            particle atom;

            long buf_ind, total_particles=0;
            cout<< "Updating particle attributes " << endl;

    		double temp_id; double temp_type;            // temporary holder for particle attributes
    		double temp_mass, temp_epot, dist_check;
    		double temp_vel_x,temp_vel_y,temp_vel_z;

            // temporary position place-holders
            double temp_pos_x,temp_pos_y,temp_pos_z;
            ivec3d cpu_fact,cell_glob_coord,cell_loc_coord, nb_cpu_gcoord;
            int loc_rank , target, cell_index, type_check;

            long particle_id;

            for (int i=0; i< (fbuffer_id.size()-1); i++){ // loop over read data

            	type_check = fbuffer_type.at(i);

            	if(type_check < mc_real_types){ // ignore spherical wall particles

            	    buf_ind = 1;                      // initialize for each buffer step
            	    // attributes assignments
            	    temp_id     = (double) fbuffer_id.at(i);
                    temp_type   = (double) fbuffer_type.at(i);
            	    temp_mass   = fbuffer_mass.at(i);
            	    // position assignments
            	    temp_pos_x  = fbuffer_pos_x.at(i);
            	    temp_pos_y  = fbuffer_pos_y.at(i);
            	    temp_pos_z  = fbuffer_pos_z.at(i);
                    // velocity
            	    temp_vel_x  = fbuffer_vel_x.at(i);
            	    temp_vel_y  = fbuffer_vel_y.at(i);
            	    temp_vel_z  = fbuffer_vel_z.at(i);
                    // potential energy
            	    temp_epot   = fbuffer_epot.at(i);

             	    // global boundary check and position update to actual box dimensions
             	    if (fbuffer_pos_x.at(i) < 0.0)               temp_pos_x = temp_pos_x + mc_simbox_x.x ;
             	    if (fbuffer_pos_x.at(i) > mc_simbox_x.x)     temp_pos_x = temp_pos_x - mc_simbox_x.x ;
             	    if (fbuffer_pos_y.at(i) < 0.0)               temp_pos_y = temp_pos_y + mc_simbox_y.y ;
             	    if (fbuffer_pos_y.at(i) > mc_simbox_y.y)     temp_pos_y = temp_pos_y - mc_simbox_y.y ;
             	    if (fbuffer_pos_z.at(i) < 0.0)               temp_pos_z = temp_pos_z + mc_simbox_z.z ;
             	    if (fbuffer_pos_z.at(i) > mc_simbox_z.z)     temp_pos_z = temp_pos_z - mc_simbox_z.z ;

             	    // global cell coordinate from particle position
             	    cell_glob_coord = cell_coordinate(temp_pos_x, temp_pos_y, temp_pos_z);

         		    // CPU global coordinate from cell global coordinate
         		    nb_cpu_gcoord = calc_cpu_coord(cell_glob_coord) ;
         		    // CPU rank
         		    loc_rank = get_cpu_neighbor(nb_cpu_gcoord.x,nb_cpu_gcoord.y,nb_cpu_gcoord.z);

         		    // particle belongs to neighbor
         		    if (loc_rank != mc_prank){
         		        //******************************************
         		        // if loc_rank != my rank fill it in target buffer
         		        //****************************

         			    // identify buffer target as per rank
         			    for(int ind=0; ind<7;ind++){
         				    if(rank_list[ind] == loc_rank) target=ind;
         			    }

         			    // buffer size counter (start filling from last index +1)
         			    buf_ind += buff_size[target];

                        // ADD TO BUFFER nb_buffer
                        *(nb_buffer+target+buf_ind++)  = temp_id;
                        *(nb_buffer+target+buf_ind++)  = temp_type;
                        *(nb_buffer+target+buf_ind++)  = temp_mass;
                        *(nb_buffer+target+buf_ind++)  = temp_pos_x;
                        *(nb_buffer+target+buf_ind++)  = temp_pos_y;
                        *(nb_buffer+target+buf_ind++)  = temp_pos_z;
                        *(nb_buffer+target+buf_ind++)  = temp_vel_x;
                        *(nb_buffer+target+buf_ind++)  = temp_vel_y;
                        *(nb_buffer+target+buf_ind++)  = temp_vel_z;
                        *(nb_buffer+target+buf_ind++)  = temp_epot;

                        buff_size[target] = buf_ind-1; // check it!!
         		    }

         		    else{ // particles from my CPU

         		        // CPU global position
         		        cpu_fact = get_cpu_gcoord(loc_rank);

         		        // get cell local coordinate
         		        cell_loc_coord = get_cell_loc_coord(cell_glob_coord,cpu_fact);

         		        // to get memory index of cell in cell list (!verify)
         		        cell_index = cell_loc_coord.x + (cell_loc_coord.y * mc_cpu_cell_dim.x) + (cell_loc_coord.z * mc_cpu_cell_dim.y * mc_cpu_cell_dim.x);

         		        particle_id = fbuffer_id.at(i);

         		        // defining particle attributes
         		        atom.set_mynumber(fbuffer_id.at(i));
         		        atom.set_mytype(fbuffer_type.at(i));
         		        atom.set_mymass(fbuffer_mass.at(i));
         		        atom.set_myposition(fbuffer_pos_x.at(i), fbuffer_pos_y.at(i), fbuffer_pos_z.at(i));
         		        atom.set_myvelocity(fbuffer_vel_x.at(i),fbuffer_vel_y.at(i),fbuffer_vel_z.at(i));
         		        atom.set_myepot(fbuffer_epot.at(i));

         		        // update particle attributes in main location
         		        update_particle(bobj,cell_index,atom);

         		    }

                }// spherical wall check
            } // for loop read data

            // assign buffer size to each neighbor CPU

            for (int r_count=0; r_count<7;){
    	      total_particles       = buff_size[r_count] * 0.1 ;  // 10 attributes for each particle
              *(nb_buffer+r_count)  = (double) total_particles;  // each buffer first location has particles count
            }

	} // acceptance condition

	else{ // rejected condition

		for(int j=0; j<7; j++){
			*(nb_buffer+j) = 0.0; // assign all neighbor buffer size to 0
		}
	}

	//****************************************************************
	//                  COMMUNICATION PART
	//                 send / receive part
	//****************************************************************

	// my cpu coordinate
	cpu_gcoord = get_cpu_gcoord(mc_prank);
	int x = cpu_gcoord.x; int y = cpu_gcoord.y; int z = cpu_gcoord.z;

	// z direction communication
	for (int ind=0;ind<4;ind++){

	     if (z%2==0){

	       // front comm
	       MPI_Send(nb_buffer[ind],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_send_phase[ind],y+y_send_phase[ind],z+z_send_phase[ind]),14,comm_name);
	       MPI_Recv(my_buffer[ind],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_recv_phase[ind],y+y_recv_phase[ind],z+z_recv_phase[ind]),15,comm_name,&status);

	     }
	     else{

	       MPI_Recv(my_buffer[ind],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_recv_phase[ind],y+y_recv_phase[ind],z+z_recv_phase[ind]),14,comm_name,&status);
	       MPI_Send(nb_buffer[ind],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_send_phase[ind],y+y_send_phase[ind],z+z_send_phase[ind]),15,comm_name);

	     }
	}

	// Next phase even-even or odd-odd communications

	// x direction communication
	if (x%2 == 0){
	         // right & left comm
	       MPI_Send(nb_buffer[4],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_send_next[0],y+y_send_next[0],z+z_send_next[0]),16,comm_name);
	       MPI_Recv(my_buffer[4],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_recv_next[0],y+y_recv_next[0],z+z_recv_next[0]),17,comm_name,&status);

	         // top-right & bottom-left comm
	       MPI_Send(nb_buffer[5],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_send_next[1],y+y_send_next[1],z+z_send_next[1]),18,comm_name);
	       MPI_Recv(my_buffer[5],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_recv_next[1],y+y_recv_next[1],z+z_recv_next[1]),19,comm_name,&status);
	}
	else{
		   // right & left comm
	       MPI_Recv(my_buffer[4],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_recv_next[0],y+y_recv_next[0],z+z_recv_next[0]),16,comm_name,&status);
	       MPI_Send(nb_buffer[4],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_send_next[0],y+y_send_next[0],z+z_send_next[0]),17,comm_name);

	       // top-right & bottom-left comm
	       MPI_Recv(my_buffer[5],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_recv_next[1],y+y_recv_next[1],z+z_recv_next[1]),18,comm_name,&status);
	       MPI_Send(nb_buffer[5],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_send_next[1],y+y_send_next[1],z+z_send_next[1]),19,comm_name);

	}

	// y direction communication

	if (y%2 == 0){
	         // top & bottom comm
	       MPI_Send(nb_buffer[6],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_send_next[2],y+y_send_next[2],z+z_send_next[2]),20,comm_name);
	       MPI_Recv(my_buffer[6],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_recv_next[2],y+y_recv_next[2],z+z_recv_next[2]),21,comm_name,&status);
	}
	else{

	       MPI_Recv(my_buffer[6],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_recv_next[2],y+y_recv_next[2],z+z_recv_next[2]),20,comm_name,&status);
	       MPI_Send(nb_buffer[6],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_send_next[2],y+y_send_next[2],z+z_send_next[2]),21,comm_name);
	}

	// synchronize communication
	MPI_Barrier(comm_name);


	//***********************************************************************
	//                if I receive updated configuration - effect changes
	//***********************************************************************

	long mybuf_ind,rec_particles=0;
	double temp_id; double temp_type;            // temporary holder for particle attributes
	double temp_mass, temp_epot, dist_check;
	double temp_vel_x,temp_vel_y,temp_vel_z;

    // temporary position place-holders
    double temp_pos_x,temp_pos_y,temp_pos_z;
    ivec3d cpu_fact,cell_glob_coord,cell_loc_coord, nb_cpu_gcoord;
    int loc_rank,cell_index;

    particle atom;

	// loop over my_received buffers and add particles to sphere cell
	for(int n_count=0;n_count<7;n_count++){

		rec_particles = (long) *(my_buffer+n_count);
		mybuf_ind =1;

        if(rec_particles !=0 ){

            atom.set_mynumber((long) *(my_buffer+mybuf_ind++));
            atom.set_mytype((int) *(my_buffer+mybuf_ind++));
            atom.set_mymass((double) *(my_buffer+mybuf_ind++));

            // shifting origin of coordinate system
            temp_pos_x = (double) *(my_buffer+mybuf_ind++);
            temp_pos_y = (double) *(my_buffer+mybuf_ind++);
            temp_pos_z = (double) *(my_buffer+mybuf_ind++);

            atom.set_myposition((double) *(my_buffer+mybuf_ind++),(double) *(my_buffer+mybuf_ind++),(double) *(my_buffer+mybuf_ind++));
            atom.set_myvelocity((double) *(my_buffer+mybuf_ind++),(double) *(my_buffer+mybuf_ind++),(double) *(my_buffer+mybuf_ind++));
            atom.set_myepot((double) *(my_buffer+mybuf_ind++));


        	// global cell coordinate from particle position
     	    cell_glob_coord = cell_coordinate(temp_pos_x, temp_pos_y, temp_pos_z);

 		    // CPU global coordinate from cell global coordinate
 		    nb_cpu_gcoord = calc_cpu_coord(cell_glob_coord) ;

 		    // CPU rank
 		    loc_rank = get_cpu_neighbor(nb_cpu_gcoord.x,nb_cpu_gcoord.y,nb_cpu_gcoord.z);


		    // to get memory index of cell in cell list (!verify)
		    cell_index = cell_loc_coord.x + (cell_loc_coord.y * mc_cpu_cell_dim.x) + (cell_loc_coord.z * mc_cpu_cell_dim.y * mc_cpu_cell_dim.x);

		    // update particle attributes in main location
		    update_particle(bobj,cell_index,atom);

        }

	}

	//********************************************************************
	//                  deallocate buffers
	//********************************************************************

    // DELETE BUFFER MEMORY ACCORDINGLY
    delete[] nb_0; delete[] nb_1;delete[] nb_2; delete[] nb_3;
    delete[] nb_4;delete[] nb_5;delete[] nb_6;

    delete[] my_0; delete[] my_1;delete[] my_2; delete[] my_3;
    delete[] my_4;delete[] my_5;delete[] my_6;


}

void construct_sphere(particle pobj, cellblock bobj, int win_id,char* filename){


	    // select neighbors as per sample window position
		int xfact = window_x[win_id]; int yfact = window_y[win_id]; int zfact = window_z[win_id];

		// send neighbors ----- odd -even process communications (in z direction)
		int x_send_phase[4] = { 0,xfact,xfact, 0}; int y_send_phase[4] = { 0, 0,yfact,yfact}; int z_send_phase[4] = {zfact,zfact,zfact,zfact};

		// Receive neighbors ----- odd -even process communications (in z direction)
		int x_recv_phase[4] = { 0,-xfact,-xfact, 0}; int y_recv_phase[4] = { 0, 0,-yfact,-yfact}; int z_recv_phase[4] = {-zfact,-zfact,-zfact,-zfact};

	    // send neighbors ----- odd-odd / even-even process communications (in x/y direction)
		int x_send_next[3] = {xfact,xfact, 0 }; int y_send_next[3] = { 0,yfact,yfact }; int z_send_next[3] = { 0, 0, 0 };

		// Receive neighbors ----- odd-odd / even-even process communications (in x/y direction)
		int x_recv_next[3] = {-xfact,-xfact, 0 }; int y_recv_next[3] = { 0,-yfact,-yfact }; int z_recv_next[3] = { 0, 0, 0 };


		// my cpu coordinate
		cpu_gcoord = get_cpu_gcoord(mc_prank);
		int x = cpu_gcoord.x; int y = cpu_gcoord.y; int z = cpu_gcoord.z;

		// randomly selected particle positions (my-my cpu and rec - neighbor cpu)
		double my_pos[3], rec_pos_0[3],rec_pos_1[3],rec_pos_2[3],rec_pos_3[3],rec_pos_4[3],rec_pos_5[3],rec_pos_6[3],rec_pos_7[3];

		double* rec_pos[8]={rec_pos_0,rec_pos_1,rec_pos_2,rec_pos_3,rec_pos_4,rec_pos_5,rec_pos_6,rec_pos_7};

		my_pos[0]=pobj.get_myposition().x; my_pos[1]=pobj.get_myposition().y; my_pos[2]=pobj.get_myposition().z;


		//***************************************************************************
		// communication part -- send/request neighbors with chosen particle
		//***************************************************************************


		// z direction communication
		for (int ind=0;ind<4;ind++){

		     if (z%2==0){

		       // front comm
		       MPI_Send(my_pos,3,MPI_DOUBLE,get_cpu_neighbor(x+x_send_phase[ind],y+y_send_phase[ind],z+z_send_phase[ind]),0,comm_name);
		       MPI_Recv(rec_pos[ind],3,MPI_DOUBLE,get_cpu_neighbor(x+x_recv_phase[ind],y+y_recv_phase[ind],z+z_recv_phase[ind]),1,comm_name,&status);

		     }
		     else{

		       MPI_Recv(rec_pos[ind],3,MPI_DOUBLE,get_cpu_neighbor(x+x_recv_phase[ind],y+y_recv_phase[ind],z+z_recv_phase[ind]),0,comm_name,&status);
		       MPI_Send(my_pos,3,MPI_DOUBLE,get_cpu_neighbor(x+x_send_phase[ind],y+y_send_phase[ind],z+z_send_phase[ind]),1,comm_name);

		     }
		}

		// Next phase even-even or odd-odd communications

		// x direction communication
		if (x%2 == 0){
		         // right & left comm
		       MPI_Send(my_pos,3,MPI_DOUBLE,get_cpu_neighbor(x+x_send_next[0],y+y_send_next[0],z+z_send_next[0]),2,comm_name);
		       MPI_Recv(rec_pos[4],3,MPI_DOUBLE,get_cpu_neighbor(x+x_recv_next[0],y+y_recv_next[0],z+z_recv_next[0]),3,comm_name,&status);

		         // top-right & bottom-left comm
		       MPI_Send(my_pos,3,MPI_DOUBLE,get_cpu_neighbor(x+x_send_next[1],y+y_send_next[1],z+z_send_next[1]),4,comm_name);
		       MPI_Recv(rec_pos[5],3,MPI_DOUBLE,get_cpu_neighbor(x+x_recv_next[1],y+y_recv_next[1],z+z_recv_next[1]),5,comm_name,&status);
		}
		else{
			   // right & left comm
		       MPI_Recv(rec_pos[4],3,MPI_DOUBLE,get_cpu_neighbor(x+x_recv_next[0],y+y_recv_next[0],z+z_recv_next[0]),2,comm_name,&status);
		       MPI_Send(my_pos,3,MPI_DOUBLE,get_cpu_neighbor(x+x_send_next[0],y+y_send_next[0],z+z_send_next[0]),3,comm_name);

		       // top-right & bottom-left comm
		       MPI_Recv(rec_pos[5],3,MPI_DOUBLE,get_cpu_neighbor(x+x_recv_next[1],y+y_recv_next[1],z+z_recv_next[1]),4,comm_name,&status);
		       MPI_Send(my_pos,3,MPI_DOUBLE,get_cpu_neighbor(x+x_send_next[1],y+y_send_next[1],z+z_send_next[1]),5,comm_name);

		}

		// y direction communication

		if (y%2 == 0){
		         // top & bottom comm
		       MPI_Send(my_pos,3,MPI_DOUBLE,get_cpu_neighbor(x+x_send_next[2],y+y_send_next[2],z+z_send_next[2]),4,comm_name);
		       MPI_Recv(rec_pos[6],3,MPI_DOUBLE,get_cpu_neighbor(x+x_recv_next[2],y+y_recv_next[2],z+z_recv_next[2]),5,comm_name,&status);
		}
		else{

		       MPI_Recv(rec_pos[6],3,MPI_DOUBLE,get_cpu_neighbor(x+x_recv_next[2],y+y_recv_next[2],z+z_recv_next[2]),4,comm_name,&status);
		       MPI_Send(my_pos,3,MPI_DOUBLE,get_cpu_neighbor(x+x_send_next[2],y+y_send_next[2],z+z_send_next[2]),5,comm_name);
		}


		// synchronize communication
		MPI_Barrier(comm_name);

		// do sweep test to construct buffers


	    //******************************************************************
	    // NOTE: HARDCODED FOR 8 CELLS IN A SAMPLER WINDOW (OWN CPU)
        // NEED TO BE GENERALIZED LATER DURING BIAS APPROACH

		// buffer allocation for receive neighbors
		int buf_count = buffer_size * 10; // buffer_size (to be allocated via suitable param file)

		// N1 N2 N4 N6 N3 N5 N7
		double* nb_buffer[7],nb_0,nb_1,nb_2,nb_3,nb_4,nb_5,nb_6; // to be allocated and send

		// allocate memory (to send to neighbors)
		nb_0 = new double [buf_count]; nb_1= new double [buf_count]; nb_2 = new double [buf_count];
		nb_3 = new double [buf_count]; nb_4= new double [buf_count]; nb_5 = new double [buf_count];
		nb_6 = new double [buf_count];

		nb_buffer[7]={nb_0,nb_1,nb_2,nb_3,nb_4,nb_5,nb_6};

		double* my_buffer [7],my_0,my_1,my_2,my_3,my_4,my_5,my_6; // to be received and filled in sphere cells

		// allocate memory (to receive for neighbors)
		my_0 = new double [buf_count]; my_1= new double [buf_count]; my_2 = new double [buf_count];
		my_3 = new double [buf_count]; my_4= new double [buf_count]; my_5 = new double [buf_count];
		my_6 = new double [buf_count];

		my_buffer={my_0,my_1,my_2,my_3,my_4,my_5,my_6};

		// neighbor cpu global coordinate
		double nb_gcoord [7][3]; // [neighbor] [gcoord_x gcoord_y gcoord_z]

		// flags for coordinate shifting
        int flag [7][3]; // [neighbor] [x y z]

        // global cpu array dimension range
      	ivec3d cpu_max = {mc_cpu_dim.x-1,mc_cpu_dim.y-1,mc_cpu_dim.z-1};            // dimension - Max
        ivec3d cpu_min = {0,0,0};                                                   // dimension - Min

        // boundary checks and flag assignments (ENSURE COORDINATES SHIFT IF REQUIRED)

        // z check  with  N1
        if( (nb_gcoord[0][2] == cpu_min.z) && (window_position[win_id].z==0) ) { flag[0][2] = -1;}
        if( (nb_gcoord[0][2] == cpu_max.z) && (window_position[win_id].z==1) ) { flag[0][2] = +1;}

        // x-z check with N2
        if( (nb_gcoord[1][0] == cpu_min.x) && (window_position[win_id].x==0) ) { flag[1][0] = -1;}
        if( (nb_gcoord[1][0] == cpu_max.x) && (window_position[win_id].x==1) ) { flag[1][0] = +1;}
        if( (nb_gcoord[1][2] == cpu_min.z) && (window_position[win_id].z==0) ) { flag[1][2] = -1;}
        if( (nb_gcoord[1][2] == cpu_max.z) && (window_position[win_id].z==1) ) { flag[1][2] = +1;}

        // x-y-z check with N4
        if( (nb_gcoord[2][0] == cpu_min.x) && (window_position[win_id].x==0) ) { flag[2][0] = -1;}
        if( (nb_gcoord[2][0] == cpu_max.x) && (window_position[win_id].x==1) ) { flag[2][0] = +1;}
        if( (nb_gcoord[2][1] == cpu_min.y) && (window_position[win_id].y==0) ) { flag[2][1] = -1;}
        if( (nb_gcoord[2][1] == cpu_max.y) && (window_position[win_id].y==1) ) { flag[2][1] = +1;}
        if( (nb_gcoord[2][2] == cpu_min.z) && (window_position[win_id].z==0) ) { flag[2][2] = -1;}
        if( (nb_gcoord[2][2] == cpu_max.z) && (window_position[win_id].z==1) ) { flag[2][2] = +1;}

        // y-z check with N6
        if( (nb_gcoord[3][1] == cpu_min.y) && (window_position[win_id].y==0) ) { flag[3][1] = -1;}
        if( (nb_gcoord[3][1] == cpu_max.y) && (window_position[win_id].y==1) ) { flag[3][1] = +1;}
        if( (nb_gcoord[3][2] == cpu_min.z) && (window_position[win_id].z==0) ) { flag[3][2] = -1;}
        if( (nb_gcoord[3][2] == cpu_max.z) && (window_position[win_id].z==1) ) { flag[3][2] = +1;}
        // x check with N3
        if( (nb_gcoord[4][0] == cpu_min.x) && (window_position[win_id].x==0) ) { flag[4][0] = -1;}
        if( (nb_gcoord[4][0] == cpu_max.x) && (window_position[win_id].x==1) ) { flag[4][0] = +1;}

        // x-y check with N5
        if( (nb_gcoord[5][0] == cpu_min.x) && (window_position[win_id].x==0) ) { flag[5][0] = -1;}
        if( (nb_gcoord[5][0] == cpu_max.x) && (window_position[win_id].x==1) ) { flag[5][0] = +1;}
        if( (nb_gcoord[5][1] == cpu_min.y) && (window_position[win_id].y==0) ) { flag[5][1] = -1;}
        if( (nb_gcoord[5][1] == cpu_max.y) && (window_position[win_id].y==1) ) { flag[5][1] = +1;}

        // y check with N7
        if( (nb_gcoord[6][1] == cpu_min.y) && (window_position[win_id].y==0) ) { flag[6][1] = -1;}
        if( (nb_gcoord[6][1] == cpu_max.y) && (window_position[win_id].y==1) ) { flag[6][1] = +1;}

		// local cell coordinates zone
	    ivec3d test, nb_test[19];

	    // window sampling cells
	    celltype sample_cells[8];

	    //window neighbor cells
        // no of each neighbor cells ( N1 N2 N4 N6 N3 N5 N7 )
        int neig_cells_size [7]={4,2,1,2,4,2,4};

        // list of neighbor receive cells
        celltype neighb[19];

        // variables for neighbor cell ranges
        int nb_xmin, nb_ymin, nb_zmin, nb_xmax, nb_ymax, nb_zmax, nb_x, nb_y, nb_z;

	    // assign neighbor cell ranges
	    nb_xmin = zone_limit[win_id].xmin; nb_xmax = zone_limit[win_id].xmax;
	    nb_ymin = zone_limit[win_id].ymin; nb_ymax = zone_limit[win_id].ymax;
	    nb_zmin = zone_limit[win_id].zmin; nb_zmax = zone_limit[win_id].zmax;

	    if(nb_xmax ==1) {nb_x = 2;}
	    else{ nb_x=0; }

	    if(nb_ymax ==1) {nb_y = 2;}
	    else{ nb_y=0; }

	    if(nb_zmax ==1) {nb_z = 2;}
	    else{ nb_z=0; }

        // prepare NEIGHBOR cell list as per sample window id
	    //Neighbor-1 (N1)
	    nb_test[0].x =nb_xmin; nb_test[0].y = nb_ymin; nb_test[0].z = nb_z; neighb[0]=bobj.cell_with_lcoord(nb_test[0]);
	    nb_test[1].x =nb_xmax; nb_test[1].y = nb_ymin; nb_test[1].z = nb_z; neighb[1]=bobj.cell_with_lcoord(nb_test[1]);
	    nb_test[2].x =nb_xmin; nb_test[2].y = nb_ymax; nb_test[2].z = nb_z; neighb[2]=bobj.cell_with_lcoord(nb_test[2]);
	    nb_test[3].x =nb_xmax; nb_test[3].y = nb_ymax; nb_test[3].z = nb_z; neighb[3]=bobj.cell_with_lcoord(nb_test[3]);

        //Neighbor-2 (N2)
        nb_test[4].x = nb_x; nb_test[4].y = nb_ymin; nb_test[4].z = nb_z; neighb[4] = bobj.cell_with_lcoord(nb_test[4]);
        nb_test[5].x = nb_x; nb_test[5].y = nb_ymax; nb_test[5].z = nb_z; neighb[5] = bobj.cell_with_lcoord(nb_test[5]);

        //Neighbor-3 (N4)
        nb_test[6].x = nb_x; nb_test[6].y = nb_y; nb_test[6].z = nb_z;    neighb[6] = bobj.cell_with_lcoord(nb_test[6]);

        //Neighbor-4 (N6)
        nb_test[7].x = nb_xmin; nb_test[7].y = nb_y; nb_test[7].z = nb_z; neighb[7] = bobj.cell_with_lcoord(nb_test[7]);
        nb_test[8].x = nb_xmax; nb_test[8].y = nb_y; nb_test[8].z = nb_z; neighb[8] = bobj.cell_with_lcoord(nb_test[8]);

        //Neighbor-5 (N3)
	    nb_test[9].x =nb_x;   nb_test[9].y = nb_ymin; nb_test[9].z = nb_zmin;  neighb[9]=bobj.cell_with_lcoord(nb_test[9]);
	    nb_test[10].x =nb_x; nb_test[10].y = nb_ymin; nb_test[10].z = nb_zmax; neighb[10]=bobj.cell_with_lcoord(nb_test[10]);
	    nb_test[11].x =nb_x; nb_test[11].y = nb_ymax; nb_test[11].z = nb_zmin; neighb[11]=bobj.cell_with_lcoord(nb_test[11]);
	    nb_test[12].x =nb_x; nb_test[12].y = nb_ymax; nb_test[12].z = nb_zmax; neighb[12]=bobj.cell_with_lcoord(nb_test[12]);

        //Neighbor-6 (N5)
        nb_test[13].x = nb_x; nb_test[13].y = nb_y; nb_test[13].z = nb_zmin; neighb[13] = bobj.cell_with_lcoord(nb_test[13]);
        nb_test[14].x = nb_x; nb_test[14].y = nb_y; nb_test[14].z = nb_zmax; neighb[14] = bobj.cell_with_lcoord(nb_test[14]);

	    //Neighbor-7 (N7)
	    nb_test[15].x =nb_xmin; nb_test[15].y = nb_y; nb_test[15].z = nb_zmin; neighb[15]=bobj.cell_with_lcoord(nb_test[15]);
	    nb_test[16].x =nb_xmin; nb_test[16].y = nb_y; nb_test[16].z = nb_zmax; neighb[16]=bobj.cell_with_lcoord(nb_test[16]);
	    nb_test[17].x =nb_xmax; nb_test[17].y = nb_y; nb_test[17].z = nb_zmin; neighb[17]=bobj.cell_with_lcoord(nb_test[17]);
	    nb_test[18].x =nb_xmax; nb_test[18].y = nb_y; nb_test[18].z = nb_zmax; neighb[18]=bobj.cell_with_lcoord(nb_test[18]);

	    //***********************************************************
	    //                  My received neighbor part - fill buffers
	    //**********************************************************

	    // sweep distance computation (with squares)
	    double r_full   = pow((mc_rsweep + mc_sphere_wall),2.0);
	    double r_sphere = pow(mc_rsweep,2.0);

	    vec3d ref_pos, curr_pos;                       // reference/current particle position

		double temp_id; double temp_type;            // temporary holder for particle attributes
		double temp_mass, temp_epot, dist_check;
		vec3d temp_pos, temp_vel;

		long total_particles=0;                       // counter for total particle no in buffers
	    int ind =0;                                   // neighbor cell counter
	    long buf_counter,part_count;

	    for(int j=0; j<7;j++ ){ // loop over all neighbor cells

	    	    buf_counter=1;

	    	    // get receive neighbor chosen particle positions
	    	    ref_pos.x = *rec_pos(j+0); ref_pos.y = *rec_pos(j+1); ref_pos.z = *rec_pos(j+2);

	    	    for (int k=0; k<neig_cells_size[j]; k++){ // loop over each neighbor sub total

	    	         part_count = neighb[ind].get_nparticles(); //get total particles



	    	         for (long val;val<part_count;val++){ // loop over particles in each cell

		    		          curr_pos =neighb[ind].get_particle(val).get_myposition();

		    		          // coordinate shifting as per flag initialization
		    		          curr_pos.x = curr_pos.x + (flag[j][0] * mc_simbox_x.x);
		    		          curr_pos.y = curr_pos.y + (flag[j][1] * mc_simbox_y.y);
		    		          curr_pos.z = curr_pos.z + (flag[j][2] * mc_simbox_z.z);

		    		          dist_check = distance_vect(ref_pos,curr_pos);
		    		          // ----- cutoff check
		    		          //
		    		          if(dist_check <= r_full ){

     		    		            // declare virtual particles on sphere boundary
	     	    		            if( dist_check > r_sphere ){

	     	    			             temp_type = (double) neighb[ind].get_particle(val).get_mytype();

	     	    			             // ignore placeholders on sphere wall
		    			                 if (temp_type != (double) 2 ) { // HC: now hardcoded for placeholders

		    			        	           temp_id     = (double) (neighb[ind].get_particle(val).get_mynumber());
		    			        	           temp_type   = (double) (temp_type + mc_real_types);  // HC: hardcoded for ntypes=3
	                                           temp_mass   =  neighb[ind].get_particle(val).get_mymass();

	                                           // assign temporary position (coordinates shifted if neighbors on cpu boundary)
	                                           temp_pos  =  curr_pos;
	                                           temp_vel  =  neighb[ind].get_particle(val).get_myvelocity();
	                                           temp_epot =  neighb[ind].get_particle(val).get_myepot();

	                                           // ADD TO BUFFER nb_buffer
	                                           *(nb_buffer+j+buf_counter++)  = temp_id;
	                                           *(nb_buffer+j+buf_counter++)  = temp_type;
	                                           *(nb_buffer+j+buf_counter++)  = temp_mass;
	                                           *(nb_buffer+j+buf_counter++)  = temp_pos.x;
	                                           *(nb_buffer+j+buf_counter++)  = temp_pos.y;
	                                           *(nb_buffer+j+buf_counter++)  = temp_pos.z;
	                                           *(nb_buffer+j+buf_counter++)  = temp_vel.x;
	                                           *(nb_buffer+j+buf_counter++)  = temp_vel.y;
	                                           *(nb_buffer+j+buf_counter++)  = temp_vel.z;
	                                           *(nb_buffer+j+buf_counter++)  = temp_epot;
		    			                 }
		    		                }
		    		                else{ // include all particles inside mc_rsweep


	 			        	            temp_id     = (double) (neighb[ind].get_particle(val).get_mynumber());
	 			        	            temp_type   = (double) (neighb[ind].get_particle(val).get_mytype());
	                                    temp_mass   =  neighb[ind].get_particle(val).get_mymass();
	                                    // assign temporary position (coordinates shifted if neighbors on cpu boundary)
	                                    temp_pos  =  curr_pos;
	                                    temp_vel  =  neighb[ind].get_particle(val).get_myvelocity();
	                                    temp_epot =  neighb[ind].get_particle(val).get_myepot();

	                                    // ADD TO BUFFER nb_nb_buffer_0
	                                    *(nb_buffer+j+buf_counter++) = temp_id;
	                                    *(nb_buffer+j+buf_counter++) = temp_type;
	                                    *(nb_buffer+j+buf_counter++) = temp_mass;
	                                    *(nb_buffer+j+buf_counter++) = temp_pos.x;
	                                    *(nb_buffer+j+buf_counter++) = temp_pos.y;
	                                    *(nb_buffer+j+buf_counter++) = temp_pos.z;
	                                    *(nb_buffer+j+buf_counter++) = temp_vel.x;
	                                    *(nb_buffer+j+buf_counter++) = temp_vel.y;
	                                    *(nb_buffer+j+buf_counter++) = temp_vel.z;
	                                    *(nb_buffer+j+buf_counter++) = temp_epot;

		             	            }

		                      }  // cut-off check


		    	     } //  end of particles loop

		    	     ind++; // neighbor cells counter

	    	    } // loop over each neighbor sub total
	    	    // no of particles in each neighbor buffer
	    	    total_particles = (buf_counter-1) * 0.1 ;     //10 attributes for each particle
                *(nb_buffer+j)  = (double) total_particles;  // each buffer first location has particles count

		} // end of Neighbor cpu loop

	    //****************************************************************
	    //          send/receive buffers to/from corresponding neighbors
	    //****************************************************************


		// synchronize communication
		MPI_Barrier(comm_name);


	    // z direction

	    for (int ind=0;ind<4;ind++){

	        if (z%2==0){

	          // front comm
	          MPI_Send(nb_buffer[ind],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_recv_phase[ind],y+y_recv_phase[ind],z+z_recv_phase[ind]),7,comm_name);
	          MPI_Recv(my_buffer[ind],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_send_phase[ind],y+y_send_phase[ind],z+z_send_phase[ind]),6,comm_name,&status);
	        }
	        else{

	          MPI_Recv(my_buffer[ind],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_send_phase[ind],y+y_send_phase[ind],z+z_send_phase[ind]),7,comm_name,&status);
	          MPI_Send(nb_buffer[ind],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_recv_phase[ind],y+y_recv_phase[ind],z+z_recv_phase[ind]),6,comm_name);
	        }
	    }

	    // Next phase even-even or odd-odd communications

	    // x direction communication
	    if (x%2 == 0){
	            // right & left comm
	          MPI_Send(nb_buffer[4],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_recv_next[0],y+y_recv_next[0],z+z_recv_next[0]),9,comm_name);
	          MPI_Recv(my_buffer[4],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_send_next[0],y+y_send_next[0],z+z_send_next[0]),8,comm_name,&status);

	          // top-right & bottom-left comm
	          MPI_Send(nb_buffer[5],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_recv_next[1],y+y_recv_next[1],z+z_recv_next[1]),11,comm_name);
	          MPI_Recv(my_buffer[5],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_send_next[1],y+y_send_next[1],z+z_send_next[1]),10,comm_name,&status);
	    }
	    else{

	          MPI_Recv(my_buffer[4],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_send_next[0],y+y_send_next[0],z+z_send_next[0]),9,comm_name,&status);
	          MPI_Send(nb_buffer[4],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_recv_next[0],y+y_recv_next[0],z+z_recv_next[0]),8,comm_name);

	          // top-right & bottom-left comm
	          MPI_Recv(my_buffer[5],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_send_next[1],y+y_send_next[1],z+z_send_next[1]),11,comm_name,&status);
	          MPI_Send(nb_buffer[5],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_recv_next[1],y+y_recv_next[1],z+z_recv_next[1]),10,comm_name);
	    }

	    // y direction communication
	    if (y%2 == 0){
	    	  // top & bottom comm
	          MPI_Send(nb_buffer[6],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_recv_next[2],y+y_recv_next[2],z+z_recv_next[2]),13,comm_name);
	          MPI_Recv(my_buffer[6],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_send_next[2],y+y_send_next[2],z+z_send_next[2]),12,comm_name,&status);
	    }
	    else{
	          MPI_Recv(my_buffer[6],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_send_next[2],y+y_send_next[2],z+z_send_next[2]),13,comm_name,&status);
	          MPI_Send(nb_buffer[6],buf_count,MPI_DOUBLE,get_cpu_neighbor(x+x_recv_next[2],y+y_recv_next[2],z+z_recv_next[2]),12,comm_name);
	    }


		// synchronize communication
		MPI_Barrier(comm_name);



	    // **************************************************************************
	    // allocate received buffers if any
	    // analyse received buffers and fill into sphere cells if necessary
	    // **************************************************************************


        celltype sphere_cell;    // contains particles in spherical domain for simulation
        particle atom,my_atom;   // particle object

	    int count=0;

	    // prepare SAMPLE cell list as per sample window id

	    for(int i=zone_limit[win_id].xmin;i<=zone_limit[win_id].xmax;i++){
	    	for(int j=zone_limit[win_id].ymin;j<=zone_limit[win_id].ymax;j++){
	    		for(int k=zone_limit[win_id].zmin;k<=zone_limit[win_id].zmax;k++){

	    			test.x = i; test.y = j; test.z = k;
	    			sample_cells[count] = bobj.cell_with_lcoord(test);

	    			// add to sample factor
	    			bobj.cell_with_lcoord(test).add_sample();

	                count++;
	    		}// k loop
	    	}// j loop
	    }// i loop


	    // Reference particle position
	    vec3d ref_position = pobj.get_myposition(); // my own selected particle

	    vec3d my_curr; // my current position

	    // loop over my sample zone and add to sphere cells with in sweep
	    // get reference particle sweep boundary
	    double ref_xmin = ref_position.x - (mc_rsweep + mc_sphere_wall);
	    double ref_ymin = ref_position.y - (mc_rsweep + mc_sphere_wall);   // minimum boundary
	    double ref_zmin = ref_position.z - (mc_rsweep + mc_sphere_wall);

	    //*************************************************************
	    //                   My neighbor part -fill sphere cells
	    //*************************************************************


	    // loop over my_received buffers and add particles to sphere cell
        // shift to origin

	    long p_count, buf_ind;
	    double temp_position_x,temp_position_y,temp_position_z;

        for(int j=0; j<7; j++){

        	p_count = (long) *(my_buffer+j); // no of particles in buffer
        	buf_ind = 1;                      // initialize for each buffer list

        	if(p_count !=0){
             	for (long i=0; i< p_count; i++){

                      atom.set_mynumber((long) *(my_buffer+j+buf_ind++));
                      atom.set_mytype((int) *(my_buffer+j+buf_ind++));
                      atom.set_mymass((double) *(my_buffer+j+buf_ind++));

                      // shifting origin of coordinate system
                      temp_position_x = (double) *(my_buffer+j+buf_ind++) - ref_xmin;
                      temp_position_y = (double) *(my_buffer+j+buf_ind++) - ref_ymin;
                      temp_position_z = (double) *(my_buffer+j+buf_ind++) - ref_zmin;

                      atom.set_myposition(temp_position_x,temp_position_y,temp_position_z);
                      atom.set_myvelocity((double) *(my_buffer+j+buf_ind++),(double) *(my_buffer+j+buf_ind++),(double) *(my_buffer+j+buf_ind++));
                      atom.set_myepot((double) *(my_buffer+j+buf_ind++));

                      // appending particle instance to sphere cell list
                      sphere_cell.add_particle(atom);
        	    }
        	}
        }

	    //*************************************************************
	    //                   My particle part -fill sphere cells
	    //*************************************************************


	    // loop over sample_cells -- perform distance check and add particles to sphere_cell
	    // shift to origin
        int my_temp_type;


        for(int j=0; j<8; j++){ // loop over sample cells

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
		                	 my_temp_type = my_temp_type + mc_real_types ;

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

        // flip its type to real Carbon
        my_atom.set_mytype(1);

        // shifting origin of coordinate system
        temp_position_x = my_atom.get_myposition().x - ref_xmin;
        temp_position_y = my_atom.get_myposition().y - ref_ymin;
        temp_position_z = my_atom.get_myposition().z - ref_zmin;

        my_atom.set_myposition(temp_position_x,temp_position_y,temp_position_z);

        sphere_cell.add_particle(my_atom);

        // *********************************************************************
        //                          file writing methods
        // *********************************************************************

        long file_number;
        int file_type;
        double file_mass, file_epot;
        vec3d file_pos, file_vel;

        // opening file stream object
        ofstream fout(filename, ios_base::out);

        cout<<" Sphere constructor writing to file : " << filename<<endl;

        //********** loop over revised sphere cell list - 2 **************
        long fp_total = sphere_cell.get_nparticles();

        for(long fp=0; fp<fp_total; fp++){

            	file_number  = sphere_cell.get_particle(fp).get_mynumber();
            	file_type    = sphere_cell.get_particle(fp).get_mytype();
            	file_mass    = sphere_cell.get_particle(fp).get_mymass();
            	file_pos     = sphere_cell.get_particle(fp).get_myposition();
            	file_vel     = sphere_cell.get_particle(fp).get_myvelocity(); // optional with #ifdef or so
            	file_epot    = sphere_cell.get_particle(fp).get_myepot();

                	  // file flush
                fout << setw(8) << file_number
                << setw(6) << file_type
                << setw(6) << setprecision(8) << file_mass
                << setw(6) << setprecision(8) << file_pos.x
                << setw(6) << setprecision(8) << file_pos.y
                << setw(6) << setprecision(8) << file_pos.z
                << setw(6) << setprecision(8) << file_vel.x
                << setw(6) << setprecision(8) << file_vel.y
                << setw(6) << setprecision(8) << file_vel.z
                << setw(6) << setprecision(8) << file_epot
                << endl;

        }

        fout.close(); // closing outfile connection

	    // DELETE BUFFER MEMORY ACCORDINGLY
	    delete[] nb_0; delete[] nb_1;delete[] nb_2; delete[] nb_3;
	    delete[] nb_4;delete[] nb_5;delete[] nb_6;

	    delete[] my_0; delete[] my_1;delete[] my_2; delete[] my_3;
	    delete[] my_4;delete[] my_5;delete[] my_6;

}

void make_mc_nblist(celltype cobj){

	//: DNOTE: later additional signature to identify zone and nbl will be constructed accordingly

	 // create neighbor list for the given cell
	 ivec3d g_max = mc_global_cell_dim;    // global cell array dimension - Max
	 ivec3d g_min = {0,0,0};               // global cell array dimension - Min

     //**************************************
	 // NOTE ** Check whether g_min = {1,1,1};
	 //****************************************

	 // my cell global coordinate
	 ivec3d ref = cobj.get_cell_glob_coord();

	 // neighbor vectors
	 ivec3d vec1,vec2,vec3,vec4,vec5,vec6,vec7,vec8,vec9;

	 int nb_con =0; // neighbor counter

	 int x_plus=0, y_plus=0, x_minus=0, y_minus=0;
	 int z_val=0;        // z value holder
	 int z_dir[3] = {0,-1,+1};

	 y_minus = --ref.y; if(ref.y == g_min.y) { y_minus = g_max.y;}
	 y_plus  = ++ref.y; if(ref.y == g_max.y) { y_plus = g_min.y;}

	 x_minus = --ref.x; if(ref.x == g_min.x) { x_minus = g_max.x;}
	 x_plus  = ++ref.x; if(ref.x == g_max.x) { x_plus = g_min.x;}

	 // loop over three possible z-directions
	 for (int i=0; i<3; i++){

		 z_val = ref.z+z_dir[i];
		 if((i==1) && (ref.z==g_min.z)) {z_val = g_max.z; }
		 if((i==2) && (ref.z==g_max.z)) {z_val = g_min.z; }

		 vec1.x = ref.x;      vec1.y  = ref.y;      vec1.z = z_val;     cobj.add_neighbor(nb_con++,vec1); // position
		 vec2.x = x_minus;    vec2.y  = ref.y;      vec2.z = z_val;     cobj.add_neighbor(nb_con++,vec2); // left
		 vec3.x = x_plus;     vec3.y  = ref.y;      vec3.z = z_val;     cobj.add_neighbor(nb_con++,vec3); // right
		 vec4.x = ref.x;      vec4.y  = y_plus;     vec4.z = z_val;     cobj.add_neighbor(nb_con++,vec4); // front
		 vec5.x = ref.x;      vec5.y  = y_minus;    vec5.z = z_val;     cobj.add_neighbor(nb_con++,vec5); // back
		 vec6.x  = x_plus;    vec6.y  = y_minus;    vec6.z = z_val;     cobj.add_neighbor(nb_con++,vec6); // east
		 vec7.x  = x_minus;   vec7.y  = y_plus;     vec7.z = z_val;     cobj.add_neighbor(nb_con++,vec7); // west
		 vec8.x = x_plus;     vec8.y  = y_plus;     vec8.z = z_val;     cobj.add_neighbor(nb_con++,vec8); // north
		 vec9.x = x_minus;    vec9.y  = y_minus;    vec9.z = z_val;     cobj.add_neighbor(nb_con++,vec9); // south
	 }

}

ivec3d get_cpu_gcoord(int myrank){

	// NEED CHANGES (can be computed from MPI Virtual topology)

	// compute cpu global coordinate based on process rank

	ivec3d cpu_array[mc_ncpus], my_coord;  // CHECK: have to be initialized via corresponding method
	int grid_coord[3];

	MPI_Cart_coords(comm_name,myrank,3,grid_coord);
	my_coord.x = grid_coord[0];
	my_coord.y = grid_coord[1];
	my_coord.z = grid_coord[2];

	//return cpu_array[myrank];
	return my_coord;
}

ivec3d get_cell_loc_coord(ivec3d glob_coord, ivec3d cpu_glob_pos){

	// compute cell local coordinates from cell global coordinate and cpu glob position

	ivec3d cell_loc_coord;

	cell_loc_coord.x  = glob_coord.x - (cpu_glob_pos.x * mc_cpu_cell_dim.x);
	cell_loc_coord.y  = glob_coord.y - (cpu_glob_pos.y * mc_cpu_cell_dim.y);
	cell_loc_coord.z  = glob_coord.z - (cpu_glob_pos.z * mc_cpu_cell_dim.z);

	return cell_loc_coord;
}

void make_particles(cellblock loc_obj){

	 // DO CHECKING

     // method for creating particle objects and fill in cell container


	 std::cout<<" make particles for process : "<< mc_prank << endl;
	 std::cout<<" total particle objects created : "<< mc_tatoms_cpu <<endl;

	 //*************************************************
	 // NOTE:before make_particles():--  cell_block()--> make_cells()--> with cell id and glob coord

	 //cellblock loc_obj;

	 //*************************************************

//	 particle* atom;

	 long m=0;
	 ivec3d cell_glob_coord,nb_cpu_gcoord;
	 int loc_rank;

	 for(long i=0;i<mc_tatoms_cpu;i++){

		 particle atom;

		 // IGNORE THIS PART IF NECESSARY
//		 atom = new particle;
//		 atom->set_mynumber(mc_atomnumber.at(i));
//		 atom->set_mytype(mc_atomtypes.at(i));
//		 atom->set_mymass(mc_atommass.at(i));
//		 atom->set_myposition(mc_positions.at(m++),mc_positions.at(m++),mc_positions.at(m++));

		 // assigning attributes - old implementation

		 atom.set_mynumber(mc_atomnumber.at(i));
 		 atom.set_mytype(mc_atomtypes.at(i));

 		 // counting particles

 		 if((mc_atomtypes.at(i)%mc_real_types) == 0) mc_real_cpu++;   // add real Fe particle
 		 if((mc_atomtypes.at(i)%mc_real_types) == 1) mc_carb_cpu++;   // add real C particle
 		 if((mc_atomtypes.at(i)%mc_real_types) == 2) mc_phold_cpu++;  // add placeholder particle


         atom.set_mymass(mc_atommass.at(i));
         // better split as position x,y & z
 		 atom.set_myposition(mc_positions.at(m++),mc_positions.at(m++),mc_positions.at(m++));
 		 atom.set_myepot(mc_epot.at(i));

 		 //**************************************************
 		 // NOTE: generate velocity
 		 //**************************************************
 		 // by now should have constructed cell with boundary

 		 // also method to return cell_id from particle position


 		 // get global cell coordinate from particle position

   		 cell_glob_coord = cell_coordinate(atom.get_myposition().x, atom.get_myposition().y, atom.get_myposition().z);

    	 // get cpu global coordinate from cell global coordinate
  		 nb_cpu_gcoord = calc_cpu_coord(cell_glob_coord) ;

  		 loc_rank = get_cpu_neighbor(nb_cpu_gcoord.x,nb_cpu_gcoord.y,nb_cpu_gcoord.z);


 		 // get cpu glob position

 		 ivec3d cpu_fact = get_cpu_gcoord(loc_rank);

 		 // get cell local coordinate

 		 ivec3d cell_loc_coord = get_cell_loc_coord(cell_glob_coord,cpu_fact);

 		 // to get memory index of cell in cell list
 		 int cell_index = cell_loc_coord.x + (cell_loc_coord.y * mc_cpu_cell_dim.x) + (cell_loc_coord.z * mc_cpu_cell_dim.x * mc_cpu_cell_dim.x);


         // add particle
 		 loc_obj.get_cell(cell_index).add_particle(atom);

 		 }

	     // clear STL container

	     mc_atomnumber.clear();
	     mc_atomtypes.clear();
         mc_atommass.clear();
         mc_positions.clear();
         mc_epot.clear();

} //make particles

void fill_mc_container(cellblock loc_obj){

	// DO CHECKING

	long par_count=0;

	long stl_counter = 0;

	// loop over cells
	for (long i=0; i<loc_obj.get_ncells();i++){

		par_count = loc_obj.get_cell(i).get_nparticles();

		for(long j=0; j<par_count;j++ ){

            // add particle attributes
			mc_atomnumber.push_back(loc_obj.get_cell(i).get_particle(j).get_mynumber());    // atom id
			mc_atomtypes.push_back(loc_obj.get_cell(i).get_particle(j).get_mytype());       // atom type
            mc_atommass.push_back(loc_obj.get_cell(i).get_particle(j).get_mymass());        // atom mass
			mc_positions.push_back(loc_obj.get_cell(i).get_particle(j).get_myposition().x);
			mc_positions.push_back(loc_obj.get_cell(i).get_particle(j).get_myposition().y); // atom positions
			mc_positions.push_back(loc_obj.get_cell(i).get_particle(j).get_myposition().z);

			mc_velocities.push_back(loc_obj.get_cell(i).get_particle(j).get_myvelocity().x);
			mc_velocities.push_back(loc_obj.get_cell(i).get_particle(j).get_myvelocity().y); // atom velocities
			mc_velocities.push_back(loc_obj.get_cell(i).get_particle(j).get_myvelocity().z);

			mc_epot.push_back(loc_obj.get_cell(i).get_particle(j).get_myepot());            // atom epot

			// delete particle
			loc_obj.get_cell(i).delete_particle(loc_obj.get_cell(i).get_particle(j),j);

		}// loop over particle

		// delete cell
		loc_obj.delete_cell(loc_obj.get_cell(i),i);

	}// loop over cells

}


void create_maxwell_velocities(cellblock loc_obj,double mc_temp, ivec3d*  mc_restriction){

	  // DO CHECKING

	  //   create and fill initial velocities for particles

      double vx,vy,vz;                      // velocity holders
      double imp_x,imp_y,imp_z;             // impulse holders
      double sum_x, sum_y, sum_z;           // sum holders
      double rtemp;                         // reduced temp variable

      int   dof_x,dof_y,dof_z;               // degrees of freedom holders
      long  tot_dof_x,tot_dof_y,tot_dof_z ;

      //cellblock loc_obj;

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
	  cout<<" initial velocities are computed for cpu id:"<< mc_prank <<endl;

}

void make_cells(cellblock loc_obj){

	// DO CHECKING

    //    method for creating cell objects and fill in cellblock container
	//    should assign cell boundary and ids correctly from the field

	std::cout<<" make cells for process : "<< mc_prank << endl;
	std::cout<<" total cell objects created : "<< calc_ncells_cpu() <<endl;

	int cell_no     = 0;

	// I: some method calls
	//  (could be hardcoded for Moving window approach as nxn )
	stack_total   = calc_cpu_stack_total();
	ncells_stack  = calc_ncells_stack(); // NOTE: verify requirement
	ncells_cpu    = calc_ncells_cpu();
	int row_total = calc_cpu_row_total();
	int col_total = calc_cpu_col_total();

    ivec3d cell_gcoord = {0,0,0}; // initialization of cell global coordinate
    ivec3d cell_lcoord = {0,0,0}; // initialization of cell local coordinate

    int cells_stack=0;
    int cell_counter = 0;

    // process/cpu global position (in MPI cartesian topology system)
    cpu_gcoord = get_cpu_gcoord(mc_prank);

    //*****************************************************************
    // create cell objects and fill into cell block container eventually
    // NOTE: hierarchy constraint: before make_particles() make_cells() should be called
    //*****************************************************************

    //cellblock loc_obj;

    // loop over stack

	for (int stack_no=0; stack_no<stack_total; stack_no++){

		for(int row_count=0; row_count<row_total; row_count++){

			for(int col_count=0; col_count<col_total; col_count++){

				  celltype cell_obj;

        		  //cell_no += (stack_no * ncells_stack) + ( ncells_cpu * mc_prank ); // cell id computation

                  // cell id computation
        		  cell_no = cell_counter+ ( ncells_cpu * mc_prank ); // cell id computation

        		  // numbering sequence: x--> y --> z

        		  // cell global coordinates
        		  cell_gcoord.x = col_count +  (col_total * cpu_gcoord.x);
        		  cell_gcoord.y = row_count +  (row_total * cpu_gcoord.y);
        		  cell_gcoord.z = stack_no  +  (stack_total * cpu_gcoord.z);


        		  // cell local coordinate
        		  cell_lcoord = get_cell_loc_coord(cell_gcoord,cpu_gcoord);

        		  // assign values
       		   	  cell_obj.set_cell_id(cell_no);

       		      cell_obj.set_glob_coord(cell_gcoord); // setting global cell coordinate

       		      cell_obj.set_loc_coord(cell_lcoord);  // setting local cell coordinate

       		      loc_obj.add_cell(cell_obj); // adding cell into cel_block

       		      cell_counter++; //  incrementing cell counter

			}//column_count

		}//row_count
		//z dim to be updated here
	}//stack loop

} // make_cells

void count_particles(cellblock loc_obj){

	// CHECK IF REQUIRED

	celltype cell_obj;
	long par_count=0;
	// loop over cells
	for (long i=0; i<loc_obj.get_ncells();i++){
		cell_obj = loc_obj.get_cell(i);
		par_count = loc_obj.get_cell(i).get_nparticles();
		for(long j=0; j<par_count;j++ ){

			if((cell_obj.get_particle(j).get_mytype()%mc_real_types) == 0){
				mc_real_cpu++;   // add real Fe particle
			}
			if((cell_obj.get_particle(j).get_mytype()%mc_real_types) == 1){
				mc_carb_cpu++;   // add real C particle
			}
			if((cell_obj.get_particle(j).get_mytype()%mc_real_types) == 2){
				mc_phold_cpu++;  // add placeholder particle
			}
		}

	}
	// counting particles

	//**** if require put back in make_particles to avoid computation overhead
//	if((mc_atomtypes.at(i)%mc_real_types) == 0) mc_real_cpu++;   // add real Fe particle
//	if((mc_atomtypes.at(i)%mc_real_types) == 1) mc_real_cpu++;   // add real C particle
//	if((mc_atomtypes.at(i)%mc_real_types) == 2) mc_phold_cpu++;  // add placeholder particle


}

void calc_mc_cpu_box(){

	// DO CHECKING
	//  computes cpu_box physical dimensions


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
  if ((0==mc_prank) && (0==mc_volume)) std::cerr("Box Edges are parallel.");

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

	// CHECK IF NECESSARY
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

vec3d calc_cell_dim(double rsample){

	//    compute MonteCarlo cell dimension

	mc_cell_dim.x = (mc_simbox_x.x) / (rsample);
	mc_cell_dim.y = (mc_simbox_y.y) / (rsample);
	mc_cell_dim.z = (mc_simbox_z.z) / (rsample);

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

	// NOT REQUIRED : (COMPUTE FROM MPI_Comm_rank())
	// method to compute cpu_rank from global cell coordinates

	ivec3d red_vec;
	int proc_rank;

	red_vec.x = cell_coord.x / mc_cpu_cell_dim.x;
	red_vec.y = cell_coord.y / mc_cpu_cell_dim.y;
	red_vec.z = cell_coord.z / mc_cpu_cell_dim.z;

	proc_rank = (cols_block * rows_block * red_vec.z ) + (cols_block * red_vec.y) + (red_vec.x);
	return proc_rank;
}

ivec3d calc_cpu_coord(ivec3d cell_coord){

	// method to compute cpu coordinate from its cell global coordinate
	ivec3d red_vec;

	red_vec.x = cell_coord.x / mc_cpu_cell_dim.x;
	red_vec.y = cell_coord.y / mc_cpu_cell_dim.y;
	red_vec.z = cell_coord.z / mc_cpu_cell_dim.z;

	return red_vec;
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

	// compute distance as sum of square
	dist = (scalar_prod(temp_v,temp_v));

	return dist;
}

// end of utility methods
