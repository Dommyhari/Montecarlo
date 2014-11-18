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

void read_update_config (char* fname,particle pobj){

	// NEED CHANGES

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
        if(f_type > mc_real_types) f_type-=mc_real_types;

        fbuffer_type.push_back(f_type);


        fin>>f_mass;      fbuffer_mass.push_back(f_mass);

        // shift back to native coordinate system
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

    //----------------------------------------------------------//
    //           2- Updating and shifting part                  //
    //----------------------------------------------------------//

    // NOTE: later could be changed with n - particles in sphere count

    cellblock loc_obj;

    particle atom;


    cout<< "Updating particle attributes " << endl;

    for (int i=0; i< (fbuffer_id.size()-1); i++){

    	// global boundary check and position update to actual box dimensions

    	if (fbuffer_pos_x.at(i) < 0.0) fbuffer_pos_x.at(i)+=mc_simbox_x.x ;

    	if (fbuffer_pos_x.at(i) > mc_simbox_x.x) fbuffer_pos_x.at(i)-=mc_simbox_x.x ;

    	if (fbuffer_pos_y.at(i) < 0.0) fbuffer_pos_y.at(i)+=mc_simbox_y.y ;

    	if (fbuffer_pos_y.at(i) > mc_simbox_y.y) fbuffer_pos_y.at(i)-=mc_simbox_y.y ;

    	if (fbuffer_pos_z.at(i) < 0.0) fbuffer_pos_z.at(i)+=mc_simbox_z.z ;

    	if (fbuffer_pos_z.at(i) > mc_simbox_z.z) fbuffer_pos_z.at(i)-=mc_simbox_z.z ;


    	// get global cell coordinate from particle position
    	ivec3d glob_coord = cell_coordinate(fbuffer_pos_x.at(i), fbuffer_pos_y.at(i), fbuffer_pos_z.at(i));

		// get cpu rank from particle physical coordinate
		int loc_rank = get_cpu_rank(glob_coord);


		if (loc_rank != mc_prank){
		//   NOTE ********************
		// if loc_rank != myrank send particle to target cpu
		//****************************
		}

		else{

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
		}

    } // for loop

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
		cpu_gcoord = get_cpu_gcoord(mc_prank); int x = cpu_gcoord.x; int y = cpu_gcoord.y; int z = cpu_gcoord.z;

		double my_pos[3], rec_pos_0[3],rec_pos_1[3],rec_pos_2[3],rec_pos_3[3],rec_pos_4[3],rec_pos_5[3],rec_pos_6[3],rec_pos_7[3];

		double* rec_pos[8]={rec_pos_0,rec_pos_1,rec_pos_2,rec_pos_3,rec_pos_4,rec_pos_5,rec_pos_6,rec_pos_7};

		my_pos[0]=pobj.get_myposition().x; my_pos[1]=pobj.get_myposition().y; my_pos[2]=pobj.get_myposition().z;


		// communication part -- request neighbors with chosen particle

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


		// do sweep test to construct buffers


	    //******************************************************************
	    // NOTE: HARDCODED FOR 8 CELLS IN A WINDOW (OWN CPU)
        // NEED TO BE GENERALIZED LATER DURING BIAS APPROACH

		// local cell coordinates zone
	    ivec3d test, nb_test[19]; ivec6d zone_limit[8];ivec6d zone_neigh[8];

	    // window sampling cells
	    celltype sample_cells[8];

	    //window neighbor cells
        celltype neighb_1[4],neighb_2[2],neighb_4[1],neighb_6[2],neighb_3[4],neighb_5[2],neighb_7[4];

        int nb_xmin,nb_ymin,nb_zmin,nb_xmax,nb_ymax,nb_zmax,nb_x,nb_y,nb_z;

	    zone_limit[0] = {0,1,0,1,0,1}; // window 0
	    zone_limit[1] = {0,1,0,1,1,2}; // window 1
	    zone_limit[2] = {1,2,0,1,0,1}; // window 2
	    zone_limit[3] = {1,2,0,1,1,2}; // window 3
	    zone_limit[4] = {0,1,1,2,0,1}; // window 4
	    zone_limit[5] = {0,1,1,2,1,2}; // window 5
	    zone_limit[6] = {1,2,1,2,0,1}; // window 6
	    zone_limit[7] = {1,2,1,2,1,2}; // window 7


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
	    nb_test[0].x =nb_xmin; nb_test[0].y = nb_ymin; nb_test[0].z = nb_z; neighb_1[0]=bobj.cell_with_lcoord(nb_test[0]);
	    nb_test[1].x =nb_xmax; nb_test[1].y = nb_ymin; nb_test[1].z = nb_z; neighb_1[1]=bobj.cell_with_lcoord(nb_test[1]);
	    nb_test[2].x =nb_xmin; nb_test[2].y = nb_ymax; nb_test[2].z = nb_z; neighb_1[2]=bobj.cell_with_lcoord(nb_test[2]);
	    nb_test[3].x =nb_xmax; nb_test[3].y = nb_ymax; nb_test[3].z = nb_z; neighb_1[3]=bobj.cell_with_lcoord(nb_test[3]);

        //Neighbor-2 (N2)
        nb_test[4].x = nb_x; nb_test[4].y = nb_ymin; nb_test[4].z = nb_z; neighb_2[0] = bobj.cell_with_lcoord(nb_test[4]);
        nb_test[5].x = nb_x; nb_test[5].y = nb_ymax; nb_test[5].z = nb_z; neighb_2[1] = bobj.cell_with_lcoord(nb_test[5]);

        //Neighbor-3 (N4)
        nb_test[6].x = nb_x; nb_test[6].y = nb_y; nb_test[6].z = nb_z; neighb_4[0] = bobj.cell_with_lcoord(nb_test[6]);

        //Neighbor-4 (N6)
        nb_test[7].x = nb_xmin; nb_test[7].y = nb_y; nb_test[7].z = nb_z; neighb_6[0] = bobj.cell_with_lcoord(nb_test[7]);
        nb_test[8].x = nb_xmax; nb_test[8].y = nb_y; nb_test[8].z = nb_z; neighb_6[1] = bobj.cell_with_lcoord(nb_test[8]);

        //Neighbor-5 (N3)
	    nb_test[9].x =nb_x;   nb_test[9].y = nb_ymin; nb_test[9].z = nb_zmin;  neighb_3[0]=bobj.cell_with_lcoord(nb_test[9]);
	    nb_test[10].x =nb_x; nb_test[10].y = nb_ymin; nb_test[10].z = nb_zmax; neighb_3[1]=bobj.cell_with_lcoord(nb_test[10]);
	    nb_test[11].x =nb_x; nb_test[11].y = nb_ymax; nb_test[11].z = nb_zmin; neighb_3[2]=bobj.cell_with_lcoord(nb_test[11]);
	    nb_test[12].x =nb_x; nb_test[12].y = nb_ymax; nb_test[12].z = nb_zmax; neighb_3[3]=bobj.cell_with_lcoord(nb_test[12]);

        //Neighbor-6 (N5)
        nb_test[13].x = nb_x; nb_test[13].y = nb_y; nb_test[13].z = nb_zmin; neighb_5[0] = bobj.cell_with_lcoord(nb_test[13]);
        nb_test[14].x = nb_x; nb_test[14].y = nb_y; nb_test[14].z = nb_zmax; neighb_5[1] = bobj.cell_with_lcoord(nb_test[14]);

	    //Neighbor-7 (N7)
	    nb_test[15].x =nb_xmin; nb_test[15].y = nb_y; nb_test[15].z = nb_zmin; neighb_7[0]=bobj.cell_with_lcoord(nb_test[15]);
	    nb_test[16].x =nb_xmin; nb_test[16].y = nb_y; nb_test[16].z = nb_zmax; neighb_7[1]=bobj.cell_with_lcoord(nb_test[16]);
	    nb_test[17].x =nb_xmax; nb_test[17].y = nb_y; nb_test[17].z = nb_zmin; neighb_7[2]=bobj.cell_with_lcoord(nb_test[17]);
	    nb_test[18].x =nb_xmax; nb_test[18].y = nb_y; nb_test[18].z = nb_zmax; neighb_7[3]=bobj.cell_with_lcoord(nb_test[18]);


	    int count=0;

	    // prepare SAMPLE cell list as per sample window id

	    for(int i=zone_limit[win_id].xmin;i<=zone_limit[win_id].xmax;i++){
	    	for(int j=zone_limit[win_id].ymin;j<=zone_limit[win_id].ymax;j++){
	    		for(int k=zone_limit[win_id].zmin;k<=zone_limit[win_id].zmax;k++){

	    			test.x = i; test.y = j; test.z = k;
	    			sample_cells[count] = bobj.cell_with_lcoord(test);

	                count++;
	    		}// k loop
	    	}// j loop
	    }// i loop




}

//void construct_sphere(particle pobj, celltype cobj, int win_id,char* filename){
//
//	// NEED CHANGES
//
//	// select neighbors as per sample window
//	int xfact = window_x[win_id];
//	int yfact = window_y[win_id];
//	int zfact = window_z[win_id];
//
//	// send neighbors ----- odd -even process communications (in z direction)
//	int x_send_phase[4] = { 0,xfact,xfact, 0};
//	int y_send_phase[4] = { 0, 0,yfact,yfact};
//	int z_send_phase[4] = {zfact,zfact,zfact,zfact};
//
//	// Receive neighbors ----- odd -even process communications (in z direction)
//	int x_recv_phase[4] = { 0,xfact,xfact, 0};
//	int y_recv_phase[4] = { 0, 0,yfact,yfact};
//	int z_recv_phase[4] = {zfact,zfact,zfact,zfact};
//
//    // send neighbors ----- odd-odd / even-even process communications (in x/y direction)
//	int x_send_next[3] = {xfact,xfact, 0 };
//	int y_send_next[3] = { 0,yfact,yfact };
//	int z_send_next[3] = { 0, 0, 0 };
//
//	// Receive neighbors ----- odd-odd / even-even process communications (in x/y direction)
//	int x_recv_next[3] = {-xfact,-xfact, 0 };
//	int y_recv_next[3] = { 0,-yfact,-yfact };
//	int z_recv_next[3] = { 0, 0, 0 };
//
//	// Construct sphere around chosen particle and write the particles with in (r_sweep = r_cut + r_sample) to configuration file
//
//	ivec3d mycell_glob;                      // chosen cell global coordinate
//	ivec3d temp_glob, temp_loc, temp_block ; // temporary holder for  coordinates
//
//	mycell_glob = cobj.get_cell_glob_coord();
//
//	vec3d temp_pos, temp_vel, vec1, vec2;
//	long temp_id; int temp_type;            // temporary holder for particle attributes
//	double temp_mass, temp_epot, dist_check;
//
//	int xflag[6]={0,0,0,0,0,0};
//	int yflag[6]={0,0,0,0,0,0};
//	int zflag[6]={0,0,0,0,0,0};     // some utility flags for boundary checks
//
//	ivec3d g_max = mc_global_cell_dim;    // global cell array dimension - Max
//	ivec3d g_min = {0,0,0};               // global cell array dimension - Min
//
//	//char * fname  = file_name ;                // file name as per simulation state from signature
//
//	// object definitions
//	particle temp_part;
//
//	celltype sphere_cell; // list of particles from neighbor cells + own cell
//	cellblock loc_obj;
//
//	// manipulator settings
//	cout << fixed << right;
//
//	// initialize reference particle
//	vec1 = pobj.get_myposition();
//
//	// get reference particle sweep boundary
//	double ref_xmin = vec1.x - (mc_rsweep + mc_sphere_wall);
//	double ref_ymin = vec1.y - (mc_rsweep + mc_sphere_wall);   // minimum boundary
//	double ref_zmin = vec1.z - (mc_rsweep + mc_sphere_wall);
//
//    // global boundary check and flag initialization
//
//	// global_xmin
//    if( temp_glob.x == g_min.x) xflag[0] =  -1;
//    // global_xmax
//    if( temp_glob.x == g_max.x) xflag[1] =   1;
//    // global_ymin
//    if( temp_glob.y == g_min.y) yflag[2] =  -1;
//    // global_ymax
//    if( temp_glob.y == g_max.y) yflag[3] =   1;
//    // global_zmin
//    if( temp_glob.z == g_min.z) zflag[4] =  -1;
//    // global_zmax
//    if( temp_glob.z == g_max.z) zflag[5] =   1;
//
//
//    //********** loop over Neighbor cells list **************
//
//    // Later neighbor cells could be constructed with SMART nbl list and looped over the same
//	for(int i=0; i<=mc_nbcells ; i++){
//
//          // getting nbl cells, cpu and cell address
//
//	      // get neighbor cell global coordinate
//		  if(i!=mc_nbcells){  temp_glob = cobj.get_nbl_id(i); }
//
//	      // include chosen particle own cell
//	      if(i==mc_nbcells){  temp_glob = mycell_glob; }
//
//	      // get cpu rank
//	      int loc_rank = get_cpu_rank(temp_glob);
//
//	      // make changes using virtual topology
//	      temp_block = get_cpu_gcoord(loc_rank);
//
//	      //  Foreign cpu cell
//	      if(mc_prank != loc_rank) {
//
//	    	   // preliminary notifier
//		       cout<<"Note from sphere construction pid: "<<mc_prank <<"\n"<<" Require process communication with : "<< loc_rank <<endl;
//
//		       // unpack datas from MPI packed buffers
//               //****  NOTE MPI receive calls here
//		       // into  buffer
//		       // after receiving and unpacking
//
//		       long p_total; // should get from MPI receive
//		       long p_c,v_c; // value holders
//
//		       // loop over n_particles received
//
//		       for( long k=0; k<p_total; k++ ){
//
//		    	  //  p_c =  k+(size*3);
//		    	  //  v_c =  k+(size*3*3);
//		    	  //  temp_type   = buffer[k];
//		    	  //  temp_id     = buffer[k+size];
//		    	  //  temp_mass   = buffer[k+(size*2)];
//		    	  //  temp_epot   = buffer[k+(size*3)]; // NOTE: need to be verified
//		    	  //  temp_pos.x  = buffer[p_c++] + (xflag[i] * mc_simbox_x.x);
//		    	  //  temp_pos.y  = buffer[p_c++] + (yflag[i] * mc_simbox_y.y);
//		    	  //  temp_pos.z  = buffer[p_c++] + (zflag[i] * mc_simbox_z.z);
//		    	  //  temp_vel.x  = buffer[v_c++];
//		    	  //  temp_vel.y  = buffer[v_c++];
//		    	  //  temp_vel.z  = buffer[v_c++];
//
//
//		          vec2 = {temp_pos.x,temp_pos.y,temp_pos.z}; // particle two position assignment
//		          dist_check = distance_vect(vec1,vec2);
//
//		          //**************************************************************
//		          //NOTE: to be updated in callee process to make dist_check
//		          //**************************************************************
//
//		          // ----- cutoff check
//
//		          if(dist_check <= mc_rsweep + mc_sphere_wall){
//
//		        	  // declare virtual particles on sphere boundary
//		        	  if( dist_check > mc_rsweep){
//
//                          // ignore placeholders on sphere wall
//		        		  if (temp_type != 2 ){ // HC: now hardcoded for placeholders
//
//				        	  temp_part.set_mytype(temp_type + mc_real_types);  // HC: hardcoded for ntypes=3
//				        	  temp_part.set_mynumber(temp_id);
//				        	  temp_part.set_mymass(temp_mass);
//
//				        	  // shift filtered particle coordinates to reference sphere
//				        	  temp_part.set_myposition((temp_pos.x - ref_xmin),(temp_pos.y - ref_ymin),(temp_pos.z - ref_zmin));
//
//				        	  temp_part.set_myvelocity(temp_vel.x,temp_vel.y,temp_vel.z);
//				        	  temp_part.set_myepot(temp_epot);
//
//				        	  sphere_cell.add_particle(temp_part);
//
//		        		  }
//		        	  }
//		        	  else{ // include all particles inside mc_rsweep
//
//			        	  temp_part.set_mytype(temp_type);
//			        	  temp_part.set_mynumber(temp_id);
//			        	  temp_part.set_mymass(temp_mass);
//
//			        	  // shift filtered particle coordinates to reference sphere
//			        	  temp_part.set_myposition((temp_pos.x - ref_xmin),(temp_pos.y - ref_ymin),(temp_pos.z - ref_zmin));
//
//			        	  temp_part.set_myvelocity(temp_vel.x,temp_vel.y,temp_vel.z);
//			        	  temp_part.set_myepot(temp_epot);
//
//			        	  sphere_cell.add_particle(temp_part);
//
//		        	  }
//
//		          } // cut-off check
//
//		       } // end of particles loop
//
//	      }//end of cpu foreign cells
//
//
//	      // cpu native cells
//	      else{
//
//		       // get cell local coordinate
//		       ivec3d cell_loc_coord = get_cell_loc_coord(temp_glob,temp_block);
//
//		       // to get memory index of cell in cell list
//		       int cell_index = cell_loc_coord.x + (cell_loc_coord.y * mc_cpu_cell_dim.x) + (cell_loc_coord.z * mc_cpu_cell_dim.x * mc_cpu_cell_dim.x);
//
//		       // get particles count
//		       long p_total = loc_obj.get_cell(cell_index).get_nparticles();
//
//		       for( long pcount; pcount<p_total; pcount++ ){
//
//		    	   // reading values
//		    	   temp_id    = loc_obj.get_cell(cell_index).get_particle(pcount).get_mynumber();
//		    	   temp_type  = loc_obj.get_cell(cell_index).get_particle(pcount).get_mytype();
//		    	   temp_mass  = loc_obj.get_cell(cell_index).get_particle(pcount).get_mymass();
//		    	   temp_pos   = loc_obj.get_cell(cell_index).get_particle(pcount).get_myposition();
//		    	   temp_vel   = loc_obj.get_cell(cell_index).get_particle(pcount).get_myvelocity();
//		    	   temp_epot  = loc_obj.get_cell(cell_index).get_particle(pcount).get_myepot();
//
//		    	   // particle 2 position assignment
//			       vec2 = temp_pos;
//
//			       dist_check = distance_vect(vec1,vec2);
//
//			          // ----- cutoff check
//
//			          //NOTE: SQUARES should be included (check distance between vectors routine)
//			          if(dist_check <= mc_rsweep + mc_sphere_wall){
//
//			        	  // declare virtual particles on sphere boundary
//			        	  if( dist_check > mc_rsweep){
//
//	                          // ignore placeholders on sphere wall
//			        		  if (temp_type != 2 ){ // HC: now hardcoded for placeholders
//                                  // introduce virtual particles
//					        	  temp_part.set_mytype(temp_type + 3);  // HC: hardcoded for ntypes=3
//					        	  temp_part.set_mynumber(temp_id);
//					        	  temp_part.set_mymass(temp_mass);
//
//					        	  // shift filtered particle coordinates to reference sphere
//					        	  temp_part.set_myposition((temp_pos.x - ref_xmin),(temp_pos.y - ref_ymin),(temp_pos.z - ref_zmin));
//
//					        	  temp_part.set_myvelocity(temp_vel.x,temp_vel.y,temp_vel.z);
//					        	  temp_part.set_myepot(temp_epot);
//
//					        	  sphere_cell.add_particle(temp_part);
//
//			        		  }
//			        	  }
//			        	  else{ // include all particles inside mc_rsweep
//
//				        	  temp_part.set_mytype(temp_type);
//				        	  temp_part.set_mynumber(temp_id);
//				        	  temp_part.set_mymass(temp_mass);
//
//				        	  // shift filtered particle coordinates to reference sphere
//				        	  temp_part.set_myposition((temp_pos.x - ref_xmin),(temp_pos.y - ref_ymin),(temp_pos.z - ref_zmin));
//
//				        	  temp_part.set_myvelocity(temp_vel.x,temp_vel.y,temp_vel.z);
//				        	  temp_part.set_myepot(temp_epot);
//
//				        	  sphere_cell.add_particle(temp_part);
//
//			        	  }
//
//			          } // cut-off check
//
//		       } // end of particles loop
//
//          }// end of native cell
//
//	}// end of sphere cell list - 1
//
//	// shift chosen particle coordinate to ref.sphere
//	particle temp_pobj = pobj;
//	vec3d temp_pos = temp_pobj.get_myposition();
//	temp_pobj.set_myposition(temp_pos.x-ref_xmin,temp_pos.y-ref_ymin,temp_pos.z-ref_zmin);
//
//	// add the randomly selected reference particle
//	sphere_cell.add_particle(temp_pobj);
//
//
//
//
//
//
//	// * THE FOLLOWING PART COULD BE USED AS SUCH
//	// file writing methods
//
//	long my_number;
//	int my_type;
//    double my_mass, my_epot;
//    vec3d my_pos, my_vel;
//
//    // opening file stream object
//    ofstream fout(filename, ios_base::out);
//
//    cout<<" Sphere constructor writing to file : " << filename<<endl;
//
//    //********** loop over revised sphere cell list - 2 **************
//
//	long fp_total = sphere_cell.get_nparticles();
//
//    for(long fp; fp<fp_total; fp++){
//
//        	  my_number  = sphere_cell.get_particle(fp).get_mynumber();
//        	  my_type    = sphere_cell.get_particle(fp).get_mytype();
//        	  my_mass    = sphere_cell.get_particle(fp).get_mymass();
//        	  my_pos     = sphere_cell.get_particle(fp).get_myposition();
//        	  my_vel     = sphere_cell.get_particle(fp).get_myvelocity(); // optional with #ifdef or so
//              my_epot    = sphere_cell.get_particle(fp).get_myepot();
//
//        	  // file flush
//        	  fout << setw(8) << my_number
//        		   << setw(6) << my_type
//        		   << setw(6) << setprecision(8) << my_mass
//        		   << setw(6) << setprecision(8) << my_pos.x
//        		   << setw(6) << setprecision(8) << my_pos.y
//        		   << setw(6) << setprecision(8) << my_pos.z
//        		   << setw(6) << setprecision(8) << my_vel.x
//        		   << setw(6) << setprecision(8) << my_vel.y
//        		   << setw(6) << setprecision(8) << my_vel.z
//        		   << setw(6) << setprecision(8) << my_epot
//        		   << endl;
//
//    }
//
//	fout.close(); // closing outfile connection
//
//}

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
