/****************************************************************************************
 *                                    Monte_classes.cpp
 *
 *                    class definitions for Montecarlo
 *****************************************************************************************/

#include<vector>
#include "Monte_classes.h"

using namespace std;

// particle class
// constructor declaration
particle::particle(void){
	mytype     = 0;
	mynumber   = 0;
	mymass     = 0.0;
	myepot     = 0.0;
	myposition = {0.0,0.0,0.0};
	myvelocity = {0.0,0.0,0.0};
}

// class method declaration
void particle::set_mytype(int type)         { mytype = type;  }
void particle::set_mynumber(long myno)     { mynumber = myno; }
void particle::set_mymass(double mass)     { mymass = mass;   }
void particle::set_myepot(double epot)     { myepot = epot;   }
void particle::set_myposition(double p_x,double p_y,double p_z){
		   myposition.x = p_x;
		   myposition.y = p_y;
		   myposition.z = p_z;
}
void particle::set_myvelocity(double v_x,double v_y,double v_z){
		   myvelocity.x = v_x;
		   myvelocity.y = v_y;
		   myvelocity.z = v_z;
}


// class celltype methods

celltype::celltype(void){
	      cell_id     = 0;
	      nparticles  = 0;
	      lcoord      = {0,0,0};
	      gcoord      = {0,0,0};
	      sample_fact = 0;
}
void celltype::set_cell_id(int uid)                   {  cell_id = uid ; }
void celltype::set_glob_coord(ivec3d cell_gcoord)     {  gcoord = cell_gcoord; }
void celltype::set_loc_coord(ivec3d cell_lcoord)      {  lcoord = cell_lcoord; }
void celltype::add_neighbor(int id, ivec3d position)  {  nbl_list[id] = position; } // nbl asssignments
void celltype::add_sample(void)                       {  sample_fact += 1;        }

//void cellblock::set_x_min(int val){      x_min = val;}
//void cellblock::set_x_max(int val){    x_max = val;  }
//void cellblock::set_y_min(int val){    y_min = val;  }
//void cellblock::set_y_max(int val){    y_max = val;  }
//void cellblock::set_z_min(int val){    y_min = val;  }
//void cellblock::set_z_max(int val){    y_max = val;  }


// class cellblock methods
//cellblock::cellblock(void){
//	  ncells = 0;
	  //mycpu  = 0;
//}

// class cellblock methods
void cellblock::set_mycpu(int cpuid){	mycpu = cpuid;}

celltype cellblock::cell_with_lcoord(ivec3d lcoord, long ncells){

	int flag=0; long index=0;
	ivec3d temp; celltype target;

	while((flag==0) && (index<ncells)){

		// avoid scope problem
		temp = cell_list[index].get_cell_loc_coord();

		// check whether object comparison works ??
		if((temp.x == lcoord.x) && (temp.y ==lcoord.y) && (temp.z == lcoord.z)) {
			target=cell_list[index];
		    flag = 1;
		}
		index++;
	}
	return target;
}
