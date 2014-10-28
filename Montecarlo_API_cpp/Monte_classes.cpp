/****************************************************************************************
 *                                    Monte_classes.cpp
 *
 *                    class definitions for Montecarlo
 *****************************************************************************************/

#include<vector>
#include<iostream>
#include "Monte_classes.h"


using namespace std;

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

void  particle::set_mytype(int type){
		   mytype = type;
}

void particle::set_mynumber(long myno){
		   mynumber = myno;
}

void particle::set_mymass(double mass){
		   mymass = mass;
}

void particle::set_myepot(double epot){
	       myepot = epot;
}
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

void celltype::set_cell_id(int uid){
		cell_id = uid ;
}

void celltype::set_glob_coord(ivec3d cell_gcoord){
	    gcoord.x = cell_gcoord.x;
	    gcoord.y = cell_gcoord.y;
	    gcoord.z = cell_gcoord.z;
}

// nbl asssignments

void celltype::set_my_left(ivec3d ip_left){

	    my_left.x = ip_left.x;
	    my_left.y = ip_left.y;
	    my_left.z = ip_left.z;

	    nbl_list[0] = my_left;
}

void celltype::set_my_right(ivec3d ip_right){

	    my_right.x = ip_right.x;
	    my_right.y = ip_right.y;
	    my_right.z = ip_right.z;

	    nbl_list[1] = my_right;
}

void celltype::set_my_front(ivec3d ip_front){

	    my_front.x = ip_front.x;
	    my_front.y = ip_front.y;
	    my_front.z = ip_front.z;

	    nbl_list[2] = my_front;
}

void celltype::set_my_back(ivec3d ip_back){

	    my_back.x = ip_back.x;
	    my_back.y = ip_back.y;
	    my_back.z = ip_back.z;

	    nbl_list[3] = my_back;
}

void celltype::set_my_top(ivec3d ip_top){

	    my_top.x = ip_top.x;
	    my_top.y = ip_top.y;
	    my_top.z = ip_top.z;

	    nbl_list[4] = my_top;
}

void celltype::set_my_bottom(ivec3d ip_bottom){

	    my_bottom.x = ip_bottom.x;
	    my_bottom.y = ip_bottom.y;
	    my_bottom.z = ip_bottom.z;

	    nbl_list[5] = my_bottom;
}

// class cellblock methods

void cellblock::set_mycpu(int cpuid){
	mycpu = cpuid;
}

void cellblock::set_x_min(int val){
    x_min = val;
}

void cellblock::set_x_max(int val){
    x_max = val;
}

void cellblock::set_y_min(int val){
    y_min = val;
}

void cellblock::set_y_max(int val){
    y_max = val;
}

void cellblock::set_z_min(int val){
    y_min = val;
}

void cellblock::set_z_max(int val){
    y_max = val;
}

//void cellblock::add_cell(celltype ct){
//	cell_list.push_back(const ct);
//}

//int main(){

//	using namespace std;

//	cout<<"yup I am working"<<endl;
//}
