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
	    gcoord = cell_gcoord;
}

// nbl asssignments

// set my layer

void celltype::set_my_left(ivec3d ip_left){
	    my_left = ip_left;
	    nbl_list[0] = my_left;
}
void celltype::set_my_right(ivec3d ip_right){
	    my_right = ip_right;
	    nbl_list[1] = my_right;
}

void celltype::set_my_back(ivec3d ip_back){
	    my_back = ip_back;
	    nbl_list[2] = my_back;
}

void celltype::set_my_front(ivec3d ip_front){
	    my_front = ip_front;
	    nbl_list[3] = my_front;
}
void celltype::set_my_east(ivec3d ip_east){
	    my_east = ip_east;
	    nbl_list[4] = my_east;
}

void celltype::set_my_west(ivec3d ip_west){
	    my_west = ip_west;
	    nbl_list[5] = my_west;
}
void celltype::set_my_south(ivec3d ip_south){
	    my_south= ip_south;
	    nbl_list[6] = my_south;
}

void celltype::set_my_north(ivec3d ip_north){
	    my_north = ip_north;
	    nbl_list[7] = my_north;
}

// set my top layer

void celltype::set_my_top(ivec3d ip_top){
	    my_top = ip_top;
	    nbl_list[8] = my_top;
}

void celltype::set_top_left(ivec3d ip_tleft){
	    top_left = ip_tleft;
	    nbl_list[9] = top_left;
}

void celltype::set_top_right(ivec3d ip_tright){
	    top_right = ip_tright;
	    nbl_list[10] = top_right;
}

void celltype::set_top_back(ivec3d ip_tback){
	    top_back = ip_tback;
	    nbl_list[11] = top_back;
}

void celltype::set_top_front(ivec3d ip_tfront){
	    top_front = ip_tfront;
	    nbl_list[12] = top_front;
}

void celltype::set_top_east(ivec3d ip_teast){
	    top_east = ip_teast;
	    nbl_list[13] = top_east;
}

void celltype::set_top_west(ivec3d ip_twest){
	    top_west = ip_twest;
	    nbl_list[14] = top_west;
}

void celltype::set_top_south(ivec3d ip_tsouth){
	    top_south = ip_tsouth;
	    nbl_list[15] = top_south;
}

void celltype::set_top_north(ivec3d ip_tnorth){
	    top_north = ip_tnorth;
	    nbl_list[16] = top_north;

}
// set my bottom layer

void celltype::set_my_bottom(ivec3d ip_bottom){
	    my_bottom = ip_bottom;
	    nbl_list[17] = my_bottom;
}

void celltype::set_bottom_left(ivec3d ip_bleft){
	    bottom_left = ip_bleft;
	    nbl_list[18] = bottom_left;
}

void celltype::set_bottom_right(ivec3d ip_bright){
	    bottom_right = ip_bright;
	    nbl_list[19] = bottom_right;
}

void celltype::set_bottom_back(ivec3d ip_bback){
	    bottom_back = ip_bback;
	    nbl_list[20] = bottom_back;
}

void celltype::set_bottom_front(ivec3d ip_bfront){
	    bottom_front = ip_bfront;
	    nbl_list[21] = bottom_front;
}

void celltype::set_bottom_east(ivec3d ip_beast){
	    bottom_east = ip_beast;
	    nbl_list[22] = bottom_east;
}

void celltype::set_bottom_west(ivec3d ip_bwest){
	    bottom_west = ip_bwest;
	    nbl_list[23]= bottom_west;
}

void celltype::set_bottom_south(ivec3d ip_bsouth){
	    bottom_south = ip_bsouth;
	    nbl_list[24] = bottom_south;
}

void celltype::set_bottom_north(ivec3d ip_bnorth){
	    bottom_north = ip_bnorth;
	    nbl_list[25] = bottom_north;

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
