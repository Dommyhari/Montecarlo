/*****************************************************************************
 *                           Monte_classes.h
 *
 *                  Headers for class definitions used in Montecarlo API
 *
 *****************************************************************************/


#ifndef MONTE_CLASSES_H_
#define MONTE_CLASSES_H_

#include<vector>
#include "Monte_globals.h"
using namespace std;


// user defined classes



// particle - fundamental user defined data type with associated attributes

class particle{

    int mytype;          // atom type
    long mynumber;       // atom id
    double mymass;       // atom mass
    double myepot;         // atom potential energy
    vec3d myposition;    // atom position
    vec3d myvelocity;    // atom velocity


    public:

	   // list of set functions definitions

       particle();

	   void set_mytype(int type);

	   void set_mynumber(long myno);

	   void set_mymass(double mass);

	   void set_myepot(double epot);

	   void set_myposition(double p_x,double p_y,double p_z);

	   void set_myvelocity(double v_x,double v_y,double v_z);


	   // list of get functions

	   int get_mytype()         { return mytype;      }

	   long get_mynumber()      { return mynumber;    }

	   double get_mymass()      { return mymass;      }

	   double get_myepot()      { return myepot;      }

	   vec3d get_myposition()   { return myposition;  }

	   vec3d get_myvelocity()   { return myvelocity;  }

};

// class celltype

class celltype{

	int cell_id;             // cell number
	long nparticles;         // total no of particles
	ivec3d gcoord;            // cell global coordinates

	// nbl list variables

	ivec3d my_left;      // left cell
	ivec3d my_right;     // right cell
	ivec3d my_front;     // front cell
	ivec3d my_back;      // back cell
	ivec3d my_top;       // top cell
	ivec3d my_bottom;    // bottom cell

	vector<particle> cell;

	//***********************************************************************
	// NOTE:  now hard coded for 6 neighbor list-- to be changed

	//ivec3d nbl_list[mc_nbcells];//
	ivec3d nbl_list[];

	//***********************************************************************

    public:

       	// set methods definition

	    void set_cell_id(int uid);

        void set_glob_coord(ivec3d cell_gcoord);

        void set_my_left(ivec3d ip_left);

        void set_my_right(ivec3d ip_right);

        void set_my_front(ivec3d ip_front);

        void set_my_back(ivec3d ip_back);

        void set_my_top(ivec3d ip_top);

        void set_my_bottom(ivec3d ip_bottom);


        // if required add other possible neighbor set methods

        void add_particle(particle p)                  { cell.push_back(p); }

        void delete_particle(particle p, int index)    { cell.erase(cell.begin()+ index);}

	    // get methods

	    int get_cell_id()                      { return cell_id;       }

	    long get_nparticles()                 { return cell.size();    }

	    ivec3d get_cell_glob_coord()          { return gcoord;         }

	    ivec3d get_my_left()                   { return my_left;        }

	    ivec3d get_my_right()                  { return my_right;       }

	    ivec3d get_my_front()                  { return my_front;       }

	    ivec3d get_my_back()                   { return my_back;        }

	    ivec3d get_my_top()                    { return my_top;         }

	    ivec3d get_my_bottom()                 { return my_bottom;      }



	    //*******************************
	    // NOTE: if required improvise it..

	    ivec3d get_nbl_id(int id)              { return nbl_list[id];   }

	    //*******************************

	    particle get_particle(long index)     { return cell[index];    }


};

// class cellblock

class cellblock{

	long ncells; // total no of cells
	int mycpu;  // my cpu id
	int x_min, x_max, y_min, y_max, z_min, z_max;

	vector<celltype> cell_list;     // container of cells


    public:

	    // set methods

	    void set_mycpu(int cpuid);

	    //*************************************************************
	    // NOTE: check for requirements of below member functions

	    void set_x_min(int val);

	    void set_y_min(int val);

	    void set_z_min(int val);

	    void set_x_max(int val);

	    void set_y_max(int val);

	    void set_z_max(int val);

	    //*************************************************************

	    void add_cell(celltype ct) { cell_list.push_back(ct);}

	    // get methods
	    long get_ncells()        { return cell_list.size();}
	    int get_mycpu()          { return mycpu; }

	    int get_x_min()          { return x_min; }
	    int get_x_max()          { return x_max; }

	    int get_y_min()          { return y_min; }
	    int get_y_max()          { return y_max; }

	    int get_z_min()          { return z_min; }
	    int get_z_max()          { return z_max; }

	    //****************************************************************

	    celltype get_cell(int index)        { return cell_list[index]; }

};

#endif /* MONTE_CLASSES_H_ */
