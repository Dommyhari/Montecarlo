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
	ivec3d gcoord;           // cell global coordinates

	// nbl list variables

	// cell layer
	ivec3d my_left;      // 1 left cell
	ivec3d my_right;     // 2 right cell
	ivec3d my_back;      // 3 back cell
	ivec3d my_front;     // 4 front cell
	ivec3d my_east;      // 5 East cell
	ivec3d my_west;      // 6 west cell
	ivec3d my_south;     // 7 south cell
	ivec3d my_north;     // 8 north cell

	ivec3d my_top;       // 9 top cell

	// cell top layer
	ivec3d top_left;     // 10  top left cell
	ivec3d top_right;    // 11 top right cell
	ivec3d top_front;    // 12 front cell
	ivec3d top_back;     // 13 back cell
	ivec3d top_east;     // 14 East cell
	ivec3d top_west;     // 15 west cell
	ivec3d top_south;    // 16 south cell
	ivec3d top_north;    // 17 north cell

	ivec3d my_bottom;    // 18 bottom cell

	// cell bottom layer
	ivec3d bottom_left;  // 19 bottom left cell
	ivec3d bottom_right; // 20 bottom right cell
	ivec3d bottom_front; // 21 bottom front cell
	ivec3d bottom_back;  // 22 bottom back cell
	ivec3d bottom_east;  // 23 bottom East cell
	ivec3d bottom_west;  // 24 bottom west cell
	ivec3d bottom_south; // 25 bottom south cell
	ivec3d bottom_north; // 26 bottom north cell


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

        // set my layer
        void set_my_left(ivec3d ip_left);
        void set_my_right(ivec3d ip_right);
        void set_my_back(ivec3d ip_back);
        void set_my_front(ivec3d ip_front);
        void set_my_east(ivec3d ip_east);
        void set_my_west(ivec3d ip_west);
        void set_my_south(ivec3d ip_south);
        void set_my_north(ivec3d ip_north);

        // set top layer
        void set_my_top(ivec3d ip_top);

        void set_top_left(ivec3d ip_tleft);
        void set_top_right(ivec3d ip_tright);
        void set_top_back(ivec3d ip_tback);
        void set_top_front(ivec3d ip_tfront);
        void set_top_east(ivec3d ip_teast);
        void set_top_west(ivec3d ip_twest);
        void set_top_south(ivec3d ip_tsouth);
        void set_top_north(ivec3d ip_tnorth);

        // set bottom layer
        void set_my_bottom(ivec3d ip_bottom);

        void set_bottom_left(ivec3d ip_bleft);
        void set_bottom_right(ivec3d ip_bright);
        void set_bottom_back(ivec3d ip_bback);
        void set_bottom_front(ivec3d ip_bfront);
        void set_bottom_east(ivec3d ip_beast);
        void set_bottom_west(ivec3d ip_bwest);
        void set_bottom_south(ivec3d ip_bsouth);
        void set_bottom_north(ivec3d ip_bnorth);


        // if required add other possible neighbor set methods

        void add_particle(particle p)                  { cell.push_back(p); }

        void delete_particle(particle p, int index)    { cell.erase(cell.begin()+ index);}

	    // get methods

	    int get_cell_id()                      { return cell_id;       }

	    long get_nparticles()                  { return cell.size();    }

	    ivec3d get_cell_glob_coord()           { return gcoord;         }

	    // get my layer
	    ivec3d get_my_left()                   { return my_left;        }
	    ivec3d get_my_right()                  { return my_right;       }
	    ivec3d get_my_front()                  { return my_front;       }
	    ivec3d get_my_back()                   { return my_back;        }
	    ivec3d get_my_east()                   { return my_east;        }
	    ivec3d get_my_west()                   { return my_west;        }
	    ivec3d get_my_south()                  { return my_south;       }
	    ivec3d get_my_north()                  { return my_north;       }

        // get my top layer
	    ivec3d get_my_top()                    { return my_top;          }
	    ivec3d get_top_left()                  { return top_left;        }
	    ivec3d get_top_right()                 { return top_right;       }
	    ivec3d get_top_front()                 { return top_front;       }
	    ivec3d get_top_back()                  { return top_back;        }
	    ivec3d get_top_east()                  { return top_east;        }
	    ivec3d get_top_west()                  { return top_west;        }
	    ivec3d get_top_south()                 { return top_south;       }
	    ivec3d get_top_north()                 { return top_north;       }

	    // get my bottom layer
	    ivec3d get_my_bottom()                 { return my_bottom;       }
	    ivec3d get_bottom_left()               { return bottom_left;     }
	    ivec3d get_bottom_right()              { return bottom_right;    }
	    ivec3d get_bottom_front()              { return bottom_front;    }
	    ivec3d get_bottom_back()               { return bottom_back;     }
	    ivec3d get_bottom_east()               { return bottom_east;     }
	    ivec3d get_bottom_west()               { return bottom_west;     }
	    ivec3d get_bottom_south()              { return bottom_south;    }
	    ivec3d get_bottom_north()              { return bottom_north;    }



	    //*******************************
	    // NOTE: if required improvise it..

	    ivec3d get_nbl_id(int id)              { return nbl_list[id];   }

	    //*******************************

	    particle get_particle(long index)     { return cell[index];     }


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

	    void delete_cell(celltype ct, int index)         { cell_list.erase(cell_list.begin() + index);}

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
