/*
 * Monte_new_classes.h
 *
 *  Created on: Mar 6, 2015
 *      Author: ganeshfx
 */

#ifndef MONTE_NEW_CLASSES_H_
#define MONTE_NEW_CLASSES_H_

#include<vector>
#include<iomanip>
#include<fstream>
#include<string>
#include<iostream>
#include<ctime>
#include<cmath>
#include<cstdlib>
#include<random>

#include "mpi.h"

using namespace std;



// user defined classes
class vec3d{  // double vector class
   public:
	 double x; double y; double z;
};
class ivec3d{ // integer vector class
   public:
	 int x; int y; int z;
};
class ivec6d{ // integer vector class (used for sample cell boundary)
   public:
	 int xmin; int xmax; int ymin; int ymax; int zmin; int zmax;
};

// particle - fundamental user defined data type with associated attributes
class particle{

    long mynumber;       // atom id
    int mytype;          // atom type
    double mymass;      // atom mass
    vec3d myposition;    // atom position
    vec3d myvelocity;    // atom velocity
    double myepot;      // atom potential energy


    public:
	   // list of set functions definitions
       particle();
	   void set_mynumber(long myno);
       void set_mytype(int type);
	   void set_mymass(double mass);
	   void set_myposition(double p_x,double p_y,double p_z);
	   void set_myvelocity(double v_x,double v_y,double v_z);
	   void set_myepot(double epot);

	   // list of get functions

	   long get_mynumber()      { return mynumber;    }
	   int get_mytype()         { return mytype;      }
	   double get_mymass()      { return mymass;      }
	   vec3d get_myposition()   { return myposition;  }
	   vec3d get_myvelocity()   { return myvelocity;  }
	   double get_myepot()      { return myepot;      }
};

// class celltype
class celltype{

	int cell_id;              // cell number
	long nparticles;          // total no of particles
	ivec3d gcoord;            // cell global coordinates
	ivec3d lcoord;            // cell local  coordinates
    int sample_fact;         // sample factor (how frequent a cell is sampled)
    vector<particle> cell;    // cell list contain particles belong to it
    ivec3d nbl_list[27];      // hard-coded for 27 neighbors (extend if required)
    long particle_counter;   // no of particles currently in list

    public:

	    celltype();
       	// set methods definition
	    void set_cell_id(int uid);
        void set_glob_coord(ivec3d cell_gcoord);
        void set_loc_coord(ivec3d cell_lcoord);

        void update_particle_counter (void);
        void add_neighbor(int id, ivec3d position);
        void add_sample(void);   // sample frequency counter

        void add_particle(particle p,long cell_ind){
        	//cell.at(ind) = p;
        	//cout << " inside add particle " << endl;
        	//cout << " cell_index  : " << cell_ind << endl;
        	//cout << " particle id : " << p.get_mynumber() << endl;

        	cell.push_back(p); }


        // for performance make manipulations
        void delete_particle(particle p, int index)    { cell.erase(cell.begin()+ index);}

	    // get methods
	    int get_cell_id()                       { return cell_id;       }  //  cell id
	    long get_nparticles()                  { return cell.size();    }  // cell total particles
	    ivec3d get_cell_glob_coord()           { return gcoord;         }  // cell global coordinate
        ivec3d get_cell_loc_coord()            { return lcoord;         }  // cell local coordinate
        int get_sample_factor()                { return sample_fact;    }  // cell sample weight
	    ivec3d get_nbl_id(int id)               { return nbl_list[id];   }  // cell neighbor
	    particle get_particle(long index)      { return cell[index];  }    // get particle with index
        long get_particle_counter()           { return particle_counter;}

	    // clear method
	    void clear_all_particles()            { cell.clear();           }   // empty container
};


// class cellblock
// cellblock - data_type analogy to CPU domain

class cellblock{

	int ncells;                                    // total no of cells
	int mycpu;                                      // my cpu id
	vector<celltype> cell_list;                     // container of cells

    public:

        //CFR
	    void set_mycpu(int cpuid);
        void set_ncells(int cell_count);

        // set methods
	    void add_cell(celltype ct)                        { cell_list.push_back(ct);}
	    void delete_cell(int index)                      { cell_list.erase(cell_list.begin() + index);}

	    // get methods
	    int get_ncells()                                  { return ncells;}      // total no cells in block
	    int get_cell_list_size()                         { return cell_list.size();}
	    int get_mycpu()                                   { return mycpu; }                // my cpu id

	    celltype get_cell(int index)                      { return cell_list[index]; }     //  get cell with index
	    celltype cell_with_lcoord(ivec3d lcoord , long ncells);                            //  cell local coordinate

	    // clear method
	    void clear_all_cells()                           { return cell_list.clear(); }
};


#endif /* MONTE_NEW_CLASSES_H_ */
