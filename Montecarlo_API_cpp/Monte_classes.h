/*****************************************************************************
 *                           Monte_classes.h
 *
 *                  Headers for class definitions used in Montecarlo API
 *
 *****************************************************************************/


#ifndef MONTE_CLASSES_H_
#define MONTE_CLASSES_H_

#include<vector>
using namespace std;


// user defined classes

class vec3d{

   public:
	 double x; double y;double z;
};

class ivec3d{

   public:
	 int x; int y; int z;
};

class ivec6d{

   public:
	 int xmin; int xmax; int ymin; int ymax; int zmin; int zmax;
};


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

	int cell_id;              // cell number
	long nparticles;          // total no of particles
	ivec3d gcoord;            // cell global coordinates
	ivec3d lcoord;            // cell local  coordinates
    int sample_fact=0;       // sample factor (how frequent a cell is sampled)

	vector<particle> cell;

	//***********************************************************************

	//ivec3d nbl_list[mc_nbcells];//
	ivec3d nbl_list[];

	//***********************************************************************

    public:

       	// set methods definition

	    void set_cell_id(int uid);
        void set_glob_coord(ivec3d cell_gcoord);
        void set_loc_coord(ivec3d cell_lcoord);
        void add_neighbor(int id, ivec3d position);
        void add_sample(void);

        // if required add other possible neighbor set methods

        void add_particle(particle p)                   { cell.push_back(p); }
        void delete_particle(particle p, int index)    { cell.erase(cell.begin()+ index);}

	    // get methods

	    int get_cell_id()                       { return cell_id;       }
	    long get_nparticles()                  { return cell.size();    }
	    ivec3d get_cell_glob_coord()           { return gcoord;         }
        ivec3d get_cell_loc_coord()            { return lcoord;         }

        int get_sample_factor()                { return sample_fact;}

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

	    celltype cell_with_lcoord(ivec3d lcoord);

};

#endif /* MONTE_CLASSES_H_ */
