/*****************************************************************************************************
 *
 *                                 Generic Montecarlo API
 *  specifications to be included
 *  Created on: Sep 4, 2014
 *      Author: ganeshfx
 ****************************************************************************************************/

/*  MC routine */

#include<iostream>
#include<ctime>
#include<cmath>
#include<cstdlib>
#include<vector>
#include "Monte_classes.h" // class definitions
#include "Monte_globals.h" // method definitions with extern

using namespace std;

// container for data's obtained from MD simulation

//**********************************************
// NOTE: later if needed,include this into particle class
//**********************************************

vector<int> mc_atomtypes;         //  vector container for atom types
vector<long> mc_atomnumber;       //  vector container for atom number
vector<double> mc_atommass;      //   vector container for atom mass
vector<double> mc_positions;     //   vector container for atom positions


// some special containers for class-objects using vector-STL

//vector<particle> cell;          // container of particles
//vector<celltype> cell_list;     // container of cells
//vector<cellblock> cellarray;    // container of cellblock


//*************************
// NOTE: add signature to get ncpus,cpu dim,MPI_COMM and other essentials list from Monte_globals.h from MD
//
// cpu_gcoord
//
//************************

extern "C" void pack_config_to_montecarlo(int md_mc_pid,long md_mc_tatoms_cpu,long *md_mc_atomnumber,int *md_mc_atomtypes,
		double  *md_mc_atommass,double  *md_mc_positions){

	    long rand_no,m=0; // random number variable

	    mc_pid = md_mc_pid;
		mc_tatoms_cpu = md_mc_tatoms_cpu; // get total particles from MD

		for(long i=0; i<mc_tatoms_cpu; i++){

			// filling data container with corresponding values
			mc_atomnumber.push_back(md_mc_atomnumber[i]);
			mc_atomtypes.push_back(md_mc_atomtypes[i]);
			mc_atommass.push_back(md_mc_atommass[i]);

			mc_positions.push_back(md_mc_positions[m++]); // holds x
			mc_positions.push_back(md_mc_positions[m++]); // holds y
			mc_positions.push_back(md_mc_positions[m++]); // holds z

		}

		// check for unallocated vector
		if(mc_atomnumber.empty()||mc_atomtypes.empty()||mc_atommass.empty()||mc_positions.empty()){
			std::cerr<<"Montecarlo: problem allocating datas from MD "<<endl;
		}
		else if((mc_atomnumber.size()==mc_tatoms_cpu)&&(mc_atomtypes.size()==mc_tatoms_cpu)&&(mc_atommass.size()==mc_tatoms_cpu)
				&&((mc_positions.size()/3)==mc_tatoms_cpu)){

			std::cout<<" ------------------------------------------- "<<endl;
			std::cout<<" Process id : "<<mc_pid<<endl;
			std::cout<<" Total particles received : " <<mc_tatoms_cpu<<endl;
			std::cout<<"Montecarlo: datas are allocated correctly " <<endl;
			std::cout<<" ------------------------------------------- "<<endl;
	    }
		else {
			std::cerr<<"Montecarlo: unknown error with allocation" <<endl;
		}


}

extern "C" void do_montecarlo(int* md_pid,long *md_tatoms_cpu,long **md_atomnumber,int **md_atomtypes,double **md_atommass,double **md_positions)
{


    // local variables
	long rand_no;

	long m=0;

	long n=0;

	srand(time(NULL));
	rand_no=(long) rand()%mc_tatoms_cpu;


	// Data mirroring procedure - Assigning pointers

	md_pid=&mc_pid;
	md_tatoms_cpu = &mc_tatoms_cpu;

	// preliminary pointer assignment

	*md_atomnumber = mc_atomnumber.data(); //error part
    *md_atomtypes  = mc_atomtypes.data();
	*md_atommass   = mc_atommass.data();  //they do have initialized pointers and memory
    *md_positions  = mc_positions.data();




    // swapping procedure

    while(mc_atomtypes[rand_no] != 2){

    	std::cout<<" chosen particle is not a placeholder "<<endl;
    	rand_no=(long) rand()%mc_tatoms_cpu;

    }

    //mc_atomtypes[rand_no] = 1;
    mc_atomtypes.at(rand_no) = 1; // do bound check instead of direct memory access


    std::cout<<" **********************************************"<<endl;
    std::cout<<" Hello from do_montecarlo with process : "<<*md_pid<<endl;
    std::cout<<" **********************************************"<<endl;
    std::cout<<" some data check "<<endl;
    std::cout<<" created random no : "<<rand_no<<endl;
    std::cout<<" md_pid  :" <<*md_pid<<" mc_pid : "<<mc_pid<<endl;
    std::cout<<" md_tatoms_cpu..  :" << *md_tatoms_cpu <<" mc_tatoms_cpu : "<<mc_tatoms_cpu<<endl;
    std::cout<<" mc_atomnumber[rand_no]  :"<<mc_atomnumber[rand_no]<<endl;
    std::cout<<" mc_atomtypes[rand_no]  :"<<mc_atomtypes[rand_no]<<endl;
    std::cout<<" mc_atommass[rand_no]  :" <<mc_atommass[rand_no]<<endl;
    std::cout<<" pos_x : mc_positions[rand_no]  :"<<mc_positions[rand_no*3]<<endl;
    std::cout<<" pos_y : mc_positions[rand_no]  :"<<mc_positions[rand_no*3+1]<<endl;
    std::cout<<" pos_z : mc_positions[rand_no]  :"<<mc_positions[rand_no*3+2]<<endl;
    std::cout<<" ------------------------------------------- "<<endl;


}

extern "C" void clean_montecarlo(){

	// clear or delete all data structures created in MC scope

	std::cout<<" **********************************************"<<endl;
	std::cout<<"  Montecarlo: clearing data structures "<<endl;
	std::cout<<" **********************************************"<<endl;
	mc_atomnumber.clear();
	mc_atomtypes.clear();
	mc_atommass.clear();
	mc_positions.clear();

	std::cout<< "   Some Post check       " <<endl;
	std::cout<<  " size : mc_atomnumber  "<<mc_atomnumber.size()<<endl;
	std::cout<<  " size : mc_atomtypes  "<<mc_atomtypes.size()<<endl;
	std::cout<<  " size : mc_atommass  "<<mc_atommass.size()<<endl;
	std::cout<<  " size : mc_positions  "<<mc_positions.size()<<endl;
	std::cout<<" ------------------------------------------- "<<endl;
}
