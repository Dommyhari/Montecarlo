#define compiler
CC = g++

#define compile options
#CFLAGS = -std=c++0x

#define path for additional headers
INCLUDE = -I/home/users/ganeshfx/workspace/Montecarlo_API_cpp/

#define archive flag
ARCHIVE = ar

#define archive option
AOPT = -rcs

#define archive source files
ASOURCE += Montecarlo_api.o	Monte_classes.o	Monte_methods.o

#all:	montecarlo_api

#Montecarlo_api.o:	Montecarlo_api.cpp
#    ${CC}	-c	Montecarlo_api.cpp
#	g++	-c	Montecarlo_api.cpp
	
Monte_classes.o:	Monte_classes.cpp  
	 g++ -std=c++0x -I/home/users/ganeshfx/workspace/Montecarlo_API_cpp/ -o Monte_classes.o Monte_classes.cpp
#    ${CC}	-std=c++0x	${INCLUDE}	-o	Monte_classes.o	Monte_classes.cpp    
#	g++	-std=gnu++0x	-c	Monte_classes.cpp		

#Monte_methods.o:	Monte_methods.cpp
#    ${CC}	${INCLUDE} -o Monte_methods.o Monte_classes.cpp
#	g++ -c 	Monte_methods.cpp
	
#libMontecarlo.a:	Montecarlo_api.o	Monte_classes.o	Monte_methods.o
#    ${ARCHIVE}	${AOPT}  Montecarlo_api.o	Monte_classes.o	Monte_methods.o
#	ar	-rcs	Montecarlo_api.o	Monte_classes.o	Monte_methods.o	 		
	
clean:	
	rm *.o 
	#montecarlo_api	libMontecarlo.a
