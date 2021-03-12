#include <stdio.h>
#include "c_dtypes.h"

// Forward declaration, in order to avoid compiler warnings, all functions below until
// init_fortran are implemented inside the file, src_for/dnamiF.for
void time_march_fortran(const int *param_int, const wp *param_float, wp *data_float);
void stored_fortran(const int *param_int, const wp *param_float, wp *data_float, int *type_st);
void filter_fortran(int *dir, const int *param_int, const wp *param_float, wp *data_float);

void pack_fortran(wp *buf, const wp *f,
		  int *ibeg, int *iend, 
		  int *jbeg, int *jend, 
		  int *kbeg, int *kend, 
		  int *sizex, int *sizey, int *sizez, 
		  int *sizenv);
void unpack_fortran(const wp *buf, wp *f,
		    int *ibeg, int *iend, 
		    int *jbeg, int *jend, 
		    int *kbeg, int *kend,
		    int *sizex,int *sizey,
		    int *sizez, 
		    int *sizenv);

void init_fortran(const int *param_int, const wp *param_float, wp *data_float);


// C-Layer functions, the following functions are called from Python 
void time_march(const int *param_int, const wp *param_float, wp *data_float){
	time_march_fortran(param_int,param_float,data_float);
}

void stored(const int *param_int, const wp *param_float, wp *data_float, int type_st){
	stored_fortran(param_int,param_float,data_float, &type_st);
}

void filter(int dir, const int *param_int, const wp *param_float, wp *data_float){
	filter_fortran(&dir, param_int,param_float,data_float);
}

void pack(wp *buf, const wp *f, int ibeg, int iend, int jbeg, int jend, int kbeg, int kend, int sizex, int sizey, int sizez, int sizenv){
	// Fortran only know pass by reference, so all parameters need to be passed by reference
	pack_fortran(buf, f, &ibeg, &iend, &jbeg, &jend, &kbeg, &kend, &sizex, &sizey, &sizez, &sizenv);
}

void unpack(const wp *buf, wp *f, int ibeg, int iend, int jbeg, int jend, int kbeg, int kend, int sizex, int sizey, int sizez, int sizenv){
	unpack_fortran(buf, f, &ibeg, &iend, &jbeg, &jend, &kbeg, &kend, &sizex, &sizey, &sizez, &sizenv);
}

void init(const int *param_int, const wp *param_float, wp *data_float){
	init_fortran(param_int, param_float, data_float);
}
