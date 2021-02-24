#include <stdio.h>

// Forward declaration, in order to avoid compiler warnings, all functions below until
// init_fortran are implemented inside the file, src_for/dnamiF.for
void time_march_fortran(const int *param_int, const double *param_float, double *data_float);
void stored_fortran(const int *param_int, const double *param_float, double *data_float, int *type_st);
void filter_fortran(int *dir, const int *param_int, const double *param_float, double *data_float);

void pack_fortran(double *buf, const double *f, 
		  int *ibeg, int *iend, 
		  int *jbeg, int *jend, 
		  int *kbeg, int *kend, 
		  int *sizex, int *sizey, int *sizez, 
		  int *sizenv);
void unpack_fortran(const double *buf, double *f, 
		    int *ibeg, int *iend, 
		    int *jbeg, int *jend, 
		    int *kbeg, int *kend,
		    int *sizex,int *sizey,
		    int *sizez, 
		    int *sizenv);

void init_fortran(const int *param_int, const double *param_float, double *data_float);


// C-Layer functions, the following functions are called from Python 
void time_march(const int *param_int, const double *param_float, double *data_float){
	time_march_fortran(param_int,param_float,data_float);
}

void stored(const int *param_int, const double *param_float, double *data_float, int type_st){
	stored_fortran(param_int,param_float,data_float, &type_st);
}

void filter(int dir, const int *param_int, const double *param_float, double *data_float){
	filter_fortran(&dir, param_int,param_float,data_float);
}

void pack(double *buf, const double *f, int ibeg, int iend, int jbeg, int jend, int kbeg, int kend, int sizex, int sizey, int sizez, int sizenv){
	// Fortran only know pass by reference, so all parameters need to be passed by reference
	pack_fortran(buf, f, &ibeg, &iend, &jbeg, &jend, &kbeg, &kend, &sizex, &sizey, &sizez, &sizenv);
}

void unpack(const double *buf, double *f, int ibeg, int iend, int jbeg, int jend, int kbeg, int kend, int sizex, int sizey, int sizez, int sizenv){
	unpack_fortran(buf, f, &ibeg, &iend, &jbeg, &jend, &kbeg, &kend, &sizex, &sizey, &sizez, &sizenv);
}

void init(const int *param_int, const double *param_float, double *data_float){
	init_fortran(param_int, param_float, data_float);
}
