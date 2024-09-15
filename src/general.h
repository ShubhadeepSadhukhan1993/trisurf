/* vim: set ts=4 sts=4 sw=4 noet : */
#ifndef _GENERAL_H
#define _GENERAL_H

#include<stdarg.h>
#include<stdio.h>
#include<gsl/gsl_complex.h>
#include <gsl/gsl_randist.h>
/* @brief This is a header file, defining general constants and structures.
  * @file header.h
  * @author Samo Penic
  * @date 5.3.2001
  * 
  * Header file for general inclusion in all the code, defining data structures
  * and general constans. All datatypes used in the code is also defined here.
  *
  * Miha: branch trisurf-polyel
  */

/* Defines */
/** @brief Return value of type bz_bool that indiceates successful function finish 
  *
  * Function usualy return some value, which are the result of certain operation. Functions that don't
  * return any parameters can return value, that indicates if the function call finished successfully.
  * In case of successful function run, the functions should return TS_SUCCESS to the caller. This define
  * is set here to get uniformity among all the functions used in program.
  *
  * Example of usage:
  *		ts_boot somefunction(ts_int param1, ....){
  *			...
  *			return TS_SUCCESS;
  *		}
  */
#define TS_SUCCESS 0

/** @brief Return value of type bz_bool that indicates unsuccessful function finish 
  *
  * Function usualy return some value, which are the result of certain operation. Functions that don't
  * return any parameters can return value, that indicates if the function call finished successfully.
  * In case of unsuccessful function run, the functions should return TS_FAIL to the caller. This define
  * is set here to get uniformity among all the functions used in program.
  *
  * Example of usage:
  *
  *		ts_boot somefunction(ts_int param1, ....){
  *			...
  *			return TS_FAIL;
  *		}
  */
#define TS_FAIL 1

/* CONSTANTS */

#define TS_ID_FILAMENT 1

/* DATA TYPES */
/** @brief Sets the default datatype for ts_double
 *
 * Requred for some functions to work, like "pow" from math.h. If ts_double is defined as
 * float program must run with "powf". Where type dependant function is used it checks this
 * define directive to decide which version to compile in. Available options
 *
 *	TS_DOUBLE_FLOAT
 *	TS_DOUBLE_DOUBLE
 *	TS_DOUBLE_LONGDOUBLE
*/
#define TS_DOUBLE_DOUBLE

/** For the purpose of greater flexibility all data types used in the program
 *  shouldn't use standard C types, but should use types defined here.
 *	ts_int (uses int)
 */
typedef int ts_int;
/** For the purpose of greater flexibility all data types used in the program
 *  shouldn't use standard C types, but should use types defined here.
 *	ts_uint (uses unsigned int)
 */
typedef unsigned int ts_uint;
/** For the purpose of greater flexibility all data types used in the program
 *  shouldn't use standard C types, but should use types defined here.
 *	ts_long (uses long)
 */
typedef long ts_long;
/** For the purpose of greater flexibility all data types used in the program
 *  shouldn't use standard C types, but should use types defined here.
 *	ts_ulong (uses unsigned long)
 */
typedef unsigned long ts_ulong;
/** For the purpose of greater flexibility all data types used in the program
 *  shouldn't use standard C types, but should use types defined here.
 *	ts_float (uses float)
 */
typedef float ts_float;
/** For the purpose of greater flexibility all data types used in the program
 *  shouldn't use standard C types, but should use types defined here.
 *	ts_double (uses double)
 */
typedef double ts_double;
/** For the purpose of greater flexibility all data types used in the program
 *  shouldn't use standard C types, but should use types defined here.
 *	ts_char (uses char)
 */
typedef char ts_char;
/** For the purpose of greater flexibility all data types used in the program
 *  shouldn't use standard C types, but should use types defined here.
 *	ts_uchar (uses unsigned char)
 */
typedef unsigned char ts_uchar;
/** For the purpose of greater flexibility all data types used in the program
 *  shouldn't use standard C types, but should use types defined here.
 *	ts_bool (uses char)
 */
typedef char ts_bool;


/* STRUCTURES */


/** @brief Data structure for keeping the coordinates in selected coordinate
 * system
 */
#define TS_COORD_CARTESIAN 0
#define TS_COORD_SPHERICAL 1
#define TS_COORD_CYLINDRICAL 2

typedef struct {
    ts_double e1;
    ts_double e2;
    ts_double e3;
    ts_uint coord_type;
} ts_coord;

/** @brief Data structure of all data connected to a vertex
 *
 *  ts_vertex holds the data for one single point (bead, vertex). To understand how to use it
 *  here is a detailed description of the fields in the data structure. */
struct ts_vertex {
        ts_uint idx;
        ts_double x; /**< The x coordinate of vertex. */
        ts_double y; /**< The y coordinate of vertex. */
        ts_double z; /**< The z coordinate of vertex. */
        ts_uint neigh_no; /**< The number of neighbours. */
        struct ts_vertex **neigh; /**< The pointer that holds neigh_no pointers to this structure. */
        ts_double *bond_length; /**< Obsolete! The bond lenght is moved to ts_bond */
        ts_double *bond_length_dual; /**< Obsolete! Bond length in dual lattice is moved to ts_bond! */
        ts_double curvature;
        ts_double energy;
        ts_double energy_h;
        ts_uint tristar_no;
        struct ts_triangle **tristar; /**< The list of triangles this vertex belongs to. This is an array of pointers to ts_triangle structure of tristar_no length */
        ts_uint bond_no;
        struct ts_bond **bond; /**< Array of pointers of lenght bond_no that stores information on bonds. */
        struct ts_cell *cell; /**< Which cell do we belong to? */
        ts_double xk;
        ts_double c;
        ts_uint id;
        ts_double projArea;
        ts_double relR;
        ts_double solAngle;


        
		struct ts_poly *grafted_poly;
		struct ts_cluster *cluster;
        /*Shubhadeep */

        ts_double Fx; 
        ts_double Fy;
        ts_double Fz; 

        
        ts_double Fshx; 
        ts_double Fshy; 
        ts_double Fshz; 

        ts_double Factx;
        ts_double Facty; 
        ts_double Factz; 

        ts_double proj;
        ts_double cx;
        

        ts_double Fbx;
        ts_double Fby;
        ts_double Fbz;

        ts_double F_vtx;
        /*Shubhadeep */
};

typedef struct ts_vertex ts_vertex;

typedef struct {
    ts_uint n;
    ts_vertex **vtx;

} ts_vertex_list;

struct ts_bond {
    	ts_uint idx;
	ts_vertex *vtx1;
	ts_vertex *vtx2;
    	ts_double bond_length;
    	ts_double bond_length_dual;
	ts_bool tainted; //TODO: remove
	ts_double energy;
	ts_double x,y,z;
};
typedef struct ts_bond ts_bond;

struct ts_bond_list {
    ts_uint n;
    ts_bond **bond;
};
typedef struct ts_bond_list ts_bond_list;

struct ts_triangle {
    ts_uint idx;
	ts_vertex *vertex[3];
	ts_uint neigh_no;
	struct ts_triangle **neigh;
	ts_double xnorm;
	ts_double ynorm;
	ts_double znorm;
    ts_double area; // firstly needed for sh.c
    ts_double volume; // firstly needed for sh.c
	ts_double energy;
};
typedef struct ts_triangle ts_triangle;

struct ts_triangle_list{
    ts_uint n;
	ts_double a0;
    ts_triangle **tria;
};
typedef struct ts_triangle_list ts_triangle_list;


typedef struct ts_cell {
    ts_uint idx;
    ts_vertex **vertex;
    ts_uint nvertex;
} ts_cell; 

typedef struct ts_cell_list{
    ts_uint ncmax[3];

    ts_uint cellno;
    ts_cell **cell;
    ts_double dcell;
    ts_double shift;
    ts_uint max_occupancy;
    ts_double c_cm[3];
	ts_double dmin_interspecies;
} ts_cell_list;


typedef struct {
    ts_uint l;
    ts_double **ulm;
    gsl_complex **ulmComplex;
    ts_double **sumUlm2;
    ts_uint N;
    ts_double **co;
    ts_double ***Ylmi;
} ts_spharm;



struct ts_poly {
	ts_vertex_list *vlist;
	ts_bond_list *blist;
	ts_vertex *grafted_vtx;
	ts_double k;
};
typedef struct ts_poly ts_poly;


struct ts_poly_list {
	ts_uint	n;
	ts_poly **poly;
};
typedef struct ts_poly_list ts_poly_list;


typedef struct{
	ts_float z_max;
	ts_float z_min;
	ts_int force_switch;
} ts_confinement_plane;


typedef struct {
	long int nshell;
	long int ncxmax;
	long int ncymax;
	long int nczmax;
	long int npoly;
	long int nmono;
	long int internal_poly;
	long int nfil;
	long int nfono;
	long int R_nucleus;
	ts_double R_nucleusX;
	ts_double R_nucleusY;
	ts_double R_nucleusZ;
	long int pswitch;
    long int constvolswitch;
    long int constareaswitch;
	long int stretchswitch;
	ts_double xkA0;
    ts_double constvolprecision;
    	char *multiprocessing;
   	long int brezveze0;
    	long int brezveze1;
    	long int brezveze2;
    	ts_double xk0;
	ts_double dmax;
	ts_double dmin_interspecies;
	ts_double stepsize;
	ts_double kspring;
	ts_double xi;
	ts_double pressure;
	long int iterations;
	long int inititer;
	long int mcsweeps;
	long int quiet;
	long int shc;
	long int number_of_vertices_with_c0;
	ts_double c0;
	ts_double w;
	ts_double F;
	long int plane_confinement_switch;
	ts_double plane_d;
	ts_double plane_F;
	long int type_of_adhesion_model;
	long int allow_xy_plane_movement;
	long int force_balance_along_z_axis;
	long int adhesion_switch;
	ts_double adhesion_cuttoff;
	ts_double adhesion_strength;
	ts_double z_adhesion;
	ts_double adhesion_radius;
    /* Shubhadeep */
    ts_double a;
    ts_double b;
    long int shear_switch;
    ts_double nx;
    ts_double ny;
    /* Shubhadeep */
    long int F_noise_switch;
    ts_double F_noise_SD;


    ts_double lamda1;
    ts_double lamda2;
    ts_double D;
    long int inhibition_switch;
    //long int inhibition_interval;
    ts_double beta;

    long int wall_switch;
    ts_double wallx;
    ts_double wally;
    ts_double wall_theta;
    ts_double cs;
    ts_double conc0;

    ts_double Fblow;
    long int blow_switch;
    ts_double wall_spring;
    ts_double adhesion_spring;

    long int box_confinement_switch;
    long int confinement_adhesion_switch;
    ts_double box_Lz;
    ts_double box_Lx1;
    ts_double box_Lx2;
    ts_double box_Ly1;
    ts_double box_Ly2;


    long int notch_switch;
    ts_double notch_angle;
    ts_double notch_range;
    ts_double notch_position;

    ts_double adhesion_strength2;
    ts_double patchx;
    ts_double patchy;
    ts_double patch_size;
    
    long int F_noise_interval;

} ts_tape;




typedef struct {
	ts_vertex_list *vlist;
	ts_bond_list *blist;
	ts_triangle_list *tlist;
	ts_cell_list *clist;
	ts_uint nshell;
	ts_double bending_rigidity;
	ts_double dmax;
	ts_double stepsize;
   	ts_double cm[3];
	ts_double volume;
	ts_spharm *sphHarmonics;
// Polymers outside the vesicle and attached to the vesicle membrane (polymer brush):
	ts_poly_list *poly_list;
// Filaments inside the vesicle (not attached to the vesicel membrane:
	ts_poly_list *filament_list;

	ts_double spring_constant;
	ts_double pressure;
	ts_int pswitch;
 	ts_tape *tape;
	ts_double R_nucleus;
	ts_double R_nucleusX;
	ts_double R_nucleusY;
	ts_double R_nucleusZ;
	ts_double nucleus_center[3];
	ts_double area;
	ts_confinement_plane confinement_plane;
	ts_double adhesion_center;
    /*Shubhadeep*/
   
    ts_double xnorm;
    ts_double ynorm;
    ts_double znorm;
    ts_double proj_min;
    ts_double proj_max;
    ts_double vmag;
    ts_double xback;
    /*Shubhadeep*/
} ts_vesicle;



struct ts_cluster{
	ts_uint nvtx;
	ts_uint idx;
	ts_vertex **vtx;
  ts_double x;
  ts_double y;
  ts_double z;
  ts_double Factx;
  ts_double Facty;
  ts_double Factz;

};

typedef struct ts_cluster ts_cluster;

typedef struct{
	ts_uint n;
	ts_cluster **cluster;
} ts_cluster_list;


/* GLOBAL VARIABLES */

extern int quiet;
extern ts_double V0;
extern ts_double A0;
extern ts_double epsvol;
extern ts_double epsarea;


extern const gsl_rng_type * TT;
extern gsl_rng * rr_noise;
/* FUNCTIONS */

/** Non-fatal error function handler:
 *      @param text is a description of an error
 *      @returns doesn't return anything
*/
void err(char *text);

/** Fatal error function handler:
 *      @param text is a description of an error
 *      @param errcode is a (non-zero) error code
 *      @returns terminates the execution of program with errcode set
*/
void fatal(char *text, ts_int errcode);

ts_uint ts_fprintf(FILE *fd, char *fmt, ...);

#define VTX(n) &(vlist->vtx[n])
#define VTX_DATA(n) vlist->vtx[n].data


/* FOR PID GENERATION ROUTINE */
#define CPF_CLOEXEC 1

int createPidFile(const char *progName, const char *pidFile, int flags);

int lockRegion(int fd, int type, int whence, int start, int len);
char *libVersion();



#endif
