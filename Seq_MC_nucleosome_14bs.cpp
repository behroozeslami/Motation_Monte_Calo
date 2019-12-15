#include <iostream>
#include <math.h>
#include <time.h>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <stdio.h>
#include <string.h>

//GSL
#include <gsl/gsl_math.h> 
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>  

//Random generator includes
#include "randomc.h"
#include "mersenne.cpp"
#include "userintf.cpp"

// make directory
#include <sys/stat.h>

using namespace std;

//Random generator object
int seed = (int)(time(0));
CRandomMersenne RanGen(seed);

// Initial Values and definitions

const double PI = 4.0*atan(1.0);

// Variables Related to DNA
const int N_bp0 = 147;
int N_bp = 147;
int length_tail1, length_tail2;

static gsl_vector* Seq_Code = gsl_vector_alloc(N_bp-1);
static gsl_vector* bp_Seq_Code = gsl_vector_alloc(N_bp);

const double a = -0.3;
const double b = 0.9;
const double c = 0.05;
const double l_0 = 0.34;
const double L = (N_bp-1)*l_0;
const double omega_0 = 1.8;
const double Tw_0 = omega_0*l_0;
const double beta = 3.0;
const double beta_seq = 3.0;

static gsl_matrix* super_Q = gsl_matrix_alloc(102,6);
static gsl_vector* super_q0 = gsl_vector_alloc(102);

static gsl_matrix* Q = gsl_matrix_alloc(6,6);
static gsl_vector* q0 = gsl_vector_alloc(6);

static gsl_matrix* d1 = gsl_matrix_alloc(3,N_bp);
static gsl_matrix* d2 = gsl_matrix_alloc(3,N_bp);
static gsl_matrix* d3 = gsl_matrix_alloc(3,N_bp);
static gsl_matrix* R = gsl_matrix_alloc(3,N_bp);

static gsl_matrix* dm1 = gsl_matrix_alloc(3,N_bp-1);
static gsl_matrix* dm2 = gsl_matrix_alloc(3,N_bp-1);
static gsl_matrix* dm3 = gsl_matrix_alloc(3,N_bp-1);
static gsl_matrix* Rm = gsl_matrix_alloc(3,N_bp-1);

static gsl_matrix* init_dm1 = gsl_matrix_alloc(3,N_bp-1);
static gsl_matrix* init_dm2 = gsl_matrix_alloc(3,N_bp-1);
static gsl_matrix* init_dm3 = gsl_matrix_alloc(3,N_bp-1);
static gsl_matrix* init_Rm = gsl_matrix_alloc(3,N_bp-1);

static gsl_matrix* d1_0 = gsl_matrix_alloc(3,N_bp0);
static gsl_matrix* d2_0 = gsl_matrix_alloc(3,N_bp0);
static gsl_matrix* d3_0 = gsl_matrix_alloc(3,N_bp0);
static gsl_matrix* R_0 = gsl_matrix_alloc(3,N_bp0);

static gsl_vector* X_Axis = gsl_vector_alloc(3);
static gsl_vector* Y_Axis = gsl_vector_alloc(3);
static gsl_vector* Z_Axis = gsl_vector_alloc(3);

static gsl_vector* X0_Axis = gsl_vector_alloc(3);
static gsl_vector* Y0_Axis = gsl_vector_alloc(3);
static gsl_vector* Z0_Axis = gsl_vector_alloc(3);

static gsl_vector* bond_array = gsl_vector_alloc(N_bp-1);

static gsl_vector* bond_steps = gsl_vector_alloc(N_bp-1);
static gsl_vector* strand_3prime_array = gsl_vector_alloc(N_bp-1);
static gsl_vector* Pos_cons = gsl_vector_alloc(N_bp-1);
static gsl_vector* Ang_cons = gsl_vector_alloc(N_bp-1);
static gsl_vector* free_wrapped_bp = gsl_vector_alloc(N_bp);
static gsl_vector* free_tail1_bp = gsl_vector_alloc(N_bp);
static gsl_vector* free_tail2_bp = gsl_vector_alloc(N_bp);

double bond_steps_number, free_wrapped_bp_number, free_tail1_bp_number, free_tail2_bp_number;

double Tot_Energy;

// Variables Related to MC
int Error_Code;

double Try_number;
double Try_number_loc;
double Try_number_pivot;
double Try_number_bond;
double Try_number_glob;
double Try_number_DisPhos;
double Try_number_Mutation;

double Accept_number;
double Accept_number_loc;
double Accept_number_pivot;
double Accept_number_bond;
double Accept_number_glob;
double Accept_number_DisPhos;
double Accept_number_Mutation;

double Accept_ratio;
double Accept_ratio_loc;
double Accept_ratio_pivot;
double Accept_ratio_bond;
double Accept_ratio_glob;
double Accept_ratio_DisPhos;
double Accept_ratio_Mutation;

int N_EQ_MCstep;
int N_Main_MCstep;

double max_rot;
double max_trans;

double Loc_Move_tail1_Percent;
double Loc_Move_wrapped_Percent;
double Loc_Move_tail2_Percent;
double Pivot_Move_tail1_Percent;
double Pivot_Move_tail2_Percent;
double Glob_Move_Percent;
double Bond_Move_Percent;
double DisPhos_Move_Percent;
double Mut_Move_Percent;
double Pivot_vs_Loc_frac;
double Bond_vs_DisPhos_frac;

// Histograms and Averages
double sum_Energy, sum_sq_Energy, ave_Energy, sigma_Energy;
int Energy_num;

double sum_EndtoEnd, sum_sq_EndtoEnd, ave_EndtoEnd, sigma_EndtoEnd;
int EndtoEnd_num;

double sum_ang_1, sum_sq_ang_1, ave_ang_1, sigma_ang_1;
int ang_1_num;

double sum_ang_2, sum_sq_ang_2, ave_ang_2, sigma_ang_2;
int ang_2_num;

double sum_ang_3, sum_sq_ang_3, ave_ang_3, sigma_ang_3;
int ang_3_num;

double min_Energy =-600.0;
double max_Energy = 1400.0;
double hist_Energy_delta = 1;
int hist_Energy_num = int((max_Energy-min_Energy)/hist_Energy_delta)+1;

double min_EndtoEnd =-1.2;
double max_EndtoEnd = 1.2;
double hist_EndtoEnd_delta = 0.001;
int hist_EndtoEnd_num = int((max_EndtoEnd-min_EndtoEnd)/hist_EndtoEnd_delta)+1;

double min_ang_1 =-1.0;
double max_ang_1 = 1.0;
double hist_ang_1_delta = 0.01;
int hist_ang_1_num = int((max_ang_1-min_ang_1)/hist_ang_1_delta)+1;

double min_ang_2 =-1.0;
double max_ang_2 = 1.0;
double hist_ang_2_delta = 0.01;
int hist_ang_2_num = int((max_ang_2-min_ang_2)/hist_ang_2_delta)+1;

double min_ang_3 =-1.0;
double max_ang_3 = 1.0;
double hist_ang_3_delta = 0.01;
int hist_ang_3_num = int((max_ang_3-min_ang_3)/hist_ang_3_delta)+1;
	
static gsl_vector* Energy_array = gsl_vector_alloc(hist_Energy_num);
static gsl_vector* hist_Energy = gsl_vector_alloc(hist_Energy_num);

static gsl_vector* EndtoEnd_array = gsl_vector_alloc(hist_EndtoEnd_num);
static gsl_vector* hist_EndtoEnd = gsl_vector_alloc(hist_EndtoEnd_num);

static gsl_vector* ang_1_array = gsl_vector_alloc(hist_ang_1_num);
static gsl_vector* hist_ang_1= gsl_vector_alloc(hist_ang_1_num);

static gsl_vector* ang_2_array = gsl_vector_alloc(hist_ang_2_num);
static gsl_vector* hist_ang_2= gsl_vector_alloc(hist_ang_2_num);

static gsl_vector* ang_3_array = gsl_vector_alloc(hist_ang_3_num);
static gsl_vector* hist_ang_3= gsl_vector_alloc(hist_ang_3_num);

static gsl_matrix* hist_Seq = gsl_matrix_alloc(N_bp-1, 16);

// Strings
string str_d1, str_d2, str_d3, str_R, str_Ax, str_En, str_Ac, str_bp_Seq, str_step;

string str_h_Energy, str_ave_Energy, str_h_EndtoEnd, str_ave_EndtoEnd, str_h_ang1, str_ave_ang1, str_h_ang2, str_ave_ang2, str_h_ang3, str_ave_ang3, str_h_Seq;

// Function dot_prod
double dot_prod(gsl_vector* V1 , gsl_vector* V2)
{
	double dot = 0.0;
	
	for (int j=0 ; j<3 ; j++)
	{
		dot=dot+gsl_vector_get(V1,j)*gsl_vector_get(V2,j);
	}
	return (dot);
}

// Function dot_prod_6
double dot_prod_6 (gsl_vector* V1 , gsl_vector* V2)
{
	double dot = 0.0;
	
	for (int j=0 ; j<6 ; j++)
	{
		dot=dot+gsl_vector_get(V1,j)*gsl_vector_get(V2,j);
	}
	return (dot);
}

// Function Rot_MAT
void Rot_MAT(double Th_1, double Th_2, double Th_3, gsl_matrix* M)
{
	gsl_matrix_set_identity(M);
	
	double Th , comp;
	
	Th = sqrt(pow(Th_1, 2)+pow(Th_2, 2)+pow(Th_3, 2));
	
	if (Th > 0.00000000001)
	{
		comp = cos(Th)+(1-cos(Th)) * pow((Th_1/Th), 2);
		gsl_matrix_set(M,0,0,comp);
	
		comp =(1-cos(Th)) * (Th_1*Th_2)/pow(Th, 2) - sin(Th) *Th_3/Th;
		gsl_matrix_set(M,0,1,comp);
	
		comp =(1-cos(Th)) * (Th_1*Th_3)/pow(Th, 2) + sin(Th) *Th_2/Th;
		gsl_matrix_set(M,0,2,comp);
	
		comp =(1-cos(Th)) * (Th_1*Th_2)/pow(Th, 2) + sin(Th) *Th_3/Th;
		gsl_matrix_set(M,1,0,comp);
	
		comp = cos(Th)+(1-cos(Th)) * pow((Th_2/Th), 2);
		gsl_matrix_set(M,1,1,comp);
	
		comp =(1-cos(Th)) * (Th_2*Th_3)/pow(Th, 2) - sin(Th) *Th_1/Th;
		gsl_matrix_set(M,1,2,comp);
	
		comp =(1-cos(Th)) * (Th_1*Th_3)/pow(Th, 2) - sin(Th) *Th_2/Th;
		gsl_matrix_set(M,2,0,comp);
	
		comp =(1-cos(Th)) * (Th_2*Th_3)/pow(Th, 2) + sin(Th) *Th_1/Th;
		gsl_matrix_set(M,2,1,comp);
	
		comp = cos(Th)+(1-cos(Th)) * pow((Th_3/Th), 2);
		gsl_matrix_set(M,2,2,comp);
	}
}

//Function Normalize_Vec
void Normalize_Vec (gsl_vector* V)
{
	double dot, comp;
	
	dot = sqrt(dot_prod(V , V));
	
	for (int i=0; i<3; i++)
	{
		comp = (1.0/dot)*gsl_vector_get(V,i);
		gsl_vector_set(V,i,comp);
	}
}

// Function Make_SO3_Matrix
void Make_SO3_Matrix (gsl_vector* V1, gsl_vector* V2, gsl_vector* V3, gsl_matrix* Mat)
{
	Normalize_Vec (V1);
	Normalize_Vec (V2);
	Normalize_Vec (V3);
	
	for (int i=0; i<3; i++)
	{	
		gsl_matrix_set(Mat,i,0,gsl_vector_get(V1,i));
		gsl_matrix_set(Mat,i,1,gsl_vector_get(V2,i));
		gsl_matrix_set(Mat,i,2,gsl_vector_get(V3,i));
		
	}
}

// Function DOF
void DOF (gsl_vector* V1_1, gsl_vector* V2_1, gsl_vector* V3_1, gsl_vector* VR_1, gsl_vector* V1_2, gsl_vector* V2_2, gsl_vector* V3_2, gsl_vector* VR_2, gsl_vector* V1_m, gsl_vector* V2_m, gsl_vector* V3_m, gsl_vector* VR_m, gsl_vector* Theta_vec, gsl_vector* Rho_vec)
{
	gsl_matrix* Rot_1 = gsl_matrix_alloc(3,3);
	gsl_matrix* Rot_2 = gsl_matrix_alloc(3,3);
	gsl_matrix* Rot_rel = gsl_matrix_alloc(3,3);
	gsl_matrix* Rot_mid = gsl_matrix_alloc(3,3);
	
	gsl_vector* half_theta = gsl_vector_alloc(3);
	gsl_vector* r_rel = gsl_vector_alloc(3);
	
	double rot_ang, frac;
	
	Make_SO3_Matrix (V1_1, V2_1, V3_1, Rot_1);
	Make_SO3_Matrix (V1_2, V2_2, V3_2, Rot_2);
	
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0, Rot_1, Rot_2, 0.0, Rot_rel);
	
	gsl_vector_set_zero(Theta_vec);
	
	rot_ang = acos(0.5*(gsl_matrix_get(Rot_rel,0,0)+gsl_matrix_get(Rot_rel,1,1)+gsl_matrix_get(Rot_rel,2,2)-1.0) );

	if (rot_ang > PI-0.00000000001)
	{
		Error_Code=1;
	}
	
	if (sin(rot_ang) > 0.00000000001)
	{
		frac = 0.5*rot_ang / sin(rot_ang);
		gsl_vector_set(Theta_vec,0, frac*(gsl_matrix_get(Rot_rel,2,1)-gsl_matrix_get(Rot_rel,1,2)) );
		gsl_vector_set(Theta_vec,1, frac*(gsl_matrix_get(Rot_rel,0,2)-gsl_matrix_get(Rot_rel,2,0)) );
		gsl_vector_set(Theta_vec,2, frac*(gsl_matrix_get(Rot_rel,1,0)-gsl_matrix_get(Rot_rel,0,1)) );
		
	}
	
	double Tilt, Roll, Twist;
    
    Tilt = gsl_vector_get(Theta_vec,0);
    Roll = gsl_vector_get(Theta_vec,1);
    Twist = gsl_vector_get(Theta_vec,2);
    
    double comp;
    
    for (int i=0; i<3; i++)
    {
		comp = 0.5*Tilt*gsl_vector_get(V1_1,i) + 0.5*Roll*gsl_vector_get(V2_1,i) + 0.5*Twist*gsl_vector_get(V3_1,i);
		gsl_vector_set(half_theta,i,comp);
	}
    
    double half_th1, half_th2, half_th3;
    
    half_th1 = gsl_vector_get(half_theta,0);
    half_th2 = gsl_vector_get(half_theta,1);
    half_th3 = gsl_vector_get(half_theta,2);
    
    Rot_MAT(half_th1, half_th2, half_th3, Rot_mid);
    
    Normalize_Vec (V1_1);
	Normalize_Vec (V2_1);
	Normalize_Vec (V3_1);
    
    gsl_blas_dgemv (CblasNoTrans, 1.0, Rot_mid, V1_1, 0.0, V1_m);
    gsl_blas_dgemv (CblasNoTrans, 1.0, Rot_mid, V2_1, 0.0, V2_m);
    gsl_blas_dgemv (CblasNoTrans, 1.0, Rot_mid, V3_1, 0.0, V3_m);
	
	for (int i=0; i<3; i++)
	{
		gsl_vector_set(VR_m,i , 0.5*(gsl_vector_get(VR_2,i)+gsl_vector_get(VR_1,i)));
		
		gsl_vector_set(r_rel,i , gsl_vector_get(VR_2,i)-gsl_vector_get(VR_1,i));
	} 
	
	double Shift, Slide, Rise;
	
	Shift = dot_prod(r_rel , V1_m);
	Slide = dot_prod(r_rel , V2_m);
	Rise = dot_prod(r_rel , V3_m);
	
	gsl_vector_set(Rho_vec, 0, Shift);
	gsl_vector_set(Rho_vec, 1, Slide);
	gsl_vector_set(Rho_vec, 2, Rise);
	
	gsl_matrix_free(Rot_1);
	gsl_matrix_free(Rot_2);
	gsl_matrix_free(Rot_rel);
	gsl_matrix_free(Rot_mid);	
	
	gsl_vector_free(r_rel);
	gsl_vector_free(half_theta);
	

}

// Function Load_Elastic_Parameters
void Load_Elastic_Parameters ()
{
	double var;
	
	ifstream stiffness_file;
	ifstream equilibrium_file;
	
	stiffness_file.open ("All_Stif_Mat_inc_hom - Curve.txt");
	
	for (int k = 0; k < 102; k++)
	{
		for (int i = 0; i < 6; i++)
		{
			stiffness_file >> var;
			gsl_matrix_set(super_Q, k, i, var);
		}
	}
	
	stiffness_file.close();
	
	equilibrium_file.open ("All_Eq_Param_inc_hom - Curve.txt");
	
	for (int k = 0; k < 102; k++)
	{
		equilibrium_file >> var;
		gsl_vector_set(super_q0, k, var);
	}
	
	equilibrium_file.close();
}

// Function init_bp_Seq_Code
void init_bp_Seq_Code ()
{
	double r;
	int bp_Code;
	
	for (int bp = 0; bp < N_bp; bp++)
	{
		r = 1 + 4*RanGen.Random();
		bp_Code = int(r);
		
		gsl_vector_set(bp_Seq_Code, bp, bp_Code);
	}
}

// Function bp_Seq_Code_to_Seq_Code
void bp_Seq_Code_to_Seq_Code (int bps_1, int bps_2, gsl_vector* bpsc, gsl_vector* sc)
{
	int bp_code_1, bp_code_2, bps_code;
	
	for (int bps = bps_1; bps < (bps_2+1); bps++)
	{
		bp_code_1 = gsl_vector_get(bpsc, bps);
		bp_code_2 = gsl_vector_get(bpsc, bps+1);
		
		bps_code = 4*(bp_code_1-1) + (bp_code_2-1);
		
		gsl_vector_set(sc, bps, bps_code);
	} 
}

// Function Read_Elastic_Parameters
void Read_Elastic_Parameters (int bp_step, gsl_vector* sc)
{
	int bp_code;
	
	bp_code = gsl_vector_get(sc, bp_step);
	
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			gsl_matrix_set(Q, i, j, gsl_matrix_get(super_Q, (i+(bp_code*6)), j));
		}
		
		gsl_vector_set(q0, i, gsl_vector_get(super_q0, (i+(bp_code*6))));
	}
}

// Function Elastic_Energy
double Elastic_Energy (gsl_matrix* d1_mat, gsl_matrix* d2_mat, gsl_matrix* d3_mat, gsl_matrix* R_mat, gsl_vector* Theta_vec, gsl_vector* Rho_vec, int bp_step, gsl_vector* sc)
{
	gsl_vector* Vec1_1 = gsl_vector_alloc(3);
	gsl_vector* Vec2_1 = gsl_vector_alloc(3);
	gsl_vector* Vec3_1 = gsl_vector_alloc(3);
	gsl_vector* RVec_1 = gsl_vector_alloc(3);
	
	gsl_vector* Vec1_2 = gsl_vector_alloc(3);
	gsl_vector* Vec2_2 = gsl_vector_alloc(3);
	gsl_vector* Vec3_2 = gsl_vector_alloc(3);
	gsl_vector* RVec_2 = gsl_vector_alloc(3);
	
	gsl_vector* Vec1_m = gsl_vector_alloc(3);
	gsl_vector* Vec2_m = gsl_vector_alloc(3);
	gsl_vector* Vec3_m = gsl_vector_alloc(3);
	gsl_vector* RVec_m = gsl_vector_alloc(3);
	
	gsl_vector* del_q = gsl_vector_alloc(6);
	gsl_vector* Q_del_q = gsl_vector_alloc(6);
	
	double E;
    
    	for (int i = 0; i < 3; i++)
	{
		gsl_vector_set(Vec1_1,i,gsl_matrix_get(d1_mat,i,bp_step));
		gsl_vector_set(Vec2_1,i,gsl_matrix_get(d2_mat,i,bp_step));
		gsl_vector_set(Vec3_1,i,gsl_matrix_get(d3_mat,i,bp_step));
		
		gsl_vector_set(RVec_1,i,gsl_matrix_get(R_mat,i,bp_step));
		
	    gsl_vector_set(Vec1_2,i,gsl_matrix_get(d1_mat,i,bp_step+1));
		gsl_vector_set(Vec2_2,i,gsl_matrix_get(d2_mat,i,bp_step+1));
		gsl_vector_set(Vec3_2,i,gsl_matrix_get(d3_mat,i,bp_step+1));
		
		gsl_vector_set(RVec_2,i,gsl_matrix_get(R_mat,i,bp_step+1));
	}
	
    DOF (Vec1_1 , Vec2_1, Vec3_1 , RVec_1, Vec1_2, Vec2_2, Vec3_2, RVec_2, Vec1_m, Vec2_m, Vec3_m, RVec_m, Theta_vec, Rho_vec);
    
    Read_Elastic_Parameters (bp_step, sc);
    
    for (int i=0; i<3; i++)
    {
		gsl_vector_set(del_q,i, gsl_vector_get(Theta_vec,i));
	}
	
	for (int i=3; i<6; i++)
    {
		gsl_vector_set(del_q,i, gsl_vector_get(Rho_vec,i-3));
	}
    
    gsl_vector_sub (del_q,q0);
    
    gsl_blas_dgemv (CblasNoTrans, 1.0,Q, del_q, 0.0, Q_del_q);
    
    E = 0.5 * dot_prod_6(del_q , Q_del_q);
    
    gsl_vector_free(Vec1_1);
	gsl_vector_free(Vec2_1);
	gsl_vector_free(Vec3_1);
	gsl_vector_free(RVec_1);
	
	gsl_vector_free(Vec1_2);
	gsl_vector_free(Vec2_2);
	gsl_vector_free(Vec3_2);
	gsl_vector_free(RVec_2);
    
    gsl_vector_free(Vec1_m);
	gsl_vector_free(Vec2_m);
	gsl_vector_free(Vec3_m);
	gsl_vector_free(RVec_m);
	
	gsl_vector_free(del_q);
	gsl_vector_free(Q_del_q);
	
    return (E);
}
// Function Spring_Energy
double Spring_Energy (gsl_matrix* d1_mid, gsl_matrix* d2_mid, gsl_matrix* d3_mid, gsl_matrix* R_mid, int bp_step, int strand_3prime, double Ks_1, double Ks_2)
{
	gsl_vector* V01_m = gsl_vector_alloc(3);
    gsl_vector* V02_m = gsl_vector_alloc(3);
    gsl_vector* V03_m = gsl_vector_alloc(3);
    gsl_vector* V0R_m = gsl_vector_alloc(3);
    
    gsl_vector* V1_m = gsl_vector_alloc(3);
    gsl_vector* V2_m = gsl_vector_alloc(3);
    gsl_vector* V3_m = gsl_vector_alloc(3);
    gsl_vector* VR_m = gsl_vector_alloc(3);
	
	gsl_matrix* Rot_1 = gsl_matrix_alloc(3,3);
	gsl_matrix* Rot_2 = gsl_matrix_alloc(3,3);
	gsl_matrix* Rot_rel = gsl_matrix_alloc(3,3);
	
	gsl_vector* r_rel = gsl_vector_alloc(3);
	
	for (int i = 0; i < 3; i++)
	{	
		gsl_vector_set(V01_m,i,gsl_matrix_get(init_dm1,i,bp_step));
		gsl_vector_set(V02_m,i,gsl_matrix_get(init_dm2,i,bp_step));
		gsl_vector_set(V03_m,i,gsl_matrix_get(init_dm3,i,bp_step));
		gsl_vector_set(V0R_m,i,gsl_matrix_get(init_Rm,i,bp_step));
		
		gsl_vector_set(V1_m,i,gsl_matrix_get(d1_mid,i,bp_step));
		gsl_vector_set(V2_m,i,gsl_matrix_get(d2_mid,i,bp_step));
		gsl_vector_set(V3_m,i,gsl_matrix_get(d3_mid,i,bp_step));
		gsl_vector_set(VR_m,i,gsl_matrix_get(R_mid,i,bp_step));
		
	}  
	
	double cos_ang, sq_Dist, E;
	int strand = (2*strand_3prime)-1;
	
	Make_SO3_Matrix (V01_m, V02_m, V03_m, Rot_1);
	Make_SO3_Matrix (V1_m, V2_m, V3_m, Rot_2);
	
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0, Rot_1, Rot_2, 0.0, Rot_rel);
	
	for (int i=0; i<3; i++)
	{
		gsl_vector_set(r_rel,i , gsl_vector_get(VR_m,i)-gsl_vector_get(V0R_m,i) + a*(gsl_vector_get(V1_m,i)-gsl_vector_get(V01_m,i)) + strand*b*(gsl_vector_get(V2_m,i)-gsl_vector_get(V02_m,i)) - strand*c*(gsl_vector_get(V3_m,i)-gsl_vector_get(V03_m,i)));
	} 
	
	cos_ang = 0.5*(gsl_matrix_get(Rot_rel,0,0)+gsl_matrix_get(Rot_rel,1,1)+gsl_matrix_get(Rot_rel,2,2)-1.0);
	
	sq_Dist = dot_prod(r_rel , r_rel);
	
	E = (0.5 * Ks_1 * sq_Dist) + (Ks_2 * (1-cos_ang));
	
	
	gsl_vector_free(V01_m);
	gsl_vector_free(V02_m);
	gsl_vector_free(V03_m);
	gsl_vector_free(V0R_m);
	
	gsl_vector_free(V1_m);
	gsl_vector_free(V2_m);
	gsl_vector_free(V3_m);
	gsl_vector_free(VR_m);
	
	gsl_matrix_free(Rot_1);
	gsl_matrix_free(Rot_2);
	gsl_matrix_free(Rot_rel);
	
	gsl_vector_free(r_rel);
}

// Function Total_Spring_Energy
double Total_Spring_Energy ()
{	
	int bps, strand_3prime;
	double Ks_1, Ks_2, E_s;
	
	E_s = 0.0;
	
	for (int n = 0; n < bond_steps_number; n++)
	{
		bps=gsl_vector_get(bond_steps,n);
		Ks_1 = gsl_vector_get(Pos_cons, n);
		Ks_2 = gsl_vector_get(Ang_cons, n);
		strand_3prime = gsl_vector_get(strand_3prime_array, n);
    
		E_s = E_s + Spring_Energy (dm1, dm2, dm3, Rm, bps, strand_3prime, Ks_1, Ks_2); 
	}
	
	return (E_s);
}

// Function Potential_Energy
double Potential_Energy(gsl_matrix* R_mat, gsl_vector* f_Vec)
{
	double E;
	
	gsl_vector* Vec = gsl_vector_alloc(3);
	
	for (int i=0;i<3;i++)
	{
		gsl_vector_set(Vec,i, (gsl_matrix_get(R_mat,i, (N_bp-1)) - gsl_matrix_get(R_mat,i, 0)) );
	}	

	E = - dot_prod (Vec , f_Vec);
	
	gsl_vector_free(Vec);
	
	return (E);
}

// Function Total_Energy
double Total_Energy(gsl_vector* f_Vec)
{
	double E_tot = 0.0;
	
	gsl_vector* theta = gsl_vector_alloc(3);
    gsl_vector* rho = gsl_vector_alloc(3);
    
    for (int bp_step = 0; bp_step < (N_bp-1); bp_step++)
    {
		E_tot = E_tot + Elastic_Energy (d1, d2, d3, R, theta, rho, bp_step, Seq_Code);
	}
	
	E_tot = E_tot + Potential_Energy(R, f_Vec);
	
	E_tot = E_tot + Total_Spring_Energy ();
	
	gsl_vector_free(theta);
	gsl_vector_free(rho);
	
    return (E_tot);
}

//Function Mid_frame
void Mid_frame (int bp_step)
{	
	gsl_vector* V1_1 = gsl_vector_alloc(3);
	gsl_vector* V2_1 = gsl_vector_alloc(3);
	gsl_vector* V3_1 = gsl_vector_alloc(3);
	gsl_vector* VR_1 = gsl_vector_alloc(3);
	
	gsl_vector* V1_2 = gsl_vector_alloc(3);
	gsl_vector* V2_2 = gsl_vector_alloc(3);
	gsl_vector* V3_2 = gsl_vector_alloc(3);
	gsl_vector* VR_2 = gsl_vector_alloc(3);
	
	gsl_vector* V1_m = gsl_vector_alloc(3);
	gsl_vector* V2_m = gsl_vector_alloc(3);
	gsl_vector* V3_m = gsl_vector_alloc(3);
	gsl_vector* VR_m = gsl_vector_alloc(3);
	
	gsl_vector* Rho_vec = gsl_vector_alloc(3);
	gsl_vector* Theta_vec = gsl_vector_alloc(3);
	
	for (int i = 0; i < 3; i++)
	{
		gsl_vector_set(V1_1,i,gsl_matrix_get(d1,i,bp_step));
		gsl_vector_set(V2_1,i,gsl_matrix_get(d2,i,bp_step));
		gsl_vector_set(V3_1,i,gsl_matrix_get(d3,i,bp_step));
		gsl_vector_set(VR_1,i,gsl_matrix_get(R,i,bp_step));
		
	    gsl_vector_set(V1_2,i,gsl_matrix_get(d1,i,bp_step+1));
		gsl_vector_set(V2_2,i,gsl_matrix_get(d2,i,bp_step+1));
		gsl_vector_set(V3_2,i,gsl_matrix_get(d3,i,bp_step+1));
		gsl_vector_set(VR_2,i,gsl_matrix_get(R,i,bp_step+1));
	}
	
	DOF (V1_1, V2_1, V3_1, VR_1, V1_2, V2_2, V3_2, VR_2, V1_m, V2_m, V3_m, VR_m, Theta_vec, Rho_vec);
	
	for (int i = 0; i < 3; i++)
	{
		gsl_matrix_set(dm1, i, bp_step, gsl_vector_get(V1_m,i));
		gsl_matrix_set(dm2, i, bp_step, gsl_vector_get(V2_m,i));
		gsl_matrix_set(dm3, i, bp_step, gsl_vector_get(V3_m,i));
		gsl_matrix_set(Rm, i, bp_step, gsl_vector_get(VR_m,i));
	}
	
	gsl_vector_free(V1_1);
	gsl_vector_free(V2_1);
	gsl_vector_free(V3_1);
	gsl_vector_free(VR_1);
	
	gsl_vector_free(V1_2);
	gsl_vector_free(V2_2);
	gsl_vector_free(V3_2);
	gsl_vector_free(VR_2);
	
	gsl_vector_free(V1_m);
	gsl_vector_free(V2_m);
	gsl_vector_free(V3_m);
	gsl_vector_free(VR_m);
	
	gsl_vector_free(Rho_vec);
	gsl_vector_free(Theta_vec);
}

// Function Load_init_config
void Load_init_config (int length_tail1, int length_tail2)
{
	gsl_vector* V1_1 = gsl_vector_alloc(3);
	gsl_vector* V2_1 = gsl_vector_alloc(3);
	gsl_vector* V3_1 = gsl_vector_alloc(3);
	gsl_vector* R_1 = gsl_vector_alloc(3);
	gsl_vector* Vm_3 = gsl_vector_alloc(3);
	gsl_matrix* Rot_rel = gsl_matrix_alloc(3,3);
	gsl_matrix* Rot_mid = gsl_matrix_alloc(3,3);
	
	double config_var, comp;
	
	ifstream config_file;
	
	config_file.open ("relx_d_1.txt");
	
		for (int k = 0; k < N_bp0; k++)
	{
		for (int i = 0; i < 3; i++)
		{
			config_file >> config_var;
			gsl_matrix_set(d1_0,i,k, config_var);
		}
	}
	
	config_file.close();
	
	config_file.open ("relx_d_2.txt");
	
		for (int k = 0; k < N_bp0; k++)
	{
		for (int i = 0; i < 3; i++)
		{
			config_file >> config_var;
			gsl_matrix_set(d2_0,i,k, config_var);
		}
	}
	
	config_file.close();
	
	config_file.open ("relx_d_3.txt");
	
		for (int k = 0; k < N_bp0; k++)
	{
		for (int i = 0; i < 3; i++)
		{
			config_file >> config_var;
			gsl_matrix_set(d3_0,i,k, config_var);
		}
	}
	
	config_file.close();
	
	config_file.open ("relx_R.txt");
	
		for (int k = 0; k < N_bp0; k++)
	{
		for (int i = 0; i < 3; i++)
		{
			config_file >> config_var;
			gsl_matrix_set(R_0,i,k, config_var);
		}
	}
	
	config_file.close();
	
	gsl_matrix_set_zero(d1);
	gsl_matrix_set_zero(d2);
	gsl_matrix_set_zero(d3);
	gsl_matrix_set_zero(R);
	
	for (int k=length_tail1; k<(length_tail1+N_bp0); k++)
	{
		for (int i=0; i<3; i++)
		{
			gsl_matrix_set(d1,i,k, gsl_matrix_get(d1_0,i,(k-length_tail1)));
			gsl_matrix_set(d2,i,k, gsl_matrix_get(d2_0,i,(k-length_tail1)));
			gsl_matrix_set(d3,i,k, gsl_matrix_get(d3_0,i,(k-length_tail1)));
			gsl_matrix_set(R,i,k, gsl_matrix_get(R_0,i,(k-length_tail1)));
		}
		
	}
	
	Rot_MAT(0.0, 0.0, Tw_0, Rot_rel);
	Rot_MAT(0.0, 0.0, 0.5*Tw_0, Rot_mid);
	
	for (int k=(length_tail1+N_bp0); k<(length_tail1+N_bp0+length_tail2); k++)
	{
		for (int i = 0; i < 3; i++)
		{
			gsl_vector_set(V1_1,i,gsl_matrix_get(d1,i,k-1));
			gsl_vector_set(V2_1,i,gsl_matrix_get(d2,i,k-1));
			gsl_vector_set(V3_1,i,gsl_matrix_get(d3,i,k-1));
			gsl_vector_set(R_1,i,gsl_matrix_get(R,i,k-1));
		}
		
		for (int j=0; j<3; j++)
		{
			comp = gsl_matrix_get(Rot_rel,0,0) * gsl_vector_get(V1_1,j)+gsl_matrix_get(Rot_rel,1,0) * gsl_vector_get(V2_1,j)+gsl_matrix_get(Rot_rel,2,0) * gsl_vector_get(V3_1,j);
			gsl_matrix_set(d1,j,k,comp);
			comp = gsl_matrix_get(Rot_rel,0,1) * gsl_vector_get(V1_1,j)+gsl_matrix_get(Rot_rel,1,1) * gsl_vector_get(V2_1,j)+gsl_matrix_get(Rot_rel,2,1) * gsl_vector_get(V3_1,j);
			gsl_matrix_set(d2,j,k,comp);
			comp = gsl_matrix_get(Rot_rel,0,2) * gsl_vector_get(V1_1,j)+gsl_matrix_get(Rot_rel,1,2) * gsl_vector_get(V2_1,j)+gsl_matrix_get(Rot_rel,2,2) * gsl_vector_get(V3_1,j);
			gsl_matrix_set(d3,j,k,comp);
			
			comp = gsl_matrix_get(Rot_mid,0,2) * gsl_vector_get(V1_1,j)+gsl_matrix_get(Rot_mid,1,2) * gsl_vector_get(V2_1,j)+gsl_matrix_get(Rot_mid,2,2) * gsl_vector_get(V3_1,j);
			gsl_vector_set(Vm_3,j,comp);
			
			comp=gsl_vector_get(R_1,j)+l_0*gsl_vector_get(Vm_3,j);
			gsl_matrix_set(R,j,k,comp);
		}
	}
	
	Rot_MAT(0.0, 0.0, -Tw_0, Rot_rel);
	Rot_MAT(0.0, 0.0, -0.5*Tw_0, Rot_mid);
	
	for (int k=(length_tail1-1); k>(-1); k--)
	{
		for (int i = 0; i < 3; i++)
		{
			gsl_vector_set(V1_1,i,gsl_matrix_get(d1,i,k+1));
			gsl_vector_set(V2_1,i,gsl_matrix_get(d2,i,k+1));
			gsl_vector_set(V3_1,i,gsl_matrix_get(d3,i,k+1));
			gsl_vector_set(R_1,i,gsl_matrix_get(R,i,k+1));
		}
		
		for (int j=0; j<3; j++)
		{
			comp = gsl_matrix_get(Rot_rel,0,0) * gsl_vector_get(V1_1,j)+gsl_matrix_get(Rot_rel,1,0) * gsl_vector_get(V2_1,j)+gsl_matrix_get(Rot_rel,2,0) * gsl_vector_get(V3_1,j);
			gsl_matrix_set(d1,j,k,comp);
			comp = gsl_matrix_get(Rot_rel,0,1) * gsl_vector_get(V1_1,j)+gsl_matrix_get(Rot_rel,1,1) * gsl_vector_get(V2_1,j)+gsl_matrix_get(Rot_rel,2,1) * gsl_vector_get(V3_1,j);
			gsl_matrix_set(d2,j,k,comp);
			comp = gsl_matrix_get(Rot_rel,0,2) * gsl_vector_get(V1_1,j)+gsl_matrix_get(Rot_rel,1,2) * gsl_vector_get(V2_1,j)+gsl_matrix_get(Rot_rel,2,2) * gsl_vector_get(V3_1,j);
			gsl_matrix_set(d3,j,k,comp);
			
			comp = gsl_matrix_get(Rot_mid,0,2) * gsl_vector_get(V1_1,j)+gsl_matrix_get(Rot_mid,1,2) * gsl_vector_get(V2_1,j)+gsl_matrix_get(Rot_mid,2,2) * gsl_vector_get(V3_1,j);
			gsl_vector_set(Vm_3,j,comp);
			
			comp=gsl_vector_get(R_1,j)-l_0*gsl_vector_get(Vm_3,j);
			gsl_matrix_set(R,j,k,comp);
		}
	}
	
	config_file.open ("init_Axes.txt");
	
	for (int i = 0; i < 3; i++)
	{
		config_file >> config_var;
		gsl_vector_set(X_Axis,i, config_var);
	}
	
	for (int i = 0; i < 3; i++)
	{
		config_file >> config_var;
		gsl_vector_set(Y_Axis,i, config_var);
	}
	
	for (int i = 0; i < 3; i++)
	{
		config_file >> config_var;
		gsl_vector_set(Z_Axis,i, config_var);
	}
	
	config_file.close();
	
	gsl_vector_set_basis(X0_Axis,0);
	gsl_vector_set_basis(Y0_Axis,1);
	gsl_vector_set_basis(Z0_Axis,2);
	
	for (int k = 0; k < (N_bp-1); k++)
	{
		Mid_frame (k);
	}
	
	gsl_matrix_memcpy(init_dm1, dm1);
	gsl_matrix_memcpy(init_dm2, dm2);
	gsl_matrix_memcpy(init_dm3, dm3);
	gsl_matrix_memcpy(init_Rm, Rm);
	
	gsl_vector_free(V1_1);
	gsl_vector_free(V2_1);
	gsl_vector_free(V3_1);
	gsl_vector_free(Vm_3);
	gsl_vector_free(R_1);
	gsl_matrix_free(Rot_rel);
	gsl_matrix_free(Rot_mid);
}

// Function Load_Binding_Sites
void Load_Binding_Sites (int bond_min, int bond_max, int length_tail1)
{
	int bond_bp_var, long_bond_bp_var1, long_bond_bp_var2, bond_var, bond_var1, bond_var2, number, number1, number2;
	double pos_var1, pos_var2, ang_var1, ang_var2;
	
	gsl_vector_set_zero(bond_array);
	
	ifstream bond_bp_file;
	ifstream Ang_cons_file;
	ifstream Pos_cons_file;
	
	bond_bp_file.open ("Bound_Phosphates_14.txt");
	Ang_cons_file.open("Ang_cons.txt");
	Pos_cons_file.open("Pos_cons.txt");
	
	bond_steps_number = 0.0;
	
	number=-1;
	
	for (int j=1; j<15; j++)
	{
		bond_bp_file >> bond_bp_var;
		long_bond_bp_var1=length_tail1+bond_bp_var-1;
		
		bond_bp_file >> bond_bp_var;
		long_bond_bp_var2=length_tail1+bond_bp_var-1;
		
		Ang_cons_file >> ang_var1;
		Ang_cons_file >> ang_var2;
		
		Pos_cons_file >> pos_var1;
		Pos_cons_file >> pos_var2;
		
		if ( (j > (bond_min-1)) && (j < (bond_max+1)) )
		{
			gsl_vector_set(bond_array,long_bond_bp_var1,j);
			gsl_vector_set(bond_array,long_bond_bp_var2,j);
			
			bond_steps_number++;
			number++;
			gsl_vector_set(bond_steps, number, long_bond_bp_var1);
			gsl_vector_set(strand_3prime_array, number, 0);
			gsl_vector_set(Pos_cons, number, pos_var1);
			gsl_vector_set(Ang_cons, number, ang_var1);
			
			bond_steps_number++;
			number++;
			gsl_vector_set(bond_steps, number, long_bond_bp_var2);
			gsl_vector_set(strand_3prime_array, number, 1);
			gsl_vector_set(Pos_cons, number, pos_var2);
			gsl_vector_set(Ang_cons, number, ang_var2);
		}
	}
	
	bond_bp_file.close();
	Ang_cons_file.close();
	Pos_cons_file.close();
	
	int bps_min=N_bp;
	int bps_max=-2;
	
	for (int i =0; i < (N_bp-1); i++)
	{
		bond_var = gsl_vector_get(bond_array,i);
		
		if (bond_var > 0)
		{
			if (i > bps_max)
			{
				bps_max = i;
			}
			
			if (i < bps_min)
			{
				bps_min = i;
			}
		}
	}
	
	free_wrapped_bp_number = 0.0; 
	free_tail1_bp_number = 0.0;
	free_tail2_bp_number = 0.0;
	
	number=-1;
	number1=-1;
	number2=-1;
	
	for (int i=0; i < N_bp; i++)
	{
		bond_var1=0;
		bond_var2=0;
		
		if (i< (N_bp-1))
		{
			bond_var1 = gsl_vector_get(bond_array,i);
		}
		
		if (i > 0)
		{
			bond_var2 = gsl_vector_get(bond_array,i-1);
		}
		
		if ((bond_var1==0) && (bond_var2==0))
		{
			if ((i > bps_min) && (i < bps_max+1))
			{
				free_wrapped_bp_number++;
				number++;
				gsl_vector_set(free_wrapped_bp,number, i);
			}
			
			if (i < bps_min)
			{
				free_tail1_bp_number++;
				number1++;
				gsl_vector_set(free_tail1_bp,number1, i);
			}
			
			if (i > bps_max+1)
			{
				free_tail2_bp_number++;
				number2++;
				gsl_vector_set(free_tail2_bp,number2, i);
			}	
		}
	}
	
	double N_tot, Sum_Glob_Mut_Move_Percents;
	
	N_tot = bond_steps_number + free_wrapped_bp_number + free_tail1_bp_number + free_tail2_bp_number;
	
	Pivot_vs_Loc_frac = 0.4;
	
	Bond_vs_DisPhos_frac = 1.0;
	
	Glob_Move_Percent = 0.0;
	
	Mut_Move_Percent = 0.3;
	
	Sum_Glob_Mut_Move_Percents = Mut_Move_Percent + Glob_Move_Percent;
	
	Bond_Move_Percent = (1-Sum_Glob_Mut_Move_Percents) * Bond_vs_DisPhos_frac * (bond_steps_number/N_tot);
	
	DisPhos_Move_Percent = (1-Sum_Glob_Mut_Move_Percents) * (1-Bond_vs_DisPhos_frac) * (bond_steps_number/N_tot);
	
	Loc_Move_wrapped_Percent = (1-Sum_Glob_Mut_Move_Percents) * (free_wrapped_bp_number/N_tot);
	
	Loc_Move_tail1_Percent = (1-Sum_Glob_Mut_Move_Percents) * (1-Pivot_vs_Loc_frac) * (free_tail1_bp_number/N_tot);
	
	Loc_Move_tail2_Percent = (1-Sum_Glob_Mut_Move_Percents) * (1-Pivot_vs_Loc_frac) * (free_tail2_bp_number/N_tot);
	
	Pivot_Move_tail1_Percent = (1-Sum_Glob_Mut_Move_Percents) * (Pivot_vs_Loc_frac) * (free_tail1_bp_number/N_tot);
	
	Pivot_Move_tail2_Percent = (1-Sum_Glob_Mut_Move_Percents) * (Pivot_vs_Loc_frac) * (free_tail2_bp_number/N_tot);
}

// Function Update_Force_vec
void Update_Force_vec (double f1, double f2, double f3, gsl_vector* Ax_1, gsl_vector* Ax_2,  gsl_vector* Ax_3,  gsl_vector* f)
{
	double comp;
	
	for (int i = 0; i < 3; i++)
	{
		comp = f1*gsl_vector_get(Ax_1,i) + f2*gsl_vector_get(Ax_2,i) +f3*gsl_vector_get(Ax_3,i);
		gsl_vector_set (f , i, comp);
	}
} 

// Function Free_MC_Move
double Free_MC_Move(double f1, double f2, double f3, int bp_i, int bp_f, gsl_vector* Rot_vec, gsl_vector* Trans_vec, gsl_matrix* Newd1, gsl_matrix* Newd2, gsl_matrix* Newd3, gsl_matrix* NewR)
{	
	gsl_matrix* Mat = gsl_matrix_alloc(3,3);
	
	gsl_vector* f_Vec = gsl_vector_alloc(3);
	    
    gsl_vector* V1 = gsl_vector_alloc(3);
    gsl_vector* V2 = gsl_vector_alloc(3);
    gsl_vector* V3 = gsl_vector_alloc(3);
    gsl_vector* VR = gsl_vector_alloc(3);
    
    gsl_vector* NewV1 = gsl_vector_alloc(3);
    gsl_vector* NewV2 = gsl_vector_alloc(3);
    gsl_vector* NewV3 = gsl_vector_alloc(3);
    gsl_vector* NewVR = gsl_vector_alloc(3);
    
    gsl_vector* Rot_Center = gsl_vector_alloc(3);
    
    gsl_vector* theta = gsl_vector_alloc(3);
    gsl_vector* rho = gsl_vector_alloc(3);
    
    Update_Force_vec (f1, f2, f3, X_Axis, Y_Axis,  Z_Axis,  f_Vec);
    
    double Th_1, Th_2, Th_3;
	
	Th_1 = gsl_vector_get(Rot_vec,0);
    Th_2 = gsl_vector_get(Rot_vec,1);
    Th_3 = gsl_vector_get(Rot_vec,2);
    
    Rot_MAT(Th_1, Th_2, Th_3, Mat);
    
    int delta_1, delta_2;
    double  fact_1, fact_2;
	
	delta_1 = (N_bp-1-bp_i)/(N_bp-1);
	delta_2 = bp_f/(N_bp-1);
    
    fact_1 = 0.5*(1-delta_1+delta_2);
    fact_2 = 0.5*(1-delta_2+delta_1);
    
	for (int i = 0; i < 3; i++)
	{
		gsl_vector_set(Rot_Center,i,(fact_2 * gsl_matrix_get(R,i,bp_f)+ fact_1 * gsl_matrix_get(R,i,bp_i)));
	}
	
	for(int i = 0; i< 3; i++)
	{
		gsl_matrix_set(NewR,i,0,gsl_matrix_get(R,i,0));
		gsl_matrix_set(NewR,i,(N_bp-1),gsl_matrix_get(R,i,(N_bp-1)));
	}
    
    for (int bp = bp_i; bp < (bp_f+1); bp++)
    {
		for (int i = 0; i < 3; i++)
		{
			gsl_vector_set(V1,i,gsl_matrix_get(d1,i,bp));
			gsl_vector_set(V2,i,gsl_matrix_get(d2,i,bp));
			gsl_vector_set(V3,i,gsl_matrix_get(d3,i,bp));
			gsl_vector_set(VR,i,(gsl_matrix_get(R,i,bp)-gsl_vector_get(Rot_Center,i)));
		}    
		
		gsl_blas_dgemv (CblasNoTrans, 1.0, Mat, V1, 0.0, NewV1);
		gsl_blas_dgemv (CblasNoTrans, 1.0, Mat, V2, 0.0, NewV2);
		gsl_blas_dgemv (CblasNoTrans, 1.0, Mat, V3, 0.0, NewV3);
		gsl_blas_dgemv (CblasNoTrans, 1.0, Mat, VR, 0.0, NewVR);
		 
		for (int i = 0; i < 3; i++)
		{
			gsl_matrix_set(Newd1,i,bp,gsl_vector_get(NewV1,i));
			gsl_matrix_set(Newd2,i,bp,gsl_vector_get(NewV2,i));
			gsl_matrix_set(Newd3,i,bp,gsl_vector_get(NewV3,i));
			gsl_matrix_set(NewR,i,bp,(gsl_vector_get(NewVR,i)+gsl_vector_get(Rot_Center,i)));
		
		}
	}
	
	for (int bp = bp_i; bp < (bp_f+1); bp++)
	{
		for (int i = 0; i < 3; i++)
		{
			gsl_vector_set(VR,i,gsl_matrix_get(NewR,i,bp));
		}
		
		gsl_vector_add(VR,Trans_vec); 
		
		for (int i = 0; i < 3; i++)
		{
			gsl_matrix_set(NewR,i,bp,gsl_vector_get(VR,i));
		}
	}

	double E_1, E_2, delta_E;
	
	E_1 = 0;
	E_2 = 0;
	
	if (bp_f < (N_bp-1))
	{
		for (int i = 0; i < 3; i++)
		{
			gsl_matrix_set(Newd1,i,bp_f+1,gsl_matrix_get(d1,i,bp_f+1));
			gsl_matrix_set(Newd2,i,bp_f+1,gsl_matrix_get(d2,i,bp_f+1));
			gsl_matrix_set(Newd3,i,bp_f+1,gsl_matrix_get(d3,i,bp_f+1));
			gsl_matrix_set(NewR,i,bp_f+1,gsl_matrix_get(R,i,bp_f+1));
		}

		E_1 = E_1 + Elastic_Energy(d1, d2, d3, R, theta, rho, bp_f, Seq_Code);
		
		E_2 = E_2 + Elastic_Energy(Newd1, Newd2, Newd3, NewR, theta, rho, bp_f, Seq_Code);
	}	
	
	if (bp_i > 0)
	{
		for (int i = 0; i < 3; i++)
		{
			gsl_matrix_set(Newd1,i,bp_i-1,gsl_matrix_get(d1,i,bp_i-1));
			gsl_matrix_set(Newd2,i,bp_i-1,gsl_matrix_get(d2,i,bp_i-1));
			gsl_matrix_set(Newd3,i,bp_i-1,gsl_matrix_get(d3,i,bp_i-1));
			gsl_matrix_set(NewR,i,bp_i-1,gsl_matrix_get(R,i,bp_i-1));
		}
		
		E_1 = E_1 + Elastic_Energy(d1, d2, d3, R, theta, rho, bp_i-1, Seq_Code);
		
		E_2 = E_2 + Elastic_Energy(Newd1, Newd2, Newd3, NewR, theta, rho, bp_i-1, Seq_Code);
	}
	
	E_1 = E_1 + Potential_Energy(R, f_Vec);
	
	E_2 = E_2 + Potential_Energy(NewR, f_Vec);
	
	delta_E = E_2 - E_1;
	
	gsl_matrix_free(Mat);
	
	gsl_vector_free(f_Vec);
	
	gsl_vector_free(V1);
	gsl_vector_free(V2);
	gsl_vector_free(V3);
	gsl_vector_free(VR);
	
	gsl_vector_free(NewV1);
	gsl_vector_free(NewV2);
	gsl_vector_free(NewV3);
	gsl_vector_free(NewVR);
	
	gsl_vector_free(Rot_Center);
	
	gsl_vector_free(rho);
	gsl_vector_free(theta);
	
	return (delta_E);
}

// Function Global_Rotation_MC_Move
double Global_Rotation_MC_Move (double f1, double f2, double f3, gsl_vector* Rot_vec, gsl_vector* NewX_Axis, gsl_vector* NewY_Axis, gsl_vector* NewZ_Axis)
{
	gsl_matrix* Mat = gsl_matrix_alloc(3,3);
	
	gsl_vector* f_Vec = gsl_vector_alloc(3);
	gsl_vector* Newf_Vec = gsl_vector_alloc(3);
	
	double Th_1, Th_2, Th_3;
	
	Th_1 = gsl_vector_get(Rot_vec,0);
    Th_2 = gsl_vector_get(Rot_vec,1);
    Th_3 = gsl_vector_get(Rot_vec,2);
	
	Rot_MAT(Th_1, Th_2, Th_3, Mat);
	
	gsl_blas_dgemv (CblasNoTrans, 1.0, Mat, X_Axis, 0.0, NewX_Axis);
	gsl_blas_dgemv (CblasNoTrans, 1.0, Mat, Y_Axis, 0.0, NewY_Axis);
	gsl_blas_dgemv (CblasNoTrans, 1.0, Mat, Z_Axis, 0.0, NewZ_Axis);
	
	Update_Force_vec (f1, f2, f3, X_Axis, Y_Axis,  Z_Axis,  f_Vec);
	
	Update_Force_vec (f1, f2, f3, NewX_Axis, NewY_Axis,  NewZ_Axis,  Newf_Vec);
	
	double E_1, E_2, delta_E;
	
	int bp_i, bp_f;
	
	E_1 = 0;
	E_2 = 0;
	
	bp_i = 0;
	bp_f = N_bp-1;
	
	E_1 = E_1 + Potential_Energy(R, f_Vec);
	
	E_2 = E_2 + Potential_Energy(R, Newf_Vec);
	
	delta_E = E_2 - E_1;
	
	gsl_matrix_free(Mat);
	
	gsl_vector_free(f_Vec);
	gsl_vector_free(Newf_Vec);
	
	return (delta_E);
}

// Function Bond_MC_Move
double Bond_MC_Move (double f1, double f2, double f3, int bps, gsl_vector* Rot_vec, gsl_vector* Trans_vec, gsl_matrix* Newd1, gsl_matrix* Newd2, gsl_matrix* Newd3, gsl_matrix* NewR)
{
	gsl_vector* f_Vec = gsl_vector_alloc(3);
	
    gsl_vector* V1_m = gsl_vector_alloc(3);
    gsl_vector* V2_m = gsl_vector_alloc(3);
    gsl_vector* V3_m = gsl_vector_alloc(3);
    gsl_vector* VR_m = gsl_vector_alloc(3);
	    
    gsl_vector* V1_1 = gsl_vector_alloc(3);
    gsl_vector* V2_1 = gsl_vector_alloc(3);
    gsl_vector* V3_1 = gsl_vector_alloc(3);
    gsl_vector* VR_1 = gsl_vector_alloc(3);
    
    gsl_vector* VR_2 = gsl_vector_alloc(3);
    
    gsl_vector* NewV1_1 = gsl_vector_alloc(3);
    gsl_vector* NewV2_1 = gsl_vector_alloc(3);
    gsl_vector* NewV3_1 = gsl_vector_alloc(3);
    
    gsl_vector* half_V1_m = gsl_vector_alloc(3);
    gsl_vector* half_V2_m = gsl_vector_alloc(3);
    gsl_vector* half_V3_m = gsl_vector_alloc(3);
    gsl_vector* half_VR_m = gsl_vector_alloc(3);
    
    gsl_vector* half_theta = gsl_vector_alloc(3);
    gsl_vector* half_rho = gsl_vector_alloc(3);
    
    gsl_vector* Rot_vec_2 = gsl_vector_alloc(3);
    
    gsl_vector* NewV1_2 = gsl_vector_alloc(3);
    gsl_vector* NewV2_2 = gsl_vector_alloc(3);
    gsl_vector* NewV3_2 = gsl_vector_alloc(3);
    
    gsl_matrix* Mat_1 = gsl_matrix_alloc(3,3);
	gsl_matrix* Mat_2 = gsl_matrix_alloc(3,3);
	gsl_matrix* Mat = gsl_matrix_alloc(3,3);
	
	gsl_vector* theta = gsl_vector_alloc(3);
    gsl_vector* rho = gsl_vector_alloc(3);
    
    Update_Force_vec (f1, f2, f3, X_Axis, Y_Axis,  Z_Axis, f_Vec);
    
    for (int i = 0; i < 3; i++)
	{
		gsl_vector_set(V1_1,i,gsl_matrix_get(d1,i,bps));
		gsl_vector_set(V2_1,i,gsl_matrix_get(d2,i,bps));
		gsl_vector_set(V3_1,i,gsl_matrix_get(d3,i,bps));
		gsl_vector_set(VR_1,i,gsl_matrix_get(R,i,bps));
		
		gsl_vector_set(VR_2,i,gsl_matrix_get(R,i,bps+1));
		
		gsl_vector_set(V1_m,i,gsl_matrix_get(dm1,i,bps));
		gsl_vector_set(V2_m,i,gsl_matrix_get(dm2,i,bps));
		gsl_vector_set(V3_m,i,gsl_matrix_get(dm3,i,bps));
		gsl_vector_set(VR_m,i,gsl_matrix_get(Rm,i,bps));
	}  
	 
	double Th_1, Th_2, Th_3;
	
	Th_1 = gsl_vector_get(Rot_vec,0);
    Th_2 = gsl_vector_get(Rot_vec,1);
    Th_3 = gsl_vector_get(Rot_vec,2);
    
    Rot_MAT(Th_1, Th_2, Th_3, Mat_1);
    
    DOF (V1_1, V2_1, V3_1, VR_1, V1_m, V2_m, V3_m, VR_m, half_V1_m, half_V2_m, half_V3_m, half_VR_m, half_theta, half_rho);
    
    double half_Th_1, half_Th_2, half_Th_3;
	
	half_Th_1 = gsl_vector_get(half_theta,0);
    half_Th_2 = gsl_vector_get(half_theta,1);
    half_Th_3 = gsl_vector_get(half_theta,2);
    
    double comp;
    
    for (int i = 0; i < 3; i++)
    {
		comp = half_Th_1*gsl_vector_get(V1_m,i) + half_Th_2*gsl_vector_get(V2_m,i) + half_Th_3*gsl_vector_get(V3_m,i);
		gsl_vector_set (Rot_vec_2,i, comp);
	}
	
	Th_1 = gsl_vector_get(Rot_vec_2,0);
    Th_2 = gsl_vector_get(Rot_vec_2,1);
    Th_3 = gsl_vector_get(Rot_vec_2,2);
    
    Rot_MAT(Th_1, Th_2, Th_3, Mat_2);
    
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0, Mat_1, Mat_2, 0.0, Mat);
    
	gsl_blas_dgemv (CblasNoTrans, 1.0, Mat, V1_m, 0.0, NewV1_1);
	gsl_blas_dgemv (CblasNoTrans, 1.0, Mat, V2_m, 0.0, NewV2_1);
	gsl_blas_dgemv (CblasNoTrans, 1.0, Mat, V3_m, 0.0, NewV3_1);
	gsl_vector_add(VR_1,Trans_vec); 
	
    gsl_blas_dgemv (CblasTrans, 1.0, Mat, V1_m, 0.0, NewV1_2);
	gsl_blas_dgemv (CblasTrans, 1.0, Mat, V2_m, 0.0, NewV2_2);
	gsl_blas_dgemv (CblasTrans, 1.0, Mat, V3_m, 0.0, NewV3_2);
	gsl_vector_sub(VR_2,Trans_vec); 
	
	for(int i = 0; i< 3; i++)
	{
		gsl_matrix_set(NewR,i,0,gsl_matrix_get(R,i,0));
		gsl_matrix_set(NewR,i,(N_bp-1),gsl_matrix_get(R,i,(N_bp-1)));
	}
	
	for (int i = 0; i < 3; i++)
	{
		gsl_matrix_set(Newd1,i,bps,gsl_vector_get(NewV1_1,i));
		gsl_matrix_set(Newd2,i,bps,gsl_vector_get(NewV2_1,i));
		gsl_matrix_set(Newd3,i,bps,gsl_vector_get(NewV3_1,i));
		gsl_matrix_set(NewR,i,bps,gsl_vector_get(VR_1,i));
		
		gsl_matrix_set(Newd1,i,bps+1,gsl_vector_get(NewV1_2,i));
		gsl_matrix_set(Newd2,i,bps+1,gsl_vector_get(NewV2_2,i));
		gsl_matrix_set(Newd3,i,bps+1,gsl_vector_get(NewV3_2,i));
		gsl_matrix_set(NewR,i,bps+1,gsl_vector_get(VR_2,i));
		
	}
	
	int bp_i, bp_f, bp_strt, bp_fnsh;
	
	bp_i=bps;
	bp_f=bp_i+1;
	
	bp_strt=bp_i;
	bp_fnsh=bp_f;
	
	if (bp_i > 0)
	{
		bp_strt=bp_i-1;
		
		for (int i = 0; i < 3; i++)
		{
			gsl_matrix_set(Newd1,i,bp_i-1,gsl_matrix_get(d1,i,bp_i-1));
			gsl_matrix_set(Newd2,i,bp_i-1,gsl_matrix_get(d2,i,bp_i-1));
			gsl_matrix_set(Newd3,i,bp_i-1,gsl_matrix_get(d3,i,bp_i-1));
			gsl_matrix_set(NewR,i,bp_i-1,gsl_matrix_get(R,i,bp_i-1));
		}
	}
	
	if (bp_f < (N_bp-1))
	{
		bp_fnsh=bp_f+1;
		
		for (int i = 0; i < 3; i++)
		{
			gsl_matrix_set(Newd1,i,bp_f+1,gsl_matrix_get(d1,i,bp_f+1));
			gsl_matrix_set(Newd2,i,bp_f+1,gsl_matrix_get(d2,i,bp_f+1));
			gsl_matrix_set(Newd3,i,bp_f+1,gsl_matrix_get(d3,i,bp_f+1));
			gsl_matrix_set(NewR,i,bp_f+1,gsl_matrix_get(R,i,bp_f+1));
		}
	}
	
	double E_1, E_2, delta_E;
	
	E_1 = 0;
	E_2 = 0;
	
	for (int bp=bp_strt; bp<bp_fnsh; bp++)
	{
		E_1 = E_1 + Elastic_Energy(d1, d2, d3, R, theta, rho, bp, Seq_Code);
						
		E_2 = E_2 + Elastic_Energy(Newd1, Newd2, Newd3, NewR, theta, rho, bp, Seq_Code);
	}
	
	E_1 = E_1 + Potential_Energy(R, f_Vec);
	
	E_2 = E_2 + Potential_Energy(NewR, f_Vec);	 
	
	delta_E = E_2 - E_1;	
	
	gsl_vector_free(f_Vec);
		
	gsl_vector_free(V1_m);
	gsl_vector_free(V2_m);
	gsl_vector_free(V3_m);
	gsl_vector_free(VR_m);
	
	gsl_vector_free(V1_1);
	gsl_vector_free(V2_1);
	gsl_vector_free(V3_1);
	gsl_vector_free(VR_1);
	
	gsl_vector_free(VR_2);

	gsl_vector_free(NewV1_1);
	gsl_vector_free(NewV2_1);
	gsl_vector_free(NewV3_1);
	
	gsl_vector_free(half_V1_m);
	gsl_vector_free(half_V2_m);
	gsl_vector_free(half_V3_m);
	gsl_vector_free(half_VR_m);
	
	gsl_vector_free(half_theta);
	gsl_vector_free(half_rho);
	
	gsl_vector_free(Rot_vec_2);
	
	gsl_vector_free(NewV1_2);
	gsl_vector_free(NewV2_2);
	gsl_vector_free(NewV3_2);
	
	gsl_vector_free(theta);
	gsl_vector_free(rho);
	
	gsl_matrix_free(Mat_1);
	gsl_matrix_free(Mat_2);
	gsl_matrix_free(Mat);
		
	return (delta_E);
}

// Function Displace_Phosphate_MC_Move
double Displace_Phosphate_MC_Move (double f1, double f2, double f3, int bps, int strand_3prime, double Ks_1, double Ks_2, gsl_vector* Rot_vec, gsl_vector* Trans_Phos_vec, gsl_matrix* Newd1, gsl_matrix* Newd2, gsl_matrix* Newd3, gsl_matrix* NewR, gsl_matrix* Newdm1, gsl_matrix* Newdm2, gsl_matrix* Newdm3, gsl_matrix* NewRm)
{
	gsl_vector* f_Vec = gsl_vector_alloc(3);
	
    gsl_vector* V1 = gsl_vector_alloc(3);
    gsl_vector* V2 = gsl_vector_alloc(3);
    gsl_vector* V3 = gsl_vector_alloc(3);
    gsl_vector* VR = gsl_vector_alloc(3);
    
    gsl_vector* V1_m = gsl_vector_alloc(3);
    gsl_vector* V2_m = gsl_vector_alloc(3);
    gsl_vector* V3_m = gsl_vector_alloc(3);
    gsl_vector* VR_m = gsl_vector_alloc(3);
    
    gsl_vector* NewV1 = gsl_vector_alloc(3);
    gsl_vector* NewV2 = gsl_vector_alloc(3);
    gsl_vector* NewV3 = gsl_vector_alloc(3);
    
    gsl_vector* NewV1_m = gsl_vector_alloc(3);
    gsl_vector* NewV2_m = gsl_vector_alloc(3);
    gsl_vector* NewV3_m = gsl_vector_alloc(3);
    
    gsl_vector* theta = gsl_vector_alloc(3);
    gsl_vector* rho = gsl_vector_alloc(3);
    
    gsl_vector* cent_to_phos_1 = gsl_vector_alloc(3);
    gsl_vector* cent_to_phos_2 = gsl_vector_alloc(3);
    gsl_vector* Trans_vec = gsl_vector_alloc(3);
    
    gsl_matrix* Mat = gsl_matrix_alloc(3,3);
    
    Update_Force_vec (f1, f2, f3, X_Axis, Y_Axis,  Z_Axis,  f_Vec);
    
    double Th_1, Th_2, Th_3;
	
	Th_1 = gsl_vector_get(Rot_vec,0);
    Th_2 = gsl_vector_get(Rot_vec,1);
    Th_3 = gsl_vector_get(Rot_vec,2);
    
    Rot_MAT(Th_1, Th_2, Th_3, Mat);
    
    for(int i = 0; i< 3; i++)
	{
		gsl_matrix_set(NewR,i,0,gsl_matrix_get(R,i,0));
		gsl_matrix_set(NewR,i,(N_bp-1),gsl_matrix_get(R,i,(N_bp-1)));
	}
    
    for (int bp = bps; bp < (bps+2); bp++)
    {
		for (int i = 0; i < 3; i++)
		{
			gsl_vector_set(V1,i,gsl_matrix_get(d1,i,bp));
			gsl_vector_set(V2,i,gsl_matrix_get(d2,i,bp));
			gsl_vector_set(V3,i,gsl_matrix_get(d3,i,bp));
		}    
		
		gsl_blas_dgemv (CblasNoTrans, 1.0, Mat, V1, 0.0, NewV1);
		gsl_blas_dgemv (CblasNoTrans, 1.0, Mat, V2, 0.0, NewV2);
		gsl_blas_dgemv (CblasNoTrans, 1.0, Mat, V3, 0.0, NewV3);
		 
		for (int i = 0; i < 3; i++)
		{
			gsl_matrix_set(Newd1,i,bp,gsl_vector_get(NewV1,i));
			gsl_matrix_set(Newd2,i,bp,gsl_vector_get(NewV2,i));
			gsl_matrix_set(Newd3,i,bp,gsl_vector_get(NewV3,i));	
		}
	}
	
	for (int i = 0; i < 3; i++)
	{
		gsl_vector_set(V1_m,i,gsl_matrix_get(dm1,i,bps));
		gsl_vector_set(V2_m,i,gsl_matrix_get(dm2,i,bps));
		gsl_vector_set(V3_m,i,gsl_matrix_get(dm3,i,bps));
	}    
		
	gsl_blas_dgemv (CblasNoTrans, 1.0, Mat, V1_m, 0.0, NewV1_m);
	gsl_blas_dgemv (CblasNoTrans, 1.0, Mat, V2_m, 0.0, NewV2_m);
	gsl_blas_dgemv (CblasNoTrans, 1.0, Mat, V3_m, 0.0, NewV3_m);
		 
	for (int i = 0; i < 3; i++)
	{
		gsl_matrix_set(Newdm1,i,bps,gsl_vector_get(NewV1_m,i));
		gsl_matrix_set(Newdm2,i,bps,gsl_vector_get(NewV2_m,i));
		gsl_matrix_set(Newdm3,i,bps,gsl_vector_get(NewV3_m,i));	
	}
	
	int strand = (2*strand_3prime)-1;
	
	for (int i=0; i<3; i++)
	{
		gsl_vector_set(cent_to_phos_1, i, a*gsl_vector_get(V1_m,i) + strand*b*gsl_vector_get(V2_m,i) - strand*c*gsl_vector_get(V3_m,i));
		
		gsl_vector_set(cent_to_phos_2, i, a*gsl_vector_get(NewV1_m,i) + strand*b*gsl_vector_get(NewV2_m,i) - strand*c*gsl_vector_get(NewV3_m,i));
	} 
	
	gsl_vector_memcpy(Trans_vec, Trans_Phos_vec);
	gsl_vector_add(Trans_vec, cent_to_phos_1);
	gsl_vector_sub(Trans_vec, cent_to_phos_2);
	
	for (int bp = bps; bp < (bps+2); bp++)
	{
		for (int i = 0; i < 3; i++)
		{
			gsl_vector_set(VR,i,gsl_matrix_get(R,i,bp));
		}
		
		gsl_vector_add(VR,Trans_vec); 
		
		for (int i = 0; i < 3; i++)
		{
			gsl_matrix_set(NewR,i,bp,gsl_vector_get(VR,i));
		}
	}
	
	for (int i = 0; i < 3; i++)
	{
		gsl_vector_set(VR_m,i,gsl_matrix_get(Rm,i,bps));
	}
		
	gsl_vector_add(VR_m,Trans_vec); 
		
	for (int i = 0; i < 3; i++)
	{
		gsl_matrix_set(NewRm,i,bps,gsl_vector_get(VR_m,i));
	}
	
	int bp_i, bp_f, bp_strt, bp_fnsh;
	
	bp_i=bps;
	bp_f=bp_i+1;
	
	bp_strt=bp_i;
	bp_fnsh=bp_f;
	
	if (bp_i > 0)
	{
		bp_strt=bp_i-1;
		
		for (int i = 0; i < 3; i++)
		{
			gsl_matrix_set(Newd1,i,bp_i-1,gsl_matrix_get(d1,i,bp_i-1));
			gsl_matrix_set(Newd2,i,bp_i-1,gsl_matrix_get(d2,i,bp_i-1));
			gsl_matrix_set(Newd3,i,bp_i-1,gsl_matrix_get(d3,i,bp_i-1));
			gsl_matrix_set(NewR,i,bp_i-1,gsl_matrix_get(R,i,bp_i-1));
		}
	}
	
	if (bp_f < (N_bp-1))
	{
		bp_fnsh=bp_f+1;
		
		for (int i = 0; i < 3; i++)
		{
			gsl_matrix_set(Newd1,i,bp_f+1,gsl_matrix_get(d1,i,bp_f+1));
			gsl_matrix_set(Newd2,i,bp_f+1,gsl_matrix_get(d2,i,bp_f+1));
			gsl_matrix_set(Newd3,i,bp_f+1,gsl_matrix_get(d3,i,bp_f+1));
			gsl_matrix_set(NewR,i,bp_f+1,gsl_matrix_get(R,i,bp_f+1));
		}
	}
	
	double E_1, E_2, delta_E;
	
	E_1 = 0;
	E_2 = 0;
	
	for (int bp=bp_strt; bp<bp_fnsh; bp++)
	{
		E_1 = E_1 + Elastic_Energy(d1, d2, d3, R, theta, rho, bp, Seq_Code);
						
		E_2 = E_2 + Elastic_Energy(Newd1, Newd2, Newd3, NewR, theta, rho, bp, Seq_Code);
	}
	
	E_1 = E_1 + Potential_Energy(R, f_Vec);
	
	E_2 = E_2 + Potential_Energy(NewR, f_Vec);	
	
	E_1 = E_1 + Spring_Energy (dm1, dm2, dm3, Rm, bps, strand_3prime, Ks_1, Ks_2); 
	
	E_2 = E_2 + Spring_Energy (Newdm1, Newdm2, Newdm3, NewRm, bps, strand_3prime, Ks_1, Ks_2); 
	
	delta_E = E_2 - E_1;
	
	gsl_vector_free(f_Vec);
	
	gsl_vector_free(V1);
	gsl_vector_free(V2);
	gsl_vector_free(V3);
	gsl_vector_free(VR);
	
	gsl_vector_free(V1_m);
	gsl_vector_free(V2_m);
	gsl_vector_free(V3_m);
	gsl_vector_free(VR_m);
	
	gsl_vector_free(NewV1);
	gsl_vector_free(NewV2);
	gsl_vector_free(NewV3);
	
	gsl_vector_free(NewV1_m);
	gsl_vector_free(NewV2_m);
	gsl_vector_free(NewV3_m);
	
	gsl_vector_free(theta);
	gsl_vector_free(rho);
	
	gsl_vector_free(cent_to_phos_1);
	gsl_vector_free(cent_to_phos_2);
	gsl_vector_free(Trans_vec);
	
	gsl_matrix_free(Mat);
	
	return (delta_E);
}

// Function Mutation_MC_Move
double Mutation_MC_Move (int bp_i, int bp_f, gsl_vector* New_bp_Seq_Code, gsl_vector* New_Seq_Code)
{
	gsl_vector* theta = gsl_vector_alloc(3);
    gsl_vector* rho = gsl_vector_alloc(3);
	
	int bps_1, bps_2;
	
	bps_1 = bp_i;
	bps_2 = bp_f-1;
	
	if (bp_i > 0)
	{
		bps_1 = bp_i-1;
		gsl_vector_set (New_bp_Seq_Code, (bp_i-1), gsl_vector_get (bp_Seq_Code, (bp_i-1)));
	}
	
	if (bp_f < (N_bp-1))
	{
		bps_2 = bp_f;
		gsl_vector_set (New_bp_Seq_Code, (bp_f+1), gsl_vector_get (bp_Seq_Code, (bp_f+1)));
	}
	
	bp_Seq_Code_to_Seq_Code (bps_1, bps_2, New_bp_Seq_Code, New_Seq_Code);
	
	double E_1, E_2, delta_E;
	
	E_1 = 0;
	E_2 = 0;
	
	for (int bps = bps_1; bps < (bps_2+1); bps++)
	{
		E_1 = E_1 + Elastic_Energy(d1, d2, d3, R, theta, rho, bps, Seq_Code);	
		
		E_2 = E_2 + Elastic_Energy(d1, d2, d3, R, theta, rho, bps, New_Seq_Code);
	}
	
	delta_E = E_2 - E_1;
	
	gsl_vector_free(theta);
	gsl_vector_free(rho);
	
	return (delta_E);
}

// Function Random_basepair
int Random_basepair (int bp_1, int bp_2)
{
	double r;
	int bp;
	r = bp_1 + (bp_2-bp_1+1)*RanGen.Random();
	bp = int(r);
	return (bp);
}

// Function Random_Rot_vector
void Random_Rot_vector (gsl_vector* Rand_vec)
{
	double a, r, b, rot_ang;
	
	a = 2*PI*RanGen.Random();
	r = RanGen.Random();
	b = acos(1-2*r);
	rot_ang = max_rot*RanGen.Random();
	
	gsl_vector_set(Rand_vec,0, rot_ang*sin(b)*cos(a));
	gsl_vector_set(Rand_vec,1, rot_ang*sin(b)*sin(a));
	gsl_vector_set(Rand_vec,2, rot_ang*cos(b));
}

// Function Random_Trans_vector
void Random_Trans_vector (gsl_vector* Rand_vec)
{	
	gsl_vector_set(Rand_vec,0, -max_trans + 2.0*max_trans*RanGen.Random());
	gsl_vector_set(Rand_vec,1, -max_trans + 2.0*max_trans*RanGen.Random());
	gsl_vector_set(Rand_vec,2, -max_trans + 2.0*max_trans*RanGen.Random());
}

// Function Accept_Loc_Move
void Accept_Loc_Move (int bp, gsl_matrix* Newd1, gsl_matrix* Newd2, gsl_matrix* Newd3, gsl_matrix* NewR, double delta_E)
{
	for (int i = 0; i < 3; i++)
	{
		gsl_matrix_set(d1,i,bp,gsl_matrix_get(Newd1,i,bp));
		gsl_matrix_set(d2,i,bp,gsl_matrix_get(Newd2,i,bp));
		gsl_matrix_set(d3,i,bp,gsl_matrix_get(Newd3,i,bp));
		gsl_matrix_set(R,i,bp,gsl_matrix_get(NewR,i,bp));
	}
	
	Accept_number_loc = Accept_number_loc + 1.0;
	
	Accept_number = Accept_number + 1.0;
		
	Tot_Energy = Tot_Energy + delta_E;
}

// Function Accept_Pivot_Move
void Accept_Pivot_Move (int bp_i, int bp_f, gsl_matrix* Newd1, gsl_matrix* Newd2, gsl_matrix* Newd3, gsl_matrix* NewR, double delta_E)
{
	for (int bp = bp_i; bp < (bp_f+1); bp++)
	{
		for (int i = 0; i < 3; i++)
		{
			gsl_matrix_set(d1,i,bp,gsl_matrix_get(Newd1,i,bp));
			gsl_matrix_set(d2,i,bp,gsl_matrix_get(Newd2,i,bp));
			gsl_matrix_set(d3,i,bp,gsl_matrix_get(Newd3,i,bp));
			gsl_matrix_set(R,i,bp,gsl_matrix_get(NewR,i,bp));
		}
		
	}
	
	Accept_number_pivot = Accept_number_pivot + 1.0;
	
	Accept_number = Accept_number + 1.0;
		
	Tot_Energy = Tot_Energy + delta_E;
}

// Function Accept_Glob_Move
void Accept_Glob_Move (gsl_vector* NewX_Axis, gsl_vector* NewY_Axis, gsl_vector* NewZ_Axis, double delta_E)
{	
	gsl_vector_memcpy(X_Axis,NewX_Axis);
	gsl_vector_memcpy(Y_Axis,NewY_Axis);
	gsl_vector_memcpy(Z_Axis,NewZ_Axis);
	
	Accept_number_glob = Accept_number_glob + 1.0;
	
	Accept_number = Accept_number + 1.0;
		
	Tot_Energy = Tot_Energy + delta_E;
}

// Function Accept_Bond_Move
void Accept_Bond_Move (int bps, gsl_matrix* Newd1, gsl_matrix* Newd2, gsl_matrix* Newd3, gsl_matrix* NewR, double delta_E)
{
	int bp_i, bp_f;
	
	bp_i = bps;
	bp_f = bps+1;
	
	for (int bp = bp_i; bp < (bp_f+1); bp++)
	{
		for (int i = 0; i < 3; i++)
		{
			gsl_matrix_set(d1,i,bp,gsl_matrix_get(Newd1,i,bp));
			gsl_matrix_set(d2,i,bp,gsl_matrix_get(Newd2,i,bp));
			gsl_matrix_set(d3,i,bp,gsl_matrix_get(Newd3,i,bp));
			gsl_matrix_set(R,i,bp,gsl_matrix_get(NewR,i,bp));
		}
		
	}

	Accept_number_bond = Accept_number_bond + 1.0;
	
	Accept_number = Accept_number + 1.0;
		
	Tot_Energy = Tot_Energy + delta_E;

}

//Function Accept_Displace_Phosphate_Move 
void Accept_Displace_Phosphate_Move (int bps, gsl_matrix* Newd1, gsl_matrix* Newd2, gsl_matrix* Newd3, gsl_matrix* NewR, gsl_matrix* Newdm1, gsl_matrix* Newdm2, gsl_matrix* Newdm3, gsl_matrix* NewRm, double delta_E)
{
	int bp_i, bp_f;
	
	bp_i = bps;
	bp_f = bps+1;
	
	for (int bp = bp_i; bp < (bp_f+1); bp++)
	{
		for (int i = 0; i < 3; i++)
		{
			gsl_matrix_set(d1,i,bp,gsl_matrix_get(Newd1,i,bp));
			gsl_matrix_set(d2,i,bp,gsl_matrix_get(Newd2,i,bp));
			gsl_matrix_set(d3,i,bp,gsl_matrix_get(Newd3,i,bp));
			gsl_matrix_set(R,i,bp,gsl_matrix_get(NewR,i,bp));
		}
		
	}
	
	for (int i = 0; i < 3; i++)
	{
		gsl_matrix_set(dm1,i,bps,gsl_matrix_get(Newdm1,i,bps));
		gsl_matrix_set(dm2,i,bps,gsl_matrix_get(Newdm2,i,bps));
		gsl_matrix_set(dm3,i,bps,gsl_matrix_get(Newdm3,i,bps));
		gsl_matrix_set(Rm,i,bps,gsl_matrix_get(NewRm,i,bps));
	}

	Accept_number_DisPhos = Accept_number_DisPhos + 1.0;
	
	Accept_number = Accept_number + 1.0;
		
	Tot_Energy = Tot_Energy + delta_E;
}

// Function Accept_Mutation_Move
void Accept_Mutation_Move (int bp_i, int bp_f, gsl_vector* New_bp_Seq_Code, gsl_vector* New_Seq_Code, double delta_E)
{
	int bps_1, bps_2;
	
	bps_1 = bp_i;
	bps_2 = bp_f-1;
	
	if (bp_i > 0)
	{
		bps_1 = bp_i-1;
	}
	
	if (bp_f < (N_bp-1))
	{
		bps_2 = bp_f;
	}
	
	for (int bp = bp_i; bp < (bp_f+1); bp++)
	{
		gsl_vector_set (bp_Seq_Code, bp, gsl_vector_get (New_bp_Seq_Code, bp));
	}
	
	for (int bps = bps_1; bps < (bps_2+1); bps++)
	{
		gsl_vector_set (Seq_Code, bps, gsl_vector_get (New_Seq_Code, bps));
	}
	
	Accept_number_Mutation = Accept_number_Mutation + 1.0;
	
	Accept_number = Accept_number + 1.0;
	
	Tot_Energy = Tot_Energy + delta_E;
}

// Function Save_Err_data
void Save_Err_data ()
{
	double config_var;
	
	ofstream config_file;
	
	config_file.open ("Err_d_1.txt");
	
	for (int k = 0; k < N_bp; k++)
	{
		for (int i = 0; i < 3; i++)
		{
			config_var = gsl_matrix_get(d1,i,k);
			config_file << config_var << "\t"; 
		}
		config_file << "\n";
	}
	
    config_file.close();

	config_file.open ("Err_d_2.txt");
	
	for (int k = 0; k < N_bp; k++)
	{
		for (int i = 0; i < 3; i++)
		{
			config_var = gsl_matrix_get(d2,i,k);
			config_file << config_var << "\t"; 
		}
		config_file << "\n";
	}
	
    config_file.close();
    
    config_file.open ("Err_d_3.txt");
	
	for (int k = 0; k < N_bp; k++)
	{
		for (int i = 0; i < 3; i++)
		{
			config_var = gsl_matrix_get(d3,i,k);
			config_file << config_var << "\t"; 
		}
		config_file << "\n";
	}
	
    config_file.close();
    
    config_file.open ("Err_R.txt");
	
	for (int k = 0; k < N_bp; k++)
	{
		for (int i = 0; i < 3; i++)
		{
			config_var = gsl_matrix_get(R,i,k);
			config_file << config_var << "\t"; 
		}
		config_file << "\n";
	}
	
    config_file.close();
    
    config_file.open ("Err_Axes.txt");
    
    for (int i = 0; i < 3; i++)
	{
		config_var = gsl_vector_get(X_Axis,i);
		config_file << config_var << "\t"; 
	}
	
	config_file << "\n";
	
	for (int i = 0; i < 3; i++)
	{
		config_var = gsl_vector_get(Y_Axis,i);
		config_file << config_var << "\t"; 
	}
	
	config_file << "\n";
	
	for (int i = 0; i < 3; i++)
	{
		config_var = gsl_vector_get(Z_Axis,i);
		config_file << config_var << "\t"; 
	}
	
	config_file.close();
	
	config_file.open ("Err_Energy.txt", ios::app);
	
	config_file << Tot_Energy << "\n"; 
	
	config_file.close();
	
	config_file.open ("Err_Accept_ratios.txt");
	
	config_file << Accept_ratio_loc << "\n"; 
	config_file << Accept_ratio_pivot << "\n"; 
	config_file << Accept_ratio_bond << "\n"; 
	config_file << Accept_ratio_glob << "\n"; 
	config_file << Accept_ratio << "\n"; 
	
	config_file.close();
}

//Function Local_Metropolis
void Local_Metropolis (double f1, double f2, double f3, gsl_matrix* Newd1, gsl_matrix* Newd2, gsl_matrix* Newd3, gsl_matrix* NewR, bool in_tail1, bool in_wrapped, bool in_tail2)
{
	gsl_vector* Rot_vec = gsl_vector_alloc(3);
	gsl_vector* Trans_vec = gsl_vector_alloc(3);
	
	int n1, n2, n3, N1, N2, N3, bp1, bp2, bp3, bp_i, bp_f, int_in_tail1, int_in_wrapped, int_in_tail2;
	double delta_E, rand_MC;
	bool cond;
	
	Try_number = Try_number + 1.0;
	Try_number_loc = Try_number_loc + 1.0;
				
	N1=int(free_tail1_bp_number)-1;
	N2=int(free_wrapped_bp_number)-1;
	N3=int(free_tail2_bp_number)-1;
	
	int_in_tail1 = int(in_tail1);
	int_in_wrapped = int(in_wrapped);
	int_in_tail2 = int(in_tail2);
	
	n1 = Random_basepair (0, N1);
	n2 = Random_basepair (0, N2);
	n3 = Random_basepair (0, N3);
	
	bp1=gsl_vector_get(free_tail1_bp, n1);
	bp2=gsl_vector_get(free_wrapped_bp, n2);
	bp3=gsl_vector_get(free_tail2_bp, n3);
	
	bp_i = int_in_tail1 * bp1 + int_in_wrapped * bp2 + int_in_tail2 * bp3;
	bp_f = bp_i;
	
	//max_rot = PI/45.0;
	//max_trans = 1/45.0;
	
	max_rot = PI/90.0;
	max_trans = 1/90.0;
	
	Random_Rot_vector (Rot_vec);
	Random_Trans_vector (Trans_vec);
	
	delta_E = Free_MC_Move(f1, f2, f3, bp_i, bp_f, Rot_vec, Trans_vec, Newd1, Newd2, Newd3, NewR);
	
	rand_MC = RanGen.Random();
	
	cond = false;
	if (delta_E <= 0)
	{cond = true;} 
	else if (rand_MC < exp(-beta*delta_E)) 
	{cond = true;}
			
	if (cond)
	{
		Accept_Loc_Move (bp_i, Newd1, Newd2, Newd3, NewR, delta_E);
	}
	
	if (Error_Code==1)
	{
		cout << "An Error Occured! MC Loop Terminated! \n";
		Save_Err_data ();
		exit(1);
	}
	
	gsl_vector_free(Rot_vec);
	gsl_vector_free(Trans_vec);
}

//Function Pivot_Metropolis
void Pivot_Metropolis (double f1, double f2, double f3, gsl_matrix* Newd1, gsl_matrix* Newd2, gsl_matrix* Newd3, gsl_matrix* NewR, bool in_tail1, bool in_tail2)
{
	gsl_vector* Rot_vec = gsl_vector_alloc(3);
	gsl_vector* Trans_vec = gsl_vector_alloc(3);
	
	int n1, n2, N1, N2, bp1, bp2, bp_i, bp_f, int_in_tail1, int_in_tail2;
	double delta_E, rand_MC;
	bool cond;
	
	Try_number = Try_number + 1.0;
	Try_number_pivot = Try_number_pivot + 1.0;
				
	N1=int(free_tail1_bp_number)-1;
	N2=int(free_tail2_bp_number)-1;
	
	int_in_tail1 = int(in_tail1);
	int_in_tail2 = int(in_tail2);
	
	
	n1 = Random_basepair (0, N1);
	n2 = Random_basepair (0, N2);
	
	bp1=gsl_vector_get(free_tail1_bp, n1);
	bp2=gsl_vector_get(free_tail2_bp, n2);
	
	bp_i = int_in_tail1 * 0 + int_in_tail2 * bp2;
	bp_f = int_in_tail1 * bp1 + int_in_tail2 * (N_bp-1);
	
	//max_rot = PI/30.0;
	//max_trans = 1/30.0;
	
	max_rot = PI/60.0;
	max_trans = 1/60.0;
	
	Random_Rot_vector (Rot_vec);
	Random_Trans_vector (Trans_vec);
	
	delta_E = Free_MC_Move(f1, f2, f3, bp_i, bp_f, Rot_vec, Trans_vec, Newd1, Newd2, Newd3, NewR);
	
	rand_MC = RanGen.Random();
	
	cond = false;
	if (delta_E <= 0)
	{cond = true;} 
	else if (rand_MC < exp(-beta*delta_E)) 
	{cond = true;}
			
	if (cond)
	{
		Accept_Pivot_Move (bp_i, bp_f, Newd1, Newd2, Newd3, NewR, delta_E);
	}
	
	if (Error_Code==1)
	{
		cout << "An Error Occured! MC Loop Terminated! \n";
		Save_Err_data ();
		exit(1);
	}
	
	gsl_vector_free(Rot_vec);
	gsl_vector_free(Trans_vec);
}

//Function Glob_Metropolis
void Glob_Metropolis (double f1, double f2, double f3, gsl_vector* NewX_Axis, gsl_vector* NewY_Axis, gsl_vector* NewZ_Axis)
{
	gsl_vector* Rot_vec = gsl_vector_alloc(3);
	
	double delta_E, rand_MC;
	bool cond;
	
	Try_number = Try_number + 1.0;
	Try_number_glob = Try_number_glob + 1.0;
	
	max_rot = PI/20.0;
	
	Random_Rot_vector (Rot_vec);
	
	delta_E = Global_Rotation_MC_Move (f1, f2, f3, Rot_vec, NewX_Axis, NewY_Axis, NewZ_Axis);
	
	rand_MC = RanGen.Random();
	
	cond = false;
	if (delta_E <= 0)
	{cond = true;} 
	else if (rand_MC < exp(-beta*delta_E)) 
	{cond = true;}
			
	if (cond)
	{
		Accept_Glob_Move (NewX_Axis, NewY_Axis, NewZ_Axis, delta_E);
	}
	
	if (Error_Code==1)
	{
		cout << "An Error Occured! MC Loop Terminated! \n";
		Save_Err_data ();
		exit(1);
	}
	
	gsl_vector_free(Rot_vec);
	
}

// Function Bond_Metropolis
void Bond_Metropolis (double f1, double f2, double f3, gsl_matrix* Newd1, gsl_matrix* Newd2, gsl_matrix* Newd3, gsl_matrix* NewR)
{
	gsl_vector* Rot_vec = gsl_vector_alloc(3);
	gsl_vector* Trans_vec = gsl_vector_alloc(3);
	
	int n, bps, N;
	double delta_E, rand_MC;
	bool cond;
	
	Try_number = Try_number + 1.0;
	Try_number_bond = Try_number_bond + 1.0;
				
	N=int(bond_steps_number)-1;
	
	n = Random_basepair (0, N);
	
	bps=gsl_vector_get(bond_steps,n);
	
	//max_rot = PI/70.0;
	//max_trans = 1/70.0;
	
	max_rot = PI/140.0;
	max_trans = 1/140.0;
	
	Random_Rot_vector (Rot_vec);
	Random_Trans_vector (Trans_vec);
	
	delta_E = Bond_MC_Move(f1, f2, f3, bps, Rot_vec, Trans_vec, Newd1, Newd2, Newd3, NewR);
	
	rand_MC = RanGen.Random();
	
	cond = false;
	if (delta_E <= 0)
	{cond = true;} 
	else if (rand_MC < exp(-beta*delta_E)) 
	{cond = true;}
			
	if (cond)
	{
		Accept_Bond_Move (bps, Newd1, Newd2, Newd3, NewR, delta_E);
	}
	
	if (Error_Code==1)
	{
		cout << "An Error Occured! MC Loop Terminated! \n";
		Save_Err_data ();
		exit(1);
	}
	
	gsl_vector_free(Rot_vec);
	gsl_vector_free(Trans_vec);

}

// Function Displace_Phosphate_Metropolis
void Displace_Phosphate_Metropolis (double f1, double f2, double f3, gsl_matrix* Newd1, gsl_matrix* Newd2, gsl_matrix* Newd3, gsl_matrix* NewR, gsl_matrix* Newdm1, gsl_matrix* Newdm2, gsl_matrix* Newdm3, gsl_matrix* NewRm)
{
	gsl_vector* Rot_vec = gsl_vector_alloc(3);
	gsl_vector* Trans_Phos_vec = gsl_vector_alloc(3);
	
	int n, bps, N, strand_3prime;
	double delta_E, rand_MC, Ks_1, Ks_2;
	bool cond;
	
	Try_number = Try_number + 1.0;
	Try_number_DisPhos = Try_number_DisPhos + 1.0;
				
	N=int(bond_steps_number)-1;
	
	n = Random_basepair (0, N);
	
	bps=gsl_vector_get(bond_steps,n);
	Ks_1 = gsl_vector_get(Pos_cons, n);
	Ks_2 = gsl_vector_get(Ang_cons, n);
	strand_3prime = gsl_vector_get(strand_3prime_array, n);
	
	max_rot = PI/45.0;;
	max_trans = 0.0;
	
	Random_Rot_vector (Rot_vec);
	Random_Trans_vector (Trans_Phos_vec);
	
	delta_E = Displace_Phosphate_MC_Move (f1, f2, f3, bps, strand_3prime, Ks_1, Ks_2, Rot_vec, Trans_Phos_vec, Newd1, Newd2, Newd3, NewR, Newdm1, Newdm2, Newdm3, NewRm);
	
	rand_MC = RanGen.Random();
	
	cond = false;
	if (delta_E <= 0)
	{cond = true;} 
	else if (rand_MC < exp(-beta*delta_E)) 
	{cond = true;}
			
	if (cond)
	{
		Accept_Displace_Phosphate_Move (bps, Newd1, Newd2, Newd3, NewR, Newdm1, Newdm2, Newdm3, NewRm, delta_E);
	}
	
	if (Error_Code==1)
	{
		cout << "An Error Occured! MC Loop Terminated! \n";
		Save_Err_data ();
		exit(1);
	}
	
	gsl_vector_free(Rot_vec);
	gsl_vector_free(Trans_Phos_vec);
}

// Function bps_Mutation_Metropolis
void bps_Mutation_Metropolis (gsl_vector* New_bp_Seq_Code, gsl_vector* New_Seq_Code)
{
	int bps, bp_1, bp_2, new_bp_Code_1, new_bp_Code_2;
	double delta_E, rand_MC;
	bool cond;
	
	Try_number = Try_number + 1.0;
	Try_number_Mutation = Try_number_Mutation + 1.0;
	
	bps = Random_basepair (0, (N_bp-2));
	
	bp_1 = bps;
	bp_2 = bps+1;
	
	new_bp_Code_1 = int(1 + 4*RanGen.Random());
	
	new_bp_Code_2 = int(1 + 4*RanGen.Random());
	
	gsl_vector_set(New_bp_Seq_Code, bp_1, new_bp_Code_1);
	gsl_vector_set(New_bp_Seq_Code, bp_2, new_bp_Code_2);
	
	delta_E = Mutation_MC_Move (bp_1, bp_2, New_bp_Seq_Code, New_Seq_Code);
	
	rand_MC = RanGen.Random();
	
	cond = false;
	if (delta_E <= 0)
	{cond = true;} 
	else if (rand_MC < exp(-beta_seq*delta_E)) 
	{cond = true;}
			
	if (cond)
	{
		Accept_Mutation_Move (bp_1, bp_2, New_bp_Seq_Code, New_Seq_Code, delta_E);
	}
	
	if (Error_Code==1)
	{
		cout << "An Error Occured! MC Loop Terminated! \n";
		Save_Err_data ();
		exit(1);
	}
}

// Function Metropolis_step
void Metropolis_step (double f1, double f2, double f3, gsl_matrix* Newd1, gsl_matrix* Newd2, gsl_matrix* Newd3, gsl_matrix* NewR, gsl_matrix* Newdm1, gsl_matrix* Newdm2, gsl_matrix* Newdm3, gsl_matrix* NewRm, gsl_vector* NewX_Axis, gsl_vector* NewY_Axis, gsl_vector* NewZ_Axis, gsl_vector* New_bp_Seq_Code, gsl_vector* New_Seq_Code)
{
	double Per_1, Per_2, Per_3, Per_4, Per_5, Per_6, Per_7, Per_8, Per_9, Per_10;
	double rand_which_move;
	
	Per_1 = 0.0;
	Per_2 = Per_1 + Loc_Move_tail1_Percent;
	Per_3 = Per_2 + Loc_Move_wrapped_Percent;
	Per_4 = Per_3 + Loc_Move_tail2_Percent;
	Per_5 = Per_4 + Pivot_Move_tail1_Percent;
	Per_6 = Per_5 + Pivot_Move_tail2_Percent;
	Per_7 = Per_6 + Glob_Move_Percent;
	Per_8 = Per_7 + Bond_Move_Percent;
	Per_9 = Per_8 + DisPhos_Move_Percent;
	Per_10 = Per_9 + Mut_Move_Percent;
	
	rand_which_move = RanGen.Random();
	
	for (int k = 0; k < N_bp; k++)
	{
		if ((Per_1 <= rand_which_move) && (rand_which_move < Per_2))
		{
			Local_Metropolis (f1, f2, f3, Newd1, Newd2, Newd3, NewR, true, false, false);
		}
		
		if ((Per_2 <= rand_which_move) && (rand_which_move < Per_3))
		{
			Local_Metropolis (f1, f2, f3, Newd1, Newd2, Newd3, NewR, false, true, false);
		}
		
		if ((Per_3 <= rand_which_move) && (rand_which_move < Per_4))
		{
			Local_Metropolis (f1, f2, f3, Newd1, Newd2, Newd3, NewR, false, false, true);
		}
		
		if ((Per_4 <= rand_which_move) && (rand_which_move < Per_5))
		{
			Pivot_Metropolis (f1, f2, f3, Newd1, Newd2, Newd3, NewR, true, false);
		}
		
		if ((Per_5 <= rand_which_move) && (rand_which_move < Per_6))
		{
			Pivot_Metropolis (f1, f2, f3, Newd1, Newd2, Newd3, NewR, false, true);
		}
		
		if ((Per_6 <= rand_which_move) && (rand_which_move < Per_7))
		{
			Glob_Metropolis (f1, f2, f3, NewX_Axis, NewY_Axis, NewZ_Axis);
		}
		
		if ((Per_7 <= rand_which_move) && (rand_which_move < Per_8))
		{
			Bond_Metropolis (f1, f2, f3, Newd1, Newd2, Newd3, NewR);
		}	
		
		if ((Per_8 <= rand_which_move) && (rand_which_move < Per_9))
		{
			Displace_Phosphate_Metropolis (f1, f2, f3, Newd1, Newd2, Newd3, NewR, Newdm1, Newdm2, Newdm3, NewRm);
		}
		
		if ((Per_9 <= rand_which_move) && (rand_which_move < Per_10))
		{
			bps_Mutation_Metropolis (New_bp_Seq_Code, New_Seq_Code);
		}		
	}
	
	Accept_ratio = Accept_number/Try_number;
    Accept_ratio_loc = Accept_number_loc/Try_number_loc;
	Accept_ratio_pivot = Accept_number_pivot/Try_number_pivot;
	Accept_ratio_bond = Accept_number_bond/Try_number_bond;
	Accept_ratio_glob = Accept_number_glob/Try_number_glob;
	Accept_ratio_DisPhos = Accept_number_DisPhos/Try_number_DisPhos;
	Accept_ratio_Mutation = Accept_number_Mutation/Try_number_Mutation;
}

// Function double2string
void double2string (double num, string& strdoub)
{
	ostringstream strs;
    strs << num;
    strdoub = strs.str();
}

// Function int2string
void int2string (int num, string& strint)
{
	ostringstream strs;
    strs << num;
    strint = strs.str();
}

// Function make_directory
void make_directory (int bond_min, int bond_max, double f1, double f2, double f3, string folder_name, int length_tail1, int length_tail2)
{
	string str1_bmin, str2_bmin, str1_bmax, str2_bmax, str_force, str_f1, str_f2, str_f3, str1_tail, str2_tail1, str2_tail2, str1_beta, str2_beta, str3_beta, str_beta, str_dash, str_slash, folder_path_0, folder_path_1, folder_path_2;
	
	str_slash = "/";
	str_dash = " - ";
	str1_bmin = "bond_min = ";
	int2string (bond_min, str2_bmin);
	
	str1_bmax = "bond_max = ";
	int2string (bond_max, str2_bmax);
	
	str_force = "force = ";
	double2string (f1, str_f1);
	double2string (f2, str_f2);
	double2string (f3, str_f3);
	
	str1_tail = "tail = ";
	int2string (length_tail1, str2_tail1);
	int2string (length_tail2, str2_tail2);
	
	str1_beta = "beta = ";
	double2string (beta, str2_beta);
	double2string (beta_seq, str3_beta);
	str_beta = str1_beta + str2_beta + str_dash + str3_beta;
	
	folder_path_0 = folder_name + str_slash + str1_tail + str2_tail1 + str_dash + str2_tail2 + str_dash + str_beta; 
	folder_path_1 = folder_path_0 + str_slash + str_force + str_f1 + str_dash + str_f2 + str_dash + str_f3;
	folder_path_2 = folder_path_1 + str_slash + str1_bmin + str2_bmin + str_dash + str1_bmax + str2_bmax; 
	
	char *str_arr;
	
	str_arr = (char*)folder_name.c_str();
	
	mkdir(str_arr, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	 
	str_arr = (char*)folder_path_0.c_str();
	
	mkdir(str_arr, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	
	str_arr = (char*)folder_path_1.c_str();
	
	mkdir(str_arr, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	
	str_arr = (char*)folder_path_2.c_str();
	
	mkdir(str_arr, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	
	
}

// Function make_string
void make_string (int bond_min, int bond_max, double f1, double f2, double f3, string folder_name, int length_tail1, int length_tail2, int Run_num)
{
	string str1_bmin, str2_bmin, str1_bmax, str2_bmax, str_force, str_f1, str_f2, str_f3, str1_tail, str2_tail1, str2_tail2, str1_beta, str2_beta, str3_beta, str_beta, str1_Run, str2_Run, str_txt, str_folder_num, str_dash, str_slash;
	string folder_path_0, folder_path_1, folder_path_2, str01, str02, str0_d1, str0_d2, str0_d3, str0_R, str0_Ax, str0_En, str0_Ac, str0_bp_Seq, str0_step; 
	string str0_h_Energy, str0_ave_Energy, str0_h_EndtoEnd, str0_ave_EndtoEnd, str0_h_ang1, str0_ave_ang1, str0_h_ang2, str0_ave_ang2, str0_h_ang3, str0_ave_ang3, str0_h_Seq;
	
	str0_d1 = "d_1";
	str0_d2 = "d_2";
	str0_d3 = "d_3";
	str0_R = "R";
	str0_Ax = "Axes";
	str0_En = "Energy";
	str0_Ac = "Accept_ratios";
	str0_Ac = "Accept_ratios";
	str0_bp_Seq = "bp_Seq_Code";
	str0_step = "MC_step.txt";
	
	str0_h_Energy = "hist_Energy";
	str0_ave_Energy = "ave_Energy";
	str0_h_EndtoEnd = "hist_EndtoEnd";
	str0_ave_EndtoEnd = "ave_EndtoEnd";
	str0_h_ang1 = "hist_Euler_Ang1";
	str0_ave_ang1 = "ave_Euler_Ang1";
	str0_h_ang2 = "hist_Euler_Ang2";
	str0_ave_ang2 = "ave_Euler_Ang2";
	str0_h_ang3 = "hist_Euler_Ang3";
	str0_ave_ang3 = "ave_Euler_Ang3";
	str0_h_Seq = "hist_Seq";
	
	str_txt = ".txt";
	str_dash = " - ";
	str_slash = "/";
	
	str1_bmin = "bond_min = ";
	int2string (bond_min, str2_bmin);
	
	str1_bmax = "bond_max = ";
	int2string (bond_max, str2_bmax);
	
	str_force = "force = ";
	double2string (f1, str_f1);
	double2string (f2, str_f2);
	double2string (f3, str_f3);
	
	str1_tail = "tail = ";
	int2string (length_tail1, str2_tail1);
	int2string (length_tail2, str2_tail2);
	
	str1_beta = "beta = ";
	double2string (beta, str2_beta);
	double2string (beta_seq, str3_beta);
	str_beta = str1_beta + str2_beta + str_dash + str3_beta;
	
	str1_Run = "Run = ";
	int2string (Run_num, str2_Run);
	
	folder_path_0 = folder_name + str_slash + str1_tail + str2_tail1 + str_dash + str2_tail2 + str_dash + str_beta; 
	folder_path_1 = folder_path_0 + str_slash + str_force + str_f1 + str_dash + str_f2 + str_dash + str_f3;
	folder_path_2 = folder_path_1 + str_slash + str1_bmin + str2_bmin + str_dash + str1_bmax + str2_bmax; 
	
	str01 = str1_bmin + str2_bmin + str_dash + str1_bmax + str2_bmax + str_dash + str_force + str_f1 + str_dash + str_f2 + str_dash + str_f3 + str_dash + str1_tail + str2_tail1 + str_dash + str2_tail2 + str_dash + str_beta + str_dash + str1_Run + str2_Run + str_txt; 
	str02 = str1_bmin + str2_bmin + str_dash + str1_bmax + str2_bmax + str_dash + str_force + str_f1 + str_dash + str_f2 + str_dash + str_f3 + str_dash + str1_tail + str2_tail1 + str_dash + str2_tail2 + str_dash + str_beta + str_txt; 
	
	str_d1 = folder_path_2 + str_slash + str0_d1 + str_dash + str02;
	str_d2 = folder_path_2 + str_slash + str0_d2 + str_dash + str02;
	str_d3 = folder_path_2 + str_slash + str0_d3 + str_dash + str02;
	str_R = folder_path_2 + str_slash + str0_R + str_dash + str02;
	str_Ax = folder_path_2 + str_slash + str0_Ax + str_dash + str02;
	str_En = folder_path_2 + str_slash + str0_En + str_dash + str02;
	str_Ac = folder_path_2 + str_slash + str0_Ac + str_dash + str01;
	str_bp_Seq = folder_path_2 + str_slash + str0_bp_Seq + str_dash + str01;
	str_step = folder_path_2 + str_slash + str0_step; 
	
	str_h_Energy = folder_path_2  + str_slash + str0_h_Energy + str_dash + str01;
	str_ave_Energy = folder_path_2 + str_slash + str0_ave_Energy + str_dash + str01;
	str_h_EndtoEnd = folder_path_2 + str_slash + str0_h_EndtoEnd + str_dash + str01;
	str_ave_EndtoEnd = folder_path_2 + str_slash + str0_ave_EndtoEnd + str_dash + str01;
	str_h_ang1 = folder_path_2 + str_slash + str0_h_ang1 + str_dash + str01;
	str_ave_ang1 = folder_path_2 + str_slash + str0_ave_ang1 + str_dash + str01;
	str_h_ang2 = folder_path_2 + str_slash + str0_h_ang2 + str_dash + str01;
	str_ave_ang2 = folder_path_2 + str_slash + str0_ave_ang2 + str_dash + str01;
	str_h_ang3 = folder_path_2 + str_slash + str0_h_ang3 + str_dash + str01;
	str_ave_ang3 = folder_path_2 + str_slash + str0_ave_ang3 + str_dash + str01;
	str_h_Seq = folder_path_2 + str_slash + str0_h_Seq + str_dash + str01;
}

// Function Save_data
void Save_data ()
{
	double config_var;
	
	ofstream config_file;
	
	char *str_arr;
	
	str_arr = (char*)str_d1.c_str();
	
	config_file.open (str_arr);
	
	for (int k = 0; k < N_bp; k++)
	{
		for (int i = 0; i < 3; i++)
		{
			config_var = gsl_matrix_get(d1,i,k);
			config_file << config_var << "\t"; 
		}
		config_file << "\n";
	}
	
    config_file.close();
    
    str_arr = (char*)str_d2.c_str();

	config_file.open (str_arr);
	
	for (int k = 0; k < N_bp; k++)
	{
		for (int i = 0; i < 3; i++)
		{
			config_var = gsl_matrix_get(d2,i,k);
			config_file << config_var << "\t"; 
		}
		config_file << "\n";
	}
	
    config_file.close();
    
    str_arr = (char*)str_d3.c_str();

	config_file.open (str_arr);
	
	for (int k = 0; k < N_bp; k++)
	{
		for (int i = 0; i < 3; i++)
		{
			config_var = gsl_matrix_get(d3,i,k);
			config_file << config_var << "\t"; 
		}
		config_file << "\n";
	}
	
    config_file.close();
    
    str_arr = (char*)str_R.c_str();

	config_file.open (str_arr);
	
	for (int k = 0; k < N_bp; k++)
	{
		for (int i = 0; i < 3; i++)
		{
			config_var = gsl_matrix_get(R,i,k);
			config_file << config_var << "\t"; 
		}
		config_file << "\n";
	}
	
    config_file.close();
    
    str_arr = (char*)str_Ax.c_str();

	config_file.open (str_arr);
    
    for (int i = 0; i < 3; i++)
	{
		config_var = gsl_vector_get(X_Axis,i);
		config_file << config_var << "\t"; 
	}
	
	config_file << "\n";
	
	for (int i = 0; i < 3; i++)
	{
		config_var = gsl_vector_get(Y_Axis,i);
		config_file << config_var << "\t"; 
	}
	
	config_file << "\n";
	
	for (int i = 0; i < 3; i++)
	{
		config_var = gsl_vector_get(Z_Axis,i);
		config_file << config_var << "\t"; 
	}
	
	config_file.close();
	
	str_arr = (char*)str_En.c_str();

	config_file.open (str_arr);
	
	config_file << Tot_Energy << "\n"; 
	
	config_file.close();
	
	str_arr = (char*)str_Ac.c_str();

	config_file.open (str_arr);
	
	config_file << Accept_ratio_loc << "\n"; 
	config_file << Accept_ratio_pivot << "\n"; 
	config_file << Accept_ratio_bond << "\n"; 
	config_file << Accept_ratio_glob << "\n"; 
	config_file << Accept_ratio_DisPhos << "\n"; 
	config_file << Accept_ratio_Mutation << "\n"; 
	config_file << Accept_ratio << "\n"; 
	
	config_file.close();
	
	str_arr = (char*)str_bp_Seq.c_str();

	config_file.open (str_arr);
    
    for (int i = 0; i < N_bp; i++)
	{
		config_var = gsl_vector_get(bp_Seq_Code,i);
		config_file << config_var << "\n"; 
	}
	
	config_file.close();
}

// Function save_MC_Step
void save_MC_Step (int b_min, int b_max, double f1, double f2, double f3, int MC_step)
{
	ofstream dat_file;
	
	char *str_arr;
	
	str_arr = (char*)str_step.c_str();
	
	dat_file.open (str_arr);
	
	dat_file << "MC step = " << MC_step << "\n";
	dat_file << "bond_min = " << b_min << "\n";
	dat_file << "bond_max = " << b_max << "\n";
	dat_file << "Force = " << f1 <<",\t" << f2 <<",\t" << f3 << "\n";

	dat_file.close();
}

// Function Calculate_Euler_Angles
void Calculate_Euler_Angles (gsl_vector* Ax_1, gsl_vector* Ax_2, gsl_vector* Ax_3, gsl_vector* Ax0_1, gsl_vector* Ax0_2, gsl_vector* Ax0_3, double& ang_1, double& ang_2, double& ang_3)
{
	double Th, cos1, cos2, cos3, sin1, sin2, sin3;
	
	Normalize_Vec (Ax_1);
	Normalize_Vec (Ax_2);
	Normalize_Vec (Ax_3);
	
	ang_1 = -2;
	ang_2 = -2;
	ang_3 = -2; 
	
	cos1 = dot_prod(Ax_3 , Ax0_3);
	Th = acos(cos1);
	
	if ((Th > 1e-7) && (PI-Th > 1e-7))
	{
		ang_2 = cos1;
		
		cos2 = -1/sin(Th) * dot_prod(Ax_2 , Ax0_3);
		
		sin2 = 1/sin(Th) * dot_prod(Ax_1 , Ax0_3);
		
		ang_1 = 1/(PI) * acos(cos2);
		
		if (sin2 < 0)
		{
			ang_1 = -ang_1;
		}
		
		cos3 = 1/sin(Th) * dot_prod(Ax_3 , Ax0_2);
		
		sin3 = 1/sin(Th) * dot_prod(Ax_3 , Ax0_1);
		
		ang_3 = 1/(PI) * acos(cos3);
		
		if (sin3 < 0)
		{
			ang_3 = -ang_3;
		}
	}
}

// Function Calculate_End_to_End
double Calculate_End_to_End (gsl_vector* Ax)
{
	double Z;
	
	gsl_vector* Vec = gsl_vector_alloc(3);
	
	for (int i=0; i<3; i++)
	{
		gsl_vector_set(Vec, i, (gsl_matrix_get(R,i,N_bp-1)- gsl_matrix_get(R,i,0)));
	}
	
	Z = 1/L*dot_prod(Vec , Ax);
	
	gsl_vector_free(Vec);
	
	return (Z);
}

// Data Analysis

// Function Initial_hist
void Initial_hist (gsl_vector* var_array, gsl_vector* hist, double hist_delta, double min_var, int hist_num, double& sum_var, double& sum_sq_var, int& dat_num)
{
	for (int j = 0; j < hist_num; j++)
	{
		gsl_vector_set(var_array,j, (min_var+(j+0.5)*hist_delta));
	}
	
	gsl_vector_set_zero(hist);
	
	sum_var = 0;
	sum_sq_var = 0;
	dat_num = 0;
}

// Function Calculate_Hist_and_Ave
void Calculate_Hist_and_Ave (gsl_vector* hist, double hist_delta, double min_var, int hist_num, double& sum_var, double& sum_sq_var, int& dat_num, double new_var)
{
	int bin_num;

	bin_num = int ((new_var-min_var)/hist_delta);

	if ((bin_num >=0) && (bin_num <= hist_num))
	{
		gsl_vector_set(hist, bin_num, (gsl_vector_get(hist, bin_num)+1) );
	
		sum_var = sum_var + new_var;
		sum_sq_var = sum_sq_var + new_var*new_var;
		dat_num++;
	}
	else
	{
		cout << "Data Out Of Range\n";
		cout << new_var << "\n";
	}
}

// Function Calculate_Seq_Hist
void Calculate_Seq_Hist ()
{
	int bps_code;
	
	for (int n = 0; n < (N_bp-1); n++)
	{
		bps_code = gsl_vector_get(Seq_Code, n);
		
		gsl_matrix_set (hist_Seq, n, bps_code, (gsl_matrix_get (hist_Seq, n, bps_code) + 1) );
	}
}

// Function Save_Hist_and_Ave
void Save_Hist_and_Ave(gsl_vector* var_array, gsl_vector* hist, int hist_num, double ave_var, double sigma_var, string str_hist, string str_ave)
{
	double dat_x, dat_y;
	
	char *str_arr;
	
	ofstream hist_file;
	
	str_arr = (char*)str_hist.c_str();
	
	hist_file.open (str_arr);
	
	for (int j = 0; j < hist_num; j++)
	{
		dat_x = gsl_vector_get(var_array,j);
		dat_y = gsl_vector_get(hist,j);
		
		hist_file << dat_x << "\t";
		hist_file << dat_y << "\n";
	}
	
	hist_file.close();
	
	str_arr = (char*)str_ave.c_str();
	
	hist_file.open (str_arr);
	
	hist_file << ave_var << "\n" << sigma_var << "\n";
	
	hist_file.close();
}

// Function Reset_All_Histograms
void Reset_All_Histograms ()
{
	Initial_hist (Energy_array, hist_Energy, hist_Energy_delta, min_Energy, hist_Energy_num, sum_Energy, sum_sq_Energy, Energy_num);
	
	/*
	
	Initial_hist (EndtoEnd_array, hist_EndtoEnd, hist_EndtoEnd_delta, min_EndtoEnd, hist_EndtoEnd_num, sum_EndtoEnd, sum_sq_EndtoEnd, EndtoEnd_num);
	
	Initial_hist (ang_1_array, hist_ang_1, hist_ang_1_delta, min_ang_1, hist_ang_1_num, sum_ang_1, sum_sq_ang_1, ang_1_num);
	
	Initial_hist (ang_2_array, hist_ang_2, hist_ang_2_delta, min_ang_2, hist_ang_2_num, sum_ang_2, sum_sq_ang_2, ang_2_num);
	
	Initial_hist (ang_3_array, hist_ang_3, hist_ang_3_delta, min_ang_3, hist_ang_3_num, sum_ang_3, sum_sq_ang_3, ang_3_num);
	*/ 
}

// Function Save_All_Hist_and_Ave
void Save_All_Hist_and_Ave ()
{
	Save_Hist_and_Ave(Energy_array, hist_Energy, hist_Energy_num, ave_Energy, sigma_Energy, str_h_Energy, str_ave_Energy);
	
	/*
	
	Save_Hist_and_Ave(EndtoEnd_array, hist_EndtoEnd, hist_EndtoEnd_num, ave_EndtoEnd, sigma_EndtoEnd, str_h_EndtoEnd, str_ave_EndtoEnd);
	
	Save_Hist_and_Ave(ang_1_array, hist_ang_1, hist_ang_1_num, ave_ang_1, sigma_ang_1, str_h_ang1, str_ave_ang1);
	
	Save_Hist_and_Ave(ang_2_array, hist_ang_2, hist_ang_2_num, ave_ang_2, sigma_ang_2, str_h_ang2, str_ave_ang2);
	
	Save_Hist_and_Ave(ang_3_array, hist_ang_3, hist_ang_3_num, ave_ang_3, sigma_ang_3, str_h_ang3, str_ave_ang3);
	*/ 
}

// Function Save_Hist_Seq
void Save_Hist_Seq ()
{
	double hist_var;
	
	char *str_arr;
	
	ofstream hist_file;
	
	str_arr = (char*)str_h_Seq.c_str();
	
	hist_file.open (str_arr);
	
	for (int n = 0; n < (N_bp-1); n++)
	{
		for (int i = 0; i < 16; i++)
		{
			hist_var = gsl_matrix_get (hist_Seq, n, i);
			
			hist_file << hist_var << "\t";
		} 
		
		hist_file << "\n";
	}
	
	hist_file.close();
}

// Function Data_Analysis
void  Data_Analysis (int count)
{
	/*
	
	double Z, ang_1, ang_2, ang_3;
	
	Z = Calculate_End_to_End (Y_Axis);
	
	Calculate_Euler_Angles (X_Axis, Z_Axis, Y_Axis, X0_Axis, Z0_Axis, Y0_Axis, ang_1, ang_2, ang_3);
	
	Calculate_Hist_and_Ave (hist_EndtoEnd, hist_EndtoEnd_delta, min_EndtoEnd, hist_EndtoEnd_num, sum_EndtoEnd, sum_sq_EndtoEnd, EndtoEnd_num, Z);
	
	Calculate_Hist_and_Ave (hist_ang_1, hist_ang_1_delta, min_ang_1, hist_ang_1_num, sum_ang_1, sum_sq_ang_1, ang_1_num, ang_1);
	
	Calculate_Hist_and_Ave (hist_ang_2, hist_ang_2_delta, min_ang_2, hist_ang_2_num, sum_ang_2, sum_sq_ang_2, ang_2_num, ang_2);
	
	Calculate_Hist_and_Ave (hist_ang_3, hist_ang_3_delta, min_ang_3, hist_ang_3_num, sum_ang_3, sum_sq_ang_3, ang_3_num, ang_3);
	*/
	
	Calculate_Hist_and_Ave (hist_Energy, hist_Energy_delta, min_Energy, hist_Energy_num, sum_Energy, sum_sq_Energy, Energy_num, Tot_Energy);
	
	ave_Energy = sum_Energy/Energy_num;
	sigma_Energy = sqrt(sum_sq_Energy/Energy_num - ave_Energy*ave_Energy);
	
	Calculate_Seq_Hist ();
	
	/*
	
	ave_EndtoEnd = sum_EndtoEnd/EndtoEnd_num;
	sigma_EndtoEnd = sqrt(sum_sq_EndtoEnd/EndtoEnd_num - ave_EndtoEnd*ave_EndtoEnd);
	
	ave_ang_1 = sum_ang_1/ang_1_num;
	sigma_ang_1 = sqrt(sum_sq_ang_1/ang_1_num - ave_ang_1*ave_ang_1);
	
	ave_ang_2 = sum_ang_2/ang_2_num;
	sigma_ang_2 = sqrt(sum_sq_ang_2/ang_2_num - ave_ang_2*ave_ang_2);
	
	ave_ang_3 = sum_ang_3/ang_3_num;
	sigma_ang_3 = sqrt(sum_sq_ang_3/ang_3_num - ave_ang_3*ave_ang_3);
	*/
	
	if (count % 100 == 0)
	{
		Save_All_Hist_and_Ave ();
		
		Save_Hist_Seq ();
	}
}

// Function Do_Monte_Carlo
void Do_Monte_Carlo (int b_min, int b_max, double f1, double f2, double f3, string folder_name, int length_tail1, int length_tail2)
{
	gsl_matrix* Newd1 = gsl_matrix_alloc(3,N_bp);
	gsl_matrix* Newd2 = gsl_matrix_alloc(3,N_bp);
	gsl_matrix* Newd3 = gsl_matrix_alloc(3,N_bp);
	gsl_matrix* NewR = gsl_matrix_alloc(3,N_bp);
	
	gsl_matrix* Newdm1 = gsl_matrix_alloc(3,N_bp-1);
	gsl_matrix* Newdm2 = gsl_matrix_alloc(3,N_bp-1);
	gsl_matrix* Newdm3 = gsl_matrix_alloc(3,N_bp-1);
	gsl_matrix* NewRm = gsl_matrix_alloc(3,N_bp-1);
	
	gsl_vector* NewX_Axis = gsl_vector_alloc(3);
	gsl_vector* NewY_Axis = gsl_vector_alloc(3);
	gsl_vector* NewZ_Axis = gsl_vector_alloc(3);
	
	gsl_vector* New_Seq_Code = gsl_vector_alloc(N_bp-1);
    gsl_vector* New_bp_Seq_Code = gsl_vector_alloc(N_bp);
	
	gsl_vector* force_Vec = gsl_vector_alloc(3);
	
	gsl_matrix_set_zero(Newd1);
	gsl_matrix_set_zero(Newd2);
	gsl_matrix_set_zero(Newd3);
	gsl_matrix_set_zero(NewR);
	
	gsl_matrix_set_zero(Newdm1);
	gsl_matrix_set_zero(Newdm2);
	gsl_matrix_set_zero(Newdm3);
	gsl_matrix_set_zero(NewRm);
	
	gsl_vector_set_zero(NewX_Axis);
	gsl_vector_set_zero(NewY_Axis);
	gsl_vector_set_zero(NewZ_Axis);
	
	gsl_vector_set_zero(New_Seq_Code);
	gsl_vector_set_zero(New_bp_Seq_Code);
	
	Try_number = 0.0;
	Try_number_loc = 0.0;
	Try_number_pivot = 0.0;
	Try_number_bond = 0.0;
	Try_number_glob = 0.0;
	Try_number_DisPhos = 0.0;
	Try_number_Mutation = 0.0;

	Accept_number = 0.0;
	Accept_number_loc = 0.0;
	Accept_number_pivot = 0.0;
	Accept_number_bond = 0.0;
	Accept_number_glob = 0.0;
	Accept_number_DisPhos = 0.0;
	Accept_number_Mutation = 0.0;
	
	Error_Code = 0;
	
	int final_count = 5000;
	
	Load_init_config (length_tail1, length_tail2);
	
	Load_Binding_Sites (b_min, b_max, length_tail1);
	
	Load_Elastic_Parameters ();
	
	Update_Force_vec (f1, f2, f3, X_Axis, Y_Axis,  Z_Axis, force_Vec);
	
	init_bp_Seq_Code ();
	
	bp_Seq_Code_to_Seq_Code (0, (N_bp-2), bp_Seq_Code, Seq_Code);
	
	Tot_Energy=Total_Energy(force_Vec);
    
	cout << "bond_min = " << b_min << "\n";
	cout << "bond_max = " << b_max << "\n";
	cout << "Force = " << f1 <<",\t" << f2 <<",\t" << f3 << "\n";
	cout<<"Total_Energy = "<< Tot_Energy <<"\n";
    
    cout<<"Equilibration MC Loop Started! \n";
    
    for (int MC_step = 0; MC_step < N_EQ_MCstep; MC_step++)
	{
		Metropolis_step (f1, f2, f3, Newd1, Newd2, Newd3, NewR, Newdm1, Newdm2, Newdm3, NewRm, NewX_Axis, NewY_Axis, NewZ_Axis, New_bp_Seq_Code, New_Seq_Code);
		
		if (MC_step % 1000 == 0)
		{
			cout<< MC_step;
			cout<<"\n";
		}
	}
	
	cout<<"Equilibration MC Loop Ended Successfully!\n";
    
	cout<<"Total Acceptance Ratio = ";
	cout<< Accept_ratio;
	cout<<"\n";
	
	cout<<"Local Acceptance Ratio = ";
	cout<< Accept_ratio_loc;
	cout<<"\n";
	
	cout<<"Pivot Acceptance Ratio = ";
	cout<< Accept_ratio_pivot;
	cout<<"\n";
	
	cout<<"Bond Acceptance Ratio = ";
	cout<< Accept_ratio_bond;
	cout<<"\n";
	
	cout<<"Global Acceptance Ratio = ";
	cout<< Accept_ratio_glob;
	cout<<"\n";
	
	cout<<"Displace Phosphate Acceptance Ratio = ";
	cout<< Accept_ratio_DisPhos;
	cout<<"\n";
	
	cout<<"Mutation Acceptance Ratio = ";
	cout<< Accept_ratio_Mutation;
	cout<<"\n";
	
	cout<<"Main MC Loop Started! \n";
	
		
	make_directory (b_min, b_max, f1, f2, f3, folder_name, length_tail1, length_tail2);
	
	int Run = 0;
	
	int count =0;
	
	make_string (b_min, b_max, f1, f2, f3, folder_name, length_tail1, length_tail2, Run);
	
	Reset_All_Histograms ();
	
	gsl_matrix_set_zero(hist_Seq);
	
	for (int MC_step = 0; MC_step < N_Main_MCstep; MC_step++)
	{
		Metropolis_step (f1, f2, f3, Newd1, Newd2, Newd3, NewR, Newdm1, Newdm2, Newdm3, NewRm, NewX_Axis, NewY_Axis, NewZ_Axis, New_bp_Seq_Code, New_Seq_Code);
		
		if (MC_step % 1000 == 0)
		{	
			cout<< MC_step <<"\n";
			
			save_MC_Step (b_min, b_max, f1, f2, f3, MC_step);
		}
		
		if (MC_step % 10 ==0)
		{
			count ++;
			
			Data_Analysis (count);
		}
		
		if (count == final_count)
		{
			Reset_All_Histograms ();
			
			gsl_matrix_set_zero(hist_Seq);
				
			Save_data ();
				
			Run++;
				
			count = 0;
				
			make_string (b_min, b_max, f1, f2, f3, folder_name, length_tail1, length_tail2, Run);
			}
	}
	
	cout<<"Main MC Loop Ended Successfully!\n";
    
	cout<<"Total Acceptance Ratio = ";
	cout<< Accept_ratio;
	cout<<"\n";
	
	cout<<"Local Acceptance Ratio = ";
	cout<< Accept_ratio_loc;
	cout<<"\n";
	
	cout<<"Pivot Acceptance Ratio = ";
	cout<< Accept_ratio_pivot;
	cout<<"\n";
	
	cout<<"Bond Acceptance Ratio = ";
	cout<< Accept_ratio_bond;
	cout<<"\n";
	
	cout<<"Global Acceptance Ratio = ";
	cout<< Accept_ratio_glob;
	cout<<"\n";
	
	cout<<"Displace Phosphate Acceptance Ratio = ";
	cout<< Accept_ratio_DisPhos;
	cout<<"\n";
	
	cout<<"Mutation Acceptance Ratio = ";
	cout<< Accept_ratio_Mutation;
	cout<<"\n";
	
	cout<<"Total_Energy = "<< Tot_Energy <<"\n";
	
	gsl_matrix_free(Newd1);
	gsl_matrix_free(Newd2);
	gsl_matrix_free(Newd3);
	gsl_matrix_free(NewR);
	
	gsl_matrix_free(Newdm1);
	gsl_matrix_free(Newdm2);
	gsl_matrix_free(Newdm3);
	gsl_matrix_free(NewRm);
	
	gsl_vector_free(NewX_Axis);
	gsl_vector_free(NewY_Axis);
	gsl_vector_free(NewZ_Axis);
	
	gsl_vector_free(New_Seq_Code);
	gsl_vector_free(New_bp_Seq_Code);
	
	gsl_vector_free(force_Vec);
}
 
// Function main
int main ()
{
	double force_1, force_2, force_3, final_Energy;
	string folder_name;
	int bond_min, bond_max;

	gsl_vector* force_Vec = gsl_vector_alloc(3);
	gsl_vector* force_2_array = gsl_vector_alloc(7);
				
	force_1 = 0.0;
	force_2 = 0.0;
	force_3 = 0.0;
	
	bond_min = 1;
	bond_max = 14;
	
	folder_name = "Data_1";
	
	length_tail1 = 0;
	length_tail2 = N_bp - N_bp0 - length_tail1;
	
	N_EQ_MCstep = 1e4;
	N_Main_MCstep = 1e7;
		
	Do_Monte_Carlo (bond_min, bond_max, force_1, force_2, force_3, folder_name, length_tail1, length_tail2);
				
	Update_Force_vec (force_1, force_2, force_3, X_Axis, Y_Axis,  Z_Axis, force_Vec);
				
	for (int k = 0; k < (N_bp-1); k++)
	{
		Mid_frame (k);
	}
	
	bp_Seq_Code_to_Seq_Code (0, (N_bp-2), bp_Seq_Code, Seq_Code);
				
	final_Energy=Total_Energy(force_Vec);
	
	cout << "final Energy = " << final_Energy<<"\n";
	
	Save_data ();
	
	return 0;
}
