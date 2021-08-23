#ifndef PHASE_FIT
#define PHASE_FIT

#include <stdio.h> 
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <math.h>

#define RELATIVE_TOLERANCE 1e-12
#define MAX_EVALS 1000
#define SCALE_OUT 1e20
#define DEBUG_ 0

typedef struct{
	size_t n;
	gsl_complex* z;
}circle_t;

typedef struct{
	size_t n;
	double * freq;
	double * phase;
}data_t;

typedef struct{
	size_t n;
	double * freq;
	gsl_complex * z;
	double fr_init;
	double Qr_init;
	double Qc_init;
	double a_mag_init;
	double fr_scale;
}data_fine_t;

void performFit(int num, double* freq, double* phase, double* theta, double* Qr, double* fr);

void refineFit(int num, double* freq, gsl_complex* z, gsl_complex* a, double* Qr, double* Qc, double* fr, double* phi0, double* delay);

void fitCircleToData(int num, gsl_complex* z, double* radius, gsl_complex* z_center);

#endif
