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
#include <nlopt.h>

#define Q_INIT 300
#define THETA_INIT 0.0
#define F_INIT 7.0
#define RELATIVE_TOLERANCE 1e-25
#define ABS_TOLERANCE_RESIDUE 1e-25
#define MAX_EVALS 1000000
#define SET_REASONABLE_LIMITS 1
#define SCALE_REASONABLE_LIMITS 4.0


typedef struct{
	size_t n;
	double * freq;
	double * phase;
}data_t;

void performFit(int num, double* freq, double* phase, double* theta, double* Qr, double* fr);

double phaseFunc_ErrorSquared(unsigned n_params, const double* params, double* grad, void* data); // Least Squares function to be minimized

#endif
