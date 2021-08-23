#include "phase_fit.h"
#define DEBUG 0

int _phaseFunc_f(const gsl_vector* params, void* data, gsl_vector* f)  // Rough Least Square Function
{
	data_t* d = (data_t*) data;
	int num = d->n;
	double* freq = d->freq;
	double* phase = d->phase;

	double theta0 = gsl_vector_get(params,0);
	double Qr = gsl_vector_get(params,1);
	double fr = gsl_vector_get(params,2);

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < num; ++i)
	{
		double phaseDiff = 2.0*atan(2.0*Qr*(1.0-(freq[i]/fr))) - theta0 - phase[i] ;
		
		while(phaseDiff > M_PI)
		{
			phaseDiff -= 2.0*M_PI;
		}

		while(phaseDiff <= -M_PI)
		{
			phaseDiff += 2.0*M_PI;
		}

		gsl_vector_set(f,i, phaseDiff);
	}

	return GSL_SUCCESS;
}

int _phaseFunc_df(const gsl_vector* params, void* data, gsl_matrix* J) // Jacobian of rough least squares function
{
	data_t* d = (data_t*) data;
	int num = d->n;
	double* freq = d->freq;
	double* phase = d->phase;

	double theta0 = gsl_vector_get(params,0);
	double Qr = gsl_vector_get(params,1);
	double fr = gsl_vector_get(params,2);

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < num; ++i)
	{
		gsl_matrix_set(J, i, 0, -1.);

		double phaseDiff = 2.0*atan(2.0*Qr*(1.0-(freq[i]/fr))) - theta0 - phase[i] ;
		
		while(phaseDiff >= 2.0*M_PI)
		{
			phaseDiff -= 2.0*M_PI;
		}

		while(phaseDiff < 0)
		{
			phaseDiff += 2.0*M_PI;
		}

	}

	return GSL_SUCCESS;
}

int _fullLorentz_f(const gsl_vector* params, void* data, gsl_vector* f)
{
	data_fine_t* d = (data_fine_t*) data;
	int num = d->n;
	double* freq = d->freq;
	gsl_complex* z = d->z;
	double fr_init = d->fr_init;
	double Qr_init = d->Qr_init;
	double Qc_init = d->Qc_init;
	double fr_scale = d->fr_scale;

	gsl_complex onePlusEpsilon = gsl_complex_polar(gsl_vector_get(params,0)*(d->a_mag_init), gsl_vector_get(params,1));
	double Qr = (gsl_vector_get(params,2) + 1.)*Qr_init;
	double Qc = (gsl_vector_get(params,3) + 1.)*Qc_init;
	double fr = fr_scale*gsl_vector_get(params,4) + fr_init;
	double phi0 = gsl_vector_get(params,5);
	double delay = gsl_vector_get(params,6)/fr_init;
	double a_arg = gsl_vector_get(params,1);

	if(a_arg > M_PI || a_arg < -M_PI)
	{
		double* p_a_arg = &(params->data[1*params->stride]);

		while(*p_a_arg < -M_PI)
		{
			*p_a_arg += 2*M_PI;
		}
		while(*p_a_arg > M_PI)
		{
			*p_a_arg -= 2*M_PI;
		}
	}

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < num; ++i)
	{
		double a = Qr/Qc;
		double b = 2.0*Qr*(freq[i]-fr)/fr;

		gsl_complex rot = gsl_complex_polar(1,-2.0*M_PI*freq[i]*delay);

		gsl_complex t = gsl_complex_polar(a,phi0);
		t = gsl_complex_div(t,gsl_complex_rect(1,b));
		t = gsl_complex_mul(rot,t);
		t = gsl_complex_sub(rot,t);
		t = gsl_complex_mul(onePlusEpsilon,t);

		t = gsl_complex_sub(t,z[i]);
		gsl_vector_set(f,2*i,GSL_REAL(t));
		gsl_vector_set(f,2*i+1,GSL_IMAG(t));
	}

	return GSL_SUCCESS;
}

int _fullLorentz_df(const gsl_vector* params, void* data, gsl_matrix* J)
{
	data_fine_t* d = (data_fine_t*) data;
	int num = d->n;
	double* freq = d->freq;
	double fr_init = d->fr_init;
	double Qr_init = d->Qr_init;
	double Qc_init = d->Qc_init;
	double a_mag_init = d->a_mag_init;
	double fr_scale = d->fr_scale;

	double a_mag = gsl_vector_get(params,0)*a_mag_init;
	double a_arg = gsl_vector_get(params,1);
	double Qr = (gsl_vector_get(params,2) + 1.)*Qr_init;
	double Qc = (gsl_vector_get(params,3) + 1.)*Qc_init;
	double fr = fr_scale*gsl_vector_get(params,4) + fr_init;
	double phi0 = gsl_vector_get(params,5);
	double delay = gsl_vector_get(params,6)/fr_init;

	if(a_arg > M_PI || a_arg < -M_PI)
	{
		double* p_a_arg = &(params->data[1*params->stride]);

		while(*p_a_arg < -M_PI)
		{
			*p_a_arg += 2*M_PI;
		}
		while(*p_a_arg > M_PI)
		{
			*p_a_arg -= 2*M_PI;
		}
	}

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < num; ++i)
	{
		double a = Qr/Qc;
		gsl_complex b = gsl_complex_rect(1,2.0*Qr*(freq[i]-fr)/fr);

		gsl_complex rot = gsl_complex_polar(a_mag, a_arg-2.0*M_PI*freq[i]*delay);

		gsl_complex t = gsl_complex_polar(a,phi0);
		t = gsl_complex_div(t,b);
		t = gsl_complex_mul(rot,t);

		gsl_matrix_set(J,2*i,2,-GSL_REAL(gsl_complex_div(t,b))*Qr_init/Qr);
		gsl_matrix_set(J,2*i+1,2,-GSL_IMAG(gsl_complex_div(t,b))*Qr_init/Qr);

		gsl_matrix_set(J,2*i,3,GSL_REAL(t)*Qc_init/Qc);
		gsl_matrix_set(J,2*i+1,3,GSL_IMAG(t)*Qc_init/Qc);

		gsl_matrix_set(J,2*i,4,GSL_IMAG(gsl_complex_div(t,b))*2.0*Qr*freq[i]*fr_scale/(fr*fr));
		gsl_matrix_set(J,2*i+1,4,-GSL_REAL(gsl_complex_div(t,b))*2.0*Qr*freq[i]*fr_scale/(fr*fr));

		gsl_matrix_set(J,2*i,5,GSL_IMAG(t));
		gsl_matrix_set(J,2*i+1,5,-GSL_REAL(t));

		t = gsl_complex_sub(rot,t);

		gsl_matrix_set(J,2*i,0,GSL_REAL(t)*a_mag_init/a_mag);
		gsl_matrix_set(J,2*i+1,0,GSL_IMAG(t)*a_mag_init/a_mag);
		gsl_matrix_set(J,2*i,1,-GSL_IMAG(t));
		gsl_matrix_set(J,2*i+1,1,GSL_REAL(t));
		gsl_matrix_set(J,2*i,6,GSL_IMAG(t)*2.0*M_PI*freq[i]/fr_init);
		gsl_matrix_set(J,2*i+1,6,-GSL_REAL(t)*2.0*M_PI*freq[i]/fr_init);
	}

	return GSL_SUCCESS;
}

int _circleFunc_f(const gsl_vector* params, void* data, gsl_vector* f)  // Circles' Least Squares function to be minimized
{
	circle_t* circ = (circle_t*) data;
	int num = circ->n;
	gsl_complex* z = circ->z;

	double x_center = gsl_vector_get(params,0);
	double y_center = gsl_vector_get(params,1);
	double radius = gsl_vector_get(params,2);

	for (int i = 0; i < num; ++i)
	{
		gsl_vector_set(f, i, (GSL_REAL(z[i]) - x_center)*(GSL_REAL(z[i]) - x_center) + (GSL_IMAG(z[i]) - y_center)*(GSL_IMAG(z[i]) - y_center) - radius*radius);
	}

	return GSL_SUCCESS;
}

int _circleFunc_df(const gsl_vector* params, void* data, gsl_matrix* J) // Jacobian of least squares function to be minimized
{
	circle_t* circ = (circle_t*) data;
	int num = circ->n;
	gsl_complex* z = circ->z;

	double x_center = gsl_vector_get(params,0);
	double y_center = gsl_vector_get(params,1);
	double radius = gsl_vector_get(params,2);

	for (int i = 0; i < num; ++i)
	{
		gsl_matrix_set(J, i, 0, 2*(x_center - GSL_REAL(z[i])));
		gsl_matrix_set(J, i, 1, 2*(y_center - GSL_IMAG(z[i])));
		gsl_matrix_set(J, i, 2, -2*radius);
	}

	return GSL_SUCCESS;
}

void _getInitValues_Kasa(int num, gsl_complex* z, double* radius, gsl_complex* z_center) // Using the fitting method as described by Kasa et al and Dale Umbach, Kerry N. Jones in "A Few Methods of Fitting Circles to Data"
{
	double a,b,c,d,e;

	double sum_xi=0, sum_yi=0, sum_xi2=0, sum_yi2=0, sum_xiyi=0, sum_xi3=0, sum_yi3=0, sum_xiyi2=0, sum_yixi2=0;

	for (int i = 0; i < num; ++i)
	{
		double x,y,x2,y2,x3,y3;
		x = z[i].dat[0];
		y = z[i].dat[1];
		x2 = x*x;
		y2 = y*y;
		x3 = x2*x;
		y3 = y2*y;

		sum_xi += x;
		sum_yi += y;
		sum_xi2 += x2;
		sum_yi2 += y2;
		sum_xiyi += x*y;
		sum_xi3 += x3;
		sum_yi3 += y3;
		sum_xiyi2 += x*y2;
		sum_yixi2 += y*x2;
	}

	a = num*sum_xi2 - sum_xi*sum_xi;
	b = num*sum_xiyi - sum_xi*sum_yi;
	c = num*sum_yi2 - sum_yi*sum_yi;
	d = 0.5*( num*sum_xiyi2 - sum_xi*sum_yi2 + num*sum_xi3 - sum_xi*sum_xi2 );
	e = 0.5*( num*sum_yixi2 - sum_yi*sum_xi2 + num*sum_yi3 - sum_yi*sum_yi2 );

	GSL_SET_COMPLEX(z_center, ((d*c - b*e)/(a*c - b*b)), ((a*e-b*d)/(a*c - b*b)));

	double r = 0;

	for (int i = 0; i < num; ++i)
	{
		r+= (z[i].dat[0] - z_center->dat[0])*(z[i].dat[0] - z_center->dat[0]) + (z[i].dat[1] - z_center->dat[1])*(z[i].dat[1] - z_center->dat[1]);
	}

	*radius = sqrt(r/num);
}

void performFit(int num, double* freq, double* phase, double* theta, double* Qr, double* fr)
{
	data_t data = {num, freq, phase};

	double init_guesses[3] = { *theta, *Qr, *fr };  /* `*`some` `initial` `guess`*` */

	double chisq, chisq0;
	int status, info;

	gsl_vector *f;
	gsl_matrix *J;
	gsl_matrix *covar = gsl_matrix_alloc (3, 3);

	const gsl_multifit_nlinear_type* T = gsl_multifit_nlinear_trust;
	gsl_multifit_nlinear_parameters params = gsl_multifit_nlinear_default_parameters();
	gsl_multifit_nlinear_workspace* w = gsl_multifit_nlinear_alloc(T, &params, num, 3);

	gsl_vector_view init_guess = gsl_vector_view_array(init_guesses,3);

	gsl_multifit_nlinear_fdf rough_phase;
	rough_phase.f = _phaseFunc_f;
	rough_phase.df = NULL;
	rough_phase.fvv = NULL;
	rough_phase.n = num;
	rough_phase.p = 3; // Using 3 parameters
	rough_phase.params = &data;

	gsl_multifit_nlinear_init(&init_guess.vector, &rough_phase, w);
	
	f = gsl_multifit_nlinear_residual(w);
	gsl_blas_ddot(f,f,&chisq0);

	status = gsl_multifit_nlinear_driver(MAX_EVALS, 1e-20, 1e-20, 0, NULL, NULL, &info, w);

	J = gsl_multifit_nlinear_jac(w);
	gsl_multifit_nlinear_covar(J, 0.0, covar);

	gsl_blas_ddot(f, f, &chisq);

	*theta = gsl_vector_get(w->x, 0);
	*Qr = gsl_vector_get(w->x, 1);
	*fr = gsl_vector_get(w->x, 2);	

	if(DEBUG_)		
	{
		fprintf(stderr, "----------------------\n");
		fprintf(stderr, "Rough Phase Fitting\n");
		fprintf(stderr, "number of iterations: %zu\n",gsl_multifit_nlinear_niter(w));
		fprintf(stderr, "function evaluations: %zu\n", rough_phase.nevalf);
		fprintf(stderr, "Jacobian evaluations: %zu\n", rough_phase.nevaldf);
		fprintf(stderr, "reason for stopping: %s\n", (info == 1) ? "small step size" : "small gradient");
		fprintf(stderr, "initial |f(x)| = %f\n", sqrt(chisq0));
		fprintf(stderr, "final |f(x)| = %f\n", sqrt(chisq));

		{
			double dof = num - 3;
			double c = GSL_MAX_DBL(1, sqrt(chisq / dof));
			fprintf(stderr, "chisq/dof = %g\n", chisq / dof);
			fprintf(stderr, "theta = %e +/- %e\n", gsl_vector_get(w->x, 0), c*sqrt(gsl_matrix_get(covar,0,0)));
			fprintf(stderr, "Qr = %lf +/- %lf\n", gsl_vector_get(w->x, 1), c*sqrt(gsl_matrix_get(covar,1,1)));
			fprintf(stderr, "fr = %lf +/- %lf\n", gsl_vector_get(w->x, 2), c*sqrt(gsl_matrix_get(covar,2,2)));
		}
		fprintf(stderr, "----------------------\n");
	}

	gsl_multifit_nlinear_free(w);
	gsl_matrix_free (covar);
}

void refineFit(int num, double* freq, gsl_complex* z, gsl_complex* a, double* Qr, double* Qc, double* fr, double* phi0, double* delay)
{
	// a would be onePlusEpsilon
	double a_mag = gsl_complex_abs(*a);
	double a_arg = gsl_complex_arg(*a);

	while(a_arg <= -M_PI)
	{
		a_arg += 2.*M_PI;
	}

	while(a_arg > M_PI)
	{
		a_arg -= 2.*M_PI;
	}

	data_fine_t data = {num, freq, z, *fr, *Qr, *Qc, a_mag, (*fr)/(*Qr)};

	double init_guesses[7] = { 1., a_arg, 0, 0, 0, *phi0, *delay};  /* `*`some` `initial` `guess`*` */

	double chisq, chisq0;
	int status, info;

	gsl_vector *f;
	gsl_matrix *J;
	gsl_matrix *covar = gsl_matrix_alloc (7, 7);

	const gsl_multifit_nlinear_type* T = gsl_multifit_nlinear_trust;
	gsl_multifit_nlinear_parameters params = gsl_multifit_nlinear_default_parameters();

	// Refining Parameters
	params.scale = gsl_multifit_nlinear_scale_more;
	params.solver = gsl_multifit_nlinear_solver_svd;
	params.fdtype = GSL_MULTIFIT_NLINEAR_CTRDIFF;

	gsl_multifit_nlinear_workspace* w = gsl_multifit_nlinear_alloc(T, &params, 2*num, 7);

	gsl_vector_view init_guess = gsl_vector_view_array(init_guesses,7);

	gsl_multifit_nlinear_fdf fine_phase;
	fine_phase.f = _fullLorentz_f;
	fine_phase.df = _fullLorentz_df;
	fine_phase.fvv = NULL;
	fine_phase.n = 2*num;
	fine_phase.p = 7; // Using 3 parameters
	fine_phase.params = &data;

	gsl_multifit_nlinear_init(&init_guess.vector, &fine_phase, w);
	
	f = gsl_multifit_nlinear_residual(w);
	gsl_blas_ddot(f,f,&chisq0);

	status = gsl_multifit_nlinear_driver(MAX_EVALS, RELATIVE_TOLERANCE, RELATIVE_TOLERANCE, 0, NULL, NULL, &info, w);

	J = gsl_multifit_nlinear_jac(w);
	gsl_multifit_nlinear_covar(J, 0.0, covar);

	gsl_blas_ddot(f, f, &chisq);

	*a = gsl_complex_polar(gsl_vector_get(w->x, 0)*a_mag, gsl_vector_get(w->x, 1));
	*Qr += (*Qr)*gsl_vector_get(w->x, 2);
	*Qc += (*Qc)*gsl_vector_get(w->x, 3);
	*fr += (*fr/(*Qr))*gsl_vector_get(w->x, 4);
	*phi0 = gsl_vector_get(w->x, 5);
	*delay = gsl_vector_get(w->x, 6)/(data.fr_init);

	if(DEBUG_)		
	{
		fprintf(stderr, "----------------------\n");
		fprintf(stderr, "Fine Phase Fitting\n");
		fprintf(stderr, "number of iterations: %zu\n",gsl_multifit_nlinear_niter(w));
		fprintf(stderr, "function evaluations: %zu\n", fine_phase.nevalf);
		fprintf(stderr, "Jacobian evaluations: %zu\n", fine_phase.nevaldf);
		fprintf(stderr, "reason for stopping: %s\n", (info == 1) ? "small step size" : "small gradient");
		fprintf(stderr, "initial |f(x)| = %f\n", sqrt(chisq0));
		fprintf(stderr, "final |f(x)| = %f\n", sqrt(chisq));

		{
			double dof = num - 7;
			double c = sqrt(chisq / dof);
			fprintf(stderr, "chisq/dof = %g\n", chisq / dof);
			fprintf(stderr, "|1+eps| = %e +/- %e\n", gsl_complex_abs(*a), c*sqrt(gsl_matrix_get(covar,0,0))*a_mag);
			fprintf(stderr, "arg(1+eps) = %e +/- %lf\n", gsl_vector_get(w->x, 1), c*sqrt(gsl_matrix_get(covar,1,1)));
			fprintf(stderr, "Qr = %lf +/- %e\n", *Qr, c*sqrt(gsl_matrix_get(covar,2,2))*(data.Qr_init));
			fprintf(stderr, "Qc = %e +/- %e\n", *Qc, c*sqrt(gsl_matrix_get(covar,3,3))*(data.Qc_init));
			fprintf(stderr, "fr = %lf +/- %e\n", *fr, c*sqrt(gsl_matrix_get(covar,4,4))*(data.fr_scale));
			fprintf(stderr, "phi0 = %lf +/- %e\n", gsl_vector_get(w->x, 5), c*sqrt(gsl_matrix_get(covar,5,5)));
			fprintf(stderr, "delay = %e +/- %e\n", gsl_vector_get(w->x, 6)/(data.fr_init), c*sqrt(gsl_matrix_get(covar,6,6))/(data.fr_init));
		}
		fprintf(stderr, "----------------------\n");
	}

	gsl_multifit_nlinear_free(w);
	gsl_matrix_free (covar);
}

void fitCircleToData(int num, gsl_complex* z, double* radius, gsl_complex* z_center) 
{
	_getInitValues_Kasa(num, z, radius, z_center);

	circle_t data = {num, z};

	double init_guesses[3] = { GSL_REAL(*z_center), GSL_IMAG(*z_center), *radius };  /* `*`some` `initial` `guess`*` */

	double chisq, chisq0;
	int status, info;

	gsl_vector *f;
	gsl_matrix *J;
	gsl_matrix *covar = gsl_matrix_alloc (3, 3);

	const gsl_multifit_nlinear_type* T = gsl_multifit_nlinear_trust;
	gsl_multifit_nlinear_parameters params = gsl_multifit_nlinear_default_parameters();
	gsl_multifit_nlinear_workspace* w = gsl_multifit_nlinear_alloc(T, &params, num, 3);

	gsl_vector_view init_guess = gsl_vector_view_array(init_guesses,3);

	gsl_multifit_nlinear_fdf circle_fitter;
	circle_fitter.f = _circleFunc_f;
	circle_fitter.df = _circleFunc_df;
	circle_fitter.fvv = NULL;
	circle_fitter.n = num;
	circle_fitter.p = 3; // Using 3 parameters
	circle_fitter.params = &data;

	gsl_multifit_nlinear_init(&init_guess.vector, &circle_fitter, w);
	
	f = gsl_multifit_nlinear_residual(w);
	gsl_blas_ddot(f,f,&chisq0);

	status = gsl_multifit_nlinear_driver(MAX_EVALS, 1e-12, 1e-12, 0, NULL, NULL, &info, w);

	J = gsl_multifit_nlinear_jac(w);
	gsl_multifit_nlinear_covar(J, 0.0, covar);

	gsl_blas_ddot(f, f, &chisq);

	*z_center = gsl_complex_rect(gsl_vector_get(w->x, 0), gsl_vector_get(w->x, 1));
	*radius = gsl_vector_get(w->x, 2);

	if(DEBUG_)		
	{
		fprintf(stderr, "----------------------\n");
		fprintf(stderr, "Circle Fitting\n");
		fprintf(stderr, "number of iterations: %zu\n",gsl_multifit_nlinear_niter(w));
		fprintf(stderr, "function evaluations: %zu\n", circle_fitter.nevalf);
		fprintf(stderr, "Jacobian evaluations: %zu\n", circle_fitter.nevaldf);
		fprintf(stderr, "reason for stopping: %s\n", (info == 1) ? "small step size" : "small gradient");
		fprintf(stderr, "initial |f(x)| = %f\n", sqrt(chisq0));
		fprintf(stderr, "final |f(x)| = %f\n", sqrt(chisq));

		{
			double dof = num - 3;
			double c = GSL_MAX_DBL(1, sqrt(chisq / dof));
			fprintf(stderr, "chisq/dof = %g\n", chisq / dof);
			fprintf(stderr, "z.x = %e +/- %e\n", gsl_vector_get(w->x, 0), c*sqrt(gsl_matrix_get(covar,0,0)));
			fprintf(stderr, "z.y = %e +/- %e\n", gsl_vector_get(w->x, 1), c*sqrt(gsl_matrix_get(covar,1,1)));
			fprintf(stderr, "r = %e +/- %e\n", gsl_vector_get(w->x, 2), c*sqrt(gsl_matrix_get(covar,2,2)));
		}
		fprintf(stderr, "----------------------\n");
	}

	gsl_multifit_nlinear_free(w);
	gsl_matrix_free (covar);
}