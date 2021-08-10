#include "phase_fit.h"
#define DEBUG 0

double phaseFunc_ErrorSquared(unsigned n_params, const double* params, double* grad, void* data) // Least Squares function to be minimized
{
	data_t* d = (data_t*) data;
	int num = d->n;
	double* freq = d->freq;
	double* phase = d->phase;

	double res = 0;

	//#pragma omp parallel for schedule(static)
	for (int i = 0; i < num; ++i)
	{
		double phaseDiff = fabs(phase[i] + params[0] - 2.0*atan(2.0*params[1]*(1.0-(freq[i]/params[2]))));
		
		while(phaseDiff > 2.0*M_PI)
		{
			phaseDiff -= 2.0*M_PI;
			//printf("Did1\n");
		}

		while(phaseDiff < -2.0*M_PI)
		{
			phaseDiff += 2.0*M_PI;
			//printf("Did2\n");
		}

		res += gsl_pow_2(fmin(phaseDiff,2*M_PI-phaseDiff));
	}

	// The cases are taken here to take care of the usage of principal domain.
	if(grad)
	{
		grad[0] = grad[1] = grad[2] = 0.0;
		for (int i = 0; i < num; ++i)
		{
			double df = params[2] - freq[i];
			double phaseDiff = 2.0*atan(2.0*params[1]*(1.0-(freq[i]/params[2]))) - (phase[i] + params[0]);

			while(phaseDiff > 2.0*M_PI)
			{
				phaseDiff -= 2.0*M_PI;
				//printf("Did1\n");
			}

			while(phaseDiff < -2.0*M_PI)
			{
				phaseDiff += 2.0*M_PI;
				//printf("Did2\n");
			}

			int caseToConsider;

			if(phaseDiff > 0)
			{
				caseToConsider = ((phaseDiff < (2*M_PI - phaseDiff)) ? 1 : 0); // Case 0 if false, 1 if true
			}
			else
			{
				caseToConsider = (((-1.0*phaseDiff) < (2*M_PI + phaseDiff)) ? 3 : 2); // Case 2 if false, 3 if true
			}

			if(caseToConsider == 0)
			{
				phaseDiff = 2*M_PI - phaseDiff;

				grad[0] -= phaseDiff;

				grad[1] -= phaseDiff*params[2]*df/(gsl_pow_2(params[2]) + gsl_pow_2(2*params[1]*df));

				grad[2] -= phaseDiff*freq[i]*params[1]/(gsl_pow_2(params[2]) + gsl_pow_2(2*params[1]*df));
			}
			else if (caseToConsider == 2)
			{
				phaseDiff = 2*M_PI + phaseDiff;

				grad[0] += phaseDiff;

				grad[1] += phaseDiff*params[2]*df/(gsl_pow_2(params[2]) + gsl_pow_2(2*params[1]*df));

				grad[2] += phaseDiff*freq[i]*params[1]/(gsl_pow_2(params[2]) + gsl_pow_2(2*params[1]*df));
			}
			else // Unchanged
			{
				grad[0] += phaseDiff;

				grad[1] += phaseDiff*params[2]*df/(gsl_pow_2(params[2]) + gsl_pow_2(2*params[1]*df));

				grad[2] += phaseDiff*freq[i]*params[1]/(gsl_pow_2(params[2]) + gsl_pow_2(2*params[1]*df));
			}
		}

		grad[0] *= -2.0;
		grad[1] *= 8.0;
		grad[2] *= 8.0;
	}

	return res;
}

void performFit(int num, double* freq, double* phase, double* theta, double* Qr, double* fr)
{
	// These are in the form {theta_0, Qr, fr}

	double f_lim_min = 0.0, f_lim_max = HUGE_VAL;
	double Qr_lim_min = 0.0, Qr_lim_max = HUGE_VAL;

	if(SET_REASONABLE_LIMITS && *fr != 0.0)
	{
		f_lim_min = *fr *(1- (1/ *Qr)*SCALE_REASONABLE_LIMITS);
		f_lim_max = *fr * (1+  (1/ *Qr)*SCALE_REASONABLE_LIMITS);

		if(DEBUG)
			printf("f_lim_min : %e ; f_lim_max : %e\n",f_lim_min,f_lim_max);
	}

	if(SET_REASONABLE_LIMITS && *Qr != 0.0)
	{
		Qr_lim_min = (*Qr)/SCALE_REASONABLE_LIMITS;
		Qr_lim_max = (*Qr)*SCALE_REASONABLE_LIMITS;

		if(DEBUG)
			printf("Qr_lim_min : %lf ; Qr_lim_max : %lf\n",Qr_lim_min,Qr_lim_max);
	}
	double lowerBound[3] = {-1.0*M_PI,Qr_lim_min,f_lim_min};
	double upperBound[3] = {M_PI,Qr_lim_max,f_lim_max};

	data_t data = {num, freq, phase};

	nlopt_opt opt;
	opt = nlopt_create(NLOPT_LD_MMA, 3); /* algorithm and dimensionality */
	nlopt_set_lower_bounds(opt, lowerBound);
	nlopt_set_upper_bounds(opt,upperBound);
	nlopt_set_min_objective(opt, phaseFunc_ErrorSquared, &data);

	nlopt_set_xtol_rel(opt, RELATIVE_TOLERANCE); // Setting Relative tolerance
	nlopt_set_xtol_abs1(opt,ABS_TOLERANCE_RESIDUE); // Setting Absolute tolerance
	nlopt_set_maxeval(opt,MAX_EVALS); // Setting the maximum number of function evaluations

	if(*fr == 0.0)
		*fr = F_INIT;

	if(*Qr == 0.0)
		*Qr = Q_INIT;

	double params[3] = { *theta, *Qr, *fr };  /* `*`some` `initial` `guess`*` */
	double minf; /* `*`the` `minimum` `objective` `value,` `upon` `return`*` */

	int success = nlopt_optimize(opt, params, &minf);
	if (success < 0) {
	    fprintf(stderr,"nlopt failed with error code %d!\n", success);
	}
	else {
	   	if(DEBUG)
	   		printf("found minimum at f(%g,%g,%g) = %0.10g with Code %d\n", params[0], params[1], params[2], minf, success);
	    if(DEBUG)
	    	printf("Number of evaluations : %d\n",nlopt_get_numevals(opt));
	}

	*theta = params[0];
	*Qr = params[1];
	*fr = params[2];


}
