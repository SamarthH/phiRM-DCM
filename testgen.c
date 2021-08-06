#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define STDEV_ANG 0.5
#define STDEV_MAG 0.00001

int main(int argc, char const *argv[]) // The inputs are cable delay(tau in ns), Qr, Qc, fr (in GHz), phi0 (in radian), freq_start (in GHz), freq_end (in GHz), step (in GHz)
{
	if(argc != 9)	
	{
		printf("Too many or too few params\n");
		return -1;
	}

	double tau=atof(argv[1]),Qr=atof(argv[2]),Qc=atof(argv[3]),fr=atof(argv[4]),phi0=atof(argv[5]),freq_start=atof(argv[6]),freq_end=atof(argv[7]),step=atof(argv[8]);

	for (double f = freq_start; f <= freq_end; f+=step)
	{
		double a = Qr/Qc;
		double b = 2.0*Qr*(f-fr)/fr;

		gsl_complex rot = gsl_complex_polar(1,-2.0*M_PI*f*tau);

		gsl_complex t = gsl_complex_polar(a,phi0);
		t = gsl_complex_div(t,gsl_complex_rect(1,b));
		t = gsl_complex_mul(rot,t);
		t = gsl_complex_sub(rot,t);
		if(1)
		{
			gsl_rng* rangen = gsl_rng_alloc(gsl_rng_ranlxd2);
			gsl_rng_set(rangen,0);
			gsl_complex noise = gsl_complex_polar(gsl_ran_gaussian(rangen,STDEV_MAG),gsl_ran_gaussian(rangen,STDEV_ANG));
			//noise = gsl_complex_add_real(noise,1.0);
			t = gsl_complex_add(t,noise);
		}

		printf("%0.11lf %0.11lf %0.11lf\n", f,  GSL_REAL(t), GSL_IMAG(t));
	}
	return 0;
}