#include <stdio.h> 
#include <stdlib.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "phase_fit.h"

#define N_MAX 300000
#define SELF_START 1
#define STDEV_OVERRIDE 100
#define STDEV_FRACTION_ALLOWED 0.01
#define DEBUG 0
#define WITH_PY 1 // Compile with WITH_PY=1 if you are going to use a python script to read the output

gsl_complex correctForEpsilon(int num, gsl_complex* z) //Rotates the initial thing and corrects for 1 + epsilon
{
	gsl_complex onePlusEpsilon;
	GSL_SET_COMPLEX(&onePlusEpsilon,z[0].dat[0],z[0].dat[1]);
	for(int i = 0; i < num; ++i)
	{
		z[i] = gsl_complex_div(z[i],onePlusEpsilon);
	}

	return onePlusEpsilon;
}

void rotate(gsl_complex* a, double phi) // Rotates 'a' by phi anticlockwise, inplace.
{
	gsl_complex z;
	GSL_SET_COMPLEX(&z,a->dat[0], a->dat[1]);
	GSL_SET_COMPLEX(a, z.dat[0]*cos(phi) - z.dat[1]*sin(phi), z.dat[1]*cos(phi) + z.dat[0]*sin(phi));
}

void removeCableDelay(int num, double* freq, gsl_complex* z, double delay)
{
	for (int i = 0; i < num; ++i)
	{
		rotate(z+i,2*M_PI*freq[i]*delay);
	}
}

void rotateAndTranslateToOrigin(int num, gsl_complex* z, gsl_complex* z_center)
{
	
	for (int i = 0; i < num; ++i)
	{
		//First translate
		z[i] = gsl_complex_sub(*z_center,z[i]);
		//Then rotate
		rotate(z+i, -1.0* gsl_complex_arg(*z_center));
	}
}

double getStdev2(gsl_complex* z)
{
	int forStdev = 0;
	while(fabs(gsl_complex_abs(z[forStdev])-gsl_complex_abs(z[0])) < STDEV_FRACTION_ALLOWED*gsl_complex_abs(z[0]))
		forStdev++;

	if(STDEV_FRACTION_ALLOWED == 0)
		forStdev = STDEV_OVERRIDE;

	double ret = 0.0;
	#pragma omp parallel for schedule(static) reduction(+:ret)
	for (int i = 0; i < forStdev; ++i)
	{
		ret += gsl_complex_abs2(gsl_complex_sub(z[i],z[i+1]));
	}

	if(DEBUG)
		printf("Num = %d, Res = %e\n",forStdev, ret);

	return ret/(2.*forStdev);
}

double getChi2(int num, gsl_complex* z, double* freq, double stdev2, double delay, double Qr, double Qc, double fr, double phi0, gsl_complex onePlusEpsilon)
{
	if(DEBUG)
		printf("stdev2 = %e\n",stdev2);

	double ret = 0.0;

	//#pragma omp parallel for schedule(static) reduction(+:ret)
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
		
		ret += (gsl_pow_2(GSL_REAL(t) - GSL_REAL(z[i])) + gsl_pow_2(GSL_IMAG(t) - GSL_IMAG(z[i])));
	}

	if(DEBUG)
		printf("ret = %e\n", ret);

	ret /= (2.*num-7);
	ret /= stdev2;

	return ret;
}

void generateExpectedValues(int num, gsl_complex* z, double* freq, double delay, double Qr, double Qc, double fr, double phi0, gsl_complex onePlusEpsilon)
{
	FILE* fp;
	fp = fopen("exp.txt","w+");

	for (int i=0; i<num; i++)
	{
		double a = Qr/Qc;
		double b = 2.0*Qr*(freq[i]-fr)/fr;

		gsl_complex rot = gsl_complex_polar(1,-2.0*M_PI*freq[i]*delay);

		gsl_complex t = gsl_complex_polar(a,phi0);
		t = gsl_complex_div(t,gsl_complex_rect(1,b));
		t = gsl_complex_mul(rot,t);
		t = gsl_complex_sub(rot,t);
		t = gsl_complex_mul(onePlusEpsilon,t);

		fprintf(fp,"%0.16e %0.16e %0.16e %0.16e %0.16e %0.16e %0.16e\n", freq[i], GSL_REAL(z[i]), GSL_IMAG(z[i]), gsl_complex_abs(z[i]), GSL_REAL(t), GSL_IMAG(t), gsl_complex_abs(t));
	}
	fclose(fp);
}

void getInitialGuesses(gsl_complex* z, double* freq, double* fr, double* Qr, int num)
{
	//Getting trial values for fr
	double s21_fr_trial = HUGE_VAL; // Some large number
	int n_fr = 0;
	for (int i = 0; i < num; ++i)
	{
		if(gsl_complex_abs(z[i]) < s21_fr_trial)
		{
			s21_fr_trial = gsl_complex_abs(z[i]);
			*fr = freq[i];
			n_fr = i;
		}
	}

	if(!SELF_START)
		*fr = 0.0;

	//Getting trial values for Qr
	double half_max = 0.7*(gsl_complex_abs(z[0]) - s21_fr_trial) + s21_fr_trial;
	*Qr = 0.0;
	int flag = 1;
	for (int i = n_fr; i >= 0 && flag; i--)
	{
		if(gsl_complex_abs(z[i+1]) < half_max && gsl_complex_abs(z[i]) >= half_max)
		{
			*Qr -= freq[i];
			flag = 0;
		}
	}
	if(flag == 0)
		flag = 1;
	for (int i = n_fr; i < num-1 && flag; i++)
	{
		if(gsl_complex_abs(z[i+1]) > half_max && gsl_complex_abs(z[i]) <= half_max)
		{
			*Qr += freq[i];
			flag = 0;
		}
	}

	if(flag == 1 || !SELF_START)
		*Qr = 0;
	else
		*Qr = 2.0*(*fr)/ *Qr;
	//Qr += 0.01;
}

void detrendInput(gsl_complex* z, double* freq, int num, int detrend_mode, int detrend_points_init, int detrend_points_fin)
{
	if (detrend_mode == 0)
	{
		return;
	}
	else if(detrend_mode != 1 && detrend_mode != 2)
	{
		fprintf(stderr, "Unsupported Detrend Mode\n");
		exit(1);
	}

	double freq_off = freq[0];

	for (int i = 0; i < num; ++i)
	{
		freq[i] -= freq_off;
	}

	double* mag = (double*) malloc(num*sizeof(double));
	double* phas = (double*) malloc(num*sizeof(double));
	for (int i = 0; i < num; ++i)
	{
		mag[i] = gsl_complex_abs(z[i]);
		phas[i] = gsl_complex_arg(z[i]);
	}

	if(detrend_mode == 1)
	{
		double slope_mag = 0.0, slope_phas = 0.0;
		
		// Extracting Slope
		if(detrend_points_init == 1 && detrend_points_fin == 1)
		{
			slope_mag = (mag[num-1] - mag[0])/(freq[num-1] - freq[0]);
			slope_phas = (phas[num-1] - phas[0])/(freq[num-1] - freq[0]);
		}
		else
		{
			double sum_x=0, sum_y_mag=0, sum_y_phas=0, sum_xx=0, sum_xy_mag=0, sum_xy_phas=0; // These will be used for the closed form solution of the least-square fit
			for (int i = 0; i < detrend_points_init; ++i)
			{
				sum_x += freq[i];
				sum_xx += freq[i]*freq[i];
				sum_xy_mag += freq[i]*mag[i];
				sum_y_mag += mag[i];
				sum_xy_phas += freq[i]*phas[i];
				sum_y_phas += phas[i];
			}
			for (int i = num-detrend_points_fin; i < num; ++i)
			{
				sum_x += freq[i];
				sum_xx += freq[i]*freq[i];
				sum_xy_mag += freq[i]*mag[i];
				sum_y_mag += mag[i];
				sum_xy_phas += freq[i]*phas[i];
				sum_y_phas += phas[i];
			}
			int n = detrend_points_init + detrend_points_fin;
			slope_mag = (n*sum_xy_mag - sum_x*sum_y_mag)/(n*sum_xx - sum_x*sum_x);
			slope_phas = (n*sum_xy_phas - sum_x*sum_y_phas)/(n*sum_xx - sum_x*sum_x);
		}
		// Obtained Slope

		// Removing Slope
		for (int i = 0; i < num; ++i)
		{
			mag[i] -= slope_mag*freq[i];
			phas[i] -= slope_phas*freq[i];
		}
	}
	else if(detrend_mode == 2)
	{

		// Fitting to ax^2 + bx + c = 0 for magnitude, and just b*x + c for phase
		double a_mag,b_mag,b_phas;

		//Obtaining components
		double sum_x = 0, sum_x2 = 0, sum_x3 = 0, sum_x4 = 0, sum_y_mag = 0, sum_y_phas = 0, sum_xy_mag = 0, sum_xy_phas = 0, sum_x2y_mag = 0;

		for (int i = 0; i < detrend_points_init; ++i)
		{
			sum_x += freq[i];
			sum_x2 += freq[i]*freq[i];
			sum_x3 += freq[i]*freq[i]*freq[i];
			sum_x4 += freq[i]*freq[i]*freq[i]*freq[i];
			sum_y_mag += mag[i];
			sum_y_phas += phas[i];
			sum_xy_mag += freq[i]*mag[i];
			sum_xy_phas += freq[i]*phas[i];
			sum_x2y_mag += freq[i]*freq[i]*mag[i];
			// sum_x2y_phas += freq[i]*freq[i];
		}
		for (int i = num-detrend_points_fin; i < num; ++i)
		{
			sum_x += freq[i];
			sum_x2 += freq[i]*freq[i];
			sum_x3 += freq[i]*freq[i]*freq[i];
			sum_x4 += freq[i]*freq[i]*freq[i]*freq[i];
			sum_y_mag += mag[i];
			sum_y_phas += phas[i];
			sum_xy_mag += freq[i]*mag[i];
			sum_xy_phas += freq[i]*phas[i];
			sum_x2y_mag += freq[i]*freq[i]*mag[i];
			// sum_x2y_phas += freq[i]*freq[i];
		}


		int n = detrend_points_init + detrend_points_fin;

		double Sxx=0, Sxx2=0, Sx2x2=0, Sx2y_mag=0, Sxy_mag=0;
		Sxx = sum_x2 - sum_x*sum_x/n;
		Sxx2 = sum_x3 - sum_x*sum_x2/n;
		Sx2x2 = sum_x4 - sum_x2*sum_x2/n;
		Sxy_mag = sum_xy_mag - sum_x*sum_y_mag/n;
		//Sxy_phas = sum_xy_phas - sum_x*sum_y_phas/n;
		Sx2y_mag = sum_x2y_mag - sum_x2*sum_y_mag/n;
		//Sx2y_phas = sum_x2y_phas - sum_x2*sum_y_phas/n;

		a_mag = (Sx2y_mag*Sxx - Sxy_mag*Sxx2)/(Sxx*Sx2x2 - Sxx2*Sxx2);
		b_mag = (Sxy_mag*Sx2x2 - Sx2y_mag*Sxx2)/(Sxx*Sx2x2 - Sxx2*Sxx2);

		b_phas = (n*sum_xy_phas - sum_x*sum_y_phas)/(n*sum_x2 - sum_x*sum_x);
		if(DEBUG)
			printf("a_mag = %lf; b_mag = %lf\n",a_mag,b_mag);
		//Obtained Components

		//Removing the components
		for (int i = 0; i < num; ++i)
		{
			mag[i] -= (a_mag*freq[i]*freq[i] + b_mag*freq[i]);
			phas[i] -=  b_phas*freq[i];
		}
	}

	double mag_max = 0.0;
	for (int i = 0; i < num; ++i)
	{
		if(mag_max< mag[i])
		{
			mag_max = mag[i];
		}
	}
	for (int i = 0; i < num; ++i)
	{
		//mag[i] /= mag_max;
		z[i] = gsl_complex_polar(mag[i],phas[i]);
	}

	for (int i = 0; i < num; ++i)
	{
		freq[i] += freq_off;
	}

	free(mag);
	free(phas);
}

int main(int argc, char const *argv[]) //Inputs are of the form <filename> <delay>
{
	//Taking inputs
	
	if(argc < 4)
	{
		printf("Not enough arguments\n");
		return -1;
	}

	FILE* input;
	input = fopen(argv[1],"r");

	int num = 0, size= 1<<6;
	double* freq = (double*) malloc(size*sizeof(double));
	gsl_complex* z = (gsl_complex*) malloc(size*sizeof(gsl_complex));

	while(fscanf(input,"%lf%lf%lf",&freq[num],&(z[num].dat[0]),&(z[num].dat[1]))==3)
	{
		//freq[num] /= 1e9;
		num++;
		if(num >= size)
		{
			size = size*2;
			freq = (double*) realloc(freq,size*sizeof(double));
			z = (gsl_complex*) realloc(z,size*sizeof(gsl_complex));
		}
	}
	if(DEBUG)
		printf("num = %d\n",num);

	freq = (double*) realloc(freq, num*sizeof(double));
	z = (gsl_complex*) realloc(z, num*sizeof(gsl_complex));

	double delay = atof(argv[2]);

	int detrend_mode = atoi(argv[3]); // 0 means no detrend, 1 means linear detrend, 2 means quadratic detrend
	int detrend_points_init = 1, detrend_points_fin = 1; // Default values
	if(argc !=6 && detrend_mode == 2)
	{
		printf("Not enough or Too many arguments for Quadratic Detrend\n");
		return -1;
	}
	else if(argc == 6)
	{
		detrend_points_init = atoi(argv[4]);
		detrend_points_fin = atoi(argv[5]);
	}

	// Done taking inputs

	// First removing Cable delay
	if(delay != 0)
		removeCableDelay(num,freq,z,delay);

	//Detrending
	detrendInput(z,freq,num,detrend_mode,detrend_points_init,detrend_points_fin);

	// Copying into array for comparision later
	gsl_complex* z_true = (gsl_complex*) malloc(num*sizeof(gsl_complex));

	for (int i = 0; i < num; ++i)
	{
		z_true[i] = z[i];
	}

	//Removing 1+epsilon
	gsl_complex onePlusEpsilon = correctForEpsilon(num,z);

	double theta,Qr,fr,Qr_init,fr_init;

	getInitialGuesses(z, freq, &fr, &Qr, num);
	fr_init = fr;
	Qr_init = Qr;

	// Fitting the circle to get some parameters
	double radius;
	gsl_complex z_center;
	fitCircleToData(num,z,&radius,&z_center);
	if(DEBUG)
	{
		//printf("RMSE of Circle Fit = %e\n",rmse);
		printf("z_center.x = %lf ; z_center.y = %lf ; radius = %lf\n",z_center.dat[0],z_center.dat[1],radius);
	}

	//Rotating and Translating to the origin
	rotateAndTranslateToOrigin(num,z,&z_center);

	//Generate the phase array for fitting
	double* phase = (double*) malloc(num*sizeof(double));

	for (int i = 0; i < num; ++i)
	{
		double arg = gsl_complex_arg(z[i]);
		phase[i] = arg;
	}

	// Fit the phase to get the other parameters
	double phi_expected = -1.0*asin((z_center.dat[1])/radius);
	if(DEBUG)
		printf("phi_exp = %lf\n", phi_expected);

	if(SELF_START && !isnan(phi_expected))
	{
		theta = phi_expected + gsl_complex_arg(z_center);
	}
	else{
		theta = 0.0;
		printf("NAN\n");
	}

	if(DEBUG)
		printf("theta_init = %lf\n",theta);

	performFit(num,freq,phase,&theta,&Qr,&fr);

	// Get other parameters

	//double Qc_PhiRM = (gsl_complex_abs(z_center)+radius)*Qr/(2.0*radius);
	
	double Qc_PhiRM = Qr/(2.0*radius);
	double phi0 = gsl_complex_arg(z_center) - theta;

	double chi2_phiRM;

	if(DEBUG)
	{
		chi2_phiRM = getChi2(num, z_true, freq, getStdev2(z_true), delay, Qr, Qc_PhiRM, fr, phi0, onePlusEpsilon);
		printf("%lf\t%lf\t%lf\t%lf\n", fr, Qr, Qc_PhiRM,chi2_phiRM);
	}

	refineFit(num, freq, z_true, &onePlusEpsilon, &Qr, &Qc_PhiRM, &fr, &phi0, &delay);

	// For DCM, introducing a correction factor
	double Qc_DCM = Qc_PhiRM/cos(phi0);

	//Finding Chi^2
	chi2_phiRM = getChi2(num, z_true, freq, getStdev2(z_true), delay, Qr, Qc_PhiRM, fr, phi0, onePlusEpsilon);

	//Printing the obtained parameters
	if(WITH_PY)
	{
		printf("fr\tQr\tQc\tQi\tchi2\n");
		printf("%lf\t%lf\t%lf\t%lf\t%lf\n", fr, Qr, Qc_PhiRM, 1.0/((1.0/Qr) - (1.0/Qc_DCM)), chi2_phiRM);
	}
	else
	{
		printf("fr: %lf (Initial Value: %lf)\n", fr, fr_init);
		printf("Qr: %lf (Initial Value: %lf)\n", Qr, Qr_init);
		printf("1/|Qc^-1|: %lf\n", Qc_PhiRM);
		printf("Qi (By PhiRM): %lf\n", 1.0/((1.0/Qr) - (1.0/Qc_PhiRM)));
		printf("Qi (By DCM): %lf\n", 1.0/((1.0/Qr) - (1.0/Qc_DCM)));
		if(DEBUG){
			printf("phi0: %lf\n", phi0);
			printf("phi_expected: %lf\n", phi_expected);
			printf("delay: %e\n", delay);
			printf("onePlusEpsilon: %e + %ej\n", GSL_REAL(onePlusEpsilon), GSL_IMAG(onePlusEpsilon));
		}
		printf("chi2_phiRM = %lf\n", chi2_phiRM);
	}
	generateExpectedValues(num, z_true, freq, delay, Qr, Qc_PhiRM, fr, phi0, onePlusEpsilon);
	free(freq);
	free(z);
	free(z_true);
	free(phase);

	return 0;
}
