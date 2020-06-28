#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>
#include <fftw3.h>

#include "sarif.h"


/**
 * Creates the range chirp vector, in frequency domain and conjugated
 */
int sarif_make_range_chirp(ERS_Raw_Parser_Params *params, double complex **out_range_chirp)
{
	int tmp_chirp_len;
	int intervals;
	int fft_size;
	double complex *ra_chirp_temp;
	double complex* range_chirp;
	double complex* fft_range_chirp;
	fftw_plan p;

	if(!params) {
		return -1;
	}

	intervals = (int)floor(params->tau*params->fs); //size_chirp_r
	printf("%g %g %g %g %d\n", params->kr, params->tau, params->fs, params->tau*params->fs, intervals);

	//Make range chirp
	range_chirp = calloc(1, sizeof(double complex) * params->ra_fft_len);
	for(int i = 0; i < intervals; i++) {
		double phase = params->ra_ph_off + M_PI*params->kr*pow(-params->tau/2+i/params->fs,2);
		range_chirp[i] = cexp(I*phase);
	}

//	printf("CHIRP\n");
//	for(int i = 0; i < tmp_chirp_len; i++) {
//		printf("[%d](%g+j%g), ", i, creal(ra_chirp_temp[i]), cimag(ra_chirp_temp[i]));
//	}
//	printf("\n");

	fft_range_chirp = calloc(1, sizeof(double complex) * params->ra_fft_len);

	p = fftw_plan_dft_1d(params->ra_fft_len, (double (*)[2])range_chirp, (double (*)[2])fft_range_chirp, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p);

//	printf("CHIRP FFT\n");
	for(int i = 0; i < params->ra_fft_len; i++) {
//		printf("[%d](%g+j%g), ", i, creal(fft_range_chirp[i]), cimag(fft_range_chirp[i]));
		fft_range_chirp[i] = conj(fft_range_chirp[i]);
	}
//	printf("\n");

	*out_range_chirp = fft_range_chirp;

	fftw_destroy_plan(p);
	free(range_chirp);

	return 0;
}

int sarif_remove_mean(ERS_Raw_Parser_Data_Patch *in)
{
	double avg_re = 0;
	double avg_im = 0;

	for(int i = 0; i < in->n_az * in->n_ra; i++) {
		avg_re += creal(in->data[i]);
		avg_im += cimag(in->data[i]);
	}
	avg_re /= (in->n_az * in->n_ra);
	avg_im /= (in->n_az * in->n_ra);
	for(int i = 0; i < in->n_az * in->n_ra; i++) {
		in->data[i] = creal(in->data[i])-avg_re + I*(cimag(in->data[i])-avg_im);
	}
	return 0;
}

int sarif_range_compression(double complex *out, ERS_Raw_Parser_Data_Patch *in, ERS_Raw_Parser_Params *params, double complex *f_conj_range_chirp)
{
	double complex *aux_line;
	double complex *f_aux_line;
	double complex *f_corr_aux_line;
	double complex *corr_aux_line;
	int range;
	fftw_plan p_forward, p_backward;

	range = in->n_ra;

	aux_line =  calloc(1, sizeof(double complex) * params->ra_fft_len);
	f_aux_line =  calloc(1, sizeof(double complex) * params->ra_fft_len);
	f_corr_aux_line =  calloc(1, sizeof(double complex) * params->ra_fft_len);
	corr_aux_line =  calloc(1, sizeof(double complex) * params->ra_fft_len);

	p_forward = fftw_plan_dft_1d(params->ra_fft_len, aux_line, f_aux_line, FFTW_FORWARD, FFTW_ESTIMATE);
	p_backward = fftw_plan_dft_1d(params->ra_fft_len, f_corr_aux_line, corr_aux_line, FFTW_BACKWARD, FFTW_ESTIMATE);

	//Range compression
	for(int i = 0; i < in->n_az; i++) {
		//Copy range values
		memcpy(aux_line, &in->data[i*range], sizeof(double complex) * range);
//		printf("IN\n"); //->ok
//		for(int j = 0; j < range; j++) {
//			printf("[%d](%g+j%g), ", j, creal(aux_line[j]), cimag(aux_line[j]));
//		}

		fftw_execute(p_forward);
		//printf("FFT\n"); //-> ok
		//for(int j = 0; j < params->ra_fft_len; j++) {
		//	printf("[%d](%g+j%g), ", j, creal(f_aux_line[j]), cimag(f_aux_line[j]));
		//}
		//printf("\n");

		//fft ok, chirp ok
		//printf("CORR\n");
		for(int j = 0; j < params->ra_fft_len; j++) {
			f_corr_aux_line[j] = f_aux_line[j] * f_conj_range_chirp[j];
		//	printf("[%d](%g+j%g), ", j, creal(f_corr_aux_line[j]), cimag(f_corr_aux_line[j]));
		}
		//printf("\n");

		fftw_execute(p_backward);
		//printf("IF\n");
		for(int j = 0; j < params->n_valid_samples; j++) {
			//printf("[%d](%g+j%g), ", j, creal(corr_aux_line[j])/params->ra_fft_len, cimag(corr_aux_line[j])/params->ra_fft_len);
			out[j + (in->az_pos + i) * params->n_valid_samples] = corr_aux_line[j]/params->ra_fft_len;
		}
		//printf("\n");
		//memcpy(out + (in->az_pos + i) * params->n_valid_samples, corr_aux_line, sizeof(double complex) * params->n_valid_samples);

		//exit(0);
		//double complex *aux = &out[i*range];
		//memcpy(aux+range/2, corr_aux_line, sizeof(double complex) * range/2);
		//memcpy(aux, corr_aux_line+range/2, sizeof(double complex) * range/2);
		//if( i == 0) {
		//	printf("SHIFTED\n");
		//	for(int j = 0; j < range; j++) {
		//		printf("[%d](%g+j%g), ", j, creal(aux[j])/params->ra_fft_len, cimag(aux[j])/params->ra_fft_len);
		//	}
		//	printf("\n");
		//}
	}
	//everything's ok!
	free(aux_line);
	free(f_aux_line);
	free(f_corr_aux_line);
	free(corr_aux_line);
	fftw_destroy_plan(p_forward);
	fftw_destroy_plan(p_backward);

	return 0;
}
