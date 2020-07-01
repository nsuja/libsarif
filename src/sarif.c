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


#define SARIF_MAX_RANGE_FFT_LEN (8192)

struct Sarif_Ctx {
	ERS_Raw_Parser_Params params;
	int has_params;
	double fc;

	double *sq_range; //squinted range

	int has_az_chirp;
	double complex *az_chirp[SARIF_MAX_RANGE_FFT_LEN]; //XXX TODO FIXME Could use a better definition (double complex *) with n_valid_samples * fft_lines, see NOTE1
	double complex *f_az_chirp[SARIF_MAX_RANGE_FFT_LEN];
	int valid_az_lines;

	int rcmc_max_offset; //RCMC matrix maximum pixel offset
	int *rcmc_offset_matrix; //Amount of "pixels" to shift left in range
};

double complex * sarif_make_azimuth_chirp_for_range(double fc, double ka, double tau, double fs, double lines);

Sarif_Ctx * sarif_ctx_alloc(ERS_Raw_Parser_Params *params)
{
	int fd_ldr = -1, fd_raw = -1;
	Sarif_Ctx *ctx;

	if(!params) {
		fprintf(stderr, "%s:: Invalid arguments (%p)", __func__, params);
		return NULL;
	}

	ctx = (Sarif_Ctx *)calloc(1, sizeof(Sarif_Ctx));
	if(!ctx) {
		fprintf(stderr, "%s:: Error calloc(): errno(%d):%s\n", __func__, errno, strerror(errno));
		return NULL;
	}

	memcpy(&ctx->params, params, sizeof(ERS_Raw_Parser_Params));
	ctx->has_params = 1;

	ctx->sq_range = calloc(1, sizeof(double) * params->n_valid_samples);
	for(int i = 0; i < params->n_valid_samples; i++) {
		ctx->sq_range[i] = params->r0 + i * (params->C/(2*params->fs));
	}

	return ctx;
}

int sarif_set_doppler_centroid(Sarif_Ctx *ctx, double fc)
{
	if(!ctx) {
		return -1;
	}
	ctx->fc = fc;
	return 0;
}

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

	//printf("CHIRP FFT\n");
	for(int i = 0; i < params->ra_fft_len; i++) {
		//printf("[%d](%g+j%g), ", i, creal(fft_range_chirp[i]), cimag(fft_range_chirp[i]));
		fft_range_chirp[i] = conj(fft_range_chirp[i]);
	}
	//printf("\n");

	*out_range_chirp = fft_range_chirp;

	fftw_destroy_plan(p);
	free(range_chirp);

	return 0;
}

double complex * sarif_make_azimuth_chirp_for_range(double fc, double ka, double tau, double fs, double lines)
{
	double complex *az_chirp;
	int intervals;

	az_chirp = calloc(1, sizeof(double complex) * lines);

	intervals = (int)floor(tau*fs);
	for(int i = 0; i < intervals; i++) {
		double phase = 2.0*M_PI*fc*(-tau/2+i/fs)+M_PI*ka*pow(-tau/2+i/fs,2);
		az_chirp[i] = cexp(I*phase);
	}

	return az_chirp;
}


int sarif_make_azimuth_chirp(Sarif_Ctx *ctx)
{
	ERS_Raw_Parser_Params *params;
	fftw_plan p;
	double complex *range, *ka;
	double complex *az_chirp_temp;
	double complex *aux_line;
	double complex *f_aux_line;

	if(!ctx || !ctx->has_params) {
		return -1;
	}

	params = &ctx->params;
	if(params->n_valid_samples > SARIF_MAX_RANGE_FFT_LEN) {
		fprintf(stderr, "%s:: More range samples (%d) than supported (%d), please increase SARIF_MAX_RANGE_FFT_LEN", __func__, params->n_valid_samples, SARIF_MAX_RANGE_FFT_LEN);
		return -1;
	}

	double ka_aux = -(2*params->velocity*params->velocity)/params->lambda;
	double aux = pow(params->lambda * ctx->fc / (2*params->velocity),2);

	range = calloc(1, sizeof(double complex) * params->n_valid_samples);
	ka = calloc(1, sizeof(double complex) * params->n_valid_samples);
	double *az_beam_width = calloc(1, sizeof(double complex) * params->n_valid_samples);
	double *az_tau = calloc(1, sizeof(double complex) * params->n_valid_samples);

	for(int i = 0; i < params->n_valid_samples; i++) {
		range[i] = ctx->sq_range[i] / sqrt(1-aux); //Squinted range
		ka[i] = ka_aux/range[i];
		az_beam_width[i] = range[i] * params->az_beam_width * 0.8; //in ground[m], 80% fixed (3db)
		az_tau[i] = az_beam_width[i]/params->velocity;
	}
	int pulses_for_full_aperture = (int)ceil(az_tau[params->n_valid_samples-1]*params->prf);
	int n_valid_lines = params->fft_lines - pulses_for_full_aperture; //fft wraparound

	//XXX TODO FIXME switch to advanced interface to avoid memcpy

	aux_line =  calloc(1, sizeof(double complex) * params->fft_lines);
	f_aux_line =  calloc(1, sizeof(double complex) * params->fft_lines);
	p = fftw_plan_dft_1d(params->fft_lines, aux_line, f_aux_line, FFTW_FORWARD, FFTW_ESTIMATE);
	//Make one azimuth chirp per range bin
	for(int i = 0; i < params->n_valid_samples; i++) {
		ctx->az_chirp[i] = sarif_make_azimuth_chirp_for_range(ctx->fc, ka[i], az_tau[i], params->prf, params->fft_lines);
		memcpy(aux_line, ctx->az_chirp[i], sizeof(double complex) * params->fft_lines);
		fftw_execute(p);
		ctx->f_az_chirp[i] = calloc(1, sizeof(double complex) * params->fft_lines);
		memcpy(ctx->f_az_chirp[i], f_aux_line, sizeof(double complex) * params->fft_lines);
	}
	ctx->has_az_chirp = 1;
	ctx->valid_az_lines = n_valid_lines;

//	//Position chirp in range vector (centered)  --> used for correlation procedure
//	double complex* azimuth_chirp = calloc(1, sizeof(double complex) * azimuth);
//	double complex* fft_azimuth_chirp = calloc(1, sizeof(double complex) * azimuth);
//	memcpy(azimuth_chirp+((azimuth-az_intervals)/2)-1, az_chirp_temp, sizeof(double complex) * az_intervals);
//	//for(int i = 0; i < azimuth; i++) {
//	//>-printf("[%d](%g+j%g), ", i, creal(azimuth_chirp[i]), cimag(azimuth_chirp[i]));
//	//}
//	//printf("\n");
//
//	p = fftw_plan_dft_1d(azimuth, azimuth_chirp, fft_azimuth_chirp, FFTW_FORWARD, FFTW_ESTIMATE);
//	fftw_execute(p);
//	for(int i = 0; i < azimuth; i++) {
//		//printf("[%d](%g+j%g), ", i, creal(fft_azimuth_chirp[i]), cimag(fft_azimuth_chirp[i]));
//		fft_azimuth_chirp[i] = conj(fft_azimuth_chirp[i]);
//	}
//	//printf("\n");
//
//	*out_azimuth_chirp = fft_azimuth_chirp;
//
	free(aux_line);
	free(f_aux_line);
	free(az_tau);
	free(az_beam_width);
	free(ka);
	free(range);
	fftw_destroy_plan(p);
//	fftw_free(azimuth_chirp);

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

int sarif_range_compression(double complex *out, ERS_Raw_Parser_Data_Patch *in, ERS_Raw_Parser_Params *params, double complex *f_conj_range_chirp, int scale_fft)
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

		//int sarif_get_rcmc_data(f_corr_aux_line, ERS_Raw_Parser_Params *params)

		fftw_execute(p_backward);

		//printf("IF\n");
		if(scale_fft) {
			for(int j = 0; j < params->n_valid_samples; j++) {
				//printf("[%d](%g+j%g), ", j, creal(corr_aux_line[j])/params->ra_fft_len, cimag(corr_aux_line[j])/params->ra_fft_len);
				out[j + (in->az_pos + i) * params->n_valid_samples] = corr_aux_line[j]/params->ra_fft_len;
			}
			//printf("\n");
		} else {
			memcpy(out + (in->az_pos + i) * params->n_valid_samples, corr_aux_line, sizeof(double complex) * params->n_valid_samples);
		}

		//double complex *aux = &out[i*params->n_valid_samples];
		//if( i == 0) {
		//	printf("SHIFTED\n");
		//	for(int j = 0; j < params->n_valid_samples; j++) {
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

double sarif_calc_doppler_centroid(double complex *in, Sarif_Doppler_Centroid_Algo algo, ERS_Raw_Parser_Params *params)
{
	double complex *sum_az;
	double *avg_phase_change;
	double avg_phase_change_mean = 0;
	double fc = NAN; //centroid
	if(!in || !params) {
		fprintf(stderr, "%s:: Invalid arguments %p %p", __func__, in, params);
		return fc;
	}
	if(algo >= SARIF_DOPPLER_CENTROID_ALGO_N) {
		fprintf(stderr, "%s:: Invalid algorithm %d", __func__, algo);
		return fc;
	}

	if(algo != SARIF_DOPPLER_CENTROID_ALGO_AVGPHASE) {
		fprintf(stderr, "%s:: Algorithm %d is not supported yet", __func__, algo);
		return fc;
	}

	sum_az = calloc(1, sizeof(double complex) * params->n_valid_samples);
	avg_phase_change = calloc(1, sizeof(double) * params->n_valid_samples);

	for(int i = 1; i < params->fft_lines; i++) {
		for(int j = 0; j < params->n_valid_samples; j++) {
			sum_az[j] += in[j + (i)*params->n_valid_samples] * conj(in[j + (i-1)*params->n_valid_samples]);
		}
	}
	for(int j = 0; j < params->n_valid_samples; j++) {
		avg_phase_change[j] = carg(sum_az[j]);
		avg_phase_change_mean += avg_phase_change[j];
	}
	avg_phase_change_mean /= (double)params->n_valid_samples;
	fc = avg_phase_change_mean*params->prf/2.0/M_PI;

	free(avg_phase_change);
	free(sum_az);

	return fc;
}

int sarif_azimuth_compression(Sarif_Ctx *ctx, double complex *out, double complex *in, int descale_range)
{
	ERS_Raw_Parser_Params *params;
	double complex *f_in;
	double complex *aux_line;
	double complex *f_aux_line;
	double complex *f_corr_aux_line;
	double complex *corr_aux_line;
	fftw_plan p_forward, p_backward;

	if(!ctx || !out) {
		fprintf(stderr, "%s:: Invalid arguments (%p, %p, %d)\n", __func__, ctx, out, descale_range);
		return -1;
	}
	if(!ctx->has_params || !ctx->has_az_chirp) {
		fprintf(stderr, "%s:: Error has_params %d has_az_chirp %d\n", __func__, ctx->has_params, ctx->has_az_chirp);
		return -1;
	}

	params = &ctx->params;

	f_in =  calloc(1, sizeof(double complex) * params->fft_lines * params->n_valid_samples);

	aux_line =  calloc(1, sizeof(double complex) * params->fft_lines);
	f_aux_line =  calloc(1, sizeof(double complex) * params->fft_lines);
	f_corr_aux_line =  calloc(1, sizeof(double complex) * params->fft_lines);
	corr_aux_line =  calloc(1, sizeof(double complex) * params->fft_lines);

	p_forward = fftw_plan_dft_1d(params->fft_lines, aux_line, f_aux_line, FFTW_FORWARD, FFTW_ESTIMATE);
	p_backward = fftw_plan_dft_1d(params->fft_lines, f_corr_aux_line, corr_aux_line, FFTW_BACKWARD, FFTW_ESTIMATE);

	//FFT in azimuth
	//XXX TODO FIXME Change azimuth and range orientation, there are more operations in az dimension and vectors have to be manually transposed or try with the advanced interface

	for(int i = 0; i < params->n_valid_samples; i++) {
		//Copy range values
		//printf("%s:: IN\n", __func__);
		for(int j = 0; j < params->fft_lines; j++) {
			if(descale_range)
				aux_line[j] = in[j*params->n_valid_samples+i]/params->ra_fft_len;
			else
				aux_line[j] = in[j*params->n_valid_samples+i];
			//printf("[%d](%g+j%g), ", j, creal(aux_line[j]), cimag(aux_line[j]));
		}
		//printf("\n");

		fftw_execute(p_forward);
		//printf("AZ FFT\n");
		for(int j = 0; j < params->fft_lines; j++) {
			//printf("[%d](%g+j%g), ", j, creal(f_aux_line[j]), cimag(f_aux_line[j]));
			f_in[j*params->n_valid_samples+i] = f_aux_line[j];
		}
		//printf("\n");
	}

	//Correlation and RCMC
	for(int i = 0; i < params->n_valid_samples; i++) {
		//printf("CORR\n");
		for(int j = 0; j < params->fft_lines; j++) {
			//Apply RCMC -> ok!
			int offset = ctx->rcmc_offset_matrix[i + j*params->n_valid_samples];
			if(offset && params->n_valid_samples - i <= offset)
				f_in[i + j*params->n_valid_samples] = 0;
			else {
				for(int k = 1; k <= ctx->rcmc_max_offset && i + k < params->n_valid_samples; k++) {
					if(ctx->rcmc_offset_matrix[i + k + j*params->n_valid_samples] == k) {
						f_in[i + j*params->n_valid_samples] = f_in[i + k + j*params->n_valid_samples];
						break;
					}
				}
			}

			//f_corr_aux_line[j] = f_in[j] * conj((*(ctx->f_az_chirp[i]+j)); //XXX NOTE1
			//printf("[%d](%g+j%g), ", j, creal(f_corr_aux_line[j]), cimag(f_corr_aux_line[j]));
		}
		//printf("\n");
	}
		exit(0);

//		fftw_execute(p_backward);
//		//if(i == 0) {
//		//printf("IF\n");
//		//for(int j = 0; j < range; j++) {
//		//printf("[%d](%g+j%g), ", j, creal(corr_aux_line[j])/azimuth, cimag(corr_aux_line[j])/azimuth);
//		//}
//		//printf("\n");
//		//}
//		//printf("Shifted\n");
//		for(int j = 0; j < azimuth; j++) {
//			int dst_j = (j + azimuth/2) % azimuth;
//			out[i+range*dst_j] = corr_aux_line[j] / azimuth;
//			//printf("[%d](%g+j%g), ", dst_j, creal(out[i+range*dst_j]), cimag(out[i+range*dst_j]));
//		}
//		//printf("\n");
	free(aux_line);
	free(f_aux_line);
	free(f_corr_aux_line);
	free(corr_aux_line);
	fftw_destroy_plan(p_forward);
	fftw_destroy_plan(p_backward);

	return 0;
}

int sarif_make_rcmc_offset_matrix(Sarif_Ctx *ctx)
{
	ERS_Raw_Parser_Params *params;
	int max_offset = -1;

	if(!ctx) {
		fprintf(stderr, "%s:: Invalid arguments %p\n", __func__, ctx);
		return -1;
	}
	if(!ctx->has_params) {
		fprintf(stderr, "%s:: Error has_params %d\n", __func__, ctx->has_params);
		return -1;
	}
	params = &ctx->params;

	//We dont need both, one array is enough. It's just to make it easier to read
	double *freq_shift = calloc(1, sizeof(double complex) * params->fft_lines);
	double *aux = calloc(1, sizeof(double complex) * params->fft_lines);

	int *offset = calloc(1, sizeof(double complex) * params->fft_lines * params->n_valid_samples);

	for(int i = 0; i < params->fft_lines; i++) {
		freq_shift[i] = ctx->fc + i * params->prf/(double)ctx->valid_az_lines;
		aux[i] = 1/sqrt(1-pow(params->lambda*freq_shift[i]/(2*params->velocity),2)) - 1;

		for(int j = 0; j < params->n_valid_samples; j++) {
			offset[j+i*params->n_valid_samples] = (int)round(aux[i] * ctx->sq_range[j] / (params->C/(2*params->fs))); //-> OK!
			if(offset[j+i*params->n_valid_samples] > max_offset)
				max_offset = offset[j+i*params->n_valid_samples];
		}
	}

	ctx->rcmc_max_offset = max_offset;
	ctx->rcmc_offset_matrix = offset;

	free(aux);
	free(freq_shift);

	return 0;
}
