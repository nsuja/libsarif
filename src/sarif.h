#pragma once
#include <math.h>
#include <complex.h>

#include <ers_raw_parser.h>

typedef enum {
	SARIF_DOPPLER_CENTROID_ALGO_AVGPHASE = 0,

	SARIF_DOPPLER_CENTROID_ALGO_N
} Sarif_Doppler_Centroid_Algo;

typedef struct Sarif_Ctx Sarif_Ctx;

Sarif_Ctx * sarif_ctx_alloc(ERS_Raw_Parser_Params *params);

int sarif_set_doppler_centroid(Sarif_Ctx *ctx, double fc);

/**
 * Creates the range chirp vector, in frequency domain and conjugated
 */
int sarif_make_range_chirp(ERS_Raw_Parser_Params *params, double complex **out_range_chirp);

/**
 * Creates the all the azimuth chirp vectors
 *
 * @param[in] fc: Doppler centroid frequency
 */
int sarif_make_azimuth_chirp(Sarif_Ctx *ctx);

int sarif_remove_mean(ERS_Raw_Parser_Data_Patch *in);

/**
 * @param[in] scale_fft: Divide the resulting FFT by the number of samples in range. This makes it run slower since it iterates over each output element
 */
int sarif_range_compression(double complex *out, ERS_Raw_Parser_Data_Patch *in, ERS_Raw_Parser_Params *params, double complex *f_conj_range_chirp, int scale_fft);

/**
 * @param[in] fc: Doppler centroid frequency
 * @param[in] descale_fft: Divide the input FFT by the number of samples in range. Set if scale_fft=0 in @see sarif_range_compression. Does not affect performance
 */
double sarif_calc_doppler_centroid(double complex *in, Sarif_Doppler_Centroid_Algo algo, ERS_Raw_Parser_Params *params);

int sarif_get_rcmc_data(double complex *in, ERS_Raw_Parser_Params *params);
