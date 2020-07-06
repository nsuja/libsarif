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

void sarif_ctx_free(Sarif_Ctx *ctx);

int sarif_set_doppler_centroid(Sarif_Ctx *ctx, double fc);

/**
 * Creates the range chirp vector, in frequency domain and conjugated
 */
int sarif_make_range_chirp(Sarif_Ctx *ctx);

/**
 * Creates the all the azimuth chirp vectors
 *
 * @param[in] fc: Doppler centroid frequency
 */
int sarif_make_azimuth_chirp(Sarif_Ctx *ctx);

int sarif_remove_mean(ERS_Raw_Parser_Data_Patch *in);

int sarif_range_compression(Sarif_Ctx *ctx, ERS_Raw_Parser_Data_Patch *in);

/**
 * @note Do not free out, data is valid until next call
 *
 * @param[out] out: Range compressed output data, values are not scaled. Output is valid until next sarif_range_compression() call, do not free
 */
int sarif_get_range_compression_out(Sarif_Ctx *ctx, double complex **out);

/**
 * Estimated doppler frequency centroid
 * @note To set the value in the context use sarif_set_doppler_centroid()
 */
double sarif_calc_doppler_centroid(double complex *in, Sarif_Doppler_Centroid_Algo algo, ERS_Raw_Parser_Params *params);

/**
 * Makes azimuth compression and gets final SLC
 * @note User must also call sarif_make_rcmc_offset_matrix() before azimuth correction
 *
 * @param[in] ctx: Valid Sarif_Ctx, must have azimuth chirp properly initialized with sarif_make_azimuth_chirp()
 */
int sarif_azimuth_compression(Sarif_Ctx *ctx);

/**
 * Gets last output SLC data patch
 * Data patch has valid values for params.n_valid_samples * az_valid_lines (@see sarif_get_az_valid_lines()). Actual allocated size is params.n_valid_samples & params.fft_lines
 * @note Do not free out. out is valid until next sarif_azimuth_compression() call
 *
 * @param[in] ctx: Valid Sarif_Ctx, must have azimuth chirp properly initialized with sarif_make_azimuth_chirp()
 * @param[out] out: Range and azimuth compressed output data, values are scaled by both fft lengths (range and azimuth)
 */
int sarif_get_slc_out(Sarif_Ctx *ctx, double complex **out);

int sarif_make_rcmc_offset_matrix(Sarif_Ctx *ctx);

int sarif_get_az_valid_lines(Sarif_Ctx *sarif_ctx);

int sarif_multilook_patch(Sarif_Ctx *ctx);

int sarif_get_multilooked_patch(Sarif_Ctx *sarif_ctx, double **multilooked, int *ml_ra, int *ml_az);
