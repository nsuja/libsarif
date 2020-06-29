#pragma once
#include <math.h>
#include <complex.h>

#include <ers_raw_parser.h>

typedef enum {
	SARIF_DOPPLER_CENTROID_ALGO_AVGPHASE = 0,

	SARIF_DOPPLER_CENTROID_ALGO_N
} Sarif_Doppler_Centroid_Algo;

int sarif_make_range_chirp(ERS_Raw_Parser_Params *params, double complex **out_range_chirp);

int sarif_remove_mean(ERS_Raw_Parser_Data_Patch *in);

int sarif_range_compression(double complex *out, ERS_Raw_Parser_Data_Patch *in, ERS_Raw_Parser_Params *params, double complex *f_conj_range_chirp);

double sarif_calc_doppler_centroid(double complex *in, Sarif_Doppler_Centroid_Algo algo, ERS_Raw_Parser_Params *params);
