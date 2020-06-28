#pragma once
#include <complex.h>

#include <ers_raw_parser.h>

int sarif_make_range_chirp(ERS_Raw_Parser_Params *params, double complex **out_range_chirp);

int sarif_remove_mean(ERS_Raw_Parser_Data_Patch *in);

int sarif_range_compression(double complex *out, ERS_Raw_Parser_Data_Patch *in, ERS_Raw_Parser_Params *params, double complex *f_conj_range_chirp);
