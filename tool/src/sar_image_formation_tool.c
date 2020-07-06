#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <fcntl.h>
#include <unistd.h>

#include <matio.h>

#include <ers_raw_parser.h>
#include <sarif.h>

#define PATH_SIZE (2000)

static int g_has_input_ldr = 0;
static char g_input_ldr[PATH_SIZE];
static int g_has_input_raw = 0;
static char g_input_raw[PATH_SIZE];
static int g_has_output = 0;
static char g_output_path[PATH_SIZE];
static int g_multilook = 0;

enum Getopt_Short{
	OPT_MULTILOOK = 1,
};

static const char *short_opts= "l:r:o:";
static struct option long_opts[] = {
	{"ldr", required_argument, NULL,'l'},
	{"raw", required_argument, NULL,'r'},
	{"output", required_argument, NULL,'o'},
	{"multilook", no_argument, NULL, OPT_MULTILOOK},
	{"help", no_argument, NULL,'h'},
	{NULL,0,NULL,0}
};

int save_to_mat_double(mat_t *matfp, double *data, int range, int azimuth)
{
	matvar_t *variable2d = NULL;
	char* fieldname2d = "data";
	size_t dim2d[2];
	double *matout_re = NULL;
	double *matout_im = NULL;

	//Create and fill Real and Imaginary arrays
	matout_re = calloc(1, sizeof(double) * azimuth * range);
	matout_im = calloc(1, sizeof(double) * azimuth * range);
	mat_complex_split_t mycomplexdouble = {matout_re, matout_im};
	for (int j = 0; j < range; j++ ) {
		for (int i = 0; i < azimuth; i++ ) {
			size_t mat_idx = azimuth*j+i;
			size_t idx = range*i+j;
			matout_re[mat_idx] = data[idx];
		}
	}

	dim2d[0] = azimuth;
	dim2d[1] = range;
	variable2d = Mat_VarCreate(fieldname2d, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dim2d, &mycomplexdouble, MAT_F_COMPLEX); //rank 2
	if(!variable2d) {
		fprintf(stderr, "%s:: Error creating variable\n", __func__);
		goto save_to_mat_error;
	}

	int dims = 1;
	if(Mat_VarWriteAppend(matfp, variable2d, MAT_COMPRESSION_NONE, dims)) {
		fprintf(stderr, "%s:: Error writing variable\n", __func__);
		goto save_to_mat_error;
	}
	Mat_VarFree(variable2d);

	free(matout_re);
	free(matout_im);

	return 0;


save_to_mat_error:
	Mat_VarFree(variable2d);
	free(matout_re);
	free(matout_im);

	return -1;
}



int save_to_mat(mat_t *matfp, complex double *data, int range, int azimuth)
{
	matvar_t *variable2d = NULL;
	char* fieldname2d = "data";
	size_t dim2d[2];
	double *matout_re = NULL;
	double *matout_im = NULL;

	//Create and fill Real and Imaginary arrays
	matout_re = calloc(1, sizeof(double) * azimuth * range);
	matout_im = calloc(1, sizeof(double) * azimuth * range);
	mat_complex_split_t mycomplexdouble = {matout_re, matout_im};
	for (int j = 0; j < range; j++ ) {
		for (int i = 0; i < azimuth; i++ ) {
			size_t mat_idx = azimuth*j+i;
			size_t idx = range*i+j;
			matout_re[mat_idx] = creal(data[idx]);
			matout_im[mat_idx] = cimag(data[idx]);
		}
	}

	dim2d[0] = azimuth;
	dim2d[1] = range;
	variable2d = Mat_VarCreate(fieldname2d, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dim2d, &mycomplexdouble, MAT_F_COMPLEX); //rank 2
	if(!variable2d) {
		fprintf(stderr, "%s:: Error creating variable\n", __func__);
		goto save_to_mat_error;
	}

	int dims = 1;
	if(Mat_VarWriteAppend(matfp, variable2d, MAT_COMPRESSION_NONE, dims)) {
		fprintf(stderr, "%s:: Error writing variable\n", __func__);
		goto save_to_mat_error;
	}
	Mat_VarFree(variable2d);

	free(matout_re);
	free(matout_im);

	return 0;

save_to_mat_error:
	Mat_VarFree(variable2d);

	free(matout_re);
	free(matout_im);

	return -1;
}

int parse_opts(int argc, char **argv)
{
	int c;

	while ((c = getopt_long(argc, argv, short_opts, long_opts, NULL)) != -1) {
		switch (c){
			case 'l':
				snprintf(g_input_ldr, PATH_SIZE, "%s", optarg);
				g_has_input_ldr = 1;
				break;
			case 'r':
				snprintf(g_input_raw, PATH_SIZE, "%s", optarg);
				g_has_input_raw = 1;
				break;
			case 'o':
				snprintf(g_output_path, PATH_SIZE, "%s", optarg);
				g_has_output = 1;
				break;
			case OPT_MULTILOOK:
				g_multilook = 1;
				break;
			case 'h':
				//print_help(argc, argv);
				return 1;
			default:
				fprintf(stderr, "Bad arguments: Try `%s --help' for further information.\n", argv[0]);
				return -1;
		}
	}

	return 0;
}

long long get_time_usec()
{
	struct timeval now;
	gettimeofday(&now, NULL);
	return now.tv_sec * 1000000 + now.tv_usec;
}

int main(int argc, char **argv)
{
	int ret;
	ERS_Raw_Parser_Ctx *ctx;
	ERS_Raw_Parser_Params params;
	ERS_Raw_Parser_Data_Patch *data;
	long long start;
	mat_t *matfp = NULL;
	Mat_LogInit("ASD");

	if(parse_opts(argc, argv)) {
		//print_help(argv[0]);
		printf("Bad arguments\n");
		return EXIT_FAILURE;
	}

	if(!g_has_input_ldr || !g_has_input_raw) {
		//print_help(argv[0]);
		printf("No input\n");
		return EXIT_FAILURE;
	}

	ctx = ers_raw_parser_alloc(g_input_ldr, g_input_raw);

	if(ers_raw_parser_get_params_from_file(ctx, &params)) {
		fprintf(stderr, "%s:: ers_raw_parser_get_params_from_file(%s)\n", __func__, g_input_ldr);
		return EXIT_FAILURE;
	}

	if(g_has_output) {
		matfp = Mat_CreateVer(g_output_path, NULL, MAT_FT_MAT73);
		if(!matfp) {
			fprintf(stderr, "%s:: Error creating MAT file: Mat_CreateVer(%s, %p, %d)\n", __func__, g_output_path, NULL, MAT_FT_MAT5);
			return EXIT_FAILURE;
		}
	}

	log_ers_params(&params);

	Sarif_Ctx *sarif_ctx;

	sarif_ctx = sarif_ctx_alloc(&params);

	sarif_make_range_chirp(sarif_ctx);

	int az_valid_lines = 0;
	for(int num_patch = 0;; num_patch++) {
		double complex *out;

		start = get_time_usec();
		fprintf(stdout, "%s:: Reading patch %d line %d\n", __func__, num_patch, az_valid_lines*num_patch);
		ret = ers_raw_parser_get_raw_data_from_file(ctx, &data, num_patch * az_valid_lines);
		if(ret == 1) {
			fprintf(stderr, "%s:: EOF?\n", __func__);
			break;
		} else if(ret < 0) {
			fprintf(stderr, "%s:: ers_raw_parser_get_raw_data_from_file(%s)\n", __func__, g_input_raw);
			break;
		}
		printf("%s:: read() took %lld us\n", __func__, get_time_usec()-start);

		start = get_time_usec();
		sarif_remove_mean(data);
		printf("%s:: remove_mean() took %lld us\n", __func__, get_time_usec()-start);

		start = get_time_usec();
		sarif_range_compression(sarif_ctx, data);
		printf("%s:: range_compress() took %lld us\n", __func__, get_time_usec()-start);

		if(num_patch == 0) {
			double fc;
			double complex *range_compressed;
			sarif_get_range_compression_out(sarif_ctx, &range_compressed);
			fc = sarif_calc_doppler_centroid(range_compressed, SARIF_DOPPLER_CENTROID_ALGO_AVGPHASE, &params);
			if(fc == NAN) {
				fprintf(stderr, "%s:: sarif_calc_doppler_centroid(%p, %d, %p) error!", __func__, data, SARIF_DOPPLER_CENTROID_ALGO_AVGPHASE, &params);
				return EXIT_FAILURE;
			} else {
				fprintf(stdout, "Doppler center: %f Hz\n", fc);
			}

			if(sarif_set_doppler_centroid(sarif_ctx, fc)) {
				fprintf(stderr, "%s:: Error: sarif_set_doppler_centroid(%p, %f)", __func__, sarif_ctx, fc);
				return EXIT_FAILURE;
			}

			if(sarif_make_azimuth_chirp(sarif_ctx)) {
				fprintf(stderr, "%s:: Error: sarif_make_azimuth_chirp(%p)", __func__, sarif_ctx);
				return EXIT_FAILURE;
			}
			az_valid_lines = sarif_get_az_valid_lines(sarif_ctx);
			sarif_make_rcmc_offset_matrix(sarif_ctx);
		}

		start = get_time_usec();
		sarif_azimuth_compression(sarif_ctx);
		printf("%s:: azimuth_compress() took %lld us\n", __func__, get_time_usec()-start);

		if(matfp) {
			if(!g_multilook) {
				sarif_get_slc_out(sarif_ctx, &out);
				save_to_mat(matfp, out, params.n_valid_samples, az_valid_lines);
			} else {
				int ml_ra, ml_az;
				double *multilooked;
				sarif_multilook_patch(sarif_ctx);
				sarif_get_multilooked_patch(sarif_ctx, &multilooked, &ml_ra, &ml_az);
				save_to_mat_double(matfp, multilooked, ml_ra, ml_az);
			}
		}

		ers_raw_parser_data_patch_free(data);
	}

	if(matfp) {
		Mat_Close(matfp);
	}

	ers_raw_parser_free(ctx);
	sarif_ctx_free(sarif_ctx);

	return 0;
}

