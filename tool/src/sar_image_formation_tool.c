#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include <ers_raw_parser.h>
#include <sarif.h>

#define INPUT_PATH_SIZE (2000)

static int g_has_input_ldr = 0;
static char g_input_ldr[INPUT_PATH_SIZE];
static int g_has_input_raw = 0;
static char g_input_raw[INPUT_PATH_SIZE];

static const char *short_opts= "l:r:";
static struct option long_opts[] = {
	{"ldr", required_argument, NULL,'l'},
	{"raw", required_argument, NULL,'r'},
	{"help", no_argument, NULL,'h'},
	{NULL,0,NULL,0}
};

int parse_opts(int argc, char **argv)
{
	int c;

	while ((c = getopt_long(argc, argv, short_opts, long_opts, NULL)) != -1) {
		switch (c){
			case 'l':
				snprintf(g_input_ldr, INPUT_PATH_SIZE, "%s", optarg);
				g_has_input_ldr = 1;
				break;
			case 'r':
				snprintf(g_input_raw, INPUT_PATH_SIZE, "%s", optarg);
				g_has_input_raw = 1;
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

int main(int argc, char **argv)
{
	int ret;
	ERS_Raw_Parser_Ctx *ctx;
	ERS_Raw_Parser_Params params;
	ERS_Raw_Parser_Data_Patch *data;

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

	log_ers_params(&params);

	if(ers_raw_parser_get_raw_data_from_file(ctx, &data)) {
		fprintf(stderr, "%s:: ers_raw_parser_get_raw_data_from_file(%s)\n", __func__, g_input_raw);
		return EXIT_FAILURE;
	}

	Sarif_Ctx *sarif_ctx;

	sarif_ctx = sarif_ctx_alloc(&params);

	double complex *f_conj_range_chirp, *f_conj_az_chirp;
	sarif_make_range_chirp(&params, &f_conj_range_chirp);

	sarif_remove_mean(data);

	double complex *processed = calloc(1, sizeof(double complex) * params.n_valid_samples * data->n_az); //1 patch
	sarif_range_compression(processed, data, &params, f_conj_range_chirp, 0);

	double fc;
	fc = sarif_calc_doppler_centroid(processed, SARIF_DOPPLER_CENTROID_ALGO_AVGPHASE, &params);
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

	double complex *out;

	sarif_make_rcmc_offset_matrix(sarif_ctx);

	sarif_azimuth_compression(sarif_ctx, &out, processed, 1);

	ers_raw_parser_data_patch_free(data);

	ers_raw_parser_free(ctx);

	return 0;
}

