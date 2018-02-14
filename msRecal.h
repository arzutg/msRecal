#ifndef __MSRECAL_H__
#define __MSRECAL_H__

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include "msRecalDefs.h"
#include "mzXMLStructures.h"
#include "pepXMLStructures.h"


// training set description structure
typedef struct training_set_type {
	char* sequence;
	double retention;	
} trainingset;


// Calibrant description type
typedef struct calibrant_type {
    double mz; /* calculated m/z */
    double peak; /* measured m/z */
    double intensity;
} calibrant;


// Sorting type
typedef struct sort_type_tp {
    long index;
    double peak;
    char charge;
    double intensity;
    int spectrum;
} sort_type;


/* functions and variables for least squares fitting */
struct data {
	double * y; /* measured f (in phony units) */
    double * mz2; /* theoretical m/z */
};


/* parameter struct */
typedef struct msrecal_param_type
{
	/* input file names */
	char* mzxml_file;
	char* pepxml_file[MAX_PEP_FILES];
	int no_pep_files;

	/* Output file names */
	char* output_mzXML_file;

	/* selection / filter input */
	int ms_start_scan;
	int ms_end_scan;

	int min_cal;				// minimum number of calibrants a spectrum should have to be recalibrated
	int crop;					// flag if output should be lossless, or if the spectra that cant be calibrated should be cropped
	float bg;					// background, or minimum intensity of peaks to be considered for recalibrating

	double lower_rel_bnd_rt;	/* lower rt window boundary */
	double upper_rel_bnd_rt;	/* upper rt window boundary */
	double recal_offset;		/* right rt window boundary */
	double mmme;				/* max mass measurement error */	

	char* score_name;			/* the name of the score parameter */
	double min_score_threshold;	/* min score filter */
	double max_score_threshold;	/* max score filter */

} msrecal_params;



/* Calibrator functions */
int calib_f(const gsl_vector *x, void *params, gsl_vector *f);
int calib_df (const gsl_vector *x, void *params, gsl_matrix *J);
int calib_fdf (const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J);

/* parameter filtering functions */
msrecal_params* read_parameters(int argc, char *argv[]);
void init_parameters(msrecal_params* params);

/* Training set filtering functions */
trainingset* build_training_set(pmsms_pipeline_analysis pepfile[], msrecal_params* params, int* pepnum);
int process_peptide(search_hit sh, msrecal_params* params);

#ifdef __cplusplus
}
#endif

#endif