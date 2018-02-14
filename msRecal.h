/* 
 * File:   msRecalFunctions.h
 * Author: Arzu Tugce Guler
 *
 * Created on 07 April 2016, 01:23
 */

#ifndef MSRECALFUNCTIONS_H
#define	MSRECALFUNCTIONS_H

#include "StringFunctions.h"
#include "msRecalDefs.h"
#include "pepXMLStructures.h"
#include "PeptideProphetDelegate.h"
#include "mzXMLStructures.h"
#include "mzXMLReader.h"

#include <gsl/gsl_vector.h> 
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

/* parameter struct */
typedef struct msrecal_param_type
{
	/* input file names */
	char* mzxml_file;               //input file to be recalibrated
	char* pepxml_file;              //input file with accurate mass

	/* Output file names */
	char* output_mzXML_file;        //output file

	/* selection / filter input */
	int ms_start_scan;              //mzxml scan start
	int ms_end_scan;                //mzxml scan end

	int min_cal;			// minimum number of calibrants a spectrum should have to be recalibrated
	int crop;			// flag if output should be lossless, or if the spectra that cant be calibrated should be cropped
	float bg;			// background, or minimum intensity of peaks to be considered for recalibrating

	double lower_rel_bnd_rt;	// lower rt window boundary
	double upper_rel_bnd_rt;	// upper rt window boundary
	double recal_offset;		// right rt window boundary
	double mmme;			// max mass measurement error 	

	char* score_name;		// the name of the score parameter 
	double min_score_threshold;	// min score filter 
	double max_score_threshold;	// max score filter 

} msrecal_params;

typedef struct peptide_set_type 
{
	char* sequence;
	double retention;	
} peptideset;

// Calibrant description type
typedef struct calibrant_type {
    double mz; /* calculated m/z */
    double peak; /* measured m/z */
    double intensity;
} calibrant;

/* functions and variables for least squares fitting */
struct data {
    double * y; /* measured f (in phony units) */
    double * mz2; /* theoretical m/z */
};

// command line
void showHelp();

// parameter filtering functions 
msrecal_params* readParameters(int argc, char *argv[]);
void initParameters(msrecal_params* params);

// reading peptides
int process_peptide(search_hit sh, msrecal_params* params);
peptideset* build_peptide_set(pmsms_pipeline_analysis pepfile, msrecal_params* params, int* pepnum);

// constructing calibrant set
void build_internal_calibrants(pmzxml_file mzXML_file, peptideset* peptide_set, int set_size, msrecal_params* params);

double calc_mass(char* sequence,int charge);

int sort_type_comp_inv_int(const void *i, const void *j);

void makeCalibrantList(int scan, pscan_peaks mzpeaks, peptideset* peptide_set, msrecal_params* params);

/* Calibrator functions */
int calib_f(const gsl_vector *x, void *params, gsl_vector *f);
int calib_df (const gsl_vector *x, void *params, gsl_matrix *J);
int calib_fdf (const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J);

double mz_recal(double peak);

int sort_type_comp_inv_err(const void *i, const void *j);

int recalibratePeaks(msrecal_params* params);

void applyCalibration(int scan, pscan_peaks mzpeaks);

#endif	/* MSRECALFUNCTIONS_H */

