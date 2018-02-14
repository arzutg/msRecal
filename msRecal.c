// msRecal.cpp : Defines the entry point for the console application.
//

/*                                                                                                                        */
/* recalibrate_using_MSMS - recalibrates 2D peaklist using Mascot MS/MS data                                              */
/* recal2 - an updated recalibrate_using_MSMS to read mzXML input instead of ASCII 2D peaklists                           */
/*                                                                                                                        */
/* (c) Magnus Palmblad, Division of Ion Physics, Uppsala University, 1999 - 2002                                          */ 
/*                                                                                                                        */
/* Usage: recal2 recal2 -m <mzXML file> -r <rt training file> -o <output> -s <scan range> -O <offset> -M <modifications>  */ 
/* where <2D peaklist> is a combined peaklist file from a 2D dataset (e.g. from mipp2D_to_peaklist) and                   */
/* <training file> is rt training file output form combine_masslist_and_Mascot_for_rt.c                                   */
/* <offset> is the offset between MGF file queries and row (spectra) numbers in 2D file                                   */
/* if using alignment (msalign), use awk '{print $1, $4}' xxx.alignment > xxx.training and the xxx.training for recal2    */
/*                                                                                                                        */
/* compile with gcc -o recal2 base64.c ramp.c recal2.c -I. -lgd -lm -lz -lgsl -lgslcblas                                  */
/*                                                                                                                        */

#include <stdio.h>
#include <stdlib.h>  
#include <ctype.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include "msRecal.h"
#include "StringFunctions.h"
#include "mzXMLStructures.h"
#include "mzXMLReader.h"
#include "mzXMLWriter.h"
#include "pepXMLReader.h"
#include "PeptideProphetDelegate.h"


static int seqlen;
static int fragment[5]; /* for elemental composition of fragment */
static double mz,Ca,Cb;
static int n_calibrants;
static calibrant calibrant_list[MAX_CALIBRANTS];

static int candidate_list[MAX_ROWS][MAX_CALIBRANTS];
static int* scan_cal_index;
static int scan_window;

static double mimass[5] = 
{
	1.0078250321,
	12,
	14.0030740052,
	15.9949146221,
	31.97207069
};

static double cyclosiloxanes[5] = 
{
	593.157605,
	667.176396,
	741.195187,
	815.213979,
	889.232770
};

static unsigned char aa[20][5] = 
{
	{5,3,1,1,0}, /* amino acid compositions */
	{12,6,4,1,0},
	{6,4,2,2,0},
	{5,4,1,3,0},
	{5,3,1,1,1}, 
	{7,5,1,3,0},
	{8,5,2,2,0},
	{3,2,1,1,0},
	{7,6,3,1,0},
	{11,6,1,1,0},
	{11,6,1,1,0},
	{12,6,2,1,0},
	{9,5,1,1,1},
	{9,9,1,1,0},
	{7,5,1,1,0},
	{5,3,1,2,0},
	{7,4,1,2,0},
	{10,11,2,1,0},
	{9,9,1,2,0},
	{9,5,1,1,0}
};
  

// Calibration function CAL2 Inverted  
int calib_f(const gsl_vector *x, void *params, gsl_vector *f)
{
	double *y = ((struct data *)params)->y;
	double *mz = ((struct data *)params)->mz2;
	double a = gsl_vector_get (x, 0);
	double b = gsl_vector_get (x, 1);
	double M;
    size_t i;
      
	for (i=0;i<n_calibrants;i++) {	  
		/* Model m = a/(f-b) (CAL2 inverted) */	 
		M = a/(y[i]-b);
		gsl_vector_set (f, i, (M-mz[i])); /* absolute or relative error? */
	}// for
      
    return GSL_SUCCESS;

}// int calib_f(const gsl_vector *x, void *params, gsl_vector *f)
  
// DF calibrator
int calib_df(const gsl_vector *x, void *params, gsl_matrix *J)
{
	double *y = ((struct data *)params)->y;
	double a = gsl_vector_get (x, 0);
	double b = gsl_vector_get (x, 1);
	size_t i;
      
	for (i=0;i<n_calibrants;i++) {
		gsl_matrix_set (J,i,0, 1/(y[i]-b) ); 
		gsl_matrix_set (J,i,1, a/((y[i]-b)*(y[i]-b)) );
	}// for
      
	return GSL_SUCCESS;

}// int calib_df (const gsl_vector *x, void *params, gsl_matrix *J)
  
// FDF Calibrator
int calib_fdf(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J)
{
	calib_f (x,params,f);
	calib_df (x,params,J);
      
	return GSL_SUCCESS;
}  

/* State printer */
void print_state(size_t iter, gsl_multifit_fdfsolver *s)
{
	printf ("iter: %3u x = %15.8f %15.8f " "|f(x)| = %g\n", iter, gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1), gsl_blas_dnrm2 (s->f));

}

/* Monisotopic mass calculator */
double calc_mass(char* sequence,int charge) 
{    
	int i,a,b;
	double mass;

	seqlen=strlen(sequence);

	/* init */
	for(i=0;i<5;i++) 
		fragment[i]=0;

	for(i=0; i<seqlen; i++) /* generate molecular formula for fragment */ {
		a=20-strlen(strchr(AMINO_ACIDS,sequence[i]));
	
		for(b=0;b<5;b++) {
			fragment[b]=fragment[b]+aa[a][b]; 
		}// for
	}// for

	fragment[0]=fragment[0]+2; 
	fragment[3]++;  /* add H2O */    

	/* calculate and return mass based on molecular formula */    
	mass=0.0;
	for(b=0;b<5;b++) 
		mass += fragment[b]*mimass[b];

	return (mass+charge*HPLUS_MASS)/charge;

}// double calc_mass

void showhelp()
{
	printf("Usage: msRexal.exe [parameters]\n" ); 
	printf("Compulsory Parameters:\n"); 
	printf("-p<string>\tpepXML file location. Multiple locations can be specified, all must be preceeded by -p.\n"); 
	printf("-m<string>\tmzXML file location.\n"); 
	printf("-o<string>\tmzXML output file location. May not be equal to input location.\n"); 
	printf("-m<float>\tMaximum mass measurement error.\n"); 
	printf("-t<string>\tScore name. Name of the score in the pepXML file, which value will be considered in comparison to bounds.\n"); 
	printf("\n"); 
	printf("Optional Parameters:\n"); 
	printf("-s<int,int>\tmzXML scan window. The first number is the beginscan, the last one the endscan.\n"); 
	printf("-r<float,float>\tRetention time window. The first number is the -delta on a peptide with respect to a scan spectrum, the second number the +delta.\n"); 
	printf("-L<float>\tLower peptide score bound.\n"); 
	printf("-U<float>\tUpper peptide score bound.\n"); 
	printf("-C<int>\tMinimum number of calibrants needed for spectrum recalibration.\n"); 
	printf("-c\tCrop flag. By specifying this, spectra that cannot be recalibrated will be emptied.\n"); 
	printf("-b\tBackground intensity. Only peaks with an intensity higher than this will be used in recalibration.\n"); 
	fflush(stdout);
}

/* reads in user-defined parameters */
msrecal_params* read_parameters(int argc, char *argv[])
{
	char temp[100];	/* temp name */
	char* p;		/* token pointer */
	int i;			/* string length */
	

	msrecal_params* params = (msrecal_params*) malloc(sizeof(msrecal_params));
	init_parameters(params);
  
	/* Assuming first parameter is the program name */
	for(i=1; i<argc; i++) {	  
		if (argv[i][0] == '-') {	
			/* filtering out the parameter string */
			p = &argv[strlen(argv[i])>2? i: i+1][strlen(argv[i])>2? 2: 0];

			/* PepXMLfile specification plus alignment output standard */
			if (argv[i][1] == 'p') {
				params->pepxml_file[params->no_pep_files] = strclone(p);
				params->no_pep_files += 1;
			}/* if */
						
			/* mzXML file input specification */
			else if (argv[i][1] == 'm') {
				params->mzxml_file = strclone(p);
			}

			/* max mass measurement error */
			else if (argv[i][1] == 'e') 
				params->mmme = atof(p);

			/* Scan range */
			else if (argv[i][1]=='s') {
				strcpy(temp, p); 
				p = strtok(temp,","); 
				params->ms_start_scan = atoi(p); 
				p = strtok('\0',","); 
				params->ms_end_scan = atoi(p);
			}// else if

			/* Retention time window */
			else if (argv[i][1]=='r') {
				strcpy(temp, p); 
				p = strtok(temp,","); 
				params->lower_rel_bnd_rt = atof(p); 
				p = strtok('\0',","); 
				params->upper_rel_bnd_rt = atof(p);
			}// else if

			/* Minimum score threshold */
			else if (argv[i][1] == 'L'){
				params->min_score_threshold = atof(p);
			}// else if

			/* Maximum score threshold */
			else if (argv[i][1] == 'U'){
				params->max_score_threshold = atof(p);
			}// else if

			/* Threshold name parameter */
			else if (argv[i][1] == 't') {
				params->score_name = strclone(p);
			}// else if

			/* Crop flag */
			else if (argv[i][1] == 'c') {
				params->crop = 1;
			}// else if

			/* Minimum Calibrant flag */
			else if (argv[i][1] == 'C') {
				params->min_cal = atoi(p);
			}// else if

			/* background / minimum intensity flag*/
			else if (argv[i][1] == 'b') {
				params->bg = atof(p);
			}// else if

			/* Output file specification */
			else if (argv[i][1]=='o') {
				params->output_mzXML_file = strclone(p);
			}// else if

			else if (argv[i][1]=='h') {
				showhelp();
				exit(-1);
			}// else if
  		}// if
	}// for
	
	// Checking integrity of the parameter file
	if (params->no_pep_files == 0 || params->mzxml_file == NULL || params->mmme < 0 || params->score_name == NULL || params->output_mzXML_file == NULL)		
		return NULL;	

	return params;

}// msrecal_params read_parameters(int argc, char *argv[])

/* Initialization routine for parameters struct */
void init_parameters(msrecal_params* params)
{
	int i;
	/* input file names */
	params->mzxml_file = NULL;

	for (i=0; i<MAX_PEP_FILES; i++)
		params->pepxml_file[i] = NULL;
	params->no_pep_files = 0;
	params->output_mzXML_file = NULL;
	params->score_name = NULL;

	/* selection / filter input */
	params->ms_start_scan = INT_MIN;
	params->ms_end_scan = INT_MAX;		

	params->crop = DEFAULT_MODE;
	params->min_cal = DEFAULT_MIN_CALIBRANTS;
	params->bg = DEFAULT_MIN_INTENSITY;

	params->lower_rel_bnd_rt = FLT_MAX;
	params->upper_rel_bnd_rt = FLT_MAX;
	params->recal_offset = 0;
	
	params->mmme = -1;				
	params->min_score_threshold = DBL_MIN;
	params->max_score_threshold = DBL_MAX;
	
}/* void init_parameters(parameters* params) */

/* Function that determines if the peptide should be processed or not */
int process_peptide(search_hit sh, msrecal_params* params)
{
	int process = 0, j;
	search_score ss;
	double* score = NULL;
	void* hookstruct;	

	for (j=0; j<strlen(sh.peptide); j++) {
		if (strchr(ANTI_ACIDS, sh.peptide[j])) {
			return 0;
		}// if
	}// if
	
	for (j=0; j<sh.search_score_count; j++) {
		ss = sh.search_score_array[j];

		if (!strstr(ss.name, params->score_name))
			continue;
		
		if (ss.value >= params->min_score_threshold && ss.value <= params->max_score_threshold) {
			process = 1;
		}/* if */
	}/* for */	

	/* None of the regular score measures applied, now looking for hooked ones from peptide prophet */
	if (!process) {
		for (j=0; j<sh.analysis_result_count; j++) {
			hookstruct = (void*) sh.analysis_result_array[j].hook;

			if (strstr(sh.analysis_result_array[j].analysis, "peptideprophet") && hookstruct)
				score = (double*) peptide_prophet_result_property(params->score_name, hookstruct);

			if (score && *score >= params->min_score_threshold && *score <= params->max_score_threshold) {
				process = 1;
			}/* if */
		}// for
	}// if

	return process;

}/* int process_peptide(search_hit sh, parameters* params) */

/* Builds a training set of peptides */
trainingset* build_training_set(pmsms_pipeline_analysis pepfile[], msrecal_params* params, int* pepnum)
{
	int i, j, filecnt, pepcnt;
	int peptide_count = 0;
	trainingset *peptide, *retval;
	spectrum_query sq;
	search_hit sh;

	// Counting the total number of peptides
	*pepnum = 0;
	pepcnt = 0;
	retval = NULL;
	for (filecnt = 0; filecnt < params->no_pep_files; filecnt++) {
		peptide_count = 0;
		for (i=0; i<pepfile[filecnt]->run_summary_count; i++) {
			peptide_count += pepfile[filecnt]->run_summary_array[i].spectrum_query_count;
		}// for

		peptide = (trainingset*) malloc(sizeof(trainingset)*peptide_count);

		/* Walking all spectrum results */
		for (i=0; i<pepfile[filecnt]->run_summary_count; i++) {
			for (j=0; j<pepfile[filecnt]->run_summary_array[i].spectrum_query_count; j++) {
				sq = pepfile[filecnt]->run_summary_array[i].spectrum_query_array[j];	/* ith search hit */
				sh = sq.search_result_array[0].search_hit_array[0];
						
				/* First we check if the peptide score is sufficient */
				if (!process_peptide(sh, params))
					continue;

				/* Found valid peptide */
				peptide[(*pepnum)-pepcnt].sequence = strclone(sh.peptide);
				peptide[(*pepnum)-pepcnt].retention = sq.retention_time_sec;
				*pepnum += 1;
			}// for
		}/* for */

		if (!retval)
			retval = (trainingset*) malloc(sizeof(trainingset) * (*pepnum));
		else
			retval = (trainingset*) realloc(retval, sizeof(trainingset) * (*pepnum));

		for (i=pepcnt; i<(*pepnum); i++) {
			retval[i] = peptide[i-pepcnt];
		}/* for */
		pepcnt = *pepnum;
		free(peptide);
	}// for

	return retval;

} // trainingset* build_training_set(pmsms_pipeline_analysis pepfile, msrecal_params* params, int* pepnum)

/* Function that builds the collection of internal calibrants */
void build_internal_calibrants(pmzxml_file mzXML_file, trainingset* training_set, int set_size, msrecal_params* params)
{
	int i, j, k, unique;
	double min_rt, max_rt;
	scan_attributes satt;

	for (i=0; i<scan_window; i++) 
		scan_cal_index[i]=0;

	for(i=0; i<set_size;i++) {
		if (params->lower_rel_bnd_rt != FLT_MAX)
			min_rt = (training_set[i].retention + params->recal_offset) - params->lower_rel_bnd_rt;
		else
			min_rt = -FLT_MAX;

		if (params->upper_rel_bnd_rt != FLT_MAX)
			max_rt = (training_set[i].retention + params->recal_offset) + params->upper_rel_bnd_rt;
		else
			max_rt = FLT_MAX;

		for(j=params->ms_start_scan; j<=params->ms_end_scan; j++) {
			satt = get_scan_attributes(mzXML_file, j);

			if (satt.retentionTime < min_rt)
				continue;
			else if (satt.retentionTime > max_rt)
				break;
			
			unique = 1;
			/* Check if the same sequence isn't already in the list for this spectrum */
			for (k=0; k<scan_cal_index[j-(params->ms_start_scan)]; k++) {
				if (strcmp(training_set[candidate_list[j-1][k]].sequence, training_set[i].sequence) == 0) {
					unique = 0;
					break;
				}// if			
			}// for

			if (unique) {
				candidate_list[j-1][scan_cal_index[j-(params->ms_start_scan)]] = i; /* copy peptide sequence to internal calibrant candidate lists for nearby MS spectra... [row][calibrant]*/				
				scan_cal_index[j-(params->ms_start_scan)]++; /* index is the calibrant number in a row (scan, spectrum) */			
			}
		}// for
	}// for
}// build_internal_calibrants(pmzxml_file mzXML_file, trainingset* training_set, int set_size, msrecal_params* params)

/* compare the integers */
int sort_type_comp_inv_int(const void *i, const void *j)
{
	if (((*(calibrant*)j).intensity - (*(calibrant*)i).intensity) > 0)
		return 1;
	else if (((*(calibrant*)j).intensity - (*(calibrant*)i).intensity) < 0)
		return - 1;
	return 0;
	
}// int comp(const void *i, const void *j)

double mz_recal(double peak)
{
	return Ca/((1/peak)-Cb);

}/* double mz_recal(double peak) */

/* compare the integers */
int sort_type_comp_inv_err(const void *i, const void *j)
{
	calibrant *ip, *jp;
	double erri, errj;

	ip = (calibrant*)i;
	jp = (calibrant*)j;

	errj = fabs((jp->mz- mz_recal(jp->peak))/jp->mz);
	erri = fabs((ip->mz - mz_recal(ip->peak))/ip->mz);

	if ((errj - erri) > 0)
		return 1;
	else if ((errj - erri) < 0)
		return -1;
	return 0;

}// int comp(const void *i, const void *j)
 
static void makeCalibrantList(int scan, pscan_peaks mzpeaks, trainingset* training_set, msrecal_params* params) 
{
	int i, j, z;
	
	n_calibrants = 0;

	// Find internal calibrants for the peaks of the spectrum 
	for (i=0; i<mzpeaks->count; i++) {
		// We do not consider peaks beneath the background theshold
		if (mzpeaks->intensities[i] < params->bg)
			continue;

		for(z=1; z<=4; z++) {
			for(j=0; j<scan_cal_index[scan-(params->ms_start_scan)]; j++) {
				mz = calc_mass(training_set[candidate_list[scan-1][j]].sequence, z);					
				if(fabs((mz - mzpeaks->mzs[i])/mz)<=EXTERNAL_CALIBRATION_TOLERANCE) {					
					calibrant_list[n_calibrants].mz = mz;
					calibrant_list[n_calibrants].peak = mzpeaks->mzs[i];
					calibrant_list[n_calibrants].intensity = mzpeaks->intensities[i];
					n_calibrants++; 
				}// if
			}// for

			/* add cyclosiloxane peaks as potential internal calibrants in row i */
			if(z==1) /* all cyclosiloxanes are singly charged, try -1 to exclude these */ {
				for(j=0; j<5; j++) {
					if(fabs((cyclosiloxanes[j] - mzpeaks->mzs[i])/ cyclosiloxanes[j]) <= EXTERNAL_CALIBRATION_TOLERANCE) {
						calibrant_list[n_calibrants].mz = cyclosiloxanes[j];
						calibrant_list[n_calibrants].peak = mzpeaks->mzs[i];
						calibrant_list[n_calibrants].intensity = mzpeaks->intensities[i];
						n_calibrants++; 
					}// if
				}// for
			}// if
		}// for
	}// for	
}// static void makeCalibrantList(int scan, pscan_peaks mzpeaks, trainingset* training_set) 


static int recalibratePeaks(msrecal_params* params)
{
	int status, SATISFIED, j;

	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	double chi;
	
	size_t iter=0;
	const size_t pp=2; /* number of free parameters in calibration function */
	double y[MAX_CALIBRANTS];
	double mz2[MAX_CALIBRANTS];
	struct data d={y,mz2};
	double x_init[2]={1.0,0.0}; /* start here, close to minimum if reasonably calibrated beforehand */

	gsl_multifit_function_fdf func;
	gsl_vector_view x=gsl_vector_view_array(x_init,pp);
	func.f = &calib_f;
	func.df = &calib_df;
	func.fdf = &calib_fdf;

	SATISFIED=0;
	while (n_calibrants >=params->min_cal && !SATISFIED) {
		/* least-squares fit first using all peaks, than removing those that don't fit */
		for (j=0;j<n_calibrants;j++) {
			d.y[j] = 1 / calibrant_list[j].peak;
			d.mz2[j] = calibrant_list[j].mz;
		}// for

		iter=0;
		T = gsl_multifit_fdfsolver_lmder;
		s = gsl_multifit_fdfsolver_alloc (T, n_calibrants, pp); /* pp = 2 parameters, Ca and Cb */
		func.n = n_calibrants;
		func.p = pp;
		func.params = &d;
		gsl_multifit_fdfsolver_set(s,&func,&x.vector);

		do {
			iter++;
			status = gsl_multifit_fdfsolver_iterate (s);

			if (status) 
				break;
			status=gsl_multifit_test_delta (s->dx, s->x, 1e-9, 1e-9);
		} while (status==GSL_CONTINUE && iter<500);
		
		Ca = gsl_vector_get(s->x,0); 
		Cb = gsl_vector_get(s->x,1);
		chi = gsl_blas_dnrm2(s->f);
		gsl_multifit_fdfsolver_free(s);
		
		/* OK, that was one internal recalibration, now lets check if all calibrants are < INTERNAL_CALIBRATION_TARGET, if not, throw these out */
		/* and recalibrate (as long as we have at least three peaks) */
		qsort(calibrant_list, n_calibrants, sizeof(calibrant), sort_type_comp_inv_err);

		for(j=n_calibrants-1; j>=0; j--)  
			if (fabs((calibrant_list[j].mz-mz_recal(calibrant_list[j].mz))/calibrant_list[j].mz)<INTERNAL_CALIBRATION_TARGET) 
				break;			
		if (j==n_calibrants-1) 
			SATISFIED=1; /* all calibrants < INTERNAL_CALIBRATION_TARGET (e.g. 2.5 ppm) */
		n_calibrants=j+1; /* remove calibrants that doesn't fit CAL2 better than e.g. 2 ppm */
	}// while

	return SATISFIED;
}


static void applyCalibration(int scan, pscan_peaks mzpeaks)
{
	int j;	

	printf("\tFinal calibration for scan %i:\n", scan);
	for(j=0;j<n_calibrants;j++) {
		printf("\t%f %f %f %.4f\n", calibrant_list[j].peak, calibrant_list[j].mz, mz_recal(calibrant_list[j].peak), 1e6*(mz_recal(calibrant_list[j].peak)-calibrant_list[j].mz)/calibrant_list[j].mz); 
		fflush(stdout);
	}// for		

	for (j=0; j<mzpeaks->count; j++) {
		mzpeaks->mzs[j] = Ca/((1/mzpeaks->mzs[j])-Cb);		
	}// for

}

int main(int argc, char *argv[]) 
{
	msrecal_params* params;
	trainingset* training_set;
	pdelegate_list dlgl;
	pdelegate_type dlg; 
	int pepnum;
 
	long scan;
	int pepfile_cnt;
	 
	int SATISFIED; 

	pmsms_pipeline_analysis pepXML_file[MAX_PEP_FILES];
	pmzxml_file mzXML_file;  
	pscan mzscan;
	scan_peaks mzpeaks;
	
	params = read_parameters(argc, argv);
	if (!params){
		showhelp();
		return -1;
	}// if

	/* reading training set file */
	dlgl = make_delegate_list();
	dlg = make_peptide_prophet_result_delegate();
	add_delegate_to_list(dlgl, dlg);

	for (pepfile_cnt=0; pepfile_cnt < params->no_pep_files; pepfile_cnt++) {
		printf("reading pepXML file %s...", params->pepxml_file[pepfile_cnt]); fflush(stdout);	
		pepXML_file[pepfile_cnt] = read_pepxml_file(params->pepxml_file[pepfile_cnt], 0, 0, dlgl);
		if (pepXML_file[pepfile_cnt] == NULL) {
			printf("error opening pepXML file (of potential internal calibrants)\n"); fflush(stdout);
			return -1;
		}// if
	}// for

	printf("done.\nFiltering peptides... "); fflush(stdout);
	
	// Filtering out peptides and making them unique
	training_set = build_training_set(pepXML_file, params, &pepnum);
	printf("done (read %i peptides in potential calibrant set)\n", pepnum); fflush(stdout);	

	// Reading mzXML file
	printf("reading mzXML dataset %s...",params->mzxml_file); fflush(stdout);
	mzXML_file = read_mzxml_file_spectrum(params->mzxml_file, 0, 0, params->ms_start_scan, params->ms_end_scan);
	if (mzXML_file == NULL) {
		printf("error opening mzXML file (of scans to be calibrated)\n"); fflush(stdout);
		return -1;
	}// if

	// Establishing borders
	if (mzXML_file->scan_num < params->ms_end_scan)
		params->ms_end_scan = mzXML_file->scan_num;
	if (params->ms_start_scan < 1)
		params->ms_start_scan = 1;		
	scan_window = (params->ms_end_scan - params->ms_start_scan)+1;
	scan_cal_index = (int*) malloc(sizeof(int) * scan_window);
	printf("done\n"); fflush(stdout);
	
	// Make internal calibrant candidate list from pepxml file
	printf("Making internal calibrant candidate list..."); fflush(stdout);
	build_internal_calibrants(mzXML_file, training_set, pepnum, params);
	printf("done. Recalibrating mzXML data...\n"); fflush(stdout);

	/* For each spectrum / scan, we load the peaks */
	for(scan=params->ms_start_scan; scan<=params->ms_end_scan; scan++) {
		mzscan = get_scan(mzXML_file, scan, 0);

		if(mzscan->attributes.msLvl != 1 ) /* MS spectrum */
			continue;
		
		// Considering the scans not valid
		if(!mzscan->peaks || mzscan->peaks->count < 0) {
			if (params->crop)
				empty_scan(mzXML_file, scan);
			continue;
		}// if 

		// Loading and filtering the peaks
		mzpeaks = load_scan_peaks(mzXML_file, scan);			

		// if none left, free the peaks and empty them
		if (mzpeaks.count == 0) {
			unload_scan_peaks(mzXML_file, scan);			
			continue;
		}// if

		// Make the list of calibrants for this spectrum  and sort them in descending order of intensity  
		makeCalibrantList(scan, &mzpeaks, training_set, params);		
		qsort(calibrant_list, n_calibrants, sizeof(calibrant), sort_type_comp_inv_int);

		// Recalibrate peaks
		SATISFIED = recalibratePeaks(params);
		
		if (SATISFIED) {
			applyCalibration(scan, &mzpeaks);	
			update_scan_peaks(mzXML_file, scan, mzpeaks.count, 32, mzpeaks.mzs, mzpeaks.intensities);
			printf("-------------------------------------------------------------\n"); fflush(stdout);
		}// if
		else {
			unload_scan_peaks(mzXML_file, scan);
			if (params->crop)
				empty_scan(mzXML_file, scan);						
		}// else										
	}// for
	
	printf("done\nWriting recalibrated mzXML file to destination %s...", params->output_mzXML_file); fflush(stdout);
	write_mzxml_file(mzXML_file, params->output_mzXML_file);
	printf("done\n"); fflush(stdout);

	return 0;
}
