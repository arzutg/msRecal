#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#include "msRecal.h"

static int candidate_list[MAX_ROWS][MAX_CALIBRANTS];
static int scan_window;
static int* numberof_cal_ms_scan;

static int seqlen;
static int fragment[5];
static calibrant calibrant_list[MAX_CALIBRANTS];

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

static int n_calibrants;
static double mz, Ca, Cb;

  
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
            }
	}

	for (j=0; j<sh.search_score_count; j++) {
            ss = sh.search_score_array[j];
            
            if (strstr(ss.name, params->score_name)== NULL) 
                continue;
            else{
                if (ss.value >= params->min_score_threshold && ss.value <= params->max_score_threshold) {
                    process = 1;
                }
                return process;
            }
            
	}
        
	/* None of the regular score measures applied, now looking for hooked ones from peptide prophet */
	if (!process) {
            
            for (j=0; j<sh.analysis_result_count; j++) {
		hookstruct = (void*) sh.analysis_result_array[j].hook;
                if (strstr(sh.analysis_result_array[j].analysis, "peptideprophet") && hookstruct){
                    score = (double*) peptide_prophet_result_property(params->score_name, hookstruct);
                }
		if (score && *score >= params->min_score_threshold && *score <= params->max_score_threshold) {
                    process = 1;
		}
            }
	}
       
	return process;

}

/*Function that builds the peptide set according to score parameters
  uses process_peptide function to check if the requirements are met*/
peptideset* build_peptide_set(pmsms_pipeline_analysis pepfile, msrecal_params* params, int* pepnum)
{    
	int i, j;
	int peptide_count = 0;
	peptideset *peptide, *retval;
	spectrum_query sq;
	search_hit sh;

	// Counting the total number of peptides
	*pepnum = 0;
	retval = NULL;
	
        peptide_count = 0;
	for (i=0; i<pepfile->run_summary_count; i++) {
            peptide_count += pepfile->run_summary_array[i].spectrum_query_count;
	}
        
	peptide = (peptideset*) malloc(sizeof(peptideset)*peptide_count);

	
	for (i=0; i<pepfile->run_summary_count; i++) {
            for (j=0; j<pepfile->run_summary_array[i].spectrum_query_count; j++) {
		sq = pepfile->run_summary_array[i].spectrum_query_array[j];	/* ith search hit */
		sh = sq.search_result_array[0].search_hit_array[0];

		//Check if the score satisfies the boundaries given as cl arg
		if (!process_peptide(sh, params)){
                    continue;
                }
         
		/* Found valid peptide */
		peptide[*pepnum].sequence = strclone(sh.peptide);
		peptide[*pepnum].retention = sq.retention_time_sec;
		*pepnum += 1;
	    }
	}

	if (!retval)
	    retval = (peptideset*) malloc(sizeof(peptideset) * (*pepnum));
	else
	    retval = (peptideset*) realloc(retval, sizeof(peptideset) * (*pepnum));

	for (i=0; i<(*pepnum); i++) {
	    retval[i] = peptide[i];
	}/* for */
	free(peptide);

	return retval;
}

/* Function that builds the collection of internal calibrants */
void build_internal_calibrants(pmzxml_file mzXML_file, peptideset* peptide_set, int set_size, msrecal_params* params)
{
	int i, j, k, unique;
	double min_rt, max_rt;
	scan_attributes satt;

    scan_window = (params->ms_end_scan - params->ms_start_scan)+1;
    numberof_cal_ms_scan = (int*) malloc(sizeof(int) * scan_window);

    //for each mzxml scan
	for (i=0; i<scan_window; i++)
		numberof_cal_ms_scan[i]=0;

	printf("\nCH 1"); fflush(stdout);
	printf("\nset size: %i", set_size); fflush(stdout);

	for(i=0; i<set_size; i++) {
		min_rt = (peptide_set[i].retention + params->recal_offset) - params->lower_rel_bnd_rt;
        max_rt = (peptide_set[i].retention + params->recal_offset) + params->upper_rel_bnd_rt;

        for(j=params->ms_start_scan; j<=params->ms_end_scan; j++) {
        	satt = get_scan_attributes(mzXML_file, j);

            if (satt.retentionTime < min_rt){
            	continue;
            }
            else if (satt.retentionTime > max_rt){
                break;
            }

            unique = 1;
            // Check if the same sequence isn't already in the list for this spectrum
            for (k=0; k<numberof_cal_ms_scan[j-(params->ms_start_scan)]; k++) {
            	if (strcmp(peptide_set[candidate_list[j-1][k]].sequence, peptide_set[i].sequence) == 0) {
                	unique = 0;
                    break;
                }
            }

            if (unique) {
            	//copy peptide index to internal calibrant candidate lists for nearby MS spectra... [mzscannumber][calibrantno]*/
                candidate_list[j-1][numberof_cal_ms_scan[j-(params->ms_start_scan)]] = i;
                numberof_cal_ms_scan[j-(params->ms_start_scan)]++; // number of calibrants for one ms scan
            }
        }
    }
        
    for(j=params->ms_start_scan; j<=params->ms_end_scan; j++) {
    	printf("\nMS scan: %i\tRT:%f\tBase peak mz:%f\n", mzXML_file->scan_id_array[j-1], get_scan_attributes(mzXML_file, j).retentionTime, get_scan_attributes(mzXML_file, j).basePeakMz); fflush(stdout);
        printf("\n%i\t%f", mzXML_file->scan_id_array[j-1], get_scan_attributes(mzXML_file, j).retentionTime); fflush(stdout);
        for(k=0; k<numberof_cal_ms_scan[j-(params->ms_start_scan)]; k++) {
        	printf("%i\t%s\tRT:%f\n", k+1, peptide_set[candidate_list[j-1][k]].sequence, peptide_set[candidate_list[j-1][k]].retention, peptide_set[candidate_list[j-1][k]]); fflush(stdout);
        }
    }
}

/* Monisotopic mass calculator */
double calc_mass(char* sequence,int charge) 
{    
	int i,a,b;
	double mass;

	seqlen=strlen(sequence);

	
	for(i=0;i<5;i++) 
            fragment[i]=0;

	for(i=0; i<seqlen; i++){ // generate molecular formula for fragment  
            a=20-strlen(strchr(AMINO_ACIDS,sequence[i]));
	
            for(b=0;b<5;b++) {
		fragment[b]=fragment[b]+aa[a][b]; 
            }
	}

	fragment[0]=fragment[0]+2; 
	fragment[3]++;  // add H2O     

	// calculate and return mass based on molecular formula    
	mass=0.0;
	for(b=0;b<5;b++) 
		mass += fragment[b]*mimass[b];

	return (mass+charge*HPLUS_MASS)/charge;

}

int sort_type_comp_inv_int(const void *i, const void *j){
	if (((*(calibrant*)j).intensity - (*(calibrant*)i).intensity) > 0)
		return 1;
	else if (((*(calibrant*)j).intensity - (*(calibrant*)i).intensity) < 0)
		return - 1;
	return 0;
	
}

void makeCalibrantList(int scan, pscan_peaks mzpeaks, peptideset* peptide_set, msrecal_params* params){
    
        int i, j, k, z;
        k = 0;
        n_calibrants = 0;
        // Find internal calibrants for the peaks of the spectrum 
        for (i=0; i<mzpeaks->count; i++) {
            //check if the intensity is above the background parameter
            if (mzpeaks->intensities[i] < params->bg){
                k++; //number of times peak < bg
                continue;
            }
            //for each charge
            for(z=1; z<=4; z++) {
                //for each potential calibrant (peptide) of that ms scan
		for(j=0; j<numberof_cal_ms_scan[scan-(params->ms_start_scan)]; j++) {
                    //calculate the theoretical mass of the peptide
                    mz = calc_mass(peptide_set[candidate_list[scan-1][j]].sequence, z);
                    if(fabs((mz - mzpeaks->mzs[i])/mz)<=params->mmme) {	
                        //printf("\n\t-Peak %i: \t\tmass= %f \tintensity= %f", i, mzpeaks->mzs[i], mzpeaks->intensities[i] ); fflush(stdout);
                        //printf("\n\t\tPeptide \tmass= %f", mz); fflush(stdout);
			calibrant_list[n_calibrants].mz = mz;
			calibrant_list[n_calibrants].peak = mzpeaks->mzs[i];
			calibrant_list[n_calibrants].intensity = mzpeaks->intensities[i];
			n_calibrants++; 
                    }
                }
                
                // add cyclosiloxane peaks as potential internal calibrants in row i 
		if(z==1) { // all cyclosiloxanes are singly charged
                    for(j=0; j<5; j++) {
			if(fabs((cyclosiloxanes[j] - mzpeaks->mzs[i])/ cyclosiloxanes[j]) <= params->mmme) {
                            //printf("\n\t-Peak %i: \t\tmass= %f \tintensity= %f", i, mzpeaks->mzs[i], mzpeaks->intensities[i] ); fflush(stdout);
                            //printf("\n\t\tCyclosiloxane \tmass= %f", cyclosiloxanes[j]); fflush(stdout);
                            calibrant_list[n_calibrants].mz = cyclosiloxanes[j];
                            calibrant_list[n_calibrants].peak = mzpeaks->mzs[i];
                            calibrant_list[n_calibrants].intensity = mzpeaks->intensities[i];
                            n_calibrants++; 
			}
                    }
		}
            }     
	}
        //printf("\n\t%i peaks are below background intensity threshold.", k); fflush(stdout);
        //sorts the calibrant_list of the spectrum w.r.t. peak intensities
        qsort(calibrant_list, n_calibrants, sizeof(calibrant), sort_type_comp_inv_int);
}
  
//Rename the function
//Add another function for orbitrap
//Calibration function CAL2 Inverted  
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

int recalibratePeaks(msrecal_params* params, int mass_analyzer){
    int status, SATISFIED, j;
        
    const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;
	double chi, chi_2;
        
    size_t iter=0;
	const size_t pp=2; /* number of free parameters in calibration function */
	const size_t opt=-3; /* degree in calibration function */
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
        if (mass_analyzer == 0){
        	while (n_calibrants >=params->min_cal && !SATISFIED) {
            /* least-squares fit first using all peaks, than removing those that don't fit */
        		for (j=0;j<n_calibrants;j++) {
        			d.y[j] = 1 / calibrant_list[j].peak;
        			d.mz2[j] = calibrant_list[j].mz;
        		}// for
            
        		iter=0;
        		T = gsl_multifit_fdfsolver_lmder;
        		s = gsl_multifit_fdfsolver_alloc (T, n_calibrants, pp, opt); /* pp = 2 parameters, Ca and Cb */
        		func.n = n_calibrants;
        		func.p = pp;
        		func.nevaldf = opt;
        		func.params = &d;
        		gsl_multifit_fdfsolver_set(s,&func,&x.vector);
            
        		do {
        			iter++;
        			status = gsl_multifit_fdfsolver_iterate (s);

        			if (status)
        				break;
        			status=gsl_multifit_test_delta (s->dx, s->x, 1e-9, 1e-9);
        		} while (status==GSL_CONTINUE && iter<500);
 
        		Ca = gsl_vector_get(s->x,1);
        		Cb = gsl_vector_get(s->x,2);
        		chi = gsl_blas_dnrm2(s->f);
        		chi_2 = gsl_blas_dnrm2(s->dx);
        		gsl_multifit_fdfsolver_free(s);
            
        		/* OK, that was one internal recalibration, now lets check if all calibrants are < INTERNAL_CALIBRATION_TARGET, if not, throw these out */
        		/* and recalibrate (as long as we have at least three peaks) */
        		qsort(calibrant_list, n_calibrants, sizeof(calibrant), sort_type_comp_inv_err);
            
        		for(j=n_calibrants-1; j>=0; j--)
        			if (fabs((calibrant_list[j].mz-mz_recal(calibrant_list[j].mz))/calibrant_list[j].mz)<params->target_mme)
        				break;
        			if (j==n_calibrants-1)
        				SATISFIED=1; /* all calibrants < INTERNAL_CALIBRATION_TARGET (e.g. 2.5 ppm) */
        			n_calibrants=j+1; /* remove calibrants that doesn't fit CAL2 better than e.g. 2 ppm */
        	}
	}

	else if (mass_analyzer == 1){
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
				status=gsl_multifit_test_delta (s->x, s->x, 1e-9, 1e-9);
		            } while (status==GSL_CONTINUE && iter<500);

		            Ca = gsl_vector_get(s->x,0);
		            Cb = gsl_vector_get(s->x,1);
		            chi = gsl_blas_dnrm2(s->f,2);
		            chi_2 = gsl_blas_dnrm2(s->f,3);
		            gsl_multifit_fdfsolver_free(s);

		            /* OK, that was one internal recalibration, now lets check if all calibrants are < INTERNAL_CALIBRATION_TARGET, if not, throw these out */
		            /* and recalibrate (as long as we have at least three peaks) */
		            qsort(calibrant_list, n_calibrants, sizeof(calibrant), sort_type_comp_inv_err);

		            for(j=n_calibrants-1; j>=0; j--)
				if (fabs((calibrant_list[j].mz-mz_recal(calibrant_list[j].mz))/calibrant_list[j].mz)<params->target_mme)
		                    break;
		            if (j==n_calibrants-1)
		                SATISFIED=1; /* all calibrants < INTERNAL_CALIBRATION_TARGET (e.g. 2.5 ppm) */
		            n_calibrants=j+1; /* remove calibrants that doesn't fit CAL2 better than e.g. 2 ppm */
		        }

	}

	else if (mass_analyzer == 2){
		//TOF function with 3 parameters

	}

	else{
		printf("Instrument type is not supported for calibration"); fflush(stdout);
		SATISFIED = -1;
	}


	return SATISFIED;
}

void applyCalibration(int scan, pscan_peaks mzpeaks)
{
	int j;	

	//printf("\tFinal calibration for scan %i:\n", scan);
	for(j=0;j<n_calibrants;j++) {
		//printf("\t%f %f %f %.4f\n", calibrant_list[j].peak, calibrant_list[j].mz, mz_recal(calibrant_list[j].peak), 1e6*(mz_recal(calibrant_list[j].peak)-calibrant_list[j].mz)/calibrant_list[j].mz);
		fflush(stdout);
	}// for		

	for (j=0; j<mzpeaks->count; j++) {
		mzpeaks->mzs[j] = Ca/((1/mzpeaks->mzs[j])-Cb);		
	}// for

}
