//13-12-2017
/* 
 * File:   main.c
 * Author: Arzu Tugce Guler
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "msRecal.h"
#include "pepXMLReader.h"
#include "PeptideProphetDelegate.h"
#include "mzXMLStructures.h"
#include "mzXMLWriter.h"


int main(int argc, char *argv[]) {
    
    msrecal_params* params;
    pdelegate_list dlgl;
    pdelegate_type dlg;
    pmsms_pipeline_analysis pepXML_file;
    pmzxml_file mzXML_file;  
    pscan mzscan;
    scan_peaks mzpeaks;
    peptideset* peptide_set;

    int pepnum, k;
    long scan;
    
    int SATISFIED;
     
    //Read and initialize parameters
    printf("\n>>Parameters being initialized..");
    params = readParameters(argc, argv);
    if (!params){
    printf("\nMissing parameters! See below..\n");
	showHelp();
        return -1;
    }
    printf("\npepXML file: %s", params->pepxml_file);
    printf("\nmzXML file: %s", params->mzxml_file);
    printf("\nOutput: %s", params->output_mzXML_file);
    printf("\nmmme: %g", params->mmme);
    printf("\nScore: %s", params->score_name);
    printf("\nMS scan interval: [%i , %i]", params->ms_start_scan, params->ms_end_scan);
    printf("\nCrop flag: %i", params->crop);
    printf("\nMinimum # of calibrants: %i", params->min_cal);
    printf("\nBackground intensity: %g", params->bg);
    printf("\nRetention time window: [-%g , +%g]", params->lower_rel_bnd_rt, params->upper_rel_bnd_rt);
    printf("\nRetention time offset: %g", params->recal_offset);
    printf("\nScore interval: [%g , %g]", params->min_score_threshold, params->max_score_threshold);
    fflush(stdout);
    printf("\ndone.\n" ); fflush(stdout);
    
    if(params->ms_start_scan > params ->ms_end_scan ){
        printf("\nStart scan cannot be after end scan, program exiting..."); fflush(stdout); 
        return 0;
    }
    
    //Get peptide prophet results
    dlgl = make_delegate_list();
    dlg = make_peptide_prophet_result_delegate();
    add_delegate_to_list(dlgl, dlg);
    
    //Read pepXML file
    printf("\n>>Reading pepXML file %s...", params->pepxml_file); fflush(stdout);
    pepXML_file = read_pepxml_file(params->pepxml_file, 0, 0, dlgl);
    if (pepXML_file == NULL) {
	printf("\nerror opening pepXML file\n"); fflush(stdout);
	return -1;
    }
    printf("\n%i MS/MS run(s)\n",pepXML_file->run_summary_count); fflush(stdout);
    printf("done.\n" ); fflush(stdout);
   
    //Filtering out peptides
    printf("\n>>Filtering peptides... "); fflush(stdout);
    peptide_set = build_peptide_set(pepXML_file, params, &pepnum);
    printf("\n%i peptides that are potential calibrants.\n", pepnum); fflush(stdout);	
    printf("done.\n" ); fflush(stdout); //CH1

    //Reading mzXML file
    printf("\n>>Reading mzXML file %s...",params->mzxml_file); fflush(stdout);
    mzXML_file = read_mzxml_file_spectrum(params->mzxml_file, 1, 0, &params->ms_start_scan, &params->ms_end_scan);
    if (mzXML_file == NULL) {
	printf("error opening mzXML file\n"); fflush(stdout);
	return -1;
    }

    printf("\nManufacturer: %s", mzXML_file->msinstrument_array->mm_value); fflush(stdout);
    printf("\nModel: %s", mzXML_file->msinstrument_array->mod_value); fflush(stdout);
    printf("\nMass analyzer: %s", mzXML_file->msinstrument_array->ma_value); fflush(stdout);
    printf("\nScan count: %i", mzXML_file->ms_scan_count); fflush(stdout);
    printf("\n%i MS scans to be calibrated within the scan window: [%i, %i]\n", params->ms_end_scan - params->ms_start_scan + 1 , mzXML_file->scan_id_array[params->ms_start_scan-1], mzXML_file->scan_id_array[params->ms_end_scan-1]); fflush(stdout);
    printf("\ndone.\n" ); fflush(stdout); //CH2

    printf("\n>>Making internal calibrant candidate list for each scan within the retention time window...\n"); fflush(stdout);

    build_internal_calibrants(mzXML_file, peptide_set, pepnum, params);

    //check the print within this function
    printf("done.\n"); fflush(stdout); //CH3

    printf("\n>>Constructing the final calibrant list for each scan and recalibrating the peaks..."); fflush(stdout);

    //for each scan
    for(scan = params->ms_start_scan; scan <= params->ms_end_scan; scan++) {
        //scan = params->ms_start_scan;
        mzscan = get_scan(mzXML_file, scan, 0);
        printf("\nMS scan: %i", mzXML_file->scan_id_array[scan-1]); fflush(stdout);
        //Check if it's MS1
        if(mzscan->attributes.msLvl != 1 ){
            printf("\tMS2 level - skipping...\n"); fflush(stdout);
            continue;
        } 
        // If crop param is 1, empty invalid scans
	if(!mzscan->peaks || mzscan->peaks->count < 0) {
            printf("\tEmpty scan..."); fflush(stdout);
            if (params->crop){
		empty_scan(mzXML_file, scan);
                printf("Cropped"); fflush(stdout);
            }
            printf("Not cropped"); fflush(stdout);
	    continue;
	}
        
        // Loading and filtering the peaks
        mzpeaks = load_scan_peaks(mzXML_file, scan);
         
        printf("Peak count: %i", mzpeaks.count); fflush(stdout);
        
	if (mzpeaks.count == 0) {
            unload_scan_peaks(mzXML_file, scan);
            //continue;
	} //CH4
        
        // Make the list of calibrants for this spectrum and sort them in descending order of intensity  
        printf("\n\t>>>Making the final calibrant list for the spectrum..."); fflush(stdout);
        makeCalibrantList(scan, &mzpeaks, peptide_set, params);
        
        printf("\n\t>>>Recalibrating peaks..."); fflush(stdout); //CH5       
        // Recalibrate peaks
	SATISFIED = recalibratePeaks(params);
       
        if (SATISFIED) {
            printf("\n\t\t-Calibration successful."); fflush(stdout); 
            //applyCalibration(scan, &mzpeaks);	//CH6
            //update_scan_peaks(mzXML_file, scan, mzpeaks.count, 32, mzpeaks.mzs, mzpeaks.intensities); //CH7
	}// if
	else {
            printf("\n\t\t-Calibration failed."); fflush(stdout); 
            unload_scan_peaks(mzXML_file, scan);
            if (params->crop)
                empty_scan(mzXML_file, scan);						
	}
    }
    
    /*
    printf("done.\n\n>>Writing recalibrated mzXML file to destination %s...", params->output_mzXML_file); fflush(stdout);
    write_mzxml_file(mzXML_file, params->output_mzXML_file);
    */
    printf("\ndone.\n"); fflush(stdout);
    
    
    
    
    return 0;	
    
}

