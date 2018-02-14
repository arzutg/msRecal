/* 
 * File:   msRecalDefs.h
 * Author: Arzu Tugce Guler
 *
 * Created on 14 April 2016, 14:13
 */

#ifndef MSRECALDEFS_H
#define	MSRECALDEFS_H

#define HPLUS_MASS 1.00727646688
#define DEFAULT_MIN_INTENSITY 3000000		//do not use peaks below MIN_INTENSITY for calibration
#define AMINO_ACIDS "ARNDCEQGHILKMFPSTWYV"
#define ANTI_ACIDS "BJOUXZ"
#define MAX_ROWS 8192
#define MAX_CALIBRANTS 80
#define DEFAULT_MODE 0				//lossless by default (0)
#define DEFAULT_MIN_CALIBRANTS 3		//minimum number of internal calibrants to recalibrate 
#define INTERNAL_CALIBRATION_TARGET 0.5e-6	//discard internal calibrants that do not fit CAL2 better than this */

#endif	/* MSRECALDEFS_H */

