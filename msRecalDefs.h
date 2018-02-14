#ifndef __MSRECALDDEFS_H__
#define __MSRECALDDEFS_H__

#define MAX_PEP_FILES 100
#define HPLUS_MASS 1.00727646688
#define DEFAULT_MIN_INTENSITY 10000			// do not use peaks below MIN_INTENSITY for calibration
#define AMINO_ACIDS "ARNDCEQGHILKMFPSTWYV"
#define ANTI_ACIDS "BJOUXZ"
#define MAX_ROWS 8192
#define MAX_CALIBRANTS 80
#define DEFAULT_MODE 0						// lossless by default
#define DEFAULT_MIN_CALIBRANTS 3			// minimum number of internal calibrants to recalibrate 
#define EXTERNAL_CALIBRATION_TOLERANCE 5e-6 // maximum error allowed to be considered calibrant 
#define INTERNAL_CALIBRATION_TARGET 0.5e-6	// discard internal calibrants that do not fit CAL2 better than this */

#endif
