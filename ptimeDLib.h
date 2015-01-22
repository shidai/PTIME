#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <fftw3.h>
#include "fitsio.h"
#include "ptimeTLib.h"
#include "ptimeLib.h"
#include "tempo2pred.h"
//#include "tempo2pred_int.h"
#include <gsl/gsl_multimin.h>

int deDM (int nphase, int npol, double *in, long double phaseShift, double *out);
int write_prof (subintegration *sub, pheader *header, double *profile);
long double phaseShiftDM (subintegration *sub, pheader *header, T2Predictor pred);
int modify_freq (subintegration *sub, pheader *header);
int createNewfile (char *input, char *output, char *ext);
