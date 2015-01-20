#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include <fftw3.h>
#include "fitsio.h"
#include "readPfits.h"
//#include "ptime.h"
#include <gsl/gsl_multimin.h>

#define NP 2048
#define K 4149.37759 // 1.0/2.41
//#define K 4148.808

typedef struct subintegration{
	char fname[128];     // name of data file
	char tname[128];     // name of template
	int mode;            // template format

	int indexSub;        // subintegration index
	int indexChn;        // channel index
	double *s_multi;
	double *p_multi;

	double *wts;
	double *freq;
	double *period;
	double Cperiod;
	double frequency;
	double offs;

	double *rms;    // rms for each profile
	double *b;      // b for each profile
	double phase, e_phase;
	long double e_dt;  
	long double t;       // TOA
} subintegration;

typedef struct params {
	int num; 
	int nchn; 
	double freqRef;
	double psrFreq;
	double dm;
	double *nfreq;
	double *rms;
	double **a_s; 
	double **a_p; 
	double **p_s; 
	double **p_p; 
} params;

int print_t2pred ( char *name );
int read_std (subintegration *sub, pheader *header);
int read_prof (subintegration *sub, pheader *header);

int dft_profiles (int N, double *in, fftw_complex *out);
int preA7 (double *s, double *p, int nphase, int nchn, params *param);
double A7 (double phase, params param);
double A7_multi (double phase, params param);
double A9 (double phase, params param);
int A9_multi (double phase, params param, double *b);
double zbrent(double (*func)(double phase, params param), double x1, double x2, double tol, params param);
double zbrent_multi(double (*func)(double phase, params param), double x1, double x2, double tol, params param);
int error (double phase, double b, double *errphase, double *errb, params param);
int error_multi (double phase, double *errphase, params param);
int cal_rms (double phase, double b, double *rms, params param);
int get_toa (double *s, double *p, subintegration *sub, pheader *header);
int get_toa_multi (subintegration *sub, pheader *header);
int form_toa_multi (subintegration *sub, pheader *header);
int form_toa (subintegration *sub, pheader *header);
int find_peak (int n0, int n, double *s, int *position);
double find_peak_value (int n, double *s);
int corr (double *s, double *p, int nphase);
int def_off_pulse (int nphase, double *in, double frac_off);
int off_pulse (int nphase, int index, double *in, double *out, double frac_off);
int remove_baseline (double *in, int index, double frac_off, int n, double *out);
int pre_diff (double *s, int nphase, int index, double frac_off, double *s_out);
int InitialGuess (double *s, double *p, int nphase, int nchn, int *chn);
int preA7_QUV (double *p, int nphase, double *real_p, double *ima_p);
int rotate (int N, double *real_p, double *real_p_rotate, double *ima_p, double *ima_p_rotate, double rot);
int align (int N, double phase, double b, double a, double *real_p, double *real_p_align, double *ima_p, double *ima_p_align, double rotate);
int inverse_dft (double *real_p, double *ima_p, int ncount, double *p_new);
int allocateMemory (params *param, int nchn, int nphase);
int deallocateMemory (params *param, int nchn);

double chiSquare (const gsl_vector *x, void *param);
int miniseNelderMead (params *param, double guess, double *phase, double *dmFit);
int getToaMultiDM (subintegration *sub, pheader *header);
double my_f (const gsl_vector *v, void *params);
int miniseNelderMeadTest (params *param, double guess, double *phase, double *dmFit);
double chiSquareTest (const gsl_vector *x, void *param);
int covariance (void *param, double phase, double dm, double *errPhase, double *errDm);
double chiSquare2 (const gsl_vector *x, void *param);
void dfChiSquare2 (const gsl_vector *x, void *param, gsl_vector *df);
void fdfChiSquare2 (const gsl_vector *x, void *params, double *f, gsl_vector *df);
int miniseD (params *param, double ini_guess, double *phase, double *dmFit);
int errInvCov (double c00, double c11, double c01, double *err0, double *err1);

void initialiseSub(subintegration *sub, pheader *header);
void demallocSub(subintegration *sub, pheader *phead);

