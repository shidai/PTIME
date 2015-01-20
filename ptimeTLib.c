// ptimeT library
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>

#include "ptimeTLib.h"
#include "ptimeLib.h"
#include "T2toolkit.h"
#include "tempo2pred.h"
#define ITMAX 100000  // Maximum allowed number of iterations.
#define EPS 1.0e-16 // Machine double floating-point precision.
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

int print_t2pred ( char *name )
{  
	fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
	int status;
  int colnum;
  long int nrows;

  status = 0;

  if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
  {
      printf( "error while openning file\n" );
  }

	fits_movnam_hdu(fptr, BINARY_TBL, (char *)"T2PREDICT",0,&status);

  if ( fits_get_num_rows(fptr, &nrows, &status) )           // get the row number
  {
      printf( "error while getting the row number\n" );
  }
  //printf ("%ld\n", nrows);
  
  //if ( fits_get_colnum(fptr, CASEINSEN, "TSUBINT", &colnum, &status) )           // get the row number
  if ( fits_get_colnum(fptr, CASEINSEN, "PREDICT", &colnum, &status) )           // get the row number
  {
		printf( "error while getting the colnum number\n" );
	}

	int frow;
  int	felem = 1;
  int nelem = 1;
  int	anynull = 0;
  char nval[]="NULL";

	char **line;
	line = (char **)malloc(sizeof(char *));
	line[0] = (char *)malloc(sizeof(char)*1024);

	// print predictor to t2pred.dat
	FILE *fp;
	char filename[]="t2pred.dat";
		    
	if ((fp = fopen(filename, "w+")) == NULL)
	{
		fprintf (stdout, "Can't open file\n");
		exit(1);
	}

	int i;
	for (i = 1; i <= nrows; i++)                             // print the results
	{
		frow = i;
		fits_read_col(fptr, TSTRING, colnum, frow, felem, nelem, nval, line, &anynull, &status);           // read the column
		fprintf (fp, "%s\n", line[0]);
	}

	if (fclose (fp) != 0)
		fprintf (stderr, "Error closing\n");


  if ( fits_close_file(fptr, &status) )
  {
		printf( " error while closing the file \n" );
  }

  return 0;
}

int read_prof (subintegration *sub, pheader *header)
{  
  fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
  int status;
  int colnum;

  status = 0;

	// open psrfits
  if ( fits_open_file(&fptr, sub->fname, READONLY, &status) )          // open the file
  {
      printf( "error while openning file\n" );
  }

	// move to subint
	fits_movnam_hdu(fptr, BINARY_TBL, (char *)"SUBINT",0,&status);

	int nbin;
  int frow;
  int felem;
  int nelem;
  int null;
  int anynull;

	/////////////////////////////////////////////////////////////////////////
	// read profile
  if ( fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status) )           // get the row number
  {
      printf( "error while getting the colnum number\n" );
	}

	nbin = header->nbin;
	frow = sub->indexSub;
	felem = 1;
	nelem = header->nbin*header->nchan*header->npol;
	null = 0;
	anynull = 0;

	fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nelem, &null, sub->p_multi, &anynull, &status);           // read the column

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// read weights
 
	if ( fits_get_colnum(fptr, CASEINSEN, "DAT_WTS", &colnum, &status) )           // get the colnum number
	{
		printf( "error while getting the colnum number\n" );
	}

	frow = sub->indexSub;
  felem = 1;
  nelem = header->nchan;
  null = 0;
  anynull = 0;

	fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nelem, &null, sub->wts, &anynull, &status);           // read the column

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// read channel frequency
 
	if ( fits_get_colnum(fptr, CASEINSEN, "DAT_FREQ", &colnum, &status) )           // get the colnum number
	{
		printf( "error while getting the colnum number\n" );
	}

	frow = sub->indexSub;
  felem = 1;
  nelem = header->nchan;
  null = 0;
  anynull = 0;

	fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nelem, &null, sub->freq, &anynull, &status);           // read the column

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// read offs
  
	if ( fits_get_colnum(fptr, CASEINSEN, "OFFS_SUB", &colnum, &status) )           // get the colnum number
	{
		printf( "error while getting the colnum number\n" );
	}

	frow = sub->indexSub;
  felem = 1;
  nelem = 1;
  null = 0;
  anynull = 0;

	fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nelem, &null, &sub->offs, &anynull, &status);           // read the column

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// close psrfits file
	if ( fits_close_file(fptr, &status) )
	{
		printf( " error while closing the file " );
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////

	double weight, frequency;
	frequency = 0.0;
	weight = 0.0;
	
	int z;
	for (z = 0; z < header->nchan; z++)
	{
		frequency += sub->freq[z]*sub->wts[z];
		weight += sub->wts[z];
	}
	frequency = frequency/weight;
	sub->frequency = frequency;
	
	// get the pulse period of this subintegration
	long double mjd0;  // the mjd of each subint
	T2Predictor pred;
	int ret;
	double period;

	// get the period
	//print_t2pred(sub->fname);   // output t2pred.dat
	T2Predictor_Init(&pred);  // prepare the predictor
	
	if (ret=T2Predictor_ReadFits(&pred,sub->fname))
	{
		printf("Error: unable to read predictor\n");
		exit(1);
	}

	// get the period at mjd0
	mjd0 = (long double)(header->imjd) + ((long double)(header->smjd) + (long double)(header->stt_offs) + (long double)(sub->offs))/86400.0L;
	//printf ("imjd is %ld \n", imjd);
	//printf ("mjd0 is %.15Lf \n", mjd0);
	sub->Cperiod = 1.0/T2Predictor_GetFrequency(&pred,mjd0,sub->frequency);

	for (z = 0; z < header->nchan; z++)
	{
		sub->period[z] = 1.0/T2Predictor_GetFrequency(&pred,mjd0,sub->freq[z]);
	}

	T2Predictor_Destroy(&pred);

	return 0;
}

int read_std (subintegration *sub, pheader *header)
{  
	int mode = sub->mode;
	int nphase = header->nbin;
	int nchn = header->nchan;

	int i,j;
	// currently, npol should be 1
	// nchn is the number of channel of the data profile
	if (mode == 0)
	{
		fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
		int status;
		int colnum;
		long int nrows;

		status = 0;

		if ( fits_open_file(&fptr, sub->tname, READONLY, &status) )          // open the file
		{
			printf( "error while openning file\n" );
		}

		fits_movnam_hdu(fptr, BINARY_TBL, (char *)"SUBINT",0,&status);

		if ( fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status) )           // get the colnum number
		{
			printf( "error while getting the colnum number\n" );
			//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
		}
		//printf ("%d\n", colnum);

		//////////////////////////////////////////////////////////////////////////
		int npol;
		if ( fits_read_key(fptr, TINT, (char *)"NPOL", &npol, NULL, &status) )           // get the pol number
		{
			printf( "error while getting the npol number\n" );
			//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
		}
		//printf ("%d\n", npol);

		int nchan;
		if ( fits_read_key(fptr, TINT, (char *)"NCHAN", &nchan, NULL, &status) )           // get the chan number
		{
			printf( "error while getting the npol number\n" );
			//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
		}
		//printf ("%d\n", nchan);
		///////////////////////////////////////////////////////////////////////////

		// read the data
		int nbin;
		int frow;
		int felem;
		int nelem;
		int null;
		int anynull;

		//double *temp;     // the array to store the data   
		//temp = ( double *)malloc( (nchan*npol*nbin) * sizeof( double ) );               // allocate space for column value

		if ( nbin != nphase )
		{
			printf ("The phase bin number of template != the phase bin number of data\n");
			exit (0);
		}
		else if (nchan != nchn && nchan != 1)
		{
			printf ("Can't not do template matching! The channel number of template should be 1 or equal to that of data.\n");
			exit (0);
		}
		else
		{
			nbin = nphase;
			frow = 1;
			felem = 1;
			nelem = nbin*nchan*npol;
			//nelem = 1024;
			null = 0;
			anynull = 0;

			//printf ("%d\n", nbin);
			fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nelem, &null, sub->s_multi, &anynull, &status);           // read the column
			//fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nelem, &null, temp, &anynull, &status);           // read the column

			if ( fits_close_file(fptr, &status) )
			{
				printf( " error while closing the file\n" );
			}

			if ( nchan == 1 )
			{
				printf ("The channel number of template != the channel number of data\n");
				for (i = 1; i < nchn; i++)
				{
					for (j = 0; j < nphase; j++)
					{
						sub->s_multi[i*nphase+j] = sub->s_multi[j];
					}
				}
			}
		}
	}
	else if (mode == 1)
	{
		tmplStruct tmpl;
		initialiseTemplate(&tmpl);
		//printf("Reading template\n");
		readTemplate(sub->tname,&tmpl);
	    //printf("Complete reading template\n");

		int i,j,k;
		double phi;
		if ( tmpl.nchan != nchn && tmpl.nchan != 1)
		{
			printf ("Can't not do template matching! The channel number of template should be 1 or equal to that of data.\n");
			exit (0);
		}
		else
		{
			for (k = 0; k < tmpl.channel[0].nstokes; k++) // must be fixed; every channel has to have the same number of stokes
			{
				for (i = 0; i < tmpl.nchan; i++)
				{
					for (j = 0; j < nphase; j++)
					{
						phi = j/(double)nphase;
						sub->s_multi[k*tmpl.nchan*nphase+i*nphase+j] = (double)evaluateTemplateChannel(&tmpl,phi,i,k,0);
					}
				}
			}

			for (k = 0; k < tmpl.channel[0].nstokes; k++)
			{
				if ( tmpl.nchan == 1 )
				{
					printf ("The channel number of template != the channel number of data\n");
					for (i = 1; i < nchn; i++)
					{
						for (j = 0; j < nphase; j++)
						{
							sub->s_multi[k*tmpl.nchan*nphase+i*nphase+j] = sub->s_multi[k*tmpl.nchan*nphase+j];
						}
					}
				}
			}
		}
	}

	return 0;
}

double A7 (double phase, params param)
//double A7 (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn)
// calculate function A7 in Taylor 92
//double A7 (int n, double *amp_s, double *amp_p, double *phi_s, double *phi_p, double phase)
{
	double A7=0;
	int i,j;

	for (i = 0; i < param.nchn; i++)
	{
		for (j = 0; j < param.num; j++)
	    {
			A7+=(j+1)*param.a_s[i][j]*param.a_p[i][j]*sin(param.p_s[i][j]-param.p_p[i][j]+(j+1)*phase);
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
	}
	
	return A7;
}

double A7_multi (double phase, params param)
//double A7_multi (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms)
// calculate function A7 in Taylor 92, for multi-frequency channel
//double A7 (int n, double *amp_s, double *amp_p, double *phi_s, double *phi_p, double phase)
{
	double A7=0;
	int i,j;
	double c1, c2, s ;

	for (i = 0; i < param.nchn; i++)
	{
		c1 = 0.0;
		c2 = 0.0;
		s = 0.0;
		for (j = 0; j < param.num; j++)
		{
			c1 += (j+1)*param.a_s[i][j]*param.a_p[i][j]*sin(param.p_s[i][j]-param.p_p[i][j]+(j+1)*phase);
			c2 += param.a_s[i][j]*param.a_p[i][j]*cos(param.p_s[i][j]-param.p_p[i][j]+(j+1)*phase);
			s += param.a_s[i][j]*param.a_s[i][j];
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
		A7 += (c1*c2)/(s*param.rms[i]*param.rms[i]);
	}
	
	return A7;
}

double A9 (double phase, params param)
//double A9 (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn)
// calculate function A9 in Taylor 92
{
	double A9=0.0, sum=0.0;
	int i,j;

	for (i = 0; i < param.nchn; i++)
	{
	    for (j = 0; j < param.num; j++)
	    {
		    A9+=param.a_s[i][j]*param.a_p[i][j]*cos(param.p_s[i][j]-param.p_p[i][j]+(j+1)*phase);
		    sum+=param.a_s[i][j]*param.a_s[i][j];
		    //printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
	}
	
	A9=A9/sum;

	return A9;
}

int A9_multi (double phase, params param, double *b)
//int A9_multi (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *b)
{
// calculate function A9 in Taylor 92, for multi-frequency channel
	double A9, sum;
	int i,j;

	for (i = 0; i < param.nchn; i++)
	{
		A9 = 0.0;
		sum = 0.0;
	  for (j = 0; j < param.num; j++)
	  {
		  A9+=param.a_s[i][j]*param.a_p[i][j]*cos(param.p_s[i][j]-param.p_p[i][j]+(j+1)*phase);
		  sum+=param.a_s[i][j]*param.a_s[i][j];
		  //printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
		b[i] = A9/sum;
	}
	
	return 0;
}

int dft_profiles (int N, double *in, fftw_complex *out)
// dft of profiles
{
	//  dft of profiles 
	///////////////////////////////////////////////////////////////////////
	
	//printf ("%lf\n", in[0]);
	int i;
	double inUse[N];
	//fftw_complex *out;
	fftw_plan p;
	
	//in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	//out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	//p = fftw_plan_dft_r2c_1d(N, inUse, out, FFTW_ESTIMATE);
	p = fftw_plan_dft_r2c_1d(N, inUse, out, FFTW_MEASURE);

	for (i = 0; i < N; i++)
	{
		inUse[i] = in[i];
	}

	fftw_execute(p);

	fftw_destroy_plan(p);
	//fftw_free(in); 
	//fftw_free(out);
  
	return 0;
}

int error (double phase, double b, double *errphase, double *errb, params param)
//int error (double phase, double b, double *errphase, double *errb, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn)
// calculate the errors of phase, a and b according to Talyor 1992  
{
	double rms,gk,s1,s2;
	int i,j,n;

	gk=0.0;
	s1=0.0;
	s2=0.0;
	n=0;

	for (i = 0; i < param.nchn; i++)
	{
	    for (j = 0; j < param.num; j++)
	    {
		    //gk+=a_p[i]*a_p[i]+b*b*a_s[i]*a_s[i]-2.0*b*a_s[i]*a_p[i]*cos(p_p[i]-p_s[i]-(i+1)*phase)+a*a*1024.0*1024.0-2.0*a*1024.0*a_p[i]*cos(p_p[i])+2.0*a*b*1024.0*a_s[i]*cos(p_s[i]+(i+1)*phase);
		    gk+=param.a_p[i][j]*param.a_p[i][j]+b*b*param.a_s[i][j]*param.a_s[i][j]-2.0*b*param.a_s[i][j]*param.a_p[i][j]*cos(param.p_p[i][j]-param.p_s[i][j]-(j+1)*phase);
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		s1+=(j+1)*(j+1)*param.a_p[i][j]*param.a_s[i][j]*cos(param.p_p[i][j]-param.p_s[i][j]-(j+1)*phase);
		s2+=param.a_s[i][j]*param.a_s[i][j];
		n++;
		}
	}
	
	rms=sqrt(gk/n);

	(*errphase)=rms/sqrt(2.0*fabs(b*s1));
	(*errb)=rms/sqrt(2.0*fabs(s2));

	return 0;
}

int error_multi (double phase, double *errphase, params param)
//int error_multi (double phase, double *errphase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms)
// calculate the errors of phase, a and b according to Talyor 1992  
{
	double s1;
	int i,j;

	s1=0.0;

	double c1, c2, c3, s;

	for (i = 0; i < param.nchn; i++)
	{
		c1 = 0.0;
		c2 = 0.0;
		c3 = 0.0;
		s = 0.0;
		for (j = 0; j < param.num; j++)
		{
			c1 += (j+1)*param.a_s[i][j]*param.a_p[i][j]*sin(param.p_s[i][j]-param.p_p[i][j]+(j+1)*phase);
			c2 += param.a_s[i][j]*param.a_p[i][j]*cos(param.p_s[i][j]-param.p_p[i][j]+(j+1)*phase);
			c3 += (j+1)*(j+1)*param.a_s[i][j]*param.a_p[i][j]*cos(param.p_s[i][j]-param.p_p[i][j]+(j+1)*phase);
			s += param.a_s[i][j]*param.a_s[i][j];
			//s2+=(a_s[i][j]*a_s[i][j])/(rms[i]*rms[i]);
		}
		s1 += (c2*c3-c1*c1)/(s*param.rms[i]*param.rms[i]);
	}
	
	(*errphase)=1.0/sqrt(2.0*fabs(s1));
	//(*errb)=1.0/sqrt(2.0*s2);

	return 0;
}

double zbrent(double (*func)(double phase, params param), double x1, double x2, double tol, params param)
//	Using Brent’s method, find the root of a function func known to lie between x1 and x2. The root, returned as zbrent, will be refined until its accuracy is tol.
{
	int iter;
	double a=x1,b=x2,c=x2,d,e,min1,min2;
	double fa=(*func)(a, param),fb=(*func)(b, param),fc,p,q,r,s,tol1,xm;
	//double fa=(*func)(a, param),fb=(*func)(b, a_s, a_p, p_s, p_p, num, nchn),fc,p,q,r,s,tol1,xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
		puts ("Root must be bracketed in zbrent\n");

	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) 
	{
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) 
		{
			c=a;   // Rename a, b, c and adjust bounding interval d.
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) 
		{
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}

		tol1=2.0*EPS*fabs(b)+0.5*tol;   // Convergence check.
		xm=0.5*(c-b);

		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) 
		{
			s=fb/fa;  // Attempt inverse quadratic interpolation.

			if (a == c) 
			{
				p=2.0*xm*s;
				q=1.0-s;
			} 
			else 
			{
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;  // Check whether in bounds.

			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);

			if (2.0*p < (min1 < min2 ? min1 : min2)) 
			{
				e=d;  // Accept interpolation.
				d=p/q;
			} 
			else 
			{
				d=xm; // Interpolation failed, use bisection.
				e=d;
			}
		} 
		else  // Bounds decreasing too slowly, use bisection.
		{
			d=xm;
			e=d;
		}
		a=b;  //  Move last best guess to a.
		fa=fb;
		if (fabs(d) > tol1)     //  Evaluate new trial root.
			b += d;
		else
			b += SIGN(tol1,xm);

		fb=(*func)(b, param);
		//fb=(*func)(b, a_s, a_p, p_s, p_p, num, nchn);
	}

	puts ("Maximum number of iterations exceeded in zbrent\n");

	return 0.0;
}

double zbrent_multi(double (*func)(double phase, params param), double x1, double x2, double tol, params param)
//double zbrent_multi(double (*func)(double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms), double x1, double x2, double tol, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms)
{
	int iter;
	double a=x1,b=x2,c=x2,d,e,min1,min2;
	double fa=(*func)(a, param),fb=(*func)(b, param),fc,p,q,r,s,tol1,xm;
	//double fa=(*func)(a, a_s, a_p, p_s, p_p, num, nchn, rms),fb=(*func)(b, a_s, a_p, p_s, p_p, num, nchn, rms),fc,p,q,r,s,tol1,xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
		puts ("Root must be bracketed in zbrent\n");

	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) 
	{
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) 
		{
			c=a;   // Rename a, b, c and adjust bounding interval d.
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) 
		{
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}

		tol1=2.0*EPS*fabs(b)+0.5*tol;   // Convergence check.
		xm=0.5*(c-b);

		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) 
		{
			s=fb/fa;  // Attempt inverse quadratic interpolation.

			if (a == c) 
			{
				p=2.0*xm*s;
				q=1.0-s;
			} 
			else 
			{
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;  // Check whether in bounds.

			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);

			if (2.0*p < (min1 < min2 ? min1 : min2)) 
			{
				e=d;  // Accept interpolation.
				d=p/q;
			} 
			else 
			{
				d=xm; // Interpolation failed, use bisection.
				e=d;
			}
		} 
		else  // Bounds decreasing too slowly, use bisection.
		{
			d=xm;
			e=d;
		}
		a=b;  //  Move last best guess to a.
		fa=fb;
		if (fabs(d) > tol1)     //  Evaluate new trial root.
			b += d;
		else
			b += SIGN(tol1,xm);

		fb=(*func)(b, param);
		//fb=(*func)(b, a_s, a_p, p_s, p_p, num, nchn, rms);
	}

	puts ("Maximum number of iterations exceeded in zbrent\n");

	return 0.0;
}

int get_toa (double *s, double *p, subintegration *sub, pheader *header)
// calculate the phase shift between template and simulated (or real) data 
// error of phase can be calculated
// initial guess of phase shift is added
// calculate the phase shift between template and simulated (or real) data 
// error of phase can be calculated
// initial guess of phase shift is added
// try to do two-dim template matching
{
	//int nphase=1024;
	int nchn=1;

	int nphase = header->nbin;
	int k;  // k=nphase/2

	//double amp_s[nchn][nphase/2],amp_p[nchn][nphase/2];  // elements for calculating A7
	//double phi_s[nchn][nphase/2],phi_p[nchn][nphase/2];
	params param;
	allocateMemory (&param, nchn, nphase);
	//double amp_s[nchn][NP],amp_p[nchn][NP];  // elements for calculating A7
	//double phi_s[nchn][NP],phi_p[nchn][NP];  // the second dim should be NP, which is large enough for different observations

	preA7(s, p, nphase, nchn, &param);
	//preA7(&k, amp_s, amp_p, phi_s, phi_p, s, p, nphase, nchn);
	//printf ("%d\n", nchn);
	
	// initial guess of the phase
	int peak_s, peak_p;	

	find_peak(0, nphase,s,&peak_s);
	find_peak(0, nphase,p,&peak_p);

	int d, chn;
	double step;
	double ini_phase,up_phase,low_phase;

	d = InitialGuess (s, p, nphase, 1, &chn);
	//d=peak_p-peak_s;
	//printf ("Initial guess: %d\n",d);
	step=2.0*M_PI/(10.0*nphase);

	if (d>=nphase/2)
	{
		if (d>0)
		{
			ini_phase=2.0*M_PI*(nphase-1-d)/nphase;
		}
		else
		{
			ini_phase=-2.0*M_PI*(nphase-1+d)/nphase;
		}
		up_phase=ini_phase+step;
		low_phase=ini_phase-step;
		while (A7(up_phase, param)*A7(low_phase, param)>0.0)
		{
		    up_phase+=step;
		    low_phase-=step;
		}
	}
	else
	{
		ini_phase=-2.0*M_PI*d/nphase;
		up_phase=ini_phase+step;
		low_phase=ini_phase-step;
		while (A7(up_phase, param)*A7(low_phase, param)>0.0)
		{
		    up_phase+=step;
		    low_phase-=step;
		}
	}

  // calculate phase shift, a and b
  double phase,b;
  phase=zbrent(A7, low_phase, up_phase, 1.0e-16, param);
  //phase=zbrent(A7, -1.0, 1.0, 1.0e-16);
  //phase=zbrent(A7, -0.005, 0.005, 1.0e-16);
  b=A9(phase, param);
  //a=A4(b);
	//(*bx) = b;
	sub->b[sub->indexChn] = b;

		
	//printf ("Phase shift: %.10lf\n", phase/(2.0*M_PI));
	//printf ("%.10lf %.10lf\n", phase, A7(phase));
	//printf ("%.10lf \n", ((phase/3.1415926)*5.75/2.0)*1.0e+3);  // microseconds
	//printf ("%.10lf \n", b);
	//printf ("%.10lf \n", a);
	//printf ("///////////////////////// \n");
		
	
	// calculate the errors of phase and b
	double errphase, errb;	

	error(phase,b,&errphase,&errb, param);
	//printf ("%.10lf %.10lf\n", ((phase/3.1415926)/(psrfreq*2.0))*1.0e+6, ((errphase/3.1415926)/(psrfreq*2.0))*1.0e+6);  // microseconds
	//printf ("%.10lf %.10lf\n", ((phase/3.1415926)*4.569651/2.0)*1.0e+3, ((errphase/3.1415926)*4.569651/2.0)*1.0e+3);  // microseconds
	//printf ("errphase %.10lf \n", ((errphase/3.1415926)*5.75/2.0)*1.0e+6);
	//printf ("errb %.10lf \n", errb);
	
	//(*phasex) = phase;
	//(*errphasex) = errphase;
	sub->phase = phase;
	sub->e_phase = errphase;

	// calculate the rms
	double rms;
	cal_rms(phase, b, &rms, param);
	sub->rms[sub->indexChn] = rms;

	deallocateMemory (&param, nchn);
	return 0;
}

int get_toa_multi (subintegration *sub, pheader *header)
{
	int nchn = header->nchan;
	int nphase = header->nbin;
	
	int k;  // k=nphase/2

	params param;
	allocateMemory (&param, nchn, nphase);
	param.rms = sub->rms;

	preA7(sub->s_multi, sub->p_multi, nphase, nchn, &param);
	
	int d, chn;
	double step;
	double ini_phase,up_phase,low_phase;

	d = InitialGuess (sub->s_multi, sub->p_multi, nphase, nchn, &chn);
	//d = InitialGuess (st, pt, nphase, nchn);
	//d=peak_p-peak_s;
	//printf ("Initial guess: %d\n",d);
	step=2.0*M_PI/(10.0*nphase);

	if (d>=nphase/2)
	{
		if (d>0)
		{
			ini_phase=2.0*M_PI*(nphase-1-d)/nphase;
		}
		else
		{
			ini_phase=-2.0*M_PI*(nphase-1+d)/nphase;
		}
		up_phase=ini_phase+step;
		low_phase=ini_phase-step;
		while (A7_multi(up_phase, param)*A7_multi(low_phase, param)>0.0)
		//while (A7_multi(up_phase, amp_s, amp_p, phi_s, phi_p, k, nchn, rms, bx)*A7_multi(low_phase, amp_s, amp_p, phi_s, phi_p, k, nchn, rms,bx)>0.0)
		{
		    up_phase+=step;
		    low_phase-=step;
		}
	}
	else
	{
		ini_phase=-2.0*M_PI*d/nphase;
		up_phase=ini_phase+step;
		low_phase=ini_phase-step;
		while (A7_multi(up_phase, param)*A7_multi(low_phase, param)>0.0)
		//while (A7_multi(up_phase, amp_s, amp_p, phi_s, phi_p, k, nchn, rms, bx)*A7_multi(low_phase, amp_s, amp_p, phi_s, phi_p, k, nchn, rms, bx)>0.0)
		{
		    up_phase+=step;
		    low_phase-=step;
		}
	}
	printf ("initial guess: %lf\n", ini_phase);

  // calculate phase shift, a and b
  double phase;
  phase=zbrent_multi(A7_multi, low_phase, up_phase, 1.0e-16, param);

	double b[nchn];
	A9_multi (phase, param, b);
		
	// calculate the errors of phase and b
  double errphase;	

	error_multi(phase, &errphase, param);
	//error_multi(phase, &errphase, amp_s, amp_p, phi_s, phi_p, k, nchn, rms, bx);
	//printf ("multi-template\n");
	//printf ("Phase shift: %.10lf+-%.10lf\n", ((phase/M_PI)/(psrfreq*2.0))*1.0e+6, ((errphase/M_PI)/(psrfreq*2.0))*1.0e+6);  // microseconds
	//printf ("%.10lf %.10lf\n", ((phase/3.1415926)*4.569651/2.0)*1.0e+3, ((errphase/3.1415926)*4.569651/2.0)*1.0e+3);  // microseconds
	//printf ("errphase %.10lf \n", ((errphase/3.1415926)*5.75/2.0)*1.0e+6);
	//printf ("errb %.10lf \n", errb);
	
	sub->phase = phase;
	sub->e_phase = errphase;

	deallocateMemory (&param, nchn);

	return 0;
}

int getToaMultiDM (subintegration *sub, pheader *header)
//int get_toa_multi (double *s, double *p, double *rms, double *bx, int nchn, double *phasex, double *errphasex, double psrfreq, int nphase)
{
	int nchn = header->nchan;
	int nphase = header->nbin;

	int k;  // k=nphase/2

	params param;
	allocateMemory (&param, nchn, nphase);
	param.nfreq = sub->freq;
	param.rms = sub->rms;
	param.psrFreq = 1.0/sub->Cperiod;
	param.dm = header->dm;
	param.freqRef = header->freq;
	//param.freqRef = 1260.0;
	printf ("obsFreq: %lf\n", param.freqRef);

	preA7(sub->s_multi, sub->p_multi, nphase, nchn, &param);
	
	int d, chn;
	double step;
	double ini_phase,up_phase,low_phase;

	d = InitialGuess (sub->s_multi, sub->p_multi, nphase, nchn, &chn);
	//printf ("chn: %d\n", chn);

	//int peak_s, peak_p;	
	//find_peak(chn, nphase,s,&peak_s);
	//find_peak(chn, nphase,p,&peak_p);
	//printf ("peak_s: %d\n", peak_s);
	//printf ("peak_p: %d\n", peak_p);
	//d=peak_p-peak_s;

	//printf ("Initial guess: %d\n",d);
	step=2.0*M_PI/(10.0*nphase);

	//printf ("initial guess: %lf\n", (double)(d)/nphase);
	//printf ("delay: %lf\n", nphase*(K*param.dm*param.psrFreq)*(1.0/(param.nfreq[chn]*param.nfreq[chn])-1.0/(param.freqRef*param.freqRef))/(2.0*3.1415926));
	//d = d - (int)(nphase*(K*param.dm*param.psrFreq)*(1.0/(param.nfreq[chn]*param.nfreq[chn])-1.0/(param.freqRef*param.freqRef)));
	//d = d - (int)(d/nphase)*nphase;
	//printf ("Initial guess: %d\n",d);
	if (fabs(d)>=nphase/2)
	{
		if (d>0)
		{
			ini_phase=2.0*M_PI*(nphase-1-d)/nphase;
		}
		else
		{
			ini_phase=-2.0*M_PI*(nphase-1+d)/nphase;
		}
	}
	else
	{
		ini_phase=-2.0*M_PI*d/nphase;
	}

	//printf ("initial guess: %lf\n", ini_phase/(2*M_PI));
  // calculate phase shift, DM, a and b
  double phase, dmFit;
	//printf ("fitDM: Initial guess %f\n",ini_phase);
	miniseNelderMead (&param, ini_phase, &phase, &dmFit);
	//miniseD (&param, ini_phase, &phase, &dmFit);

	// calculate the errors of phase and DM
  double errphase, errDm;	
	covariance (&param, phase, dmFit, &errphase, &errDm);

	// test
	//miniseNelderMeadTest (&param, ini_phase, &phase, &dmFit);
	//errphase = 0.01;
	//errDm = 0.01;

	printf ("multi-template\n");
	//printf ("Phase shift: %.10lf+-%.10lf\n", phase, errphase);  // microseconds
	//printf ("Phase shift: %.10lf+-%.10lf\n", ((phase/3.1415926)/(psrfreq*2.0))*1.0e+6, ((errphase/3.1415926)/(psrfreq*2.0))*1.0e+6);  // microseconds
	printf ("DM: %.10lf   %.10lf\n", dmFit, errDm);
	//printf ("%.10lf %.10lf\n", ((phase/3.1415926)*4.569651/2.0)*1.0e+3, ((errphase/3.1415926)*4.569651/2.0)*1.0e+3);  // microseconds
	//printf ("errphase %.10lf \n", ((errphase/3.1415926)*5.75/2.0)*1.0e+6);
	//printf ("errb %.10lf \n", errb);
	
	sub->phase = phase;
	sub->e_phase = errphase;

	deallocateMemory (&param, nchn);

	return 0;
}

int preA7 (double *s, double *p, int nphase, int nchn, params *param)
//int preA7 (int *k, double amp_s[][NP], double amp_p[][NP], double phi_s[][NP], double phi_p[][NP], double *s, double *p, int nphase, int nchn)
// preparation for calculating A7 of Talyor 1992  
{
	// nphase is the dimention of one profile, nchn is number of profiles
	// k is the dimention of amp of one profile 
	param->nchn = nchn;
	int i,j;
	
	/////////////////////////////////////////////////////////////////////////////////

	fftw_complex *out_s;
	fftw_complex *out_p;
	
	out_s = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nphase);
	out_p = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nphase);
	
	double s_temp[nphase];  // store one template and profile
	double p_temp[nphase];  

	int n;
	double r_s[nphase/2],im_s[nphase/2];
	double r_p[nphase/2],im_p[nphase/2];
	for (i = 0; i < nchn; i++)
	{
	    for (j=0;j<nphase;j++)
	    {
		    s_temp[j]=s[i*nphase + j];
		    p_temp[j]=p[i*nphase + j];
	    }

	    dft_profiles(nphase,s_temp,out_s);
	    //printf ("%lf %lf\n", out_s[1][0], out_s[1][1]);

	    dft_profiles(nphase,p_temp,out_p);

	    //double amp_s[N/2],phi_s[N/2];
	    //double amp_p[N/2],phi_p[N/2];

		n = 0;
	    for (j = 0; j <= nphase/2-1; j++)
	    {
		    r_s[j]=out_s[j+1][0];
		    im_s[j]=out_s[j+1][1];
		    r_p[j]=out_p[j+1][0];
		    im_p[j]=out_p[j+1][1];
		    //printf ("%lf %lf\n", r_p[i], im_p[i]);
		    //printf ("%lf %lf\n", out_s[i][0], out_s[i][1]);
		    n++;
	    }
	    //printf ("%d\n", n);
	    //printf ("%d %d\n", nphase, nchn);

	    for (j = 0; j < n; j++)
	    {
		    param->a_s[i][j]=sqrt(r_s[j]*r_s[j]+im_s[j]*im_s[j]);
		    param->a_p[i][j]=sqrt(r_p[j]*r_p[j]+im_p[j]*im_p[j]);
		    param->p_s[i][j]=atan2(im_s[j],r_s[j]);
		    param->p_p[i][j]=atan2(im_p[j],r_p[j]);
		    //printf ("%lf %lf %lf\n", r_s[i], im_s[i], amp_s[i]);
		    //printf ("%lf %lf %lf\n", r_p[i], im_p[i], amp_p[i]);
		    //printf ("%lf\n", amp_s[i]);
		    //printf ("%lf\n", amp_p[i]);
	    }
	}
	//(*k)=n;
	param->num = n;

	fftw_free(out_s); 
	fftw_free(out_p); 
	//fftw_free(out_t); 

	return 0;
}

int cal_rms (double phase, double b, double *rms, params param)
// calculate the rms of each subchannel  
{
	double gk;
	int i,j,n;

	gk=0.0;
	n=0;

	for (i = 0; i < param.nchn; i++)
	{
	    for (j = 0; j < param.num; j++)
	    {
		    //gk+=a_p[i]*a_p[i]+b*b*a_s[i]*a_s[i]-2.0*b*a_s[i]*a_p[i]*cos(p_p[i]-p_s[i]-(i+1)*phase)+a*a*1024.0*1024.0-2.0*a*1024.0*a_p[i]*cos(p_p[i])+2.0*a*b*1024.0*a_s[i]*cos(p_s[i]+(i+1)*phase);
		    gk+=param.a_p[i][j]*param.a_p[i][j]+b*b*param.a_s[i][j]*param.a_s[i][j]-2.0*b*param.a_s[i][j]*param.a_p[i][j]*cos(param.p_p[i][j]-param.p_s[i][j]-(j+1)*phase);
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		n++;
		}
	}
	
	(*rms)=sqrt(gk/n);

	return 0;
}

int form_toa_multi (subintegration *sub, pheader *header)
{
	long double dt;  
	// transform phase shift to TOAs

	// transform phase shift to time shift
	//dt = (phase/M_PI)*period/2.0;
	//e_dt = (e_phase/M_PI)*period/2.0;
	dt = ((long double)(sub->phase)/M_PI)*((long double)(sub->Cperiod))/2.0L;
	sub->e_dt = ((long double)(sub->e_phase)/M_PI)*((long double)(sub->Cperiod))/2.0L;
	//printf ("dt is %.10Lf +/- %.10Lf\n", dt, e_dt);

	// calculate the TOA
  sub->t = (long double)(header->imjd) + ((long double)(header->smjd) + (long double)(header->stt_offs) - (long double)(dt) + (long double)(sub->offs))/86400.0L;
		
	//printf ("offset is %lf\n", offset);
	//fprintf (fp, "%s  %lf  %.15Lf  %Lf  7\n", fname, frequency, t, e_dt*1e+6);

	return 0;
}

int form_toa (subintegration *sub, pheader *header)
// chn is the channel to form toa
// nchn is the total number of subchn
{
	long double dt;  
	
	// transform phase shift to time shift
	//dt = (phase/M_PI)*period/2.0;
	//e_dt = (e_phase/M_PI)*period/2.0;
	dt = ((long double)(sub->phase)/M_PI)*((long double)(sub->period[sub->indexChn]))/2.0L;
	sub->e_dt = ((long double)(sub->e_phase)/M_PI)*((long double)(sub->period[sub->indexChn]))/2.0L;

	// calculate the TOA
	sub->t = (long double)(header->imjd) + ((long double)(header->smjd) + (long double)(header->stt_offs) - (long double)(dt) + (long double)(sub->offs))/86400.0L;
		
	return 0;
}

int find_peak (int n0, int n, double *s, int *position)
{
	int i;
	double temp[n];
	double peak;

	for (i = n*n0; i < n*n0+n; i++)
	{
		temp[i-n*n0] = s[i];
	}

	double a, b, c;
	for (i = 0; i < n-1; i++)
	{
		//a = fabs(temp[i]);
		//b = fabs(temp[i+1]);
		a = temp[i];
		b = temp[i+1];
		c = (a >= b ? a : b);

		temp[i+1] = c;
	}
	peak = temp[n-1];

	for (i = n*n0; i < n*n0+n; i++)
	{
		if (fabs(peak-s[i]) < 1.0e-3)
		{
			(*position) = i;
		}
	}

	return 0;
}

double find_peak_value (int n, double *s)
{
	int i;
	double temp[n];

	for (i = 0; i < n; i++)
	{
		temp[i] = s[i];
	}

	double a, b, c;
	for (i = 0; i < n-1; i++)
	{
		a = temp[i];
		b = temp[i+1];
		//a = fabs(temp[i]);
		//b = fabs(temp[i+1]);
		//c = (fabs(a) >= fabs(b) ? a : b);
		c = (a >= b ? a : b);

		temp[i+1] = c;
	}

	return temp[n-1];
}

int corr (double *s, double *p, int nphase)
{
	/*
	FILE *fp1, *fp2;

	if ((fp1 = fopen(argv[1], "r")) == NULL)
	{
		fprintf (stdout, "Can't open file\n");
		exit(1);
	}

	if ((fp2 = fopen(argv[2], "r")) == NULL)
	{
		fprintf (stdout, "Can't open file\n");
		exit(1);
	}

	float x1[1024], x2[1024];
	int i = 0;
	while (fscanf (fp1, "%f", &x1[i]) == 1)
	{
		i++;
	}

	i = 0;
	while (fscanf (fp2, "%f", &x2[i]) == 1)
	{
		i++;
	}
	*/

	double r[nphase];
	int i, j;
	for (j = 0; j < nphase; j++)
	{
		r[j] = 0.0;
		for (i = 0; i < nphase; i++)
		{
			if ((i+j) > (nphase-1))
			{
				r[j] += p[i]*s[i+j-(nphase-1)];
			}
			else
			{
				r[j] += p[i]*s[i+j];
			}
			//printf ("%f %f\n", x1[i], x2[i]);
		}
	}

	int shift;
	find_peak (0, nphase, r,  &shift);
	/*
	for (j = 0; j < 1024; j++)
	{
		printf ("%f\n", r[j]);
	}
	*/

	return -shift;
}

int def_off_pulse (int nphase, double *in, double frac_off)
// define the off pulse region based on I, return the starting index of off pulse region
// using frac_off to calculate the off pulse region
{
	int n = nphase;
	int num_off = (int)(n*frac_off);
	int i,j;
	double small;
	double temp;
	int index = 0;

	for (i = 0; i < n; i++)
	{
		if (i == 0)
		{
			small = 0.0;
			for(j = 0; j < num_off; j++)
			{
				small += (in[j]+30000.0)*(in[j]+30000.0);  // make all numbers positive
			}
			small = sqrt(small/num_off);
		}
			
		temp = 0.0;
		for(j = 0; j < num_off; j++)
		{
			if ((i+j) > n-1)
			{
				temp += (in[(i+j)-(n-1)]+30000.0)*(in[(i+j)-(n-1)]+30000.0);
			}
			else 
			{
				temp += (in[i+j]+30000.0)*(in[i+j]+30000.0);
			}
		}
		temp = sqrt(temp/num_off);

		small = (temp <= small ? temp : small);
		index = (temp <= small ? i : index);
		//printf ("%d %lf %lf\n", index, small, ave);
	}

	return index;
}

int off_pulse (int nphase, int index, double *in, double *out, double frac_off)
// get the off_pulse region
{
	int n = nphase;
	int num_off = (int)(n*frac_off);
	int i;

	for (i = 0; i < num_off; i++)
	{
		if ((index+i) > n-1)
		{
			out[i] = in[(index+i)-(n-1)];
		}
		else 
		{
			out[i] = in[index+i];
		}
	}

	return 0;
}

int remove_baseline (double *in, int index, double frac_off, int n, double *out)
{
	// define the off_pulse range, frac_off is the fraction of the phase
	// index is the starting point of the off_pulse range
	int num_off = (int)(n*frac_off);
	double off_0[num_off];

	off_pulse (n, index, in, off_0, frac_off);

	int i;
	double baseline = 0.0;
    for (i = 0; i < num_off; i++)
    {
        baseline += off_0[i];
        //average_s += s_off[i];
    }
	baseline = baseline/num_off;

    //printf ("the baseline of std is: %lf \n", baseline);
    //printf ("average is: %lf %lf\n", average, average_s);

	for (i = 0; i < n; i++)
	{
		out[i] = (in[i]-baseline);
		//s_norm[i] = (s[i]-baseline)/(s_peak-baseline);
	}
	
	return 0;
}

int pre_diff (double *s, int nphase, int index, double frac_off, double *s_out)
{
	int n = nphase;
	
	// remove the baseline
	remove_baseline (s, index, frac_off, n, s_out);

    return 0;
}


int InitialGuess (double *s, double *p, int nphase, int nchn, int *chn) 
{
	int index;
	double frac_off = 0.05;  // set to be 0.05

	double ptemp[nphase];
	double temp[nchn];
	int i, h;
	for (i = 0; i < nchn; i++)
	{
		for (h = 0; h < nphase; h++)
		{
			ptemp[h] = p[i*nphase+h];
		}
		int x;
		x = def_off_pulse (nphase, ptemp, frac_off);

		double ptemp_out[nphase];
		pre_diff (ptemp, nphase, x, frac_off, ptemp_out);

		temp[i] = find_peak_value(nphase,ptemp_out);
	}

	int peak;
	find_peak (0, nchn, temp, &peak);
	//printf ("%d\n",peak);
	(*chn) = peak;

	double p_use[nphase];
	double s_use[nphase];
	for (h = 0; h < nphase; h++)
	{
		p_use[h] = p[peak*nphase+h];
		s_use[h] = s[peak*nphase+h];
	}

	// remove the baseline of template
	index = def_off_pulse (nphase, s_use, frac_off);

	double s_out[nphase];
	pre_diff (s_use, nphase, index, frac_off, s_out);

	// remove the baseline of profile
	index = def_off_pulse (nphase, p_use, frac_off);

	double p_out[nphase];
	pre_diff (p_use, nphase, index, frac_off, p_out);

	// Guess the phase shift
	int d;
	d = corr (s_out, p_out, nphase);

	return d;
}

int preA7_QUV (double *p, int nphase, double *real_p, double *ima_p)
// preparation for calculating A7 of Talyor 1992  
{
	// nphase is the dimention of one profile, nchn is number of profiles
	// k is the dimention of amp of one profile 
	int i,j;
	
	/////////////////////////////////////////////////////////////////////////////////

	fftw_complex *out_p;
	
	out_p = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nphase);
	
	double p_temp[nphase];  // store one template and profile

	for (j=0;j<nphase;j++)
	{
	  p_temp[j]=p[i*nphase + j];
	}

	dft_profiles(nphase,p_temp,out_p);

	//double amp_s[N/2],phi_s[N/2];
	//double amp_p[N/2],phi_p[N/2];

	for (j = 0; j < nphase/2+1; j++)                                                  
	{                                                                      
		real_p[j]=out_p[j][0];                                             
		ima_p[j]=out_p[j][1];                                              
	}
										
	fftw_free(out_p); 
	//fftw_free(out_t); 

	return 0;
}

int rotate (int N, double *real_p, double *real_p_rotate, double *ima_p, double *ima_p_rotate, double rot)
{
	// k is the dimention of amp, N is the dimention of s
	int i;

	// for substraction 
	double amp,cosina,sina;
	for (i=0;i<N/2+1;i++)
	{
		// calculate the sin(phi) and cos(phi) of the profile
		amp=sqrt(real_p[i]*real_p[i]+ima_p[i]*ima_p[i]);
		cosina=real_p[i]/amp;
		sina=ima_p[i]/amp;

		// rotate profile
		real_p_rotate[i]=amp*(cosina*cos(-i*rot*M_PI)-sina*sin(-i*rot*M_PI));
		ima_p_rotate[i]=amp*(sina*cos(-i*rot*M_PI)+cosina*sin(-i*rot*M_PI));
		//real_p_rotate[i]=amp*(cosina*cos(-i*M_PI)-sina*sin(-i*M_PI));
		//ima_p_rotate[i]=amp*(sina*cos(-i*M_PI)+cosina*sin(-i*M_PI));
		
	}

	return 0;
}

int align (int N, double phase, double b, double a, double *real_p, double *real_p_align, double *ima_p, double *ima_p_align, double rotate)
{
	// k is the dimention of amp, N is the dimention of s
	int i;

	// for substraction 
	double amp,cosina,sina;
	for (i=0;i<N/2+1;i++)
	{
		// calculate the sin(phi) and cos(phi) of the profile
		amp=sqrt(real_p[i]*real_p[i]+ima_p[i]*ima_p[i]);
		cosina=real_p[i]/amp;
		sina=ima_p[i]/amp;

		// add phase shift to the profile, phase
		//real_p_align[i]=amp*(cosina)/b;
		//ima_p_align[i]=amp*(sina)/b;
		//real_p_align[i]=amp*(cosina*cos(-i*phase)-sina*sin(-i*phase));
		//ima_p_align[i]=amp*(sina*cos(-i*phase)+cosina*sin(-i*phase));
		//real_p_align[i]=amp*(cosina*cos(-i*phase)-sina*sin(-i*phase))/b;
		//ima_p_align[i]=amp*(sina*cos(-i*phase)+cosina*sin(-i*phase))/b;
		real_p_align[i]=(amp*(cosina*cos(-i*(phase+rotate*M_PI))-sina*sin(-i*(phase+rotate*M_PI))))/b;
		ima_p_align[i]=(amp*(sina*cos(-i*(phase+rotate*M_PI))+cosina*sin(-i*(phase+rotate*M_PI))))/b;
		
	}

	return 0;
}

int inverse_dft (double *real_p, double *ima_p, int ncount, double *p_new)
{
	double *dp;
    fftw_plan plan;
	fftw_complex *cp;

    dp = (double *)malloc(sizeof (double) * ncount);
	cp = (fftw_complex *)fftw_malloc(sizeof (fftw_complex) * ncount);
	memset(dp, 0, sizeof (double) * ncount);
	memset(cp, 0, sizeof (fftw_complex) * ncount);

	// initialize the dft...
	double *dp_t;
    fftw_plan plan_t;
	fftw_complex *cp_t;

    dp_t = (double *)malloc(sizeof (double) * ncount);
	cp_t = (fftw_complex *)fftw_malloc(sizeof (fftw_complex) * ncount);
	memset(dp_t, 0, sizeof (double) * ncount);
	memset(cp_t, 0, sizeof (fftw_complex) * ncount);

	int i;
    double real,ima,amp,cosina,sina;

	for (i = 0; i < ncount; i++)
	{
		if (i < ncount/2+1)
		{
            real = real_p[i];
            ima = ima_p[i];
			amp = sqrt(real*real+ima*ima);
			cosina = real/amp;
			sina = ima/amp;

			cp[i][0] = amp*(cosina);
			cp[i][1] = amp*(sina);
			//cp[i][0] = amp*(cosina*cos(-i*3.1415926)-sina*sin(-i*3.1415926));
			//cp[i][1] = amp*(sina*cos(-i*3.1415926)+cosina*sin(-i*3.1415926));
			//cp[i][0]=real_s[i]-real_p[i];
			//cp[i][1]=ima_s[i]-ima_p[i];
			//cp[i][0]=-real_s[i]+real_p[i];
			//cp[i][1]=-ima_s[i]+ima_p[i];
			cp_t[i][0] = real_p[i];
			cp_t[i][1] = ima_p[i];
			//cp[i][0]=real_p[i];
			//cp[i][1]=ima_p[i];
		}
		else
		{
			cp[i][0]=0.0;
			cp[i][1]=0.0;
			cp_t[i][0]=0.0;
			cp_t[i][1]=0.0;
		}
	}

    plan_t = fftw_plan_dft_c2r_1d(ncount, cp_t, dp_t, FFTW_MEASURE);

    fftw_execute(plan_t);

    fftw_destroy_plan(plan_t);

	/////////////////////////////////////////////////////////////////

    plan = fftw_plan_dft_c2r_1d(ncount, cp, dp, FFTW_MEASURE);

    fftw_execute(plan);

    fftw_destroy_plan(plan);

	for (i = 0; i < ncount; i++)
	{
		p_new[i] = dp[i]/ncount;  // normalized by the ncount
		//printf ("%lf\n", p_new[i]);
	}

	return 0;
}

int allocateMemory (params *param, int nchn, int nphase)
{
	param->a_s = (double **)malloc(sizeof(double *)*nchn);
	param->a_p = (double **)malloc(sizeof(double *)*nchn);
	param->p_s = (double **)malloc(sizeof(double *)*nchn);
	param->p_p = (double **)malloc(sizeof(double *)*nchn);
	//printf ("test\n");

	int i;
	for (i = 0; i < nchn; i++)
	{
		param->a_s[i] = (double *)malloc(sizeof(double)*nphase);
		param->a_p[i] = (double *)malloc(sizeof(double)*nphase);
		param->p_s[i] = (double *)malloc(sizeof(double)*nphase);
		param->p_p[i] = (double *)malloc(sizeof(double)*nphase);
	}
}

int deallocateMemory (params *param, int nchn)
{
	int i;
	for (i = 0; i < nchn; i++)
	{
		free(param->a_s[i]);
		free(param->a_p[i]);
		free(param->p_s[i]);
		free(param->p_p[i]);
	}
		
	free(param->a_s);
	free(param->a_p);
	free(param->p_s);
	free(param->p_p);
}

// fit DM functions
double my_f (const gsl_vector *v, void *params)
{
	double x, y;
	double *p = (double *)params;
	x = gsl_vector_get(v, 0);
	y = gsl_vector_get(v, 1);
	return p[2] * (x - p[0]) * (x - p[0]) +
		p[3] * (y - p[1]) * (y - p[1]) + p[4];
}

double chiSquare (const gsl_vector *x, void *param)
{
	double phase = gsl_vector_get (x,0);
	double dm = gsl_vector_get (x,1);

	int nchn = ((params *)param)->nchn;
	int num = ((params *)param)->num;
	double psrFreq = ((params *)param)->psrFreq;
	double *nfreq = ((params *)param)->nfreq;
	double *rms = ((params *)param)->rms;
	double **a_s = ((params *)param)->a_s;
	double **a_p = ((params *)param)->a_p;
	double **p_s = ((params *)param)->p_s;
	double **p_p = ((params *)param)->p_p;
	double freqRef;
	freqRef = ((params *)param)->freqRef;

	int i,j;
	double chi2;

	double phaseNchn;

	/*
	double s, c00, nu0;
	double c001, c002, c003;
	for (i = 0; i < nchn; i++)
	{
		s = 0.0;
		c001 = 0.0;
		c002 = 0.0;
		c003 = 0.0;
		//printf ("nchn freq: %lf\n",nfreq[i]);
		//phaseNchn = phase + (K*dm*psrFreq)*(1.0/(nfreq[i]*nfreq[i]));
		phaseNchn = phase + (K*dm*psrFreq)*(1.0/(nfreq[i]*nfreq[i]));
		for (j = 0; j < num; j++)
		{
			s += a_s[i][j]*a_s[i][j];
			c001 += a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
			c002 += (j+1)*a_s[i][j]*a_p[i][j]*sin(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
			c003 += (j+1)*(j+1)*a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
		c00 += 2.0*((-c002*c002+c001*c003)/s)/(rms[i]*rms[i]);
		nu0 += (1.0/(nfreq[i]*nfreq[i]))*((-c002*c002+c001*c003)/s)/(rms[i]*rms[i]);
	}

	freqRef = sqrt(c00/(2.0*nu0));
	*/

	chi2 = 0.0;
	double P, PS, S;
	for (i = 0; i < nchn; i++)
	{
		P = 0.0;
		PS = 0.0;
		S = 0.0;
		//printf ("nchn freq: %lf\n",nfreq[i]);
		phaseNchn = phase - (2.0*M_PI)*(K*dm*psrFreq)*(1.0/(nfreq[i]*nfreq[i])-1.0/(freqRef*freqRef));
		//phaseNchn = phase - (2.0*M_PI)*(K*dm*psrFreq)*(1.0/(nfreq[i]*nfreq[i]));
		for (j = 0; j < num; j++)
		{
			PS += a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
			S += a_s[i][j]*a_s[i][j];
			P += a_p[i][j]*a_p[i][j];
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
		chi2 += (P-PS*PS/S)/(rms[i]*rms[i]);
	}
	
	return chi2;
}

int miniseNelderMead (params *param, double ini_guess, double *phase, double *dmFit)
{
	double psrFreq = ((params *)param)->psrFreq;

	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;

	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;

	gsl_multimin_function minex_func;

	size_t iter = 0;
	int status;
	double size;

	/* Starting point */
	double guess, dmGuess;
	//guess = 0.88;
	guess = ini_guess;
	//dmGuess = param->dm;
	dmGuess = 0.0;
	//printf ("phase guess: %.10lf(%.10lf); DM guess: %.5lf\n", (guess/(3.1415926*2.0)), guess, dmGuess);	
	//printf ("phase guess: %.10lf(%.10lf); DM guess: %.5lf\n", ((guess/3.1415926)/(psrFreq*2.0))*1.0e+6, guess, dmGuess);	

	x = gsl_vector_alloc (2);
	gsl_vector_set (x, 0, guess);
	gsl_vector_set (x, 1, dmGuess);

	// test function
	//gsl_vector_set (x, 0, 5.0);
	//gsl_vector_set (x, 1, 7.0);

	/* Set initial step sizes to 1 */
	
	ss = gsl_vector_alloc (2);
	gsl_vector_set_all (ss, 0.0001);

	/* Initialize method and iterate */
	minex_func.n = 2;
	minex_func.f = chiSquare;
	minex_func.params = param;

	// test function
	//minex_func.f = my_f;
	//double par[5] = {1.0, 2.0, 10.0, 20.0, 30.0};
	//minex_func.params = par;

	s = gsl_multimin_fminimizer_alloc (T, 2);
	printf ("Fit for phase shift and DM.\n");
	gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		if (status)
			break;

		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 1e-6);

		if (status == GSL_SUCCESS)
		{
			printf ("converged to minimum at\n");
		}

		//printf ("%5d %10.3e %10.3e f() = %7.3f size = %.6f\n", iter, gsl_vector_get (s->x, 0)/(3.1415926*2.0), gsl_vector_get (s->x, 1), s->fval, size);
		(*phase) = gsl_vector_get (s->x, 0);
		(*dmFit) = gsl_vector_get (s->x, 1);
	}
	while (status == GSL_CONTINUE && iter < 1000);

	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);

	return status;
}

int covariance (void *param, double phase, double dm, double *errPhase, double *errDm)
{
	int nchn = ((params *)param)->nchn;
	int num = ((params *)param)->num;
	double psrFreq = ((params *)param)->psrFreq;
	double *nfreq = ((params *)param)->nfreq;
	double *rms = ((params *)param)->rms;
	double **a_s = ((params *)param)->a_s;
	double **a_p = ((params *)param)->a_p;
	double **p_s = ((params *)param)->p_s;
	double **p_p = ((params *)param)->p_p;
	double freqRef;

	freqRef = ((params *)param)->freqRef;

	int i,j;

	double s;
	double c00,c11,c01,nu0;
	double c001,c002,c003;
	double c111,c112;
	double c121;
	double phaseNchn;

	/*
	c00 = 0.0;
	nu0 = 0.0;
	for (i = 0; i < nchn; i++)
	{
		s = 0.0;
		c001 = 0.0;
		c002 = 0.0;
		c003 = 0.0;
		//printf ("nchn freq: %lf\n",nfreq[i]);
		//phaseNchn = phase + (K*dm*psrFreq)*(1.0/(nfreq[i]*nfreq[i]));
		phaseNchn = phase + (K*dm*psrFreq)*(1.0/(nfreq[i]*nfreq[i]));
		for (j = 0; j < num; j++)
		{
			s += a_s[i][j]*a_s[i][j];
			c001 += a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
			c002 += (j+1)*a_s[i][j]*a_p[i][j]*sin(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
			c003 += (j+1)*(j+1)*a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
		c00 += 2.0*((-c002*c002+c001*c003)/s)/(rms[i]*rms[i]);
		nu0 += (1.0/(nfreq[i]*nfreq[i]))*((-c002*c002+c001*c003)/s)/(rms[i]*rms[i]);
	}
	freqRef = sqrt(c00/(2.0*nu0));
	*/

	double A = 2.0*M_PI*K*psrFreq;
	//printf ("A: %lf\n", A);
	c00 = 0.0;
	c11 = 0.0;
	c01 = 0.0;
	nu0 = 0.0;
	for (i = 0; i < nchn; i++)
	{
		s = 0.0;
		c001 = 0.0;
		c002 = 0.0;
		c003 = 0.0;
		c111 = 0.0;
		c112 = 0.0;
		c121 = 0.0;
		//printf ("nchn freq: %lf\n",nfreq[i]);
		phaseNchn = phase - (A*dm)*(1.0/(nfreq[i]*nfreq[i])-1.0/(freqRef*freqRef));
		//phaseNchn = phase - (2.0*3.1415926)*(K*dm*psrFreq)*(1.0/(nfreq[i]*nfreq[i])-1.0/(freqRef*freqRef));
		//printf ("phaseNchn: %lf %lf\n", nfreq[i], phaseNchn);
		//printf ("phaseNchn: %lf %.10lf\n", nfreq[i], (1.0/(nfreq[i]*nfreq[i])-1.0/(freqRef*freqRef)));
		for (j = 0; j < num; j++)
		{
			s += a_s[i][j]*a_s[i][j];
			c001 += a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
			c002 += (j+1)*a_s[i][j]*a_p[i][j]*sin(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
			c003 += (j+1)*(j+1)*a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);

			c111 += ((j+1)*A*(1.0/(nfreq[i]*nfreq[i])-1.0/(freqRef*freqRef)))*a_s[i][j]*a_p[i][j]*sin(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
			//c111 += ((j+1)*K*psrFreq*(1.0/(nfreq[i]*nfreq[i])-1.0/(freqRef*freqRef)))*a_s[i][j]*a_p[i][j]*sin(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
			c112 += pow(((j+1)*A*(1.0/(nfreq[i]*nfreq[i])-1.0/(freqRef*freqRef))),2.0)*a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
			//c112 += pow(((j+1)*K*psrFreq*(1.0/(nfreq[i]*nfreq[i])-1.0/(freqRef*freqRef))),2.0)*a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);

			c121 += ((j+1)*(j+1)*A*(1.0/(nfreq[i]*nfreq[i])-1.0/(freqRef*freqRef)))*a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
			//c121 += ((j+1)*(j+1)*K*psrFreq*(1.0/(nfreq[i]*nfreq[i])-1.0/(freqRef*freqRef)))*a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
		c00 += 2.0*((-c002*c002+c001*c003)/s)/(rms[i]*rms[i]);
		c11 += 2.0*((-c111*c111+c001*c112)/s)/(rms[i]*rms[i]);
		c01 += 2.0*((-c111*c002+c001*c121)/s)/(rms[i]*rms[i]);
		//nu0 += (1.0/(nfreq[i]*nfreq[i]))*((-c002*c002+c001*c003)/s)/(rms[i]*rms[i]);
	}

	printf ("c00: %lf; c11: %lf; c12: %lf\n", c00, c11, c01);
	double err0, err1;
	errInvCov (c00, c11, c01, &err0, &err1);
	(*errPhase) = err0;
	(*errDm) = err1;
	//(*errPhase) = sqrt(2.0/fabs(c00));
	//(*errDm) = sqrt(2.0/fabs(c11));

	//printf ("phase error: %lf; DM error: %lf\n",  ((err0/3.1415926)/(psrFreq*2.0))*1.0e+6, err1);
	//printf ("phase error: %lf; DM error: %lf\n",  ((sqrt(2.0/fabs(c00))/3.1415926)/(psrFreq*2.0))*1.0e+6, sqrt(2.0/fabs(c11)));
	//printf ("freqRef: %lf\n", freqRef);
	
	return 0;
}

double chiSquareTest (const gsl_vector *x, void *param)
{
	double phase = gsl_vector_get (x,0);
	double dm = ((params *)param)->dm;
	//double dm = 1.19761;

	int nchn = ((params *)param)->nchn;
	int num = ((params *)param)->num;
	int psrFreq = ((params *)param)->psrFreq;
	double *nfreq = ((params *)param)->nfreq;
	double *rms = ((params *)param)->rms;
	double **a_s = ((params *)param)->a_s;
	double **a_p = ((params *)param)->a_p;
	double **p_s = ((params *)param)->p_s;
	double **p_p = ((params *)param)->p_p;
	double freqRef = ((params *)param)->freqRef;

	int i,j;
	double chi2;

	chi2 = 0.0;
	double P, PS, S;
	double phaseNchn;
	for (i = 0; i < nchn; i++)
	{
		P = 0.0;
		PS = 0.0;
		S = 0.0;
		//phaseNchn = phase;
		//phaseNchn = phase + (K*dm*psrFreq)*(1.0/(nfreq[i]*nfreq[i])-1.0/(nfreq[nchn/2]*nfreq[nchn/2]));
		phaseNchn = phase - (2.0*M_PI)*(K*dm*psrFreq)*(1.0/(nfreq[i]*nfreq[i])-1.0/(freqRef*freqRef));
		for (j = 0; j < num; j++)
		{
			PS += a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
			S += a_s[i][j]*a_s[i][j];
			P += a_p[i][j]*a_p[i][j];
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
		chi2 += (P-PS*PS/S)/(rms[i]*rms[i]);
	}
	
	return chi2;
}

int miniseNelderMeadTest (params *param, double ini_guess, double *phase, double *dmFit)
{
	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;

	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;

	gsl_multimin_function minex_func;

	size_t iter = 0;
	int status;
	double size;

	/* Starting point */
	double guess;
	//guess = 0.490874;
	guess = ini_guess;
	x = gsl_vector_alloc (1);
	gsl_vector_set (x, 0, guess);
	//printf ("initial guess: %lf\n", ((guess/3.1415926)/(param->psrFreq*2.0))*1.0e+6);
	//gsl_vector_set (x, 1, param->dm);

	// test function
	//gsl_vector_set (x, 0, 5.0);
	//gsl_vector_set (x, 1, 7.0);

	/* Set initial step sizes to 1 */
	
	ss = gsl_vector_alloc (1);
	gsl_vector_set_all (ss, 0.01);

	/* Initialize method and iterate */
	minex_func.n = 1;
	minex_func.f = chiSquareTest;
	minex_func.params = param;

	// test function
	//minex_func.f = my_f;
	//double par[5] = {1.0, 2.0, 10.0, 20.0, 30.0};
	//minex_func.params = par;

	s = gsl_multimin_fminimizer_alloc (T, 1);
	printf ("Fit for phase shift and DM (Test).\n");
	gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		if (status)
			break;

		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 1e-6);

		if (status == GSL_SUCCESS)
		{
			printf ("converged to minimum at\n");
		}

		printf ("%5d %10.3e f() = %7.3f size = %.3f\n", iter, gsl_vector_get (s->x, 0), s->fval, size);
		(*phase) = gsl_vector_get (s->x, 0);
		(*dmFit) = 15.9898;
	}
	while (status == GSL_CONTINUE && iter < 100);

	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);

	return status;
}

double chiSquare2 (const gsl_vector *x, void *param)
{
	double phase = gsl_vector_get (x,0);
	double dm = gsl_vector_get (x,1);

	int nchn = ((params *)param)->nchn;
	int num = ((params *)param)->num;
	double psrFreq = ((params *)param)->psrFreq;
	double *nfreq = ((params *)param)->nfreq;
	double *rms = ((params *)param)->rms;
	double **a_s = ((params *)param)->a_s;
	double **a_p = ((params *)param)->a_p;
	double **p_s = ((params *)param)->p_s;
	double **p_p = ((params *)param)->p_p;
	double freqRef;
	freqRef = ((params *)param)->freqRef;

	int i,j;
	double chi2;

	double phaseNchn;

	chi2 = 0.0;
	double P, PS, S;
	for (i = 0; i < nchn; i++)
	{
		P = 0.0;
		PS = 0.0;
		S = 0.0;
		//printf ("nchn freq: %lf\n",nfreq[i]);
		//phaseNchn = phase + (K*dm*psrFreq)*(1.0/(nfreq[i]*nfreq[i]));
		phaseNchn = phase - (2.0*M_PI)*(K*dm*psrFreq)*(1.0/(nfreq[i]*nfreq[i])-1.0/(freqRef*freqRef));
		for (j = 0; j < num; j++)
		{
			PS += a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
			S += a_s[i][j]*a_s[i][j];
			P += a_p[i][j]*a_p[i][j];
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
		chi2 += (P-PS*PS/S)/(rms[i]*rms[i]);
	}
	
	return chi2;
}

void dfChiSquare2 (const gsl_vector *x, void *param, gsl_vector *df)
{
	double phase = gsl_vector_get (x,0);
	double dm = gsl_vector_get (x,1);

	int nchn = ((params *)param)->nchn;
	int num = ((params *)param)->num;
	double psrFreq = ((params *)param)->psrFreq;
	double *nfreq = ((params *)param)->nfreq;
	double *rms = ((params *)param)->rms;
	double **a_s = ((params *)param)->a_s;
	double **a_p = ((params *)param)->a_p;
	double **p_s = ((params *)param)->p_s;
	double **p_p = ((params *)param)->p_p;
	double freqRef;

	freqRef = ((params *)param)->freqRef;

	int i,j;

	double s;
	double c0,c1;
	double c01,c02;
	double c11;
	double phaseNchn;

	c0 = 0.0;
	c1 = 0.0;
	for (i = 0; i < nchn; i++)
	{
		s = 0.0;
		c01 = 0.0;
		c02 = 0.0;
		c11 = 0.0;
		//printf ("nchn freq: %lf\n",nfreq[i]);
		//phaseNchn = phase + (K*dm*psrFreq)*(1.0/(nfreq[i]*nfreq[i]));
		phaseNchn = phase - (2.0*M_PI)*(K*dm*psrFreq)*(1.0/(nfreq[i]*nfreq[i])-1.0/(freqRef*freqRef));
		for (j = 0; j < num; j++)
		{
			s += a_s[i][j]*a_s[i][j];
			c01 += a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
			c02 += (j+1)*a_s[i][j]*a_p[i][j]*sin(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);

			c11 += ((j+1)*K*psrFreq*(1.0/(nfreq[i]*nfreq[i])-1.0/(freqRef*freqRef)))*a_s[i][j]*a_p[i][j]*sin(p_s[i][j]-p_p[i][j]+(j+1)*phaseNchn);
		}
		c0 += 2.0*((c01*c02)/s)/(rms[i]*rms[i]);
		c1 += 2.0*((c01*c11)/s)/(rms[i]*rms[i]);
		//nu0 += (1.0/(nfreq[i]*nfreq[i]))*((-c002*c002+c001*c003)/s)/(rms[i]*rms[i]);
	}

	//printf ("phase error: %lf; DM error: %lf\n",  ((sqrt(2.0/fabs(c00))/3.1415926)/(psrFreq*2.0))*1.0e+6, sqrt(2.0/fabs(c11)));
	//printf ("c00: %lf; c11: %lf; c12: %lf\n", c00, c11, c01);
	//printf ("freqRef: %lf\n", freqRef);
	
	gsl_vector_set(df, 0, c0);
	gsl_vector_set(df, 1, c1);
}

void fdfChiSquare2 (const gsl_vector *x, void *params, double *f, gsl_vector *df)
{
	*f = chiSquare2(x, params);
	dfChiSquare2(x, params, df);
}

int miniseD (params *param, double ini_guess, double *phase, double *dmFit)
{
	double psrFreq = ((params *)param)->psrFreq;

	size_t iter = 0;
	int status;

	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;

	gsl_vector *x;

	gsl_multimin_function_fdf my_func;

	/* Initialize method and iterate */
	my_func.n = 2;
	my_func.f = chiSquare2;
	my_func.df = dfChiSquare2;
	my_func.fdf = fdfChiSquare2;
	my_func.params = param;

	/* Starting point */
	double guess, dmGuess;
	//guess = 0.497010;
	guess = ini_guess;
	dmGuess = param->dm;
	//dmGuess = 0.0;
	//printf ("phase guess: %.10lf(%.10lf); DM guess: %.5lf\n", ((guess/3.1415926)/(psrFreq*2.0))*1.0e+6, guess, dmGuess);	

	x = gsl_vector_alloc (2);
	gsl_vector_set (x, 0, guess);
	gsl_vector_set (x, 1, dmGuess);

	
	T = gsl_multimin_fdfminimizer_conjugate_fr;
	s = gsl_multimin_fdfminimizer_alloc (T, 2);
	gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, 1e-3);

	printf ("Fit for phase shift and DM.\n");

	do
	{
		iter++;
		status = gsl_multimin_fdfminimizer_iterate (s);
		printf ("Status: %d\n", status);
		if (status)
			break;
		status = gsl_multimin_test_gradient (s->gradient, 1e-3);
		if (status == GSL_SUCCESS)
			printf ("Minimum found at:\n");
		printf ("%5d %.5f %.5f %10.5f\n", iter, gsl_vector_get (s->x, 0)/(3.1415926*2.0), gsl_vector_get (s->x, 1), s->f);
		(*phase) = gsl_vector_get (s->x, 0);
		(*dmFit) = gsl_vector_get (s->x, 1);
	}
	while (status == GSL_CONTINUE && iter < 100);

	gsl_multimin_fdfminimizer_free (s);
	gsl_vector_free(x);

	return status;
}

int errInvCov (double c00, double c11, double c01, double *err0, double *err1)
{
	double cov[] = {c00, c01, c01, c11};

	gsl_matrix_view m = gsl_matrix_view_array (cov, 2, 2);

	int s;
	gsl_permutation * p = gsl_permutation_alloc (2);
	gsl_matrix * inverse = gsl_matrix_alloc (2, 2);

	gsl_linalg_LU_decomp (&m.matrix, p, &s);

	gsl_linalg_LU_invert (&m.matrix, p, inverse);

	printf ("Covariance Matrix --> c00: %.15lf; c11: %.15lf\n", gsl_matrix_get(inverse, 0, 0), gsl_matrix_get(inverse, 1, 1));
	(*err0) = sqrt(2*fabs(gsl_matrix_get(inverse, 0, 0)));
	(*err1) = sqrt(2*fabs(gsl_matrix_get(inverse, 1, 1)));

	return 0;
}

void initialiseSub(subintegration *sub, pheader *header)
{
	int nphase;
	int nchn;
	int nsub;
	int npol;
	int mode = sub->mode;
	
	nchn = header->nchan; 
	npol = header->npol; 
	nsub = header->nsub; 
	nphase = header->nbin; 	

	sub->rms = (double *)malloc(sizeof(double)*nchn);
	sub->b = (double *)malloc(sizeof(double)*nchn);

	sub->wts = (double *)malloc(sizeof(double)*nchn);
	sub->freq = (double *)malloc(sizeof(double)*nchn);
	sub->period = (double *)malloc(sizeof(double)*nchn);

	sub->s_multi = (double *)malloc(sizeof(double)*nchn*npol*nphase);
	sub->p_multi = (double *)malloc(sizeof(double)*nchn*npol*nphase);
}

void demallocSub(subintegration *sub, pheader *phead)
{
	free(sub->rms);
	free(sub->b);

	free(sub->wts);
	free(sub->freq);
	free(sub->period);

	free(sub->s_multi);
	free(sub->p_multi);
	free(sub);
}
