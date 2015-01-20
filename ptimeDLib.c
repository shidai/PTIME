#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"
#include "ptimeDLib.h"

int write_prof (subintegration *sub, pheader *header, double *profile)
{  
  fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
  int status;
  int colnum;

  status = 0;

  if ( fits_open_file(&fptr, sub->fname, READWRITE, &status) )          // open the file
  {
		printf( "error while openning file\n" );
  }

	fits_movnam_hdu(fptr, BINARY_TBL, (char *)"SUBINT",0,&status);
  
  if ( fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status) )           // get the row number
  {
		printf( "error while getting the colnum number\n" );
	}

	///////////////////////////////////////////////////////////////////////////

  int frow;
  int felem;
  int nelem;
	int null;
  //int anynull;

  frow = sub->indexSub;
  felem = 1;
  nelem = header->nbin*header->nchan*header->npol;
  null = 0;

  //fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nelem, &null, profile, &anynull, &status);           // read the column
	fits_write_col(fptr, TDOUBLE, colnum, frow, felem, nelem, profile, &status);           // read the column

  if ( fits_close_file(fptr, &status) )
  {
      printf("error while closing the file\n");
  }

  return 0;
}

int modify_freq (subintegration *sub, pheader *header)
//int main (int argc, char *argv[] )
{  
  fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
  int status;
  int colnum;

  status = 0;

  if ( fits_open_file(&fptr, sub->fname, READWRITE, &status) )          // open the file
  //if ( fits_open_file(&fptr, argv[1], READONLY, &status) )          // open the file
  {
      printf( "error while openning file\n" );
  }

	fits_movnam_hdu(fptr, BINARY_TBL, (char *)"SUBINT",0,&status);

  if ( fits_get_colnum(fptr, CASEINSEN, "DAT_FREQ", &colnum, &status) )           // get the colnum number
  {
      printf( "error while getting the colnum number\n" );
	//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
  //printf ("%d\n", colnum);

  int frow;
  int	felem = 1;
  int nelem = header->nchan;
  //int	anynull = 0;

	//int subint = 1;
	//int nchan = 8;
	//double freq[nchan];
	frow = sub->indexSub;

	double freq[header->nchan];
	int i;
	for (i = 0; i < header->nchan; i++)
	{
		freq[i] = header->freq;
	}

  //fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nelem, &null, freq, &anynull, &status);           // read the column
	fits_write_col(fptr, TDOUBLE, colnum, frow, felem, nelem, freq, &status);           // read the column

  if ( fits_close_file(fptr, &status) )
  {
      printf( " error while closing the file \n" );
  }

  return 0;
}

int createNewfile (char *input, char *output, char *ext)
{
	int i, n;
	int nchar = strlen(input);

	for (i = 0; i < nchar-1; i++)
	{
		if (input[i] == '.')
			n = i;
	}

	char temp[n+2];
	for (i = 0; i < n+1; i++)
	{
		temp[i] = input[i];
	}
	temp[n+1] = '\0';

	//strcat(temp,ext);
	strcpy(output,temp);
	strcat(output,ext);

	/////////////////////////////////////////////////////////////////
  fitsfile *fptrNew, *fptr;       // pointer to the FITS file, defined in fitsio.h 
	int status = 0;

  if ( fits_create_file(&fptrNew, output, &status) )  // open the file
	{
		printf( "error while creating file\n" );
	}

  if ( fits_open_file(&fptr, input, READONLY, &status) ) // open the file
	{
		printf( "error while openning file\n" );
	}

	fits_copy_file(fptr, fptrNew, 1, 1, 1, &status);

  if ( fits_close_file(fptrNew, &status) )
  {
      printf( " error while closing the file \n" );
  }

  if ( fits_close_file(fptr, &status) )
  {
      printf( " error while closing the file \n" );
  }

	return 0;
}

double phaseShiftDM (subintegration *sub, pheader *header, T2Predictor pred)
{
	double dm = header->dm;
	double freq = sub->freq[sub->indexChn];
	double psrFreq = 1.0/sub->Cperiod;
	double freqRef = header->freq;

	long double mjd0;
	mjd0 = (long double)(header->imjd) + ((long double)(header->smjd) + (long double)(header->stt_offs) + (long double)(sub->offs))/86400.0L;

	double phase;
	double phaseShift;

	if (sub->mode == 1)
	{
		phase = (K*dm*psrFreq)*(1.0/(freq*freq)-1.0/(freqRef*freqRef));
		phaseShift = -(2.0*M_PI)*(phase - floor(phase));
		//phaseShift = (K*dm*psrFreq)*(1.0/(freq*freq)-1.0/(freqRef*freqRef));
	
		//phaseShift = -(2.0*M_PI)*(K*dm*psrFreq)*(1.0/(freq*freq)-1.0/(freqRef*freqRef));
		//printf ("%lf %lf\n", freq, phaseShift);
	}
	else 
	{
		phase = T2Predictor_GetPhase(&pred, mjd0, freq);
		phaseShift = (2.0*M_PI)*(phase - floor(phase));
		//phaseShift = -T2Predictor_GetPhase(&pred, mjd, freq);
	}

	//printf ("Predictor: %lf %Lf %lf\n", freq, mjd0, phase);

	return phaseShift;
}

int deDM (int nphase, int npol, double *in, double phaseShift, double *out)
// de-disperse 
{
	int i, j;
	
	double I_in[nphase], Q_in[nphase], U_in[nphase], V_in[nphase];
	double I_out[nphase], Q_out[nphase], U_out[nphase], V_out[nphase];
	
	for (i = 0; i < npol; i++)
	{
		for (j = 0; j < nphase; j++)
		{
			if (i == 0)
			{
				I_in[j] = in[i*nphase+j];
			}
			else if (i == 1)
			{
				Q_in[j] = in[i*nphase+j];
			}
			else if (i == 2)
			{
				U_in[j] = in[i*nphase+j];
			}
			else if (i == 3)
			{
				V_in[j] = in[i*nphase+j];
			}
			else
			{
				printf ("Wrong npol!\n");
			}
		}
	}

	double I_in_real[nphase/2], I_in_ima[nphase/2];
	double Q_in_real[nphase/2], Q_in_ima[nphase/2];
	double U_in_real[nphase/2], U_in_ima[nphase/2];
	double V_in_real[nphase/2], V_in_ima[nphase/2];

	preA7_QUV (I_in, nphase, I_in_real, I_in_ima);

	if (npol == 4)
	{
		preA7_QUV (Q_in, nphase, Q_in_real, Q_in_ima);
  	preA7_QUV (U_in, nphase, U_in_real, U_in_ima);
  	preA7_QUV (V_in, nphase, V_in_real, V_in_ima);
	}

	double I_out_real[nphase/2], I_out_ima[nphase/2];
	double Q_out_real[nphase/2], Q_out_ima[nphase/2];
	double U_out_real[nphase/2], U_out_ima[nphase/2];
	double V_out_real[nphase/2], V_out_ima[nphase/2];
	rotate (nphase, I_in_real, I_out_real, I_in_ima, I_out_ima, phaseShift);
	if (npol == 4)
	{
		rotate (nphase, Q_in_real, Q_out_real, Q_in_ima, Q_out_ima, phaseShift);
		rotate (nphase, U_in_real, U_out_real, U_in_ima, U_out_ima, phaseShift);
		rotate (nphase, V_in_real, V_out_real, V_in_ima, V_out_ima, phaseShift);
	}

	inverse_dft (I_out_real, I_out_ima, nphase, I_out);
	if (npol == 4)
	{
		inverse_dft (Q_out_real, Q_out_ima, nphase, Q_out);
		inverse_dft (U_out_real, U_out_ima, nphase, U_out);
		inverse_dft (V_out_real, V_out_ima, nphase, V_out);
	}

	for (i = 0; i < npol; i++)
	{
		for (j = 0; j < nphase; j++)
		{
			if (i == 0)
			{
				out[i*nphase+j] = I_out[j];
				//printf ("%d %lf %lf\n", j, I_in[j], I_out[j]);
			}
			else if (i == 1)
			{
				out[i*nphase+j] = Q_out[j];
			}
			else if (i == 2)
			{
				out[i*nphase+j] = U_out[j];
			}
			else if (i == 3)
			{
				out[i*nphase+j] = V_out[j];
			}
			else
			{
				printf ("Wrong npol!\n");
			}
		}
	}

	return 0;
}

