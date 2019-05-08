// calculate the phase shift between template and simulated (or real) data 
// error of phase can be calculated
// initial guess of phase shift is added
// try to do two-dim template matching
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "ptimeDLib.h"
#include "T2toolkit.h"
//#include "tempo2pred.h"

int main (int argc, char *argv[])
{
	int h,i,j,k;

	pheader *header;
	header = (pheader *)malloc(sizeof(pheader));

	fitsfile *fp;

	subintegration *sub;
	sub = (subintegration *)malloc(sizeof(subintegration));

	//////////////////////////////////////////////////////
	char inName[128];   // name of input data file
	char predName[128];   // name of tempo2 predictor file
	char ext[128];   // extension of new data file
	char ext0[]="D";   // default extension of new data file
	int nstokes;
	int mode = 0;  // default: creat new file ".D"
	int pmode = 0;  // default: use predictor; 1: don't use any predictor; 2: use a specified predictor

	int index, n;
	for (i=0;i<argc;i++)
	{
		if (strcmp(argv[i],"-f") == 0)
		{
			index = i + 1;
			n = 0;
			while ( (index + n) < argc && strcmp(argv[index+n],"-e") != 0 && strcmp(argv[index+n],"-np") != 0 && strcmp(argv[index+n],"-pred") != 0)
			{
				n++;
			}
		}
		else if (strcmp(argv[i],"-e") == 0)
		{
			strcpy(ext,argv[++i]);
			mode = 1;  // creat new file with new extension
		}
		else if (strcmp(argv[i],"-np") == 0)  // not use predictor
		{
			pmode = 1;  
		}
		else if (strcmp(argv[i],"-pred") == 0)  // read predictor from a specified file
		{
			strcpy(predName,argv[++i]);
			pmode = 2;  
		}
	}

	sub->mode = pmode;
		
	T2Predictor pred;

	T2Predictor_Init(&pred);  // prepare the predictor
				
	/////////////////////////////////////////////////////////////////////////////////
	// start to deal with different data file
	for (k = index; k < index + n; k++)
	{
		// get the data file name
		if (mode == 0)
		{
			strcpy(inName,argv[k]);
			createNewfile(inName, sub->fname, ext0);
			//printf ("%s\n", sub->fname);
		}
		else
		{
			strcpy(inName,argv[k]);
			createNewfile(inName, sub->fname, ext);
			//printf ("%s\n", sub->fname);
		}

		// open psrfits file
		fp = openFitsFile(sub->fname);
				
		// read header info
		loadPrimaryHeader(fp,header);
		closeFitsFile(fp);

		////////////////////////////////////////////////////
	
		//printf ("DM0: %.4lf\n", header->dm);
	
		////////////////////////////////////////////////
		int nphase;
		int nchn;
		int nsub;
		int npol;
	
		nchn   = header->nchan;
		npol   = header->npol; 
		nsub   = header->nsub; 
		nphase = header->nbin; 

		//printf ("%d\n", nchn);
		////////////////////////////////////////////////

    initialiseSub(sub, header);

		if (sub->mode == 0)
		{
			printf("Reading predictor from the header of %s\n", sub->fname);
			if (T2Predictor_ReadFits(&pred,sub->fname))
			{
				printf("Error: unable to read predictor\n");
				exit(1);
			}
		}
		else if (sub->mode == 2)
		{
			printf("Reading predictor from %s\n", predName);
			if (T2Predictor_Read(&pred, predName))
			{
				printf("Error: unable to read predictor\n");
				exit(1);
			}
		}

		////////////////////////////////////////////////////////////////////////////////

		double *p_multi_deDM;
		p_multi_deDM = (double *)malloc(sizeof(double)*nchn*npol*nphase);

		double *p_temp, *p_temp_deDM;
		p_temp = (double *)malloc(sizeof(double)*npol*nphase);
		p_temp_deDM = (double *)malloc(sizeof(double)*npol*nphase);

		long double phaseShift;

		int n;
		// start to derive toa from different subint
		for (h = 1; h <= nsub; h++)
		{
			sub->indexSub = h;
			      
			// read profiles from data file
			read_prof(sub, header);

			// start to derive toas for different channels
			for (i = 0; i < nchn; i++)
			{
				sub->indexChn = i;
				//printf ("Chan%d\n", sub->indexChn);

				n = 0;
				for (nstokes = 0; nstokes < npol; nstokes++)
				{
					for (j = 0; j < nphase; j++)
					{
						p_temp[n] = sub->p_multi[nstokes*nchn*nphase + i*nphase + j];
						//printf ("%d %lf\n", n, p_temp[n]);
						n++;
					}
				}

				// dedisperse
				phaseShift = phaseShiftDM (sub,header,pred);
				deDM (nphase, npol, p_temp, phaseShift, p_temp_deDM);

				for (nstokes = 0; nstokes < npol; nstokes++)
				{
					for (j = 0; j < nphase; j++)
					{
						p_multi_deDM[nstokes*nchn*nphase + i*nphase + j] = p_temp_deDM[j];
						//printf ("%d %lf\n", j, p_temp_deDM[j]);
					}
				}
			}
			write_prof (sub, header, p_multi_deDM);
			//modify_freq (sub, header);
		}

		free(p_multi_deDM);
		free(p_temp);
		free(p_temp_deDM);
	}

	free(header);
	demallocSub(sub, header);

	T2Predictor_Destroy(&pred);

	return 0;
}

//long double phaseShiftDM (subintegration *sub, pheader *header, T2Predictor pred)
//{
//	//double dm = header->dm;
//	double dm = sub->dmFit;
//	double freq = sub->freq[sub->indexChn];
//	double psrFreq = 1.0/sub->Cperiod;
//	double freqRef = header->freq;
//
//	long double mjd0;
//	//mjd0 = (long double)(header->imjd) + ((long double)(header->smjd) + (long double)(header->stt_offs))/86400.0L;
//	mjd0 = (long double)(header->imjd) + ((long double)(header->smjd) + (long double)(header->stt_offs) + (long double)(sub->offs))/86400.0L;
//
//	long double phase;
//	long double phaseShift;
//
//	phase = T2Predictor_GetPhase(&pred, mjd0, freq);
//	phaseShift = (2.0*M_PI)*(phase - floor(phase));
//
//	return phaseShift;
//}
//
//int deDM (int nphase, double *in, long double phaseShift, double *out)
//// de-disperse 
//{
//	int i, j;
//	
//	double I_in[nphase];
//	double I_out[nphase];
//	
//	for (j = 0; j < nphase; j++)
//	{
//		I_in[j] = in[j];
//	}
//
//	double I_in_real[nphase/2], I_in_ima[nphase/2];
//
//	preA7_QUV (I_in, nphase, I_in_real, I_in_ima);
//
//	double I_out_real[nphase/2], I_out_ima[nphase/2];
//
//	rotate (nphase, I_in_real, I_out_real, I_in_ima, I_out_ima, phaseShift);
//	inverse_dft (I_out_real, I_out_ima, nphase, I_out);
//
//	for (j = 0; j < nphase; j++)
//	{
//		out[i*nphase+j] = I_out[j];
//	}
//
//	return 0;
//}

