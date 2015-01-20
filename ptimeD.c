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
	char ext[128];   // extension of new data file
	char ext0[]="D";   // default extension of new data file
	int nstokes;
	int mode = 0;  // default: creat new file ".D"
	int pmode = 0;  // default: use predictor

	int index, n;
	for (i=0;i<argc;i++)
	{
		if (strcmp(argv[i],"-f") == 0)
		{
			index = i + 1;
			n = 0;
			while ( (index + n) < argc && strcmp(argv[index+n],"-e") != 0 && strcmp(argv[index+n],"-np") != 0 )
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
	}

	sub->mode = pmode;
		
	T2Predictor pred;
	int ret;

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
			printf ("%s\n", sub->fname);
		}
		else
		{
			strcpy(inName,argv[k]);
			createNewfile(inName, sub->fname, ext);
			printf ("%s\n", sub->fname);
		}

		// open psrfits file
		fp = openFitsFile(sub->fname);
				
		// read header info
		loadPrimaryHeader(fp,header);
		closeFitsFile(fp);

		////////////////////////////////////////////////////
	
		printf ("DM0: %.4lf\n", header->dm);
	
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

		if (ret=T2Predictor_ReadFits(&pred,sub->fname))
		{
			printf("Error: unable to read predictor\n");
			exit(1);
		}

		////////////////////////////////////////////////////////////////////////////////

		double *p_multi_deDM;
		p_multi_deDM = (double *)malloc(sizeof(double)*nchn*npol*nphase);

		double *p_temp, *p_temp_deDM;
		p_temp = (double *)malloc(sizeof(double)*npol*nphase);
		p_temp_deDM = (double *)malloc(sizeof(double)*npol*nphase);

		double phaseShift;

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

				n = 0;
				for (nstokes = 0; nstokes < npol; nstokes++)
				{
					for (j = 0; j < nphase; j++)
					{
						p_multi_deDM[nstokes*nchn*nphase + i*nphase + j] = p_temp_deDM[j];
						//printf ("%d %lf\n", j, p_temp_deDM[j]);
						n++;
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

	T2Predictor_Destroy(&pred);

	return 0;
}
