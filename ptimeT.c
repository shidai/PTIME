// calculate the phase shift between template and simulated (or real) data 
// error of phase can be calculated
// initial guess of phase shift is added
// try to do two-dim template matching
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "ptimeTLib.h"
//#include "T2toolkit.h"
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
	char oname[128];   // name of output .tim
	int tmode;
	int fitDM = 0;         // fit DM or not; fitDM = 0, don't fit; fitDM = 1, fit
	int nstokes;

	double DM;
	int DMstatus = 0;

	int index, n;
	for (i=0;i<argc;i++)
	{
		if (strcmp(argv[i],"-f") == 0)
		{
			index = i + 1;
			n = 0;
			while ( (index + n) < argc && strcmp(argv[index+n],"-std") != 0 && strcmp(argv[index+n],"-pt") != 0 && strcmp(argv[index+n],"-o") != 0 && strcmp(argv[index+n],"-single") != 0 && strcmp(argv[index+n],"-fitDM") != 0 && strcmp(argv[index+n],"-I") != 0 && strcmp(argv[index+n],"-Q") != 0 && strcmp(argv[index+n],"-U") != 0 && strcmp(argv[index+n],"-V") != 0 && strcmp(argv[index+n],"-DM") != 0)
			{
				n++;
			}
		}
		else if (strcmp(argv[i],"-std")==0)
		{
			strcpy(sub->tname,argv[++i]);
			sub->mode = 0; // standard template format
			printf ("standard template format\n");
		}
		else if (strcmp(argv[i],"-pt")==0)
		{
			strcpy(sub->tname,argv[++i]);
			sub->mode = 1; // ptime template
			printf ("ptime template format\n");
		}
		else if (strcmp(argv[i],"-o")==0)
		{
			strcpy(oname,argv[++i]);
		}
		else if (strcmp(argv[i],"-single")==0)
		{
			tmode = 0; // do freq-dependent matching, and get one TOA
		}
		else if (strcmp(argv[i],"-multi")==0)
		{
			tmode = 1; // do freq-dependent matching, and get TOA for each channel
		}
		else if (strcmp(argv[i],"-fitDM")==0)
		{
			fitDM = 1; // do freq-dependent matching, and get TOA for each channel
		}
		else if (strcmp(argv[i],"-I")==0)
		{
			nstokes = 0; // do freq-dependent matching, and get TOA for each channel
		}
		else if (strcmp(argv[i],"-V")==0)
		{
			nstokes = 3; // do freq-dependent matching, and get TOA for each channel
		}
		else if (strcmp(argv[i],"-Q")==0)
		{
			nstokes = 1; // do freq-dependent matching, and get TOA for each channel
		}
		else if (strcmp(argv[i],"-U")==0)
		{
			nstokes = 2; // do freq-dependent matching, and get TOA for each channel
		}
		else if (strcmp(argv[i],"-DM")==0)
		{
			DM = atof(argv[++i]);
			DMstatus = 1; // use input DM to de-dispersion
		}
	}

	/////////////////////////////////////////////////////////////////////////////////
	// open file to write toa 
	FILE *fpt;
	if ((fpt = fopen(oname, "w+")) == NULL)
	{
		fprintf (stdout, "Can't open file\n");
		exit(1);
	}
	
	fprintf (fpt, "FORMAT 1\n");

	/////////////////////////////////////////////////////////////////////////////////
	// start to deal with different data file
	for (k = index; k < index + n; k++)
	{
		// get the data file name
		strcpy(sub->fname,argv[k]);
		printf ("%s\n", sub->fname);

		// open psrfits file
		fp = openFitsFile(sub->fname);

		// read header info
		loadPrimaryHeader(fp,header);
		closeFitsFile(fp);
		////////////////////////////////////////////////////

		if (DMstatus == 1)
		{
			header->dm = DM;
		}
		printf ("DM0: %.4lf\n", header->dm);
	
		////////////////////////////////////////////////
		int nphase;
		int nchn;
		int nsub;
		int npol;
	
		nchn = header->nchan; 
		npol = header->npol; 
		nsub = header->nsub; 
		nphase = header->nbin; 	

		//printf ("%d\n", nchn);
		////////////////////////////////////////////////

		double s_temp[nphase];
		double p_temp[nphase];

		initialiseSub(sub, header);
		read_std (sub, header);

		// start to derive toa from different subint
		for (h = 1; h <= nsub; h++)
		{
			sub->indexSub = h;
			// read profiles from data file
			read_prof(sub, header);

			// get pulse period of this subintegration

			// start to derive toas for different channels
			for (i = 0; i < nchn; i++)
			{
				sub->indexChn = i;
				for (j = 0; j < nphase; j++)
				{
					//printf ("%lf %lf\n", p_multi[j], s[j]);
					//s_multi[i*nphase + j] = s[j];
					if (nstokes == 0)
					{
						p_temp[j] = sub->p_multi[nstokes*nchn*nphase + i*nphase + j];
						s_temp[j] = sub->s_multi[nstokes*nchn*nphase + i*nphase + j];
					}
					else if (nstokes == 1)
					{
						p_temp[j] = sub->p_multi[nstokes*nchn*nphase + i*nphase + j];
						s_temp[j] = sub->s_multi[nstokes*nchn*nphase + i*nphase + j];
					}
					else if (nstokes == 2)
					{
						p_temp[j] = sub->p_multi[nstokes*nchn*nphase + i*nphase + j];
						s_temp[j] = sub->s_multi[nstokes*nchn*nphase + i*nphase + j];
					}
					else if (nstokes == 3)
					{
						p_temp[j] = sub->p_multi[nstokes*nchn*nphase + i*nphase + j];
						s_temp[j] = sub->s_multi[nstokes*nchn*nphase + i*nphase + j];
					}
					//s_temp[j] = s_multi[i*nphase + j];
				}

				// calculate toa, rms for each channel
				get_toa (s_temp, p_temp, sub, header);

				// if tmode == 1, get TOA for each channel, and transform phase shifts to MJD TOAs
				if (tmode == 1)
				{
					form_toa (sub, header);
					fprintf (fpt, "%s  %lf  %.15Lf  %Lf  7 -f c%d\n", sub->fname, sub->freq[sub->indexChn], sub->t, sub->e_dt*1e+6, i+1);
				}
			}

			// if tmode == 0, do freq-dependent template matching, get one phase shift
			if (tmode == 0)
			{
				if (fitDM == 0 )
				{
					get_toa_multi (sub, header);
				}
				else
				{
					getToaMultiDM (sub, header);
				}

				// transform phase shifts to MJD TOAs
				form_toa_multi (sub, header);

				fprintf (fpt, "%s  %lf  %.15Lf  %Lf  7\n", sub->fname, sub->frequency, sub->t, sub->e_dt*1e+6);
			}
		}
	}

	if (fclose (fpt) != 0)
		fprintf (stderr, "Error closing\n");

	free(header);
	demallocSub(sub, header);
	
	return 0;
}

