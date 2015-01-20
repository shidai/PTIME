#!/bin/sh

gcc -lm -lgsl -lgslcblas -lfftw3 -L/usr/local/lib/cfitsio -I/usr/include/cfitsio/ -lcfitsio ptimeT.c ptimeLib.c ptimeTLib.c readPfits.c T2toolkit.c tempo2pred.c cheby2d.c t1polyco.c -o ptimeT 

gcc -lm -lgsl -lgslcblas -lfftw3 -L/usr/local/lib/cfitsio -I/usr/include/cfitsio/ -lcfitsio ptimeD.c ptimeLib.c ptimeTLib.c ptimeDLib.c readPfits.c T2toolkit.c tempo2pred.c cheby2d.c t1polyco.c -o ptimeD 

gcc -lm -lpgplot -lcpgplot ptime_plot.c  ptimeLib.c -o ptime_plot 

gcc -lm -lpgplot -lcpgplot -L/usr/local/lib/cfitsio -I/usr/include/cfitsio/ -lcfitsio -lgsl -lgslcblas ptime_create.c ptimeLib.c readPfits.c -o ptime_create   

gcc -lm -lpgplot -lcpgplot -L/usr/local/lib/cfitsio -I/usr/include/cfitsio/ -lcfitsio ptime_modify.c ptimeLib.c readPfits.c -o ptime_modify   

