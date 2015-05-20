#ifndef BFGS_H
#define BFGS_H

/** 
   \file bfgs.h
   Header file for interface function used in BFGS optimization - setULbounds()
   
   This implementation was extracted from HMMld 0.801.
   The code was apparently translated from fortran to C by f2c (version 19991025),
   so it is extremely messy.
   Documentation exists in the C file
   
   Original disclaimer:
   "  Copyright (C) 2008 Thorfinn Sand Korneliussen  thorfinn@binf.ku.dk
      
      HMMld is free software: you can redistribute it and/or modify
      it under the terms of the GNU General Public License as published by
      the Free Software Foundation, either version 3 of the License, or
      (at your option) any later version.                                    "
*/
int setULbounds(int *n, int *m, double *x, double *l, double *u, 
			     int *nbd, double *f, double *g, double *factr, 
			     double *pgtol, double *wa, int *iwa, 
			     char *task, int *iprint, char *csave, int *lsave,
			     int *isave, double *dsave);
				 

#endif

