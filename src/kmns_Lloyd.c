/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2004   The R Development Core Team.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */


//**************************************************************************
//     Implementation R function Kmeans (Lloyd agorithm)
//     (R source code is stored in the same directory for referece.)
//     Author: Xin Zhou
//     Date  : Aug. 19, 2010
//**************************************************************************



#include <stdlib.h>
// #include <stdio.h>
// #include <time.h>
#include <memory.h>
// #include <string.h>

#include "util.h"
#include "rand.h"
#include "kmns_Lloyd.h"


//**************************************************************************
//   Main function
//     Parameters:
//        x, an n by p matrix storing the sample points, one row as one sample point
//        n & p, dimensions of x
//        k, the number of clusters
//        maxiter, the maximum number of iterations allowed
//        nstart, the same meaning with the one in R function
//        (*piter), number of iterations during the clustering
//        cl, a vector with size n storing the cluster id for each sample point
//        centers, a k by p matrix storing the centers of clusters
//        wss, a vector with size k,  each element is the within-cluster sum of squares for the corresponding cluster
//        nc, a vector of size k,  the number of points in each cluster
//     Return values:
//        1 for clustering success,     
//        2 for not converge in maxiter iterations, 
//        3 for existing empty cluster, may try to set a large nstart,  
//        4 for number of unique points less than the number of clusters
//          piter, cl, centers, wss, and nc are also return values
//     Notes:
//        The parameter "centers" is only used as return values. It is not the initial set of centers
//        The program random selects sample points as the initial centers. 
//        By setting a large number to nstart, different sets of initial centers may be compared to get the better clustering
//         
//****************************************************************************
int kmns_lloyd(double *x, int n, int p, int k, int maxiter, int nstart, int* piter, int* cl, double *centers, double *wss, int *nc)
{
	double *cn = NULL;
	int *sampling = NULL;  // sampling stores the indices of sample points as the initial centers
	register int i, j;
	int flagDuplicateCenters = 0;  // indicator for duplicated centers
	int nunique = n;  // number of unique points
	int flagStatus = 1;  // return value

	int* cl1 = NULL, *nc1 = NULL;
	int iter=*piter, iter1; // record the iterations for kmeans
	double* wss1 = NULL;
	double* centers1 = NULL;
	double bestsumwss = 0.0, sumwss;

	// srand (time(NULL));   //it can be initialized in the parent function
	sampling = (int*)malloc(k*sizeof(int));
	if (nstart==1) {
		sampling_by_index(n, sampling, k, 0);
		for(i=0; i<k; i++)  {
			memcpy((void*)(centers+i*p), (void*)(x+sampling[i]*p), p*sizeof(double));
		}
		flagDuplicateCenters = duplicated(centers, k, p);
	}

	if (nstart>1 || flagDuplicateCenters)  {
		cn = unique(x, n, p, &nunique);
		if (nunique < k) {
			flagStatus = 4;  // number of unique points less than the number of clusters
			goto stop;
		}
		sampling_by_index(nunique, sampling, k, 0);
		for(i=0; i<k; i++)  {
			memcpy((void*)(centers+i*p), (void*)(cn+sampling[i]*p), p*sizeof(double));
		}
	}
	kmeans_Lloyd(x, n, p, centers, k, cl, maxiter, &iter, nc, wss);

	if (nstart>1 && (cn!=NULL))  {
		for(i=0; i<k; i++)  bestsumwss += wss[i];   // save the sum of wss for comparison
		
		centers1 = (double*)malloc(k*p*sizeof(double));
		cl1 = (int*)malloc(n*sizeof(int));
		nc1 = (int*)malloc(k*sizeof(int));
		wss1 = (double*)malloc(k*sizeof(double));

		for(j=1; j<nstart; j++)  {
			sampling_by_index(nunique, sampling, k, 0);
			for(i=0; i<k; i++)  {
				memcpy((void*)(centers1+i*p), (void*)(cn+sampling[i]*p), p*sizeof(double));
			}
			kmeans_Lloyd(x, n, p, centers1, k, cl1, maxiter, &iter1, nc1, wss1);
			sumwss = 0.0;
			for(i=0; i<k; i++)  sumwss += wss1[i];
			if( sumwss < bestsumwss )  {
				bestsumwss = sumwss;
				memcpy((void*)centers, (void*)centers1, k*p*sizeof(double));
				memcpy((void*)cl, (void*)cl1, n*sizeof(int));
				memcpy((void*)nc, (void*)nc1, k*sizeof(int));
				memcpy((void*)wss, (void*)wss1, k*sizeof(double));
				iter = iter1;
			}
		}
	}
	*piter = iter;
	if (iter > maxiter) flagStatus = 2;  // not converge
	for (i=0; i<k; i++) { if (nc[i]==0) { flagStatus = 3; break;} } // empty cluster
	free(centers1);
	free(cl1);
	free(nc1);
	free(wss1);
stop:
	if(cn != NULL) free(cn);
	free(sampling);
	return flagStatus;
}



//**************************************************************************
//   Modified from main function by only returning centers, and discarding piter, nc, cl, wss
//     Parameters:
//        x, an n by p matrix storing the sample points, one row as one sample point
//        n & p, dimensions of x
//        k, the number of clusters
//        maxiter, the maximum number of iterations allowed
//        nstart, the same meaning with the one in R function
//        centers, a k by p matrix storing the centers of clusters
//     Return values:
//        1 for clustering success,     
//        2 for not converge in maxiter iterations, 
//        3 for existing empty cluster, may try to set a large nstart,  
//        4 for number of unique points less than the number of clusters
//          centers is also a return value
//     Notes:
//        The parameter "centers" is only used as return values. It is not the initial set of centers
//        The program random selects sample points as the initial centers. 
//        By setting a large (or proper) number to nstart, different sets of initial centers may be compared to get the better clustering
//         
//****************************************************************************


int kmns_lloyd_compact_ret(double *x, int n, int p, int k, int maxiter, int nstart, double *centers, double *sumwss)
{
	double *cn = NULL;
	int *sampling = NULL;  // sampling stores the indices of sample points as the initial centers
	register int i, j;
	int flagDuplicateCenters = 0;  // indicator for duplicated centers
	int nunique = n;  // number of unique points
	int flagStatus = 1;  // return value

	int* cl1 = NULL, *nc1 = NULL;
	int* cl = NULL, *nc=NULL;
	int iter, iter1; // record the iterations for kmeans
	double* wss=NULL, *wss1 = NULL;
	double* centers1 = NULL;
	double bestsumwss = 0.0, sumwss1;

	// srand (time(NULL));    // may be called in parent function or here
	sampling = (int*)malloc(k*sizeof(int));
	if (nstart==1) {
		sampling_by_index(n, sampling, k, 0);
		for(i=0; i<k; i++)  {
			memcpy((void*)(centers+i*p), (void*)(x+sampling[i]*p), p*sizeof(double));
		}
		flagDuplicateCenters = duplicated(centers, k, p);
	}

	if (nstart>1 || flagDuplicateCenters)  {
		cn = unique(x, n, p, &nunique);
		if (nunique < k) {
			flagStatus = 4;  // number of unique points less than the number of clusters
			goto stop;
		}
		sampling_by_index(nunique, sampling, k, 0);
		for(i=0; i<k; i++)  {
			memcpy((void*)(centers+i*p), (void*)(cn+sampling[i]*p), p*sizeof(double));
		}
	}
	cl = (int*)malloc(n*sizeof(int));
	nc = (int*)malloc(k*sizeof(int));
	wss = (double*)malloc(k*sizeof(double));
	kmeans_Lloyd(x, n, p, centers, k, cl, maxiter, &iter, nc, wss);

	if (nstart>1 && (cn!=NULL))  {
		bestsumwss = 0.0;
		for(i=0; i<k; i++)  bestsumwss += wss[i];   // save the sum of wss for comparison
		
		centers1 = (double*)malloc(k*p*sizeof(double));
		cl1 = (int*)malloc(n*sizeof(int));
		nc1 = (int*)malloc(k*sizeof(int));
		wss1 = (double*)malloc(k*sizeof(double));

		for(j=1; j<nstart; j++)  {
			sampling_by_index(nunique, sampling, k, 0);
			for(i=0; i<k; i++)  {
				memcpy((void*)(centers1+i*p), (void*)(cn+sampling[i]*p), p*sizeof(double));
			}
			kmeans_Lloyd(x, n, p, centers1, k, cl1, maxiter, &iter1, nc1, wss1);
			sumwss1 = 0.0;
			for(i=0; i<k; i++)  sumwss1 += wss1[i];
			if( sumwss1 < bestsumwss )  {
				bestsumwss = sumwss1;
				memcpy((void*)centers, (void*)centers1, k*p*sizeof(double));
				memcpy((void*)cl, (void*)cl1, n*sizeof(int));
				memcpy((void*)nc, (void*)nc1, k*sizeof(int));
				memcpy((void*)wss, (void*)wss1, k*sizeof(double));
				iter = iter1;
			}
		}
	}
	if (iter > maxiter) flagStatus = 2;  // not converge
	for (i=0; i<k; i++) { if (nc[i]==0) { flagStatus = 3; break;} } // empty cluster
	*sumwss = bestsumwss;
	free(cl);
	free(nc);
	free(wss);
	free(centers1);
	free(cl1);
	free(nc1);
	free(wss1);
stop:
	if(cn != NULL) free(cn);
	free(sampling);
	return flagStatus;
}


//**************************************************************************
//**  kmeans_Lloyd code is from R Team
//**  Modification: 
//**     1) Remove "R.h" and corresponding definitions (such as Rboolean, R_PosInf)
//**     2) Change the some parameters (n, p, k and maxiter) from pointers to non-pointers
//**        Add a pointer parameter piter to return the number of iterations (avoid changing parameter maxiter
//**     3) Change matrices x & center to c convention (stored by rows in the memory)
//**************************************************************************

void kmeans_Lloyd(double *x, int n, int p, double *cen, int k, int *cl,
		  int maxiter, int *piter, int *nc, double *wss)
{
    int iter, i, j, c, it, inew = 0;
    double best, dd, tmp;
    int updated;  // 1 for TRUE, 0 for FALSE

    for(i = 0; i < n; i++) cl[i] = -1;
    for(iter = 0; iter < maxiter; iter++) {
		updated = 0;
		for(i = 0; i < n; i++) {
			/* find nearest centre for each point */
			best = PosInf;
			for(j = 0; j < k; j++) {
				dd = 0.0;
				for(c = 0; c < p; c++) {
					tmp = x[i*p+c] - cen[j*p+c];
					dd += tmp * tmp;
				}
				if(dd < best) {
					best = dd;
					inew = j+1;
				}
			}
			if(cl[i] != inew) {
				updated = 1;
				cl[i] = inew;
		    }
		}
		if(!updated) break;
		/* update each centre */
		for(j = 0; j < k*p; j++) cen[j] = 0.0;
		for(j = 0; j < k; j++) nc[j] = 0;
		for(i = 0; i < n; i++) {
			it = cl[i] - 1;
			nc[it]++;
			for(c = 0; c < p; c++) cen[it*p+c] += x[i*p+c];
		}
		for(j=0; j<k; j++)
			for(i=0; i<p; i++)
				cen[j*p+i] /= nc[j];
    }

    *piter = iter + 1;
    for(j = 0; j < k; j++) wss[j] = 0.0;
    for(i = 0; i < n; i++) {
		it = cl[i] - 1;
		for(c = 0; c < p; c++) {
			tmp = x[i*p+c] - cen[it*p+c];
			wss[it] += tmp * tmp;
		}
    }
}
