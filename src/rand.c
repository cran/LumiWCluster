#include <stdlib.h>
#include <time.h>
// #include <math.h>
// #include <string.h>
#include <memory.h>
// #include "util.h"

#include "rand.h"


int sampling_int(int* pool, int nPool, int* sample, int nSize, int replace)
{
	// Select nSize samples from the pool with or without replacement 
	// Parameters:
	//   pool: the pool for sampling, a vector of size nPool
	//   nPool: the size of the pool
	//   sample: the selected sample, a vector of nSize, it is the return value
	//   replace: indicator of replacement. replace = 1, with replacement, replace = 0 without replacement
	// Return values:
	//   1 for success 0 for failure
	//   sample with nSize is the return value. Be careful of freeing the memory in the parent function

	int *tmppool=NULL, nSampleIndex;
	register int i;
	// srand (time(NULL));
	if (replace) {
		for(i=0; i<nSize; i++) {
			nSampleIndex = rand() % nPool;
			sample[i] = pool[nSampleIndex];
		}
	} 
	else {
		if (nSize>nPool) {
			// printf("The size of sample should be less than the size of pool!");
			return 0;
		}
		tmppool = (int*)malloc(nPool*sizeof(int));
		memcpy((void*)tmppool, (void*)pool, nPool*sizeof(int));
		for(i=0; i<nSize; i++)  {
			nSampleIndex = rand() % nPool;
			sample[i] = pool[nSampleIndex];
			if (nPool-nSampleIndex > 1) memcpy((void*)(pool+nSampleIndex), (void*)(pool+nSampleIndex+1), (nPool-nSampleIndex-1)*sizeof(int));
			nPool--;
		}
	}
	free(tmppool);
	return 1;
}


int sampling_double(double* pool, int nPool, double* sample, int nSize, int replace)
{
	// Select nSize samples from the pool with or without replacement 
	// Parameters:
	//   pool: the pool for sampling, a vector of size nPool
	//   nPool: the size of the pool
	//   sample: the selected sample, a vector of nSize, it is the return value
	//   replace: indicator of replacement. replace = 1, with replacement, replace = 0 without replacement
	// Return values:
	//   1 for success 0 for failure
	//   sample with nSize is the return value. Be careful of freeing the memory in the parent function

	double *tmppool=NULL;
	int nSampleIndex;
	register int i;
	srand (time(NULL));
	if (replace) {
		for(i=0; i<nSize; i++) {
			nSampleIndex = rand() % nPool;
			sample[i] = pool[nSampleIndex];
		}
	} 
	else {
		if (nSize>nPool) {
			// printf("The size of sample should be less than the size of pool!");
			return 0;
		}
		tmppool = (double*)malloc(nPool*sizeof(double));
		memcpy((void*)tmppool, (void*)pool, nPool*sizeof(double));
		for(i=0; i<nSize; i++)  {
			nSampleIndex = rand() % nPool;
			sample[i] = pool[nSampleIndex];
			if (nPool-nSampleIndex > 1) memcpy((void*)(pool+nSampleIndex), (void*)(pool+nSampleIndex+1), (nPool-nSampleIndex-1)*sizeof(double));
			nPool--;
		}
	}
	free(tmppool);
	return 1;
}



int sampling_by_index(int nPool, int* sample, int nSize, int replace)
{
	// Select index of nSize samples from the pool (1:nPool) with or without replacement 
	// Parameters:
	//   nPool: the size of the pool, the index of pool is 0:nPool-1
	//   sample: the selected sample index, a vector of nSize, it is the return value
	//   replace: indicator of replacement. replace = 1, with replacement, replace = 0 without replacement
	// Return values:
	//   1 for success 0 for failure
	//   sample with nSize is the return value. Be careful of freeing the memory in the parent function

	int *pool=NULL;
	int flag;
	register int i;
	pool = (int*)malloc(nPool*sizeof(int));
	for(i=0; i<nPool; i++) pool[i] = i;
	flag = sampling_int(pool, nPool, sample, nSize, replace);
	free(pool);
	return flag;
}


/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000-2002 the R Development Core Team
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
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *    double exp_rand(void);
 *
 *  DESCRIPTION
 *
 *    Random variates from the standard exponential distribution.
 *
 *  REFERENCE
 *
 *    Ahrens, J.H. and Dieter, U. (1972).
 *    Computer methods for sampling from the exponential and
 *    normal distributions.
 *    Comm. ACM, 15, 873-882.
 */

double exp_rand(void)
{
    /* q[k-1] = sum(log(2)^k / k!)  k=1,..,n, */
    /* The highest n (here 8) is determined by q[n-1] = 1.0 */
    /* within standard precision */
    const static double q[] =
    {
	0.6931471805599453,
	0.9333736875190459,
	0.9888777961838675,
	0.9984959252914960,
	0.9998292811061389,
	0.9999833164100727,
	0.9999985691438767,
	0.9999998906925558,
	0.9999999924734159,
	0.9999999995283275,
	0.9999999999728814,
	0.9999999999985598,
	0.9999999999999289,
	0.9999999999999968,
	0.9999999999999999,
	1.0000000000000000
    };
    double a, u, ustar, umin;
    int i;

    a = 0.;
    /* precaution if u = 0 is ever returned */
    u = unif_rand();
    while(u <= 0.0 || u >= 1.0) u = unif_rand();
    for (;;) {
		u += u;
		if (u > 1.0)
			break;
		a += q[0];
    }
    u -= 1.;

    if (u <= q[0])
	return a + u;

    i = 0;
    ustar = unif_rand();
    umin = ustar;
    do {
		ustar = unif_rand();
		if (ustar < umin)
			umin = ustar;
		i++;
    } while (u > q[i]);
    return a + umin * q[0];
}


/* A version of Marsaglia-MultiCarry */

static unsigned int I1=1234, I2=5678;

void set_seed(unsigned int i1, unsigned int i2)
{
    I1 = i1; I2 = i2;
}

void get_seed(unsigned int *i1, unsigned int *i2)
{
    *i1 = I1; *i2 = I2;
}


double unif_rand(void)
{
    I1= 36969*(I1 & 0177777) + (I1>>16);
    I2= 18000*(I2 & 0177777) + (I2>>16);
    return ((I1 << 16)^(I2 & 0177777)) * 2.328306437080797e-10; /* in [0,1) */
}

