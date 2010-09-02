#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "util.h"

void mtranspose(double* A, double* B, int n, int p)
{
	// B is transpose of A. A is an n by p matrix, and B is a p by n matrix
	
	register int i, j;
	double * C = NULL;
	int flag = 0;
	if (A==B) {
		// for the calling A = t(A)
		C = (double*)malloc(n*p*sizeof(double));
		flag = 1;
	} 
	else {
		C = B;
	}
	for (i = 0; i < n; i++) 
		for (j = 0; j < p; j++)
			C[j * n + i] = A[i * p + j];
	
	if(flag) {
		memcpy((void*)B, (void*)C, n*p*sizeof(double));
		free(C);
	}
}

void mmul(double* A, double* B, double* C, int n, int p, int q)
{
	// C = A*B. A is a n by p matrix, B is a p by q matrix, C is a n by q matrix
	// avoid any calling like A=A*B or B=A*B in the parent function
	register int i, j, k;
	for (i = 0; i < n; i++)
		for (j = 0; j < q; j++)
			C[i * q + j] = 0;

	for (i = 0; i < n; i++)
		for (j = 0; j < q; j++)
			for (k = 0; k < p; k++)
				C[i * q + j] += A[i * p + k] * B[k * q + j];
}

double norm(double* A, int n)
{
	// norm of vector A with dimension n
	double normvalue;
	mmul(A, A, &normvalue, 1, n, 1);
	return sqrt(normvalue);
}

void madd(double* A, double* B, double* C, int n)
{
	// C = A + B
	register int i;
	for (i = 0; i < n; i++)
		C[i] = A[i] + B[i];
}

void msub(double* A, double* B, double* C, int n)
{
	// C = A - B
	register int i;
	for (i = 0; i < n; i++)
		C[i] = A[i] - B[i];;
}

void kmul(double k, double* A, double* B, int n)
{
	// scalar * matrix
	// B = k * A
	register int i;
	for (i = 0; i < n; i++)
		B[i] = A[i] * k;
}

void dmul(double* A, double* B, double* C, int n)
{
	// C = A.*B, that is C(i,j) = A(i,j)*B(i,j);
	register int i;
	for (i = 0; i < n; i++)
		C[i] = A[i] * B[i];
}

void ludcmp(double *a, int n, int *indx, double *d)
{
	int i,imax=0,j,k;
	double big,dum,sum,temp;
	double *vv;

	vv=(double *)malloc(n*sizeof(double));
	*d=1.0;
	for (i=0;i<n;i++) 
	{
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a[i*n+j])) > big) big=temp;
		if (big == 0.0)
		{
			// printf("Singular matrix in routine ludcmp\n");
			return;
		}
		vv[i]=1.0/big;
	}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a[i*n+j];
			for (k=0;k<i;k++) sum -= a[i*n+k]*a[k*n+j];
			a[i*n+j]=sum;
		}
		big=0.0;
		for (i=j;i<n;i++) {
			sum=a[i*n+j];
			for (k=0;k<j;k++)
				sum -= a[i*n+k]*a[k*n+j];
			a[i*n+j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=0;k<n;k++) {
				dum=a[imax*n+k];
				a[imax*n+k]=a[j*n+k];
				a[j*n+k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j*n+j] == 0.0) a[j*n+j]=TINY;
		if (j != n) {
			dum=1.0/(a[j*n+j]);
			for (i=j+1;i<n;i++) a[i*n+j] *= dum;
		}
	}
	free(vv);
}

void lubksb(double *a, int n, int *indx, double b[])
{
	int i,ii=-1,ip,j;
	double sum;

	for (i=0;i<n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii+1)
			for (j=ii;j<=i-1;j++) sum -= a[i*n+j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n-1;i>=0;i--) {
		sum=b[i];
		for (j=i+1;j<n;j++) sum -= a[i*n+j]*b[j];
		b[i]=sum/a[i*n+i];
	}
}

void inverse(double *a, int n, double *y)
{
	double d, *col=NULL;
	int i,j,*indx=NULL;
	double *aa=NULL;
	aa = (double*)malloc(n*n*sizeof(double));
	memcpy((void*)aa, (void*)a, n*n * sizeof(double));
	indx = (int*)malloc(n*sizeof(int));
	col = (double*)malloc(n*sizeof(double));

	ludcmp(aa,n,indx,&d);
	for(j=0;j<n;j++)
	{
		for(i=0;i<n;i++) col[i]=0.0;
		col[j]=1.0;
		lubksb(aa,n,indx,col);
		for(i=0;i<n;i++) y[i*n+j]=col[i];
	}

	free(aa);
	free(indx);
	free(col);
}

void mmean(double *a, int n, int p, double *b, int type)
{
	// matrix a is n by p, b is the mean vector
	// type is used to guide the direction of the average;
	// type = 2, average along the column, b has p elements
	// type = 1, average along the row, b has n elements
	int nn=0, pp=0, i;
	switch(type)
	{
		case 2:
			nn = n;
			pp = p;
			break;
		case 1:
			nn = p;
			pp = n;
			break;
	}
	msum(a,n,p,b,type);
	for(i=0;i<pp;i++)
	{
		b[i] /= nn;
	}
}

void msum(double *a, int n, int p, double *b, int type)
{
	// matrix a is n by p, b is the mean vector
	// type is used to guide the direction of sum;
	// type = 2, sum along the column, b has p elements
	// type = 1, sum along the row, b has n elements
	double *c=NULL;
	int nn=0, pp=0;
	register int i,j;
	c = (double*)malloc(n*p*sizeof(double));
	switch(type)
	{
		case 2:
			memcpy((void*)c, (void*)a, n*p * sizeof(double));
			nn = n;
			pp = p;
			break;
		case 1:
			mtranspose(a, c, n, p);
			nn = p;
			pp = n;
			break;
	}
	for(i=0;i<pp;i++)
		b[i] = 0;
	for(i=0;i<pp;i++)
	{
		for(j=0;j<nn;j++)
			b[i] += c[j*pp+i];
	}
	free(c);
}

void covariance(double *a, int n, int p, double *cov)
{
	double *St=NULL;
	double *m=NULL;
	double *mt=NULL;
	double k;
	register int i,j;

	mt = (double*)malloc(p*sizeof(double));
	m = (double*)malloc(p*sizeof(double));
	St = (double*)malloc(p*p*sizeof(double));

	k = 1.0/(n-1); 

	mmean(a, n, p, m, 1);

	for (i=0; i<p*p; i++) cov[i] = 0;

	for (i=0; i<n; i++)
	{
		for (j=0; j<p; j++)
			mt[j] = a[i*p+j] - m[j];
		mmul(mt, mt, St, p, 1, p);
		for (j=0; j<p*p; j++)
			cov[j] += k*St[j];
	}

	free(m);
	free(mt);
	free(St);
}



void kron(double *A, int rowA, int colA, double *B, int rowB, int colB, double *C)
{
	// C = A(*)B. Kronecker tensor product of A and B
	// A is a rowA-by-colA matrix
	// B is a rowB-by-colB matrix
	// C is a (rowA*rowB)-by-(colA*colB) matrix
	int i,j,k,l;
	int alpha, beta;
	int rowC, colC;

	rowC = rowA * rowB;
	colC = colA * colB;

	for(i=0;i<rowA;i++)
		for(j=0;j<colA;j++)
			for(k=0;k<rowB;k++)
				for(l=0;l<colB;l++)
				{
					alpha = rowB*i+k;
					beta = colB*j+l;
					C[alpha*colC+beta] = A[i*colA+j]*B[k*colB+l];
				}
}

double ncdf(double x)
{

	//	   e = ncdf(double x);
	//     double e -> normal CDF (0<=e<=1)
	//     double x -> real value

	//   - normal CDF (translation of the MatLab
	//     implementation of an algorithm by 
	//     W. J. Cody; see his paper "Rational Chebyshev
	//     approximations for the error function",
	//     Math. Comp., 1969, pp. 631-638)

	const double a[] = {3.16112374387056560e00,
						1.13864154151050156e02,
						3.77485237685302021e02,
						3.20937758913846947e03,
						1.85777706184603153e-1};
	const double b[] = {2.36012909523441209e01,
						2.44024637934444173e02,
						1.28261652607737228e03,
						2.84423683343917062e03};
	const double c[] = {5.64188496988670089e-1,
						8.88314979438837594e00,
						6.61191906371416295e01,
						2.98635138197400131e02,
						8.81952221241769090e02,
						1.71204761263407058e03,
						2.05107837782607147e03,
						1.23033935479799725e03,
						2.15311535474403846e-8};
	const double d[] = {1.57449261107098347e01,
						1.17693950891312499e02,
						5.37181101862009858e02,
						1.62138957456669019e03,
						3.29079923573345963e03,
						4.36261909014324716e03,
						3.43936767414372164e03,
						1.23033935480374942e03};
	const double p[] = {3.05326634961232344e-1,
						3.60344899949804439e-1,
						1.25781726111229246e-1,
						1.60837851487422766e-2,
						6.58749161529837803e-4,
						1.63153871373020978e-2};
	const double q[] = {2.56852019228982242e00,
						1.87295284992346047e00,
						5.27905102951428412e-1,
						6.05183413124413191e-2,
						2.33520497626869185e-3};
	const double xbreak = 0.46875;
	const double pi = acos(-1);
	double y,z,xnum,xden,r,del,e;
	int i;

	x = x/sqrt(2);
	y = fabs(x);

	/* evaluate erf for |x|<=0.46875 */

	if (fabs(x) <= xbreak)
	{
		z = y*y;
		xnum = a[4]*z;
		xden = z;
		for (i=0;i<3;i++)
		{
			xnum = (xnum+a[i])*z;
			xden = (xden+b[i])*z;
		}
		r = x*(xnum+a[3])/(xden+b[3]);
	}
	/* evaluate erf for 0.46875<=|x|<=4.0 */
	else if ((fabs(x) > xbreak) && (fabs(x) <= 4.0))
	{
		xnum = c[8]*y;
		xden = y;
		for (i=0; i<7; i++)
		{
			 xnum = (xnum+c[i])*y;
			 xden = (xden+d[i])*y;
		}
		r = (xnum+c[7])/(xden+d[7]);
		z = floor(y*16.0)/16.0;
		del = (y-z)*(y+z);
		r = exp(-z*z)*exp(-del)*r;
		if (x>0)
			r = 1-r;
		else
			r = r-1;
	}
	/* evaluate erf for |x|>4.0 */
	else
	{
		z = 1.0/(y*y);
		xnum = p[5]*z;
		xden = z;
		for (i=0; i<4; i++)
		{
			xnum = (xnum+p[i])*z;
			xden = (xden+q[i])*z;
		}
		r = z*(xnum+p[4])/(xden+q[4]);
		r = ((1/sqrt(pi))-r)/y;
		z = floor(y*16.0)/16.0;
		del = (y-z)*(y+z);
		r = exp(-z*z)*exp(-del)*r;
		if (x>0)
			r = 1-r;
		else
			r = r-1;
	}

	e = 0.5*(1+r);
	if (e>1) e=1;

	return e;
}

double normdiff(double* A, double* B, int n)
{
	// norm of vector A-B with dimension n
	double* C=NULL;
	double normvalue;
	C = (double*)malloc(n*sizeof(double));
	msub(A, B, C, n);
	mmul(C, C, &normvalue, 1, n, 1);
	free(C);
	return sqrt(normvalue);
}

int duplicated(double* x, int n, int p)
{
	//  to check whether there are duplicated rows in the matrix x
	//  Parameters:
    //    x, an n by p matrix
	//    n & p, the dimensions of x
	//  Return values:
	//    1 for yes, 0 for no

	double * y=x, * z=x;
	register int i, j;

	for(i=0;i<n-1;i++)
	{
		z = x+(i+1)*p;
		for(j=i+1;j<n;j++)
		{
			if (!memcmp((void*)y, (void*)z, p*sizeof(double))) return 1;
			z += p;
		}
		y += p;
	}
	return 0;
}

double* unique(double* x, int n, int p, int *pnunique)
{
	// to extract the unique rows of x
	// Parameters:
	//   n & p, the dimensions of x
	//   x, an n by p matrix
	//   *pnunique, the number of unique rows of x
	//   
	// Return Values:
	//   an (*pnunique) by p matrix, the rows storing the unique rows of x
	//   pnunique is also a return value.
	//   ************* REMEMBER to FREE the memory after using it in the parent function

	double * y=x;
	double *xunique=NULL, *z = NULL, *z1=NULL;
	int nunique = *pnunique;
	register int i, j;

	z = (double*)malloc(n*p*sizeof(double));
	z1 = z;
	memcpy((void*)z1, (void*)y, p * sizeof(double));
	nunique = 1;
	y += p;
	for(i=1;i<n;i++)  {
		z1 = z;
		for(j=0; j<nunique; j++)  {
			if (!memcmp((void*)y, (void*)z1, p*sizeof(double))) break;
			z1 += p;
		}
		if (j==nunique) {
			memcpy((void*)z1, (void*)y, p * sizeof(double));
			nunique++;
		}
		y += p;
	}
	xunique = (double*)malloc(nunique*p*sizeof(double));
	memcpy((void*)xunique, (void*)z, nunique*p*sizeof(double));
	free(z);
	(*pnunique) = nunique;
	return xunique;
}

/*
void unique(double* x, int *pn, int *pp, double* xunique, int *pnunique)
{
	// to extract the unique rows of x
	// Parameters:
	//   n & p, the dimensions of x
	//   x, an n by p matrix
	//   *pnunique, the number of unique rows of x
	//   xunique, an n by p matrix, the first (*pnunique) rows storing the unique rows of x, the others are set to zero
	// Return Values:
	//   1 as all rows are unique, 0 otherwise ( 0 for duplicated rows existing )  
	//   xunique and pnunique are also return values.

	int n=*pn, p=*pp;
	double * y=x, *z = xunique;
	int nunique = *pnunique;
	register int i, j;

	memeset((void*)z, '\0', n*p*sizeof(double));
	memcpy((void*)z, (void*)y, p * sizeof(double));
	nunique = 1;
	y += p;
	for(i=1;i<n;i++)  {
		z = xunique;
		for(j=0; j<nunique; j++)  {
			if (!memcmp((void*)y, (void*)z, p*sizeof(double))) break;
			z += p;
		}
		if (j==nunique) {
			memcpy((void*)z, (void*)y, p * sizeof(double));
			nunique++;
		}
		y += p;
	}
	(*pnunique) = nunique;
}
*/

void mmax(double *x, int n, int p, double* y, int type)
{
	// maximal elements of a matrix
	// matrix x is n by p, y is the vector of maximal elements
	// type is used to guide the direction of maximum;
	// type = 1, max along the row, b has n elements
	//           R code: y <- apply(x, 1, max)
	// type = 2, max along the column, b has p elements
	//           R code: y <- apply(x, 2, max)

	register int i, j;
	double *c=NULL;
	int nn=0, pp=0;
	c = (double*)malloc(n*p*sizeof(double));
	switch(type)
	{
		case 2:
			memcpy((void*)c, (void*)x, n*p * sizeof(double));
			nn = n;
			pp = p;
			break;
		case 1:
			mtranspose(x, c, n, p);
			nn = p;
			pp = n;
			break;
	}
	memcpy((void*)y, (void*)c, pp*sizeof(double));
	for(i=0; i<pp; i++)  {
		for (j=1; j<nn; j++) {
			if (c[j*pp+i]>y[i]) y[i] = c[j*pp+i];
		}
	}
	free(c);
}


void msubv(double *x, int n, int p, double *y, double *z, int type)
{
	// matrix minus vector in R fashion
	// x is an n by p matrix
	// y is a vector, its length depends on the parameter 'type'
	// z is the return value. z is also an n by p matrix
	//      it is allowed here that z is equal to x
	// type is used to guide the direction of the vector
	// type = 1, y is a row vector with size p
	//           R code: z <- t(t(x)-y)
	// type = 2, y is a column vector with size n
	//           R code: z <- (x-y)

	double * c = NULL;
	register int i;
	int nn=n, pp=p;

	if(type==2) {nn=p; pp=n;}

	c = malloc(n*p*sizeof(double));
	for(i = 0; i < nn; i++)  memcpy((void*)(c+i*pp), (void*)y, pp*sizeof(double));
	if (type==2) mtranspose(c, c, nn, pp);
	msub(x, c, c, n*p);
	memcpy((void*)z, (void*)c, n*p*sizeof(double));
	free(c);
}

void maddv(double *x, int n, int p, double *y, double *z, int type)
{
	// matrix plus vector in R fashion
	// x is an n by p matrix
	// y is a vector, its length depends on the parameter 'type'
	// z is the return value. z is also an n by p matrix
	//      it is allowed here that z is equal to x
	// type is used to guide the direction of the vector
	// type = 1, y is a row vector with size p
	//           R code: z <- t(t(x)+y)
	// type = 2, y is a column vector with size n
	//           R code: z <- (x+y)

	double * c = NULL;
	register int i;
	int nn=n, pp=p;

	if(type==2) {nn=p; pp=n;}

	c = malloc(n*p*sizeof(double));
	for(i = 0; i < nn; i++)  memcpy((void*)(c+i*pp), (void*)y, pp*sizeof(double));
	if (type==2) mtranspose(c, c, nn, pp);
	madd(x, c, c, n*p);
	memcpy((void*)z, (void*)c, n*p*sizeof(double));
	free(c);
}

void mmulv(double *x, int n, int p, double *y, double *z, int type)
{
	// matrix multiply vector in R fashion
	// x is an n by p matrix
	// y is a vector, its length depends on the parameter 'type'
	// z is the return value. z is also an n by p matrix
	//      it is allowed here that z is equal to x
	// type is used to guide the direction of the vector
	// type = 1, y is a row vector with size p
	//           R code: z <- t(t(x)*y)
	// type = 2, y is a column vector with size n
	//           R code: z <- (x*y)

	double * c = NULL;
	register int i;
	int nn=n, pp=p;

	if(type==2) {nn=p; pp=n;}

	c = malloc(n*p*sizeof(double));
	for(i = 0; i < nn; i++)  memcpy((void*)(c+i*pp), (void*)y, pp*sizeof(double));
	if (type==2) mtranspose(c, c, nn, pp);
	dmul(x, c, c, n*p);
	memcpy((void*)z, (void*)c, n*p*sizeof(double));
	free(c);
}

void mdivv(double *x, int n, int p, double *y, double *z, int type)
{
	// divide a matrix by a vector in R fashion
	// x is an n by p matrix
	// y is a vector, its length depends on the parameter 'type'
	// z is the return value. z is also an n by p matrix
	//      it is allowed here that z is equal to x
	// type is used to guide the direction of the vector
	// type = 1, y is a row vector with size p
	//           R code: z <- t(t(x)/y)
	// type = 2, y is a column vector with size n
	//           R code: z <- (x/y)

	double * c = NULL;
	register int i;
	int nn=n, pp=p;

	if(type==2) {nn=p; pp=n;}

	c = malloc(n*p*sizeof(double));
	for(i = 0; i < nn; i++)  memcpy((void*)(c+i*pp), (void*)y, pp*sizeof(double));
	for(i=0; i<n*p; i++)  {
		if (c[i]==0) { c[i] = TINY; }  // in case of dividing by zero
		c[i] = 1.0 / c[i];
	}

	if (type==2) mtranspose(c, c, nn, pp);
	dmul(x, c, c, n*p);
	memcpy((void*)z, (void*)c, n*p*sizeof(double));
	free(c);
}

void vdivm(double *y, double *x, int n, int p, double *z, int type)
{
	// divide a vector by a matrix in R fashion
	// y is a vector, its length depends on the parameter 'type'
	// x is an n by p matrix
	// z is the return value. z is also an n by p matrix
	//      it is allowed here that z is equal to x
	// type is used to guide the direction of the vector
	// type = 1, y is a row vector with size p
	//           R code: z <- t(y/t(x))
	// type = 2, y is a column vector with size n
	//           R code: z <- (y/x)

	double * c = NULL;
	register int i;
	int nn=n, pp=p;

	if(type==2) {nn=p; pp=n;}

	c = malloc(n*p*sizeof(double));
	for(i = 0; i < nn; i++)  memcpy((void*)(c+i*pp), (void*)y, pp*sizeof(double));
	if (type==2) mtranspose(c, c, nn, pp);

	for(i=0; i<n*p; i++)  {
		if (x[i]==0) {
			c[i] /= TINY;  // in case of dividing by zero
		}
		else c[i] /= x[i];
	}

	memcpy((void*)z, (void*)c, n*p*sizeof(double));
	free(c);
}


void vsubm(double *y, double *x, int n, int p, double *z, int type)
{
	// vector minus matrix in R fashion
	// y is a vector, its length depends on the parameter 'type'
	// x is an n by p matrix
	// z is the return value. z is also an n by p matrix
	//      it is allowed here that z is equal to x
	// type is used to guide the direction of the vector
	// type = 1, y is a row vector with size p
	//           R code: z <- t(y-t(x))
	// type = 2, y is a column vector with size n
	//           R code: z <- (y-x)


	register int i;
	msubv(x, n, p, y, z, type);
	for(i=0; i<n*p; i++) *(z+i) = -(*(z+i));
}

