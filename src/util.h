#define TINY 1.0e-10
#define PosInf 1.0e20

void mtranspose(double*, double*, int, int);
void mmul(double*, double*, double*, int, int, int);
void madd(double*, double*, double*, int);
void msub(double*, double*, double*, int);
void kmul(double, double*, double*, int);
void dmul(double*, double*, double*, int);
void ludcmp(double *, int, int *, double *);
void lubksb(double *, int, int *, double *);
void inverse(double *, int, double *);
void mmean(double *, int, int, double *, int);
void msum(double *, int, int, double *, int);
void covariance(double*, int, int, double*);
void kron(double *, int, int, double *, int, int, double *);
double normsquare(double*, int);
double ncdf(double);
double normdiff(double*, double*, int);

int duplicated(double*, int, int);
//int unique(double*, int, int, double*, int *);
double* unique(double*, int, int, int *);
void mmax(double *, int, int, double*, int);

// operation between matrix and vector in R fashion
void msubv(double *, int, int, double *, double *, int);
void maddv(double *, int, int, double *, double *, int);
void mmulv(double *, int, int, double *, double *, int);
void mdivv(double *, int, int, double *, double *, int);
// x is an n by p matrix
// y is a vector, its length depends on the parameter 'type'
// z is the return value. z is also an n by p matrix
//      it is allowed here that z is equal to x
// type is used to guide the direction of the vector
// type = 1, y is a row vector with size p
//           R code: z <- t(t(x)*y)
// type = 2, y is a column vector with size n
//           R code: z <- (x*y)

void vdivm(double *, double *, int, int, double *, int);
void vsubm(double *, double *, int, int, double *, int);

