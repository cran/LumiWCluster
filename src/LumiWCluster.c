#define MAXATTEMPT 10
#define MAXKMNSATTEMPT 10
// #define PI 3.14159265358979323846264338327956288

#include <stdlib.h>
#include <time.h>
#include <memory.h>
#include <string.h>
#include <R.h>

#include "util.h"
#include "rand.h"
#include "kmns_Lloyd.h"

//**************************************************************************
//   Main function
//     Parameters:
//        orig_y, the original y, an n by p matrix storing the original sample points (logitbeta)
//        y, an n by p matrix storing the normalized sample points (logitbeta), one row as one sample point
//        n & p, dimensions of y
//        mean_y, a vector with size p. all 0 or mean of original y
//        norm_y, a vector with size p. all 1 or sd of original y
//                that is, original y = y*norm_y+mean_y
//        w, a vector with size n. weights for each sample point, any element of w != 0
//        g, a vector with size p. weights for each probe, any element of g != 0
//        vK, a vector with size nK. The number of clusters
//        nK, the length of K
//        maxK, max value of K
//        vLambda, a vector with size nLambda. Penalty parameter
//        nLambda, the length of lambda
//        adaptive, indictor of adaptive weighting. 1 for yes, 0 no
//        epsilon, a small number (e.g. 1e-6)
//        maxiter, maximum iterations allowed in EM
//        trace, print EM iterations or not. 1 for yes, 0 no
// 
//           the following parameters are return values
//
//        BIC, an nK by nLambda matrix storing the BIC for each combination of K and lambda
//        optBIC, the smallest BIC
//        optK, the K to achieve optBIC 
//        optLambda, the lambda to achieve optBIC
//        ClusterID, a vector with size n. Clustering result for each observation
//        converge, indicator of convergence for EM. 1 for yes, 0 no.
//        pi_vec, a vector with size maxK. Only need first optK elements
//        mu, a maxK by p matrix. Only need first optK rows
//        mu_centering, a maxK by p matrix. Only need first optK rows
//        init_mu, a maxK by p matrix. Only need first optK rows
//        sigma2_vec, a vector with size p
//        iter, the EM iterations
//
//     Return values:
//        None
//         
//     Notes:
//        No parameter checking in the C code
//        The parameters passed from R should be checked carefully there to avoid any unexpected breakdown here.
//****************************************************************************



void LumiWCluster(double *orig_y, double *y, int *pn, int *pp, double *mean_y, double *norm_y, 
				 double *w, double *g, int *vK, int *pnK, int *pmaxK, double *vLambda,
				 int *pnLambda, int *padaptive, double *pepsilon, int *pmaxiter, int *ptrace,
				 double *BIC, double *poptBIC, int *poptK, double *poptLambda, int *ClusterID,
				 int *pconverge, double *pi_vec, double *mu, double *mu_centering,
				 double *init_mu, double *sigma2_vec, int *piter)
{
	int n=*pn, p=*pp, nK=*pnK, maxK = *pmaxK, nLambda=*pnLambda;
	int adaptive=*padaptive, maxiter=*pmaxiter, trace=*ptrace;
	double epsilon=*pepsilon, optBIC = PosInf;
	int converge = 0, converge1 = 0;
	int optK = 1;
	double optLambda = 0.0;
	register int i, j, r, s, t; // counter
	int K = 2;
	double lambda = 0.0;
	int npara;
	double lik1, lik2;

	// pi_vec1, mu1, mu_centering1, init_mu1, sigma2_vec1, iter1 record the results for each K & lambda combination
	// pi_vec, mu, mu_centering, init_mu, sigma2_vec, iter keep the best ones
	double *pi_vec1 = NULL, *mu1=NULL, *mu_centering1=NULL, *init_mu1=NULL, *sigma2_vec1=NULL;
	int iter = 0, iter1 = 0;
	double tmp;
	double *z=NULL;  // a maxK by n matrix
	double *theta=NULL;  // a maxK by p matrix
	double *weight1 = NULL;     // a vector with size p
	double *weight2 = NULL;     // a maxK by p matrix
	double *lambda_theta = NULL;  // a maxK by p matrix
	double *mu_old = NULL;   // a maxK by p matrix
	double *kern_mat = NULL;  // a maxK by n matrix
	double *tmp_y = NULL;  // an n by p matrix, temporary storage
	double *tmp_kern_mat = NULL;  // a maxK by n matrix, temporary storage for kern_mat
	double *sigma2_mat = NULL;   // a maxK by p matrix
	double *mu0 = NULL; // a maxK by p matrix
	double *gamma = NULL, *tmp_gamma=NULL; // a vector with size p
	double *tmp_theta = NULL;  // a maxK by p matrix, temporary storage for theta
	double *tmp_mu = NULL;   // a maxK by p matrix
	double *max_kern_mat = NULL;  // a vector with size n

	int nAttempt, nKmnsAttempt;
	double diff = 1.0;
	double * best_kmns_mu = NULL;
	int kmns_status;
	double bestsumwss, sumwss1;

//	srand (time(NULL));   //it can be initialized in the parent function
	srand(2222);
	pi_vec1 = (double*) malloc( maxK * sizeof(double) );
	mu1 = (double*)malloc( maxK * p * sizeof(double) );
	mu_centering1 = (double*)malloc( maxK * p * sizeof(double) );
	init_mu1 = (double*)malloc( maxK * p * sizeof(double) );
	sigma2_vec1 = (double*) malloc( p*sizeof(double) );
	z = (double*)malloc(maxK*n*sizeof(double));
	gamma = (double*)malloc(p*sizeof(double));
	theta = (double*)malloc(maxK*p*sizeof(double));
	weight1 = (double*)malloc(p*sizeof(double));
	weight2 = (double*)malloc(maxK*p*sizeof(double));
	lambda_theta = (double*)malloc(maxK*p*sizeof(double));
	mu_old = (double*)malloc(maxK*p*sizeof(double));
	kern_mat = (double*)malloc(maxK*n*sizeof(double));
	tmp_y = (double*)malloc(n*p*sizeof(double));
	tmp_kern_mat = (double*)malloc(maxK*n*sizeof(double));
	sigma2_mat = (double*)malloc(maxK*p*sizeof(double));
	mu0 = (double*)malloc(maxK*p*sizeof(double));
	tmp_theta = (double*)malloc(maxK*p*sizeof(double));
	tmp_gamma = (double*)malloc(p*sizeof(double));
	tmp_mu = (double*)malloc(maxK*p*sizeof(double));
	
	best_kmns_mu = (double*)malloc(maxK*p*sizeof(double));
	max_kern_mat = (double*)malloc(n*sizeof(double));

	for (i=0; i<nK; i++) {
		K = vK[i];
		Rprintf("For K = %d :\n", K);
		kmns_status = kmns_lloyd_compact_ret(y, n, p, K, 10*n, 1, best_kmns_mu, &bestsumwss);  			// maxiter for kmeans is 10n, and nstart = n
		if (kmns_status == 4) {
			Rprintf("    K = %d is greater than the number of unique samples. Please try smaller K.\n", K);
			continue;
		}
		if (kmns_status == 3)  {
			for(nKmnsAttempt = 0; nKmnsAttempt< MAXKMNSATTEMPT; nKmnsAttempt++)  {
				if(kmns_lloyd_compact_ret(y, n, p, K, 10*n, 1, best_kmns_mu, &bestsumwss)!=3) break;
			}
			if (nKmnsAttempt==MAXKMNSATTEMPT) {
				Rprintf("    Empty clusters exist. Please try smaller K.\n");
				continue;
			}
		}
		if (kmns_status == 2)  {
			Rprintf("    KMeans clustering does not converge. \n");
		}

		for(r=0; r<p; r++) weight1[r] = 1.0/p/g[r];
		for(r=0; r<K*p; r++) weight2[r] = 1.0;

		for (j=0; j<nLambda; j++) {
			lambda = vLambda[j];
			Rprintf("    For lambda = %f :\n", lambda);

			// ***  initialization  ***
			//	pi.vec <- rep(1/K,K)
			tmp = 1.0 / K;
			for(r=0; r<K; r++) pi_vec1[r] = tmp;
			
			//	z <- matrix(0,nrow=K,ncol=n)
			memset(z, 0, maxK*n*sizeof(double));

			//	sigma2.vec <- rep(1,p)
			for(r=0; r<p; r++) sigma2_vec1[r] = 1.0;

			//	if (missing(init.mu)|length(init.mu)==0) {
			//		mu <- kmeans(y,K)$centers
			//		init.mu <- mu
			//	} else mu <- init.mu
			memcpy((void*)mu1, (void*)best_kmns_mu, K*p*sizeof(double));
			memcpy((void*)init_mu1, (void*)mu1, K*p*sizeof(double));
			// return value of kmns??????
			
			
			//	gamma <- rep(1,p)
			//  theta <- mu				
			memcpy((void*)gamma, (void*)sigma2_vec1, p*sizeof(double));
			memcpy((void*)theta, (void*)mu1, K*p*sizeof(double));

			//	if(adaptive){
			//		weight1 <- 1/(g*p*apply(abs(mu),2,max))
			//    	weight2 <- 1/abs(mu)
			//  	lambda.theta <- lambda*weight2
			//	} else{
			//		weight1 <- 1/(g*p)
			//		weight2 <- matrix(1,nrow=K,ncol=p)
			//		lambda.theta <- lambda*weight2
			//	}
			if(adaptive)  {
				for(r=0; r<K*p; r++) weight2[r] = fabs(mu1[r])+TINY;  // plus TINY to avoid Inf by division
				mmax(weight2, K, p, weight1, 2);
				for(r=0; r<p; r++) weight1[r] = 1.0/p/g[r]/weight1[r];
				for(r=0; r<K*p; r++) weight2[r] = 1.0/weight2[r];
			} 
			for(r=0; r<K*p; r++) lambda_theta[r] = lambda * weight2[r];
			
			// ***  EM algorithm ***
			//	attempt_ID <- 0
			//	diff <- 1
			//	iter <- 0
			nAttempt = 0;
			diff = 1.0;
			iter1 = 0;
			converge1 = 0;

			//	while((diff>=epsilon) & (iter<=max.iter)) {
			//while(diff>=epsilon && iter1<=maxiter) {
			//while (nAttempt > MAXATTEMPT) {
			while (1) {
				//		iter <- iter + 1
				//		mu.old <- mu
				iter1++;
				memcpy((void*)mu_old, (void*)mu1, K*p*sizeof(double));

				// *** E-step ***
				//		kern.mat <-  matrix(0,nrow=K,ncol=n)
				memset((void*)kern_mat, 0, K*n*sizeof(double));

				//		for (k in 1:K) kern.mat[k,] <- -0.5*apply((t(y)-mu[k,])^2/sigma2.vec, 2, sum)
				for(r=0; r<K; r++)  {
					msubv(y, n, p, mu1+r*p, tmp_y, 1);
					dmul(tmp_y, tmp_y, tmp_y, n*p);
					mdivv(tmp_y, n, p, sigma2_vec1, tmp_y, 1);
					msum(tmp_y, n, p, kern_mat+r*n, 1);
				}
				kmul(-0.5, kern_mat, kern_mat, K*n);

				//   	     	for (k in 1:K) {
				//	         		tmp <- exp(t(t(kern.mat)-kern.mat[k,]))*pi.vec
				//		      		z[k,] <- pi.vec[k]/apply(tmp,2,sum)
				//	    		}
				for (r=0; r<K; r++)  {
					msubv(kern_mat, K, n, kern_mat+r*n, tmp_kern_mat, 1);
					for (s=0; s<K*n; s++) tmp_kern_mat[s] = exp(tmp_kern_mat[s]);
					mmulv(tmp_kern_mat, K, n, pi_vec1, tmp_kern_mat, 2);
					msum(tmp_kern_mat, K, n, z+r*n, 2);
				}
				for (r=0; r<K*n; r++) {
					if (z[r] < TINY) 
						z[r] = TINY;
				}
				vdivm(pi_vec1, z, K, n, z, 2);

				// *** M-step ***
				//		pi.vec <- apply(w*t(z),2,sum)
				mmulv(z, K, n, w, tmp_kern_mat, 1);
				msum(tmp_kern_mat, K, n, pi_vec1, 1);
				for (r=0; r<K; r++) {
					if ( pi_vec1[r]<TINY ) 
						pi_vec1[r] = TINY;
				}

				//		sigma2.mat <- matrix(0,nrow=K,ncol=p)
				memset((void*)sigma2_mat, 0, K*p*sizeof(double));
//LOC1			
				memset((void*)mu0, 0, K*p*sizeof(double));
				
				//		for (k in 1:K) sigma2.mat[k,] <- apply(t(t(y)-mu[k,])^2*w*z[k,],2,sum)
				for(r=0; r<K; r++)  {
					msubv(y, n, p, mu1+r*p, tmp_y, 1);
					dmul(tmp_y, tmp_y, tmp_y, n*p);
					dmul(w, z+r*n, tmp_kern_mat, n);
					mmulv(tmp_y, n, p, tmp_kern_mat, tmp_y, 2);
					msum(tmp_y, n, p, sigma2_mat+r*p, 2);
//LOC2
					for(s=0; s<n; s++) tmp_kern_mat[s] /= pi_vec1[r];
					mmulv(y, n, p, tmp_kern_mat, tmp_y, 2);
					msum(tmp_y, n, p, mu0+r*p, 2);
				}
				
				//		sigma2.vec <- apply(sigma2.mat,2,sum)
				msum(sigma2_mat, K, p, sigma2_vec1, 2); 

				//		mu0 <- matrix(0,nrow=K,ncol=p)
				// implemented at LOC1
				
				//		for (k in 1:K) mu0[k,] <- apply(y*w*z[k,]/sum(w*z[k,]),2,sum)
				// implemented at LOC2

				//		gamma <- (apply(theta*mu0*n*apply(w*t(z),2,sum),2,sum)-weight1*sigma2.vec)/apply(theta^2*n*apply(w*t(z),2,sum),2,sum)
				dmul(theta, mu0, tmp_theta, K*p);
				kmul(n, tmp_theta, tmp_theta, K*p);
				mmulv(tmp_theta, K, p, pi_vec1, tmp_theta, 2);
				msum(tmp_theta, K, p, gamma, 2);
				dmul(weight1, sigma2_vec1, tmp_gamma, p);
				msub(gamma, tmp_gamma, gamma, p);
				kmul(1.0/n, gamma, gamma, p);
				dmul(theta, theta, tmp_theta, K*p);
				mmulv(tmp_theta, K, p, pi_vec1, tmp_theta, 2);
				msum(tmp_theta, K, p, tmp_gamma, 2);
				mdivv(gamma, p, 1, tmp_gamma, gamma, 2);
				
				//		gamma[gamma<1e-10] <- 1e-10
				for (r=0; r<p; r++) {
					if (gamma[r]<TINY) 
						gamma[r]=0;
				}
				for (r=0; r<p; r++)  {
					for (s=0; s<K; s++)  {
						t = s*p + r;
						if (theta[t]>TINY || theta[t]<-TINY) { 
							break;
						}
					}
					if (s==K) gamma[r] = 0;
				}

				//    		theta <- t(t(abs(mu0))/gamma)-t(t(lambda.theta)*sigma2.vec/gamma^2)/(n*apply(w*t(z),2,sum))
				for(r=0; r<K*p; r++) tmp_mu[r] = fabs(mu0[r]);
				mdivv(tmp_mu, K, p, gamma, tmp_mu, 1);
				dmul(gamma, gamma, tmp_gamma, p);
				mdivv(sigma2_vec1, p, 1, tmp_gamma, tmp_gamma, 2);
				mmulv(lambda_theta, K, p, tmp_gamma, tmp_theta, 1);
				kmul(1.0/n, tmp_theta, tmp_theta, K*p);
				mdivv(tmp_theta, K, p, pi_vec1, tmp_theta, 2);
				msub(tmp_mu, tmp_theta, theta, K*p);

				//    		theta[theta<1e-10] <- 1e-10
				for (r=0; r<K*p; r++) if (theta[r]<TINY) theta[r]=TINY;
				for(r=0; r<p; r++) {
					if (gamma[r]<TINY) for(s=0; s<K; s++) theta[s*p+r]=0;
				}

				//    		theta <- theta*sign(mu0)
				for(r=0; r<K*p; r++) {
					tmp = mu0[r];
					if (tmp==0) theta[r] = 0;
					else {
						if (tmp<0) theta[r] = -theta[r];
					}
				}

				//    		mu <- as.matrix(t(t(theta)*gamma))
				mmulv(theta, K, p, gamma, mu1, 1);
				//			diff <- max(abs(mu-mu.old))
				msub(mu1, mu_old, tmp_mu, K*p);
				for(r=0; r<K*p; r++) tmp_mu[r] = fabs(tmp_mu[r]);
				mmax(tmp_mu, 1, K*p, &diff, 1);

				if (trace) Rprintf("            Iteration: %d  ; Difference: %E\n", iter1, diff);
				if (diff < epsilon) {
					converge1 = 1;
					break;
				}
				if (iter1 > maxiter)  {
					nAttempt++;
					if (nAttempt > MAXATTEMPT) {
						break;
					}
					Rprintf("            Restart attempt %d\n", nAttempt);
					kmns_status = kmns_lloyd_compact_ret(y, n, p, K, 10*n, 1, mu1, &sumwss1);  			// maxiter for kmeans is 10n, and nstart = 4n
					if (kmns_status == 3)  {
						for(nKmnsAttempt = 0; nKmnsAttempt< MAXKMNSATTEMPT; nKmnsAttempt++)  {
							if(kmns_lloyd_compact_ret(y, n, p, K, 10*n, 1, mu1, &sumwss1)!=3) break;
						}
						if (nKmnsAttempt==MAXKMNSATTEMPT) {
							memcpy((void*)mu1, (void*)best_kmns_mu, K*p*sizeof(double));
						} else {
							if(sumwss1<bestsumwss) {
								bestsumwss = sumwss1;
								memcpy((void*)best_kmns_mu, (void*)mu1, K*p*sizeof(double));
							}
						}
					}

					tmp = 0.0;
					for(r=0; r<K; r++) { pi_vec1[r] = 10*exp_rand(); tmp +=pi_vec1[r]; }
					for(r=0; r<K; r++) pi_vec1[r] /= tmp;
				
					//	z <- matrix(0,nrow=K,ncol=n)
					memset(z, 0, maxK*n*sizeof(double));

					//	sigma2.vec <- rep(1,p)
					for(r=0; r<p; r++) sigma2_vec1[r] = 1.0;

					memcpy((void*)init_mu1, (void*)mu1, K*p*sizeof(double));
					
					
					//	gamma <- rep(1,p)
					//  theta <- mu				
					memcpy((void*)gamma, (void*)sigma2_vec1, p*sizeof(double));
					memcpy((void*)theta, (void*)mu1, K*p*sizeof(double));

					//	if(adaptive){
					//		weight1 <- 1/(g*p*apply(abs(mu),2,max))
					//    	weight2 <- 1/abs(mu)
					//  	lambda.theta <- lambda*weight2
					//	} else{
					//		weight1 <- 1/(g*p)
					//		weight2 <- matrix(1,nrow=K,ncol=p)
					//		lambda.theta <- lambda*weight2
					//	}
					if(adaptive)  {
						for(r=0; r<K*p; r++) weight2[r] = fabs(mu1[r])+TINY;  // plus TINY to avoid Inf by division
						mmax(weight2, K, p, weight1, 2);
						for(r=0; r<p; r++) weight1[r] = 1.0/p/g[r]/weight1[r];
						for(r=0; r<K*p; r++) weight2[r] = 1.0/weight2[r];
					} 
					for(r=0; r<K*p; r++) lambda_theta[r] = lambda * weight2[r];
					
					// ***  EM algorithm ***
					//	attempt_ID <- 0
					//	diff <- 1
					//	iter <- 0
					diff = 1.0;
					iter1 = 0;
				}
			
			}
			if(converge1) { Rprintf("            The algorithm converges.\n"); }
			else { Rprintf("            The algorithm does not converge.\n"); }

			npara = K-1 + p;
			for(r=0; r<K*p; r++) {
				if(mu1[r]<=epsilon && mu1[r]>=-epsilon) {
					mu1[r]=0;
				} else npara++;
			}

			Rprintf("            npara = %d.\n", npara); 

			mmulv(mu1, K, p, norm_y, mu_centering1, 1);

			maddv(mu_centering1, K, p, mean_y, mu1, 1);
			for(r=0; r<p; r++) sigma2_vec1[r] = sigma2_vec1[r]*norm_y[r]*norm_y[r];

			lik1 = 0.0;
			for(r=0; r<p; r++)  lik1 += log(sigma2_vec1[r]);
			lik1 = -0.5*p*log(2*PI) - 0.5*lik1;
//test
			Rprintf("            lik1 = %f.\n", lik1); 
			

				for(r=0; r<K; r++)  {
					msubv(y, n, p, mu1+r*p, tmp_y, 1);
					dmul(tmp_y, tmp_y, tmp_y, n*p);
					mdivv(tmp_y, n, p, sigma2_vec1, tmp_y, 1);
					msum(tmp_y, n, p, kern_mat+r*n, 1);
				}
				kmul(-0.5, kern_mat, kern_mat, K*n);
			
			
			
			memset((void*)kern_mat, 0, K*n*sizeof(double));
			for(r=0; r<K; r++)  {
				msubv(orig_y, n, p, mu1+r*p, tmp_y, 1);
				dmul(tmp_y, tmp_y, tmp_y, n*p);
				mdivv(tmp_y, n, p, sigma2_vec1, tmp_y, 1);
				msum(tmp_y, n, p, kern_mat+r*n, 1);
			}
			kmul(-0.5, kern_mat, kern_mat, K*n);
			mmax(kern_mat, K, n, max_kern_mat, 2);
			msubv(kern_mat, K, n, max_kern_mat, kern_mat, 1);
			for(r=0; r<K*n; r++) kern_mat[r] = exp(kern_mat[r]);
			mmulv(kern_mat, K, n, pi_vec1, kern_mat, 2);
			msum(kern_mat, K, n, kern_mat, 2);
			lik2 = 0.0;
			for(r=0; r<n; r++) lik2 += w[r]*(log(kern_mat[r])+max_kern_mat[r]);
//test
			Rprintf("            lik2 = %f.\n", lik2); 
			tmp = -2*n*(lik1+lik2)+log(n)*npara;
			BIC[i*nLambda+j] =  tmp;
//test
			Rprintf("            BIC = %f.\n", tmp); 
			if (tmp < optBIC)  {
				optBIC = tmp;
				converge = converge1;
				optK = K;
				optLambda = lambda;
				memcpy((void*)pi_vec, (void*)pi_vec1, K*sizeof(double));
				memcpy((void*)mu, (void*)mu1, K*p*sizeof(double));
				memcpy((void*)mu_centering, (void*)mu_centering1, K*p*sizeof(double));
				memcpy((void*)init_mu, (void*)init_mu1, K*p*sizeof(double));
				memcpy((void*)sigma2_vec, (void*)sigma2_vec1, p*sizeof(double));
				iter = iter1;
				for(r=0; r<n; r++) {
					ClusterID[r] = 1;
					for(s=1; s<K; s++)  {
						if (z[s*n+r]>z[r]) {
							ClusterID[r] = s+1;
							z[r] = z[s*n+r];
						}
					}
				}
			}
		}
	}

	*piter = iter;
	*poptBIC = optBIC;
	*poptK = optK;
	*poptLambda = optLambda;
	*pconverge = converge;

	if (pi_vec1!= NULL) free(pi_vec1);
	if (mu1 != NULL) free(mu1);
	if (mu_centering1 != NULL) free(mu_centering1);
	if (init_mu1 != NULL) free(init_mu1);
	if (sigma2_vec1 != NULL) free(sigma2_vec1);
	if (z != NULL) free(z);
	if (gamma != NULL) free(gamma);
	if (theta != NULL) free(theta);
	if (weight1 != NULL) free(weight1);
	if (weight2 != NULL) free(weight2);
	if (lambda_theta != NULL) free(lambda_theta);
	if (mu_old != NULL) free(mu_old);
	if (kern_mat != NULL) free(kern_mat);
	if (tmp_y != NULL) free(tmp_y);
	if (tmp_kern_mat != NULL) free(tmp_kern_mat);
	if (sigma2_mat != NULL) free(sigma2_mat);
	if (mu0 != NULL) free(mu0);
	if (tmp_theta != NULL) free(tmp_theta);
	if (tmp_gamma != NULL) free(tmp_gamma);
	if (tmp_mu != NULL) free(tmp_mu);
	if (best_kmns_mu!=NULL) free(best_kmns_mu);
	if (max_kern_mat!=NULL) free(max_kern_mat);
}
