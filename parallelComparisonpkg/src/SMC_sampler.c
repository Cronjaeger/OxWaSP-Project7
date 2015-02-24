#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//#include "square.c"

/* Begin constants */
const double SIGMA = 0.55;
const double MIXTURE_WEIGHTS[4] = {0.25,0.25,0.25,0.25};
const double MU_TRUE[4] = {-3,0,3,6};
const int N_MIXTURE_COMPONENTS = 4; // should be no larger than an unsigned char
const int M = 200;
const int MH_STEPS = 10;


/* begin function-signatures */
int test(void);
double dnorm(const double x,const double mu);
double log_likelihood(const double * y, const int N_y, const double* mu);
_Bool ESS_is_largeEnough(const double * W,const int N_W);
double pi_updateStep(const double* y,const int N_y, const double* x, const int n);
void smc_sampler(const double* y_obs,const int N_y_obs ,const int n_particles , double** particles, double* weights);
void randParticle(double * x);
//void resample(double* q_weights, double** particles);
void resample(double ** particles, double * weights,const int n_particles, const double w0, gsl_rng * r, gsl_ran_discrete_t * wSampler);
void propagate_particle(double * x,const int n, const double* y_obs,const int N_y_obs,gsl_rng * r);

int main(){
  return test();
}

/* Main function for running SMC sampler
 y_obs is an a vector of N_y_obs observations

 n_particles is the desired number of particles to be simulated

 particles should be a double[n_particles][N_MIXTURE_COMPONENTS] array.
 It will be used for storing results and doing computaations.

 weights should be a double[n_particles] array.
 It will be used for storing results and doing computations.
*/
void smc_sampler(
  const double* y_obs,
  const int N_y_obs,
  const int n_particles,
  double** particles,
  double* weights)
  {

  // /* An auxiliary array for containing the quantile of the weight-vector */
  // double* q_weights;
  // q_weights = (double *)malloc(sizeof(double)*n_particles);

  size_t i; //size_t j; //summation indices

  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_ran_discrete_t * wSampler;

  /* Initialize  weights and particles*/
  double w0 = 1 / (double) n_particles; // default weight
  for(i=0; i < n_particles ; ++i){
    weights[i] = w0;
    randParticle(particles[i]);
  }

  /* Iterate from n=1 to n=M */
  int n = 1;
  while(1){

    /* Resample and reset weights, if nessecary */
    if(!(ESS_is_largeEnough(weights,n_particles))){

      //Preprocessing for sampling from a discrete distribution
      wSampler = gsl_ran_discrete_preproc (n_particles,weights);

      //Resample (AND reset) weights
      resample(particles,weights,n_particles,w0,r,wSampler);
    }


    /* Update particle weights */
    for(i = 0 ; i<n_particles ; ++i){
      weights[i] *= pi_updateStep(y_obs, N_y_obs, particles[i] ,n);
    }


    /* break out when n == M */
    if(n==M) break;


    /* Propagate particles */
    for(i=0; i<n_particles ; ++i){
      propagate_particle(particles[i], n, y_obs, N_y_obs,r);
    }

    /*increment counter*/
    ++n;
  }
  // free(q_weights);
  gsl_rng_free(r);
  gsl_ran_discrete_free(wSampler);
}






/**********************
* Auxiliary Functions *
**********************/
double dnorm(const double x, const double mu){
  return  gsl_sf_erf_Z((x - mu)/SIGMA)/SIGMA;
}

/* returns log( p( y | mu , sigma , weights ) )
 y should be a vector of doubles corresponding to observations
 N_y corresponds to the length of y
 mu coresponds to a vector of modes of length N_MIXTURE_COMPONENTS */
double log_likelihood(const double* y, const int N_y,const double* mu){
  int i; int j; // initialize summation-indices

  /* Verify that all modes are in [-10,10] */
  for(i=0; i<N_MIXTURE_COMPONENTS; ++i){
    //INFINITY is a macro defined in math.h
    if( mu[i] > 10.0 || mu[i] < -10.0 ){
      return -INFINITY;
      printf("WOAH!\n");
    }
  }
  /* Compute the log-likelihood iteratively */
  double sum; // inner summation (over mu)
  double Sum = 0; //outer summation (over y)
  for(i = 0;i < N_y; ++i){
    sum = 0;
    for(j = 0 ; j < N_MIXTURE_COMPONENTS ; ++j){
      sum += MIXTURE_WEIGHTS[j] * dnorm(y[i], mu[j]);
//      printf("%5sum=%.15f\n",(double) sum);
    }
    Sum += log(sum);
//    printf("Sum=%.15f\nsum=%.15f\n",(double) Sum, (double) sum);
  }
  return Sum;
}



/* Generates a ranodm particle ie. 4 ind. samples uniformly from -10 to 10 */
/* Results are stored in x */
void randParticle(double* x){
  unsigned char i;
  for(i=0 ; i<N_MIXTURE_COMPONENTS; ++i){
    x[i] = (double) -10 + (double) 20 * (double) rand() / (double) RAND_MAX;
  }
}

_Bool ESS_is_largeEnough(const double * W,const int N_W){
  double sumSq = 0;
  int i;
  for(i = 0 ; i < N_W ; ++i){
    sumSq += W[i] * W[i];
  }
  return (_Bool) 1 >= sumSq * N_W / 2.0;
}

void resample(
  double ** particles,
  double * weights,
  const int n_particles,
  const double w0,
  gsl_rng * r,
  gsl_ran_discrete_t * wSampler)
  {
  size_t i; size_t j;
  for(i=0; i< n_particles; ++i){
    //iterate over particles

    size_t i_new = gsl_ran_discrete(r,wSampler);
    // copy the particle in question
    for(j = 0 ; j < N_MIXTURE_COMPONENTS ; ++j){
      particles[i][j] = particles[i_new][j];
    }
  }
  //reset particle weights
  for(i=0; i < n_particles ; ++i){
    weights[i] = w0;
  }
}

// an auxiliary function for computing pi_n(x_(n-1)) / pi_(n-1)(x_(n-1))
// \eqn{\propto} p( x_(n-1) | y ) ^ ((2n -1)/M^2). Used when updating weights.
// @param y Observed data
// @param x Our proposed value for Mu (a vector of modes)
double pi_updateStep(const double* y,const int N_y, const double* x,const int n){
  return exp( ( ((float) 2*n -1)/((float) (M*M)) ) * log_likelihood(y,N_y,x));
}

void propagate_particle(double * x,const int n, const double* y_obs,const int N_y_obs,gsl_rng * r){
  size_t i = 0; size_t j;
  double xNew[N_MIXTURE_COMPONENTS]; //requires -std=c99 or later. otherwise uncomment the code below.
  // double* xNew;
  // xNew = (double *)malloc(sizeof(double) * N_MIXTURE_COMPONENTS);

  for(i = 0 ; i < MH_STEPS ; ++i){

    // Generate new observation
    for(j = 0 ; j < N_MIXTURE_COMPONENTS ; ++j){
      xNew[j] = x[j] + gsl_ran_gaussian(r,SIGMA);
    }

    double alpha
    = fmin(
        1,
        exp(
          (n/M)*(n/M)
          *
            (log_likelihood(y_obs,N_y_obs,xNew)
             - log_likelihood(y_obs,N_y_obs,x))
          )
        );

    if( ((double) rand() / (double) RAND_MAX) < alpha){
      for(j = 0 ; j < N_MIXTURE_COMPONENTS ; ++j){
        x[j] = xNew[j];
      }
    }
  }
}

/* A generic function to test during development
Return 0 if succcessfull*/
int test(void){
/*
  double x = 1; double y;
  double mu = 0; double sigma = 5;
  y = dnorm(x,mu,sigma);
  printf("Hello World!\nNormal(1;0,5) = %.10f\n",y);
*/
/*
  double y[2] = {0,2};
  double mu[4] = {-3,0,3,6};
  double x = log_likelihood(y,2,mu);
  printf("Hello World!\nx = %.10f\n",x);
*/
/*  double mu[4] = {-3,0,3,6};
  if( mu[2] < 10 ) printf("Hello World!\n");
*/

/*
  int N_W = 6;
  double W1[6] = {1,0,0,0,0,0};
  double W2[6] = {0.2,0.2,0.2,0.2,0.2,0};

  if(!(ESS_is_largeEnough(W1,N_W))) printf("ESS1 is too small!\n");
  if(ESS_is_largeEnough(W2,N_W)) printf("ESS2 on the other hand is groovy!\n");
*/

/*
  double y[2] = {0,2};
  double mu[4] = {-3,0,3,6};
  double x = pi_updateStep(y,2,mu,180);
  printf("Hello World!\nx = %.10f\n",x);
*/

  int i; int j;
  double particles[10][4];
  for(i = 0 ; i<10; ++i){
    randParticle(particles[i]);
  }

  for(i=0 ; i<10 ; ++i){
    printf("rand particles = (");
    for(j=0 ; j< 4 ; ++j){
      printf(" %.2f ",particles[i][j]);
    }
    printf( ")\n");
  }

  return 0;
}
