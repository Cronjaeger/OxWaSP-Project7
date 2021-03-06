#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* Begin constants */
#define VERBOSE 0 //when set to 1, running smc_sampler(...) will generate A LOT of output.
#define M 200
#define MH_STEPS 10
#define N_MIXTURE_COMPONENTS 4
#define SIGMA 0.55
#define HLINE printf("\n--------------------------------------------------------------------------------\n")

const double MIXTURE_WEIGHTS[N_MIXTURE_COMPONENTS] = {0.25,0.25,0.25,0.25};
static double* xNew; //an array used for swapping

/* begin function-signatures */
int test(void);
void allocateParticle_swap(void);
double dnorm(const double x,const double mu);
double log_likelihood(const double * y, const int N_y, const double* mu);
_Bool ESS_is_largeEnough(const double * W,const int N_W);
double pi_updateStep(const double* y,const int N_y, const double* x, const int n);
void smc_sampler(const double* yObs,const int N_yObs ,const int n_particles , double** particles, double* weights);
void randParticle(double * x);
void resample(double ** particles, double * weights,const int n_particles, const double w0, gsl_rng * r, gsl_ran_discrete_t * wSampler);
void propagate_particle(double * x,const int n, const double* yObs,const int N_yObs,gsl_rng * r);
void smc_sampler_for_R( double* yObs, int* N_yObs, int* N_particles, double* X_vec, double* W, double* time);
void freeParticle_swap(void);


int main(){
  allocateParticle_swap();
  return test();
  freeParticle_swap();
}

/* a wrapper function for calling smc_sampler from R */
void smc_sampler_for_R(
  double* yObs,
  int* N_yObs,
  int* N_particles,
  double* X_vec,
  double* W,
  double* t_SMC )
  {

  // allocate memmory for swapping
  allocateParticle_swap();

  // Generate an (N_particles x N_MIXTURE_COMPONENTS) array used to
  // internally represent particles
  double** X;
  X = (double**) malloc(sizeof(double *) * (*N_particles));
  size_t i = 0;
  for(i=0 ; i < *N_particles ; ++i){
    X[i] = (double*) malloc(sizeof(double)*N_MIXTURE_COMPONENTS);
  }

  //start timer
  time_t t0 = time(NULL);

  //run sampler
  smc_sampler(yObs,*N_yObs,*N_particles,X,W);

  //stop timer
  time_t t1 = time(NULL);

  //return runtime in secconds
  *t_SMC = difftime(t1,t0);

  //write to output vector
  size_t j;
  for( i = 0 ; i < (*N_particles) ; ++i){
    for( j = 0 ; j < N_MIXTURE_COMPONENTS ;  ++j){
      X_vec[N_MIXTURE_COMPONENTS * i + j] = X[i][j];
    }
  }

  //free up memory;
  for(i=0 ; i < *N_particles ; ++i){
    free(X[i]);
  }
  free(X);

  freeParticle_swap();

  /* //for debugging.
  int i;
  for(i=0 ; i < *N_particles ; ++i){
    X_vec[i] = (double) i;
  }
  */
}

/* Main function for running SMC sampler
 yObs is an a vector of N_yObs observations

 n_particles is the desired number of particles to be simulated

 particles should be a double[n_particles][N_MIXTURE_COMPONENTS] array.
 It will be used for storing results and doing computaations.

 weights should be a double[n_particles] array.
 It will be used for storing results and doing computations.
*/
void smc_sampler(
  const double* yObs,
  const int N_yObs,
  const int n_particles,
  double** particles,
  double* weights)
  {

  time_t t0 = time(NULL);
  time_t t1;
  size_t i; //summation indices

  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_ran_discrete_t * wSampler;

  if(VERBOSE) printf("Initializinng...");
  /* Initialize  weights and particles*/
  double w0 = 1 / (double) n_particles; // default weight
  for(i=0; i < n_particles ; ++i){
    weights[i] = w0;
    randParticle(particles[i]);
  }
  if(VERBOSE) printf("DONE!\n");

  if(VERBOSE) time(&t1);
  if(VERBOSE) printf("Elapsed time = %.0f sec\n",difftime(t1,t0));

  if(VERBOSE) printf("Running sampling steps from 1 to 200...\n");

  /* Iterate from n=1 to n=M */
  int n = 1;
  unsigned int resampleCounter = 0;
  while(1){

    // for diagnosing running time and output
    if((n-1)%10 == 0){
      if(VERBOSE) time(&t1);
      if(VERBOSE) printf("Running steps %i...%i; \t",n,n+9);
      if(VERBOSE) printf("elapsed time so far = %.0f sec\n",difftime(t1,t0));

      if(VERBOSE){
        printf("Printing the first 10 particles to check validity of output:\n");
        size_t j;
        for(i=0;i<10;++i){
          printf("  x[%i]_(%i) = ",(int) i,n);
          for(j=0; j<N_MIXTURE_COMPONENTS ; ++j){
            if(particles[i][j] >= 0) printf(" ");
            if(fabs(particles[i][j]) < 10) printf(" ");
            printf("%.3f\t",particles[i][j]);
          }
          printf(" \tW[%i](%i) = %f",(int) i , n, weights[i]);
          printf("\n");
        }
      }

      double sumW = 0;
      for(i = 0 ; i < n_particles ; ++i){
        sumW += weights[i];
      }
      if(VERBOSE) printf("sum of weights = %f",sumW);
      if(VERBOSE) HLINE;

     }

    /* Resample and reset weights, if nessecary */
    if(!(ESS_is_largeEnough(weights,n_particles))){
      if(VERBOSE) ++resampleCounter;
      if(VERBOSE) HLINE;
      if(VERBOSE) printf("resampling in step n = %i\n",n);
      if(VERBOSE) HLINE;
      //Preprocessing for sampling from a discrete distribution
      wSampler = gsl_ran_discrete_preproc (n_particles,weights);

      //Resample (AND reset) weights
      resample(particles,weights,n_particles,w0,r,wSampler);
    }

    /* Update particle weights */
    double sumW = 0;
    for(i = 0 ; i<n_particles ; ++i){
      weights[i] *= pi_updateStep(yObs, N_yObs, particles[i] ,n);
      sumW += weights[i];
    }
    double Z = (double) 1 / sumW;
    for(i = 0 ; i<n_particles ; ++i){
      weights[i] *= Z;
    }



    /* break out when n == M */
    if(n==M) break;


    /* Propagate particles */
    for(i=0; i<n_particles ; ++i){
      propagate_particle(particles[i], n, yObs, N_yObs,r);
    }

    /*increment counter*/
    ++n;
  }
  if(VERBOSE) printf("total number of resampling-steps = %i\n",resampleCounter);
  gsl_rng_free(r);
  // gsl_ran_discrete_free(wSampler);
}






/**********************
* Auxiliary Functions *
**********************/

void allocateParticle_swap(void){
  //modify in case of parralelization
  xNew = (double*) malloc(sizeof(double)*N_MIXTURE_COMPONENTS);
}

void freeParticle_swap(void){
  //modify in case of parralelization
  free(xNew);
}

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
      if(VERBOSE) printf("WOAH! particle out of range\n");
      return -INFINITY;
    }
  }
  /* Compute the log-likelihood iteratively */
  double sum; // inner summation (over mu)
  double Sum = 0; //outer summation (over y)
  for(i = 0;i < N_y; ++i){
    sum = 0;
    for(j = 0 ; j < N_MIXTURE_COMPONENTS ; ++j){
      sum += MIXTURE_WEIGHTS[j] * dnorm(y[i], mu[j]);
    }
    Sum += log(sum);
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

/*
A n auxiliary function for computing pi_n(x_(n-1)) / pi_(n-1)(x_(n-1))
\eqn{\propto} p( x_(n-1) | y ) ^ ((2n -1)/M^2). Used when updating weights.
@param y Observed data
@param x Our proposed value for Mu (a vector of modes)
*/
double pi_updateStep(const double* y,const int N_y, const double* x,const int n){
  return exp( ( ((float) 2*n -1)/((float) (M*M)) ) * log_likelihood(y,N_y,x));
}

void propagate_particle(double * x,const int n, const double* yObs,const int N_yObs,gsl_rng * r){

  //NOTE: xNew MUST have been allocated before calling this method!

  size_t i = 0; size_t j;
  for(i = 0 ; i < MH_STEPS ; ++i){

    // Generate new observation
    for(j = 0 ; j < N_MIXTURE_COMPONENTS ; ++j){
      xNew[j] = x[j] + gsl_ran_gaussian(r,SIGMA);
    }

/*
    // REMOVE WHEN COMPILING FOR SPEED!
    if(VERBOSE){
      // printf("Acceptance! alpha = %0.3f\n",alpha);
      printf("  x\t= ");
      for(j=0; j<N_MIXTURE_COMPONENTS ; ++j){
        if(x[j] >= 0) printf(" ");
        if(fabs(x[j]) < 10) printf(" ");
        printf("%.3f\t",x[j]);
      }
      printf("\n  x_new\t=");
      for(j=0; j<N_MIXTURE_COMPONENTS ; ++j){
        if(xNew[j] >= 0) printf(" ");
        if(fabs(xNew[j]) < 10) printf(" ");
        printf("%.3f\t",xNew[j]);
      }
      printf("\n");
    }
*/

    // double logL_new = log_likelihood(yObs, N_yObs, xNew);
    // double logL_old = log_likelihood(yObs, N_yObs, x);
    double logL_Diff = log_likelihood(yObs, N_yObs, xNew) - log_likelihood(yObs, N_yObs, x);
    double alpha = 1;
    _Bool autoAccept = 1;

    // Iff logL_Diff < 0 true, the sign of the argument fo exp() will be <1.
    if(logL_Diff < 0){
      double exponent = (double) n / (double) M;
      exponent *= exponent; // square the exponent
      alpha = exp(exponent * logL_Diff );
      autoAccept = 0;
      if(VERBOSE) printf("alpha = %.4f",alpha);
    }

    if( autoAccept || ((double) rand() / (double) RAND_MAX) < alpha){
      if( VERBOSE && !autoAccept ) printf("\t new point accepted anyway\n");
      for(j = 0 ; j < N_MIXTURE_COMPONENTS ; ++j){
        x[j] = xNew[j];
      }
    } else {
      if(VERBOSE) printf("\n");
    }
  }
  if(VERBOSE) printf("\n");
}

/* A generic function to test during development
Returns 0 if test succcessfull*/
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

/*
  double mu[4] = {-3,0,3,6};
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

/*
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
*/
/*
  int N = 1000;

  double* W;
  W = (double*) malloc(sizeof(double)*N);

  double** X;
  X = (double**) malloc(sizeof(double *) * N);
  size_t i = 0;
  for(i=0 ; i < N ; ++i){
    X[i] = (double*) malloc(sizeof(double*)*N_MIXTURE_COMPONENTS);
  }

  // double W;
  // double* X[N][N_MIXTURE_COMPONENTS];

  // data pasted in from R
  // generated using sampleMM(100) from mixtureModel_SMC.R
  int N_yObs = 100;
  double yObs[100] =
   {0.270469630282116 , 0.0948544266135598 , 4.82353957890236 ,
    2.71176569013062 , -2.10140287048951 , -2.4639397242585 ,
    -2.47096320311181 , 5.59163596897931 , -3.50775981783882 ,
    3.26972797997964 , 3.32256384988976 , 5.73616849148524 ,
    5.74816296903159 , -2.97620992270215 , -0.211030444261816 ,
    5.56902016608288 , -3.69397287391297 , 5.71718272922272 ,
    -0.375056085230512 , 0.290727561775271 , 3.53116678441137 ,
    -2.92396872511173 , 5.44623183851942 , 5.60306594746839 ,
    5.56855197639855 , 2.28105552916777 , 2.71638272346431 ,
    3.7321825583144 , 6.10978497180474 , 6.7089442115652 ,
    3.14182326021059 , -2.72341774114777 , -0.0919451805310823 ,
    -3.05334875104901 , 6.00697089830724 , -2.81827687988051 ,
    6.03028599825704 , 3.77081092101111 , -0.17554472859167 ,
    5.57845673340545 , 6.25654424410659 , 6.40046681340197 ,
    -2.39322730700275 , 0.937990298556008 , -3.61930611127869 ,
    2.99419558349795 , -2.95427533963002 , 2.53999917081115 ,
    -4.35170347149833 , 2.56389630336686 , 2.78219272857881 ,
    -2.53796084626858 , 2.75783249644654 , 0.265653681397089 ,
    -4.08489661583957 , 0.141285949971235 , -3.09744117607388 ,
    -3.29683116823359 , -0.56348650880506 , 2.13298104427159 ,
    5.81313774243827 , 2.90131666078861 , -0.870942485382193 ,
    -0.653182465687133 , 3.24697479240955 , 3.46323186816872 ,
    3.70853665171654 , -3.79582094987125 , 6.55169236217858 ,
    -3.63836564940765 , 0.132152004366257 , 5.24827155186988 ,
    2.53806422373366 , 5.6550753862474 , 2.42737960470922 ,
    5.88968129179486 , 6.83133504027115 , -0.473386824297909 ,
    -0.883286661476473 , 6.94044430889499 , 2.73913195305667 ,
    3.02781182874733 , 6.53392928375413 , 6.53623960926176 ,
    6.44586265929949 , 0.737757409838858 , -2.81327223680057 ,
    0.534431297323472 , -3.28504638650391 , -0.0839400339092673 ,
    -3.30271081947717 , 5.97247453668681 , -3.55468774001656 ,
    -0.211152666301356 , 6.12181217531004 , 2.56643856048767 ,
    0.200557127297898 , -3.97176332345347 , 0.0176492920069539 ,
    2.87381823065611};

  printf("Running SMC-sampler with %i particles, based on %i observations.\n",N,N_yObs);
  time_t t1 = time(NULL);
  smc_sampler(yObs,N_yObs,N,X,W);
  time_t t2 = time(NULL);
  printf("DONE!\nTotal elapsed time = %.0f\n",difftime(t2,t1));

  HLINE;

  printf("Printing the first 10 particles to check validity of output:\n");
  size_t j;
  for(i=0;i<10;++i){
    printf("  x[%i] = ",(int) i);
    for(j=0; j<N_MIXTURE_COMPONENTS ; ++j){
      if(X[i][j] >= 0) printf(" ");
      if(fabs(X[i][j]) < 10) printf(" ");
      printf("%.3f\t",X[i][j]);
    }
    printf("\n");
  }

  free(W);
  // free(X)
  for(i=0 ; i<N; ++i){
    free(X[i]);
  }
  free(X);
*/

/*
  int N_yObs = 100;

  //DATA from R, so that the S-method and the .c method can be directly compated
  double yObs[100] = {6.42614409154078,0.260646235868524,6.20077626783781,-3.54124388162526,3.16810352365003,-3.02631343546788,6.22743478871312,6.8131742684717,6.24776924281844,3.65204896173014,6.01102502511306,-2.78204404408537,6.54496067301763,2.86661732139087,6.41368106980899,-2.54549938986914,6.08430821908108,2.65819528014716,-4.10369954335742,-3.70836228135375,3.55973960494948,-3.99408019214758,-3.06751044524848,-0.560058224387682,0.625592793169344,2.58052593365431,0.689049384542289,0.494930560119956,0.390397430673082,2.87213201595585,6.08140477180053,-3.63073724259346,1.76368358083185,0.0978387885170348,6.49576360025151,6.65592881975105,-2.31116539097971,6.69822824187681,0.0197128170022586,-3.6584094557839,6.37219590110012,5.98608504754753,5.76477302651333,-2.49473884182791,3.68861876383596,6.71753596416937,-2.86065603329243,2.36234293889529,-2.13229209410826,-3.1821797074668,3.67640450999392,2.35200323180728,-2.90422204445312,3.85237110998024,5.73406945312462,3.61994874049274,-4.0561902744938,5.66696946993138,5.96245114627653,6.59816990256299,0.231879318383607,2.34472949459441,2.72937582033292,-0.444334698629388,2.22482121394123,3.15359838488669,5.79310254628412,2.95627482456672,0.0754282164405023,-3.27128946837834,-2.77302236660845,7.25710145940199,6.2324724323204,5.9296066592542,5.60999425487975,-0.286465744070193,0.587073138580418,-2.67791031284446,6.98129381109518,3.36137498170339,-3.16490919686681,-3.62156309636938,5.41982668735859,-3.07991914749654,-2.89255426727392,5.61010574600101,-3.10612256078391,3.48991878214866,2.58896632171936,0.703862729104843,4.72001741709978,-0.244802119787893,6.45395925612269,0.253461730596016,2.69213945052578,-2.90338536462303,3.28208239426602,-2.57167632631676,0.920364746339414,0.993277872559861};

  int N_X = 5;
  double X[5][4] = {
    {0,3,-3,6},
    {-3,0,3,6},
    {0,0,0,0},
    {9,9,9,9},
    {0,6,3,15}
  };

  size_t i;
  for(i = 0 ; i<N_X ; ++i){
    double LL = log_likelihood(yObs,N_yObs,X[i]);
    printf("log_L(X[%i]) = %f\n",(int) i,LL);
  }
*/
  return 0;
}
