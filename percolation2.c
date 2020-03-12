/*=========================================================/
/  C code (OpenMP) for PX425 assignment 3 2019             /
/  Generates random graphs and invokes analysis routines   /
/  to determine the largest cluster of connected vertices  /
/  and the overall number of clusters.                     /
/                                                          /
/  Original code created by D.Quigley - November 2014      /
/  Updated by Nicholas Hine - October 2019                 /
/=========================================================*/
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <time.h>
#include "mt19937ar.h"

/* Function prototypes for cluster analysis routines */
void find_clusters_recursive(long Nvert,long Maxcon,long *Ncon,long *Lcon,long *lclus, long *nclus);
void find_clusters_eqclass  (long Nvert,long Maxcon,long *Ncon,long *Lcon,long *lclus, long *nclus);

int main () {

  /* Number of vertices */
  const long Nvert = 5000;

  /* Numnber of graphs to average over at each Pcon */
  long Ngraphs = 20;

  /* Number of Pcon values to test */
  long Np = 8;
  
  /* Step incrementing probability of connecting any two vertices */
  const double Pcon_step = 0.0004;

  /* Current value of probability of connecting any two vertices */
  /* Range is from Pcon_step to Pcon_step*Np */
  double Pcon;

  /* Maximum number of edges per vertex */
  const long Maxcon = 80;

  /* Array holding the number of edges at each vertex */
  long *Ncon;

  /* Array holding a list of connections from each vertex ie if 
     Lcon[i*Maxcon+k] = j then the kth connection from vertex
     i is to vertex j and we count from the zeroth connection */
  long *Lcon;

  /* Size of largest cluster, and number of clusters */
  long lclus,nclus,avlclus,avnclus;

  /* Initial and final times */
  struct timespec ti,tf;
  double t_start,t_end;

  /* Random number */
  double xi;

  /* Time taken to generate statistics   */
  double gentime;

  /* Loop counters */
  long i,j,igraph,ip;

  /*------------------------------------/
  / Initialise random number generator  /
  /------------------------------------*/
  unsigned long seed=2383791387;
  unsigned long myseed;
#pragma omp parallel default(shared) private(myseed)
  {
    myseed = seed;
#ifdef _OPENMP    
    /* Modify seed to be unique on each thread */
    myseed = myseed*(omp_get_thread_num()+1);
#endif
    init_genrand(myseed);  
  }

  /*------------/
  / Start timer /
  /------------*/
#ifndef _OPENMP    
  /* Use clock_gettime */
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&ti);
  t_start = (double)ti.tv_sec + ti.tv_nsec/1.0E9;
#endif
#ifdef _OPENMP    
  /* Use OpenMP clock */
  t_start = omp_get_wtime();
#endif

 

  /*--------------------------------/
  / Loop over Pcon values           /
  /--------------------------------*/
#pragma omp parallel for  private(i, ip, igraph, j, Ncon, Lcon, Pcon, xi, avnclus, avlclus) shared(nclus, lclus) ordered
for (ip=0;ip<Np;ip++) {

    /* Compute Pcon from ip */
    Pcon = (double)(ip+1)*Pcon_step;

    /*----------------------/
    / Initialise averages   /
    /----------------------*/
    avlclus = 0.0;
    avnclus = 0.0;

    /*------------------/
    / Loop over graphs  /
    /------------------*/
    for (igraph=0;igraph<Ngraphs;igraph++) {

      /*--------------------------------------------/
      / Allocate memory to hold graph connectivity  /
      /--------------------------------------------*/
      Ncon = (long *)malloc(Nvert*sizeof(long));
      if (Ncon==NULL) { printf("Error allocating Ncon array\n") ; exit(EXIT_FAILURE); }
      
      Lcon = (long *)malloc(Nvert*Maxcon*sizeof(long));
      if (Lcon==NULL) { printf("Error allocating Lcon array\n") ; exit(EXIT_FAILURE); }
      
      /*-----------------------------/
      / Generate a new random graph  /
      /-----------------------------*/
      for (i=0;i<Nvert;i++)        { Ncon[i] =  0; } /* Initialise num. connections */
      for (i=0;i<Nvert*Maxcon;i++) { Lcon[i] = -1; } /* Initialise connection list  */


      for (i=0;i<Nvert;i++) {         /* Loop over vertices i */
	    for (j=i+1;j<Nvert;j++) {         /* Loop over other vertices j */
	        xi = genrand();                  /* Generate a random number xi */
	        if ( xi < Pcon ) {               /* Randomly choose to form an edge */

	        Ncon[i] = Ncon[i] + 1;          /* Increment edges involving i */
	        Ncon[j] = Ncon[j] + 1;          /* Increment edges involving j */

	        /* Check that we will not overrun the Lcon array */
	        if ( ( Ncon[i] > Maxcon-1 ) || ( Ncon[j] > Maxcon-1 ) ) {
		        printf("Error generating random graph.\n");
		        printf("Maximum number of edges per vertex exceeded!\n");
		        exit(EXIT_FAILURE);
	        }

	        Lcon[Maxcon*i+Ncon[i]-1] = j;               /* j is connected to i */
	        Lcon[Maxcon*j+Ncon[j]-1] = i;               /* i is connected to j */

	        } /* if  */
	    } /* j   */
      } /* i   */

      

      /*--------------------------------------------------------/
      / Identify the clusters with a recusive search over edges /
      /--------------------------------------------------------*/
      /* find_clusters_recursive(Nvert,Maxcon,Ncon,Lcon,&lclus,&nclus);  */

      /*----------------------------------------------------------/
      / Identify clusters via a equivalence class algorithm       /
      /----------------------------------------------------------*/
      find_clusters_eqclass(Nvert,Maxcon,Ncon,Lcon,&lclus,&nclus);  

      /* Accumulate averages */
      avlclus += (double)lclus;
      avnclus += (double)nclus;

      /*----------------/
      / Release memory  /
      /----------------*/
      free(Lcon);
      free(Ncon);


    } /* igraph */
    #pragma omp ordered
    printf("Pcon = %12.4f Av. Num. Clusters. = %12.4f Av. Largest Cluster = %12.4f\n",
	   Pcon,avnclus/(double)Ngraphs,avlclus/(double)Ngraphs);

  } /* end loop over ip */

  
  /*------------------------/
  / Check timer             /
  /------------------------*/
#ifndef _OPENMP    
  /* Use clock_gettime */
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&tf);
  t_end =  (double)tf.tv_sec + tf.tv_nsec/1.0E9;
#endif
#ifdef _OPENMP    
  /* Use OpenMP clock */
  t_end = omp_get_wtime();
#endif

  gentime = t_end-t_start;
  printf("Time elapsed : %12.4f seconds\n",gentime);
  
  exit(EXIT_SUCCESS);

}
