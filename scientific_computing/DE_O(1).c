#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void read_input(double *C, double *L, int *nx, double *t_F,double *t_out);

double S_func(double x,double L) {
  if(x<L/4||x>L/2) {
    return 0.0;
  } else {
    return 1.0;
  }
}

int main(void) {
  // **********
  // Parameters
  // **********

  // Number of grid points                  
  int nx;
  // Length of domain
  double L;
  // Equation coefficients
  double C;
  //  Length of time to run simulation.
  double t_F;
  // How frequently in time to output.
  double output_timestep;
  // Read in from file; 
  read_input(&C, &L, &nx, &t_F, &output_timestep);
  // Grid spacing
  double dx = L/(nx-1);		      
  // Time step
  double dt = 0.5*dx/C;
  // ************
  // Grid Storage 
  // ************
  double *y, *y_next;  //y at current and next timestep

  /* Allocate memory according to size of nx */
  y       = malloc(nx*sizeof(double));
  y_next  = malloc(nx*sizeof(double));
  if (y==NULL||y_next==NULL) {
    printf("Memory allocation failed\n");
    return 1;
  }
  
  int j;
  double x;

  // **************
  // initialisation 
  // **************
  for(j=0;j<nx;j++) {
    x = j*dx;
    y[j]  = 0.0;
  }
  
  // Output at start of simulation.
  double ctime = 0.0;
  for (j=0; j<nx; j++ ) {
    x = j*dx;
    printf("%g %g %g\n",ctime,x,y[j]);
  }
  double next_output_time = output_timestep;


  //loop over timesteps 
  while (ctime<t_F){
    double dt0 = dt;
    int output = 0;
    // If we would go past the next output step, reduce the timestep.
    if (ctime+dt0>next_output_time) {
      dt0 = next_output_time - ctime;
      output = 1;
    }
    
    //loop over points 
    for (j=0; j<nx; j++) {
      x = j*dx;
      int jm = (j-1+(nx-1))%(nx-1);
      // Need upwinding for stability.
      double dydx = (y[j] - y[jm])/dx;
      y_next[j] = y[j] + dt0*(S_func(x,L) -C*dydx);
    }
	   
    // Efficiently copy next values at timestep to y array.
    double *tmp;
    tmp = y_next;
    y_next = y; 
    y = tmp;

    // Increment time.   
    ctime += dt0;
    if (output) {
      for (j=0; j<nx; j++ ) {
	x = j*dx;
	printf("%g %g %g \n",ctime,x,y[j]);
      }
      next_output_time += output_timestep;
    }
  }

  free(y);
  free(y_next);
  return 0;
}

void read_input(double *C, double *L, int *nx, double *t_F, double *t_out) {
   FILE *infile;
   if(!(infile=fopen("input.txt","r"))) {
       printf("Error opening file\n");
       exit(1);
   }
   if(5!=fscanf(infile,"%lf %lf %d %lf %lf",C,L,nx,t_F,t_out)) {
       printf("Error reading parameters from file\n");
       exit(1);
   }
   fclose(infile);
}
