#include <stdlib.h>
#include <stdio.h>
int main(void) {
  double *T, *T_n;  /* Array to store current and next time step values */
  double length = 1.0;
  long nx = 10;
  /* Need to include both boundary points at x_0 and x_N, so
     that domain is length 1.0 */
  double dx = length/(nx-1);
  double time_final = 1.0; 
  double dt = 0.25*(dx*dx);  /* Time step set by stability limit dt/dx^2 < 0.5 */
  long nsteps = time_final/dt;
  /* In order to finish _exactly_ at time_final, increase the number of steps by one and 
  decrease the time step slightly: stability still satisfied. */
  nsteps+=1;
  dt = time_final/nsteps;   
  T   = malloc(sizeof(double)*nx); 
  T_n = malloc(sizeof(double)*nx); 
  long j; long k; /* Space and time grid indices respectively */
  /* Initialisation: delta function at mid-radius */
  for(j=0;j<nx;j++) {
    T[j] = 0.0;
  }
  T[nx/2] = 1.0;

  /* Print first time step */
  double t0 = 0.0;
  for(j=0;j<nx;j++) {
    printf("%g %g %g \n",t0,j*dx,T[j]);
  }
  FILE *fp = fopen("outparams.txt","w"); /* WRITES OUTPUT PARAMETERS TO FILE */
  fprintf(fp,"%ld %ld %g %g \n", nx, nsteps, dx, dt);
  fclose(fp);

  for(k=0;k<nsteps;k++) {
    /* Loop over interior points */
    for(j=1;j<(nx-1);j++) {
      double deriv = (T[j-1] + T[j+1] - 2*T[j])/(dx*dx);
      T_n[j] = T[j] + deriv * dt;
    }
    /* Copy temperature at next step back into current step */
    for(j=0;j<nx;j++) {
      T[j] = T_n[j];
      printf("%g %g %g \n",k*dt,j*dx,T[j]);
    }
  }
  free(T);
  free(T_n);
  return 0;
}
