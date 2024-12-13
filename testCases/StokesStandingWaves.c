/**
 * @file  StokesStandingWaves.c
 * @brief Simulation of standing waves. 
 * @author Vatsal Sanjay
 * @version 0.1
 * @date Dec 13, 2024

## Description:
 * This code simulates the dynamics of a standing Stokes wave. 
 * 
 * Usage:
 * ./program maxLevel We Ohd tf Ohs De Ec tmax
 * where:
 *   maxLevel: Maximum refinement level for adaptive mesh
 *   We: Weber number (ratio of kinetic energy to surface energy)
 *   Ohd: Ohnesorge number for the drop (ratio of viscous to inertial-capillary forces)
 *   tf: thickness of the shell
 *   Ohs: Ohnesorge number for the shell (ratio of viscous to inertial-capillary forces) -- note that this is the solvent viscosity of the shell. In case of purely elastic (rather Kelvin–Voigt) shell, it is the background viscosity responsible for dissipation in solids.
 *   De: Deborah number (ratio of relaxation time to flow time). For Elastic/Kelvin–Voigt solids, keep this at 1e30 (hard coded to be infinity).
 *   Ec: Elasto-capillary number (ratio of elastic to surface tension forces) of the shell. 
 *   tmax: Maximum simulation time

## Some assumptions to start with:
 * No gravity (easy to change).
 * No contact between the compound drop and the substrate (not planned to be added in the future).
 * Shell - drop density ratio is 1 (very easy to change).
 * Air-water viscosity ratio (2e-2) is used for ambient air - drop (very easy to change).
 * Drop-shell and shell-air surface tensions are the same (both equal to 1). In the future, might be instructive to test with different surface tensions. Also, in that case, we must decide which surface tension will be 1.

## Reating variables:
 * Radius of the drop (without the shell). 
 * Surface tension (currently both shell-drop and shell-air are the same).
 * Density of the drop (currently, both shell and drop densities are the same).
*/

#include "navier-stokes/centered.h"
#define FILTERED // Smear density and viscosity jumps
#include "two-phase.h"
#include "tension.h"
#include "reduced.h"
#include "distance.h"

#define tsnap (1e-2) // 0.001 only for some cases. 
// Error tolerancs
#define fErr (1e-3)                                 // error tolerance in f1 VOF
#define KErr (1e-6)                                 // error tolerance in VoF curvature calculated using heigh function method (see adapt event)
#define VelErr (1e-3)                               // error tolerances in velocity -- Use 1e-2 for low Oh and 1e-3 to 5e-3 for high Oh/moderate to high J

// Numbers!
#define Ldomain 2.0


int MAXlevel;
// Ga -> Gallileo number 
// Bo -> Bond number

double Ga, Bo, A0, tmax;
int ORDER;
char nameOut[80], dumpFile[80], logFile[80];

int  main(int argc, char const *argv[]) {

  L0 = Ldomain;
  origin (-L0/2., -L0/2);
  
  /*
  Values taken from the terminal. Here we use some representative values. In production run, you can pass it from the command line.
  */
  MAXlevel = 8; //atoi(argv[1]);
  Ga = 1e6;
  Bo = 1e1;
  A0 = 0.1;
  ORDER = -1;
  tmax = 1e1; //atof(argv[6]);

  init_grid (1 << 5);
  // Create a folder named intermediate where all the simulation snapshots are stored.
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  // Name of the restart file. See writingFiles event.
  sprintf (dumpFile, "restart");
  sprintf(logFile, "logData.dat");

/**
 * Physical Properties:
 * rho1, rho2: Density of liquid and gas phase
 * mu1, mu2: Dynamic viscosity of liquid and gas phase
*/
  rho1 = 1., rho2 = 1e-3;
  mu1 = 1e0/sqrt(Ga), mu2 = 1e-2/sqrt(Ga);
  G.y = -1e0;
  f.sigma = 1.0/Bo;

  // //  BEWARE of these for stability issues. 
  // dtmax = 1e-5; 
  // TOLERANCE=1e-4;
  // CFL = 1e-1;

  run();
}

double stokes_coefficient (int n) {
  if (n == 1) return 1.0;
  int k = 2*n - 1;
  double double_factorial = 1.0;
  for (int i = 1; i <= k; i += 2)
    double_factorial *= i;
  return double_factorial / (pow(2,(n-1)) * tgamma(n+1));
}

double eta_stokes (double x, double a, double k, int order) {
  double eta = 0.0;
  double epsilon = k*a;
  for (int n = 1; n <= order; n++) {
    double c = stokes_coefficient(n);
    eta += c * pow(epsilon, n)*cos(n*k*x);
  }
  return eta * a;
}

event init (t = 0) {
  if (!restore (file = dumpFile)) {
    refine (y > -1.25*A0 && y < 1.25*A0 && level < MAXlevel);
    if (ORDER > 0){
    fraction (f, -(y - eta_stokes(x, A0, 2*pi, ORDER)));
    } else {
      fprintf (ferr, "Using best fit using order 1 and 6.\n");
      char filename[60];
      sprintf(filename,"expBestFit.dat");
      FILE * fp = fopen(filename,"rb");
        if (fp == NULL){
          fprintf(ferr, "There is no file named %s\n", filename);
          // try in folder one level up
          sprintf(filename,"../../initialConditionTest/expBestFit.dat");
          fp = fopen(filename,"rb");
          if (fp == NULL){
            fprintf(ferr, "There is no file named %s\n", filename);
            return 1;
          }
        }
      coord* InitialShape;
      InitialShape = input_xy(fp);
      fclose (fp);
      scalar d[];
      distance (d, InitialShape);

      while (adapt_wavelet ((scalar *){f, d}, (double[]){1e-8, 1e-8}, MAXlevel).nf);
      
    // The distance function is defined at the center of each cell, we have
    // to calculate the value of this function at each vertex. 
      vertex scalar phi[];
      foreach_vertex(){
        phi[] = -(d[] + d[-1] + d[0,-1] + d[-1,-1])/4.;
      }
      
    // We can now initialize the volume fraction of the domain. 
      fractions (phi, f);
    }
  }
}

/**
## Adaptive Mesh Refinement
*/
event adapt(i++){
  scalar KAPPA[];
  curvature(f, KAPPA);

   adapt_wavelet ((scalar *){f, u.x, u.y, KAPPA},
      (double[]){fErr, VelErr, VelErr, KErr},
      MAXlevel, MAXlevel-6);
}

/**
## Dumping snapshots
*/
event writingFiles (t = 0; t += tsnap; t <= tmax) {
  dump (file = dumpFile);
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}

/**
## Ending Simulation
*/
event end (t = end) {
  if (pid() == 0)
    fprintf(ferr, "Level %d, Ga %2.1e, Bo %2.1e, A0 %2.1e\n", MAXlevel, Ga, Bo, A0);
}

/**
## Log writing
*/
event logWriting (i++) {
  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += (0.5*rho(f[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
  }

  if (pid() == 0) {
    static FILE * fp;
    if (i == 0) {
      fprintf(ferr, "Level %d, Ga %2.1e, Bo %2.1e, A0 %2.1e\n", MAXlevel, Ga, Bo, A0);
      fprintf (ferr, "i dt t ke\n");
      fp = fopen (logFile, "w");
      fprintf(fp, "Level %d, Ga %2.1e, Bo %2.1e, A0 %2.1e\n", MAXlevel, Ga, Bo, A0);
      fprintf (fp, "i dt t ke\n");
      fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
      fclose(fp);
    } else {
      fp = fopen (logFile, "a");
      fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
      fclose(fp);
    }
    fprintf (ferr, "%d %g %g %g\n", i, dt, t, ke);

    assert(ke > -1e-10);

    if (ke > 1e2 && i > 1e2) {
      fprintf(ferr, "The kinetic energy blew up. Stopping simulation\n");
      fp = fopen (logFile, "a");
      fprintf(fp, "The kinetic energy blew up. Stopping simulation\n");
      fclose(fp);
      dump(file=dumpFile);
      return 1;
    }
    
    if (ke < 1e-6 && i > 1e2) {
      fprintf(ferr, "kinetic energy too small now! Stopping!\n");
      dump(file=dumpFile);
      fp = fopen (logFile, "a");
      fprintf(fp, "kinetic energy too small now! Stopping!\n");
      fclose(fp);
      return 1;
    }
  }
}
