/**
 * @file  StokesStandingWaves.c
 * @brief Simulation of standing waves using the Basilisk flow solver. 
 * @author Vatsal Sanjay
 * @version 0.1
 * @date Dec 24, 2024
 *
## Description:
 * This code simulates the dynamics of standing waves using adaptive mesh refinement. 
 * The simulation uses a two-phase flow solver with surface tension and implements
 * analytical Stokes wave solutions or experimental best-fit initial conditions.
 * 
 * Usage:
 * ./StokesStandingWaves maxLevel Ga Bo A0 ORDER tmax
 * where:
 *   maxLevel: Maximum refinement level for adaptive mesh
 *   Ga: Gallileo number (ratio of gravitational to viscous forces)
 *   Bo: Bond number (ratio of gravitational to surface tension forces)
 *   A0: Amplitude of the standing wave
 *   ORDER: Order of the initial condition: 0,1,2...6,upto 8.. or if you use -1, it uses the best fit. see [../initialConditionTest](../initialConditionTest)
 *   tmax: Maximum simulation time
 *
## Key Features:
 * Two-phase flow with surface tension
 * Adaptive mesh refinement for interface tracking
 * Choice between analytical Stokes wave or experimental best-fit initial conditions
 * Volume of Fluid (VoF) method for interface tracking
 *
## Physical Parameters:
 * Density ratio (rho1/rho2): 1000 (water-air like)
 * Viscosity ratio (mu2/mu1): 0.01 (water-air like)
 * Surface tension coefficient: 1.0/Bo
 * Gravity: -1.0 (non-dimensionalized)
 *
## Numerical Parameters:
 * Time interval between snapshots: tsnap = 1e-2
 * Error tolerances:
   - VOF (f): fErr = 1e-3
   - Curvature (K): KErr = 1e-6
   - Velocity: VelErr = 1e-3 (adjust based on Oh number)
 * Domain size: 2.0 × 2.0
 *
## Reference Variables:
 * Length scale: Wavelength of the standing wave
 * Time scale: Based on gravitational acceleration
 * All quantities are dimensionless
*/

#include "navier-stokes/centered.h"
#define FILTERED // Smear density and viscosity jumps
#include "two-phase.h"
#include "tension.h"
#include "reduced.h"
#include "distance.h"

/**
 * @brief Error tolerance parameters for numerical stability and accuracy
 */
#define tsnap (1e-2) // Time interval between snapshots. Use 0.001 only for specific cases requiring higher temporal resolution.
// Error tolerances
#define fErr (1e-3)                                 // error tolerance in f1 VOF (Volume of Fluid)
#define KErr (1e-6)                                 // error tolerance in VOF curvature calculated using height function method
#define VelErr (1e-3)                               // error tolerances in velocity -- Use 1e-2 for low Oh and 1e-3 to 5e-3 for high Oh/moderate to high J

// Domain size
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
  MAXlevel = 7; //atoi(argv[1]);
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

/**
 * @brief Calculates the Stokes coefficient for nth order wave
 * @param n Order of the wave
 * @return Coefficient value based on double factorial calculation
 */
double stokes_coefficient (int n) {
  if (n == 1) return 1.0;
  int k = 2*n - 1;
  double double_factorial = 1.0;
  for (int i = 1; i <= k; i += 2)
    double_factorial *= i;
  return double_factorial / (pow(2,(n-1)) * tgamma(n+1));
}

/**
 * @brief Calculates the surface elevation for Stokes wave
 * @param x Horizontal position
 * @param a Wave amplitude
 * @param k Wave number
 * @param order Order of Stokes expansion
 * @return Surface elevation at position x
 */
double eta_stokes (double x, double a, double k, int order) {
  double eta = 0.0;
  double epsilon = k*a;
  for (int n = 1; n <= order; n++) {
    double c = stokes_coefficient(n);
    eta += c * pow(epsilon, n)*cos(n*k*x);
  }
  return eta * a;
}

/**
 * @brief Initializes the simulation domain and wave profile
 * Either uses analytical Stokes wave solution or reads experimental best-fit data
 */
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
        phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1])/4.;
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
 * @brief Logs simulation data and implements stopping criteria
 * Records:
 * - Kinetic energy
 * - Maximum and minimum surface positions
 * - Simulation parameters
 * 
 * Stopping conditions:
 * - Kinetic energy exceeds 1e2 (blow-up)
 * - Kinetic energy falls below 1e-6 (stagnation)
 */
event logWriting (i++) {
  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += (0.5*rho(f[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
  }

  scalar pos[];
  position (f, pos, {0,1});
  double max = statsf(pos).max;
  double min = statsf(pos).min;

  if (pid() == 0) {
    static FILE * fp;
    if (i == 0) {
      fprintf(ferr, "Level %d, Ga %2.1e, Bo %2.1e, A0 %2.1e\n", MAXlevel, Ga, Bo, A0);
      fprintf (ferr, "i dt t ke at\n");
      fp = fopen (logFile, "w");
      fprintf(fp, "Level %d, Ga %2.1e, Bo %2.1e, A0 %2.1e\n", MAXlevel, Ga, Bo, A0);
      fprintf (fp, "i dt t ke at\n");
      fprintf (fp, "%d %g %g %g %g\n", i, dt, t, ke, max-min);
      fclose(fp);
    } else {
      fp = fopen (logFile, "a");
      fprintf (fp, "%d %g %g %g %g\n", i, dt, t, ke, max-min);
      fclose(fp);
    }
    fprintf (ferr, "%d %g %g %g %g\n", i, dt, t, ke, max-min);

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
