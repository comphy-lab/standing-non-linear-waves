/**
 * @file coatedDrops.c
 * @brief Simulation of falling drops coated with an elastic shell using the Basilisk framework
 * @author Vatsal Sanjay
 * @version 0.1
 * @date Dec 13, 2024

## Description:
 * This code simulates the dynamics of falling drops coated with an elastic shell.
 * The simulation uses a three-phase flow model with viscoelastic properties,
 * implemented using the log-conformation approach for numerical stability.
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

#include "axi.h"
#include "navier-stokes/centered.h"

/*
see: V. Sanjay, Zenodo, DOI: 10.5281/zenodo.14210635 (2024) for details
*/
#define _SCALAR // uncomment to use the scalar version of the viscoelastic code, use scalar for easier extenson to 3D.
#ifndef _SCALAR
#include "log-conform-viscoelastic.h" 
#else 
#include "log-conform-viscoelastic-scalar-2D.h"
#endif

/**
 * Simulation Parameters:
 * FILTERED: Enable density and viscosity jump smoothing
 * tsnap: Time interval between snapshots (default: 1e-2)
 * fErr: Error tolerance for volume fraction (1e-3)
 * KErr: Error tolerance for curvature calculation (1e-6)
 * VelErr: Error tolerance for velocity field (1e-3)
 * AErr: Error tolerance for conformation tensor (1e-3)
 * Ldomain: Domain size in characteristic lengths (8)
*/
#define FILTERED // Smear density and viscosity jumps
#include "three-phase-viscoelastic.h"
#include "tension.h"

#define tsnap (1e-2) // 0.001 only for some cases. 
// Error tolerancs
#define fErr (1e-3)                                 // error tolerance in f1 VOF
#define KErr (1e-6)                                 // error tolerance in VoF curvature calculated using heigh function method (see adapt event)
#define VelErr (1e-3)                               // error tolerances in velocity -- Use 1e-2 for low Oh and 1e-3 to 5e-3 for high Oh/moderate to high J
#define AErr (1e-3)                             // error tolerances in conformation inside the liquid

// Numbers!
#define Ldomain 8

// boundary conditions - outflow on the right boundary
u.n[right] = neumann(0.);
p[right] = dirichlet(0.);

// left wall
u.t[left] = dirichlet(0.);
f1[left] = dirichlet(0.);
f2[left] = dirichlet(0.);


int MAXlevel;
// Oh -> Solvent Ohnesorge number
// Oha -> air Ohnesorge number
// De -> Deborah number
// Ec -> Elasto-capillary number

double We, Ohd, tf, Ohs, De, Ec, Oha, tmax;
char nameOut[80], dumpFile[80], logFile[80];

int  main(int argc, char const *argv[]) {

  L0 = Ldomain;
  origin (0., 0.);
  
  /*
  Values taken from the terminal. Here we use some representative values. In production run, you can pass it from the command line.
  */
  MAXlevel = 10; //atoi(argv[1]);
  We = 4e0; //atof(argv[1]);
  Ohd = 1e-2; //atof(argv[2]);
  tf = 5e-2; // atof(argv[3]);
  Ohs = 1e-2; //atof(argv[4]); -- Elastic implies Ohs \to 0
  De = 1e30; //atof(argv[2]); // Use a value of 1e30 to simulate the De \to \infty limit. 
  Ec = 1e0; //atof(argv[3]);
  tmax = 3e0; //atof(argv[6]);

  // Ensure that all the variables were transferred properly from the terminal or job script.
  // if (argc < 8){
  //   fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments\n", 8-argc);
  //   return 1;
  // }
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
 * rho1, rho2, rho3: Density of shell, drop and gas phases
 * mu1, mu2, mu3: Dynamic viscosity of shell, liquid and gas phases
 * lambda1, lambda2, lambda3: Relaxation times
 * G1, G2, G3: Elastic moduli
*/
  rho1 = 1., rho2 = 1e0; rho3 = 1e-3;
  Oha = 2e-2 * Ohd; // air-water viscosity ratio
  mu1 = Ohs/sqrt(We), mu2 = Ohd/sqrt(We), mu3 = Oha/sqrt(We);
  lambda1 = De*sqrt(We); lambda2 = 0., lambda3 = 0.;
  G1 = Ec/We; G2 = 0., G3 = 0.;

  f1.sigma = 1.0/We; f2.sigma = 1.0/We;

  //  BEWARE of these for stability issues. 
  dtmax = 1e-5; 
  TOLERANCE=1e-4;
  CFL = 1e-1;

  run();
}

// initial shape of the drop and shell
#define xDist 0.02
#define R2(x,y,z,TH) (sq(x-(1e0+xDist+TH)) + sq(y))

event init (t = 0) {
  if (!restore (file = dumpFile)){
    refine((R2(x,y,z,tf) < 1.21*sq(1e0+tf)) && (level < MAXlevel)); 
    fraction (f2, 1. - R2(x,y,z,tf));
    fraction (f1, sq(1e0+tf) - R2(x,y,z,tf));
    foreach(){
      u.x[] = -1.0*f1[];
    }
  }
}

/**
## Adaptive Mesh Refinement
*/
event adapt(i++){
  scalar KAPPA1[], KAPPA2[];
  curvature(f1, KAPPA1);
  curvature(f2, KAPPA2);

  #ifndef _SCALAR
   adapt_wavelet ((scalar *){f1, f2, u.x, u.y, conform_p.x.x, conform_p.y.y, conform_p.y.x, conform_qq, KAPPA1, KAPPA2},
      (double[]){fErr, fErr, VelErr, VelErr, AErr, AErr, AErr, AErr, KErr, KErr},
      MAXlevel, MAXlevel-6);
  #else
   adapt_wavelet ((scalar *){f1, f2, u.x, u.y, A11, A22, A12, AThTh, KAPPA1, KAPPA2},
      (double[]){fErr, fErr, VelErr, VelErr, AErr, AErr, AErr, AErr, KErr, KErr},
      MAXlevel, MAXlevel-6);
  #endif
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
    fprintf(ferr, "Level %d, We %2.1e, Ohd %2.1e, tf %2.1e, Ohs %2.1e, De %2.1e, Ec %2.1e, Oha %2.1e,\n", MAXlevel, We, Ohd, tf, Ohs, De, Ec, Oha);
}

/**
## Log writing
*/
event logWriting (i++) {
  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += (2*pi*y)*(0.5*rho(f1[], f2[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
  }
  if (pid() == 0) {
    static FILE * fp;
    if (i == 0) {
      fprintf(ferr, "Level %d, We %2.1e, Ohd %2.1e, tf %2.1e, Ohs %2.1e, De %2.1e, Ec %2.1e, Oha %2.1e,\n", MAXlevel, We, Ohd, tf, Ohs, De, Ec, Oha);
      fprintf (ferr, "i dt t ke\n");
      fp = fopen (logFile, "w");
      fprintf(fp, "Level %d, We %2.1e, Ohd %2.1e, tf %2.1e, Ohs %2.1e, De %2.1e, Ec %2.1e, Oha %2.1e,\n", MAXlevel, We, Ohd, tf, Ohs, De, Ec, Oha);
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

  if (ke > 1e2 && i > 1e1){
    if (pid() == 0){
      fprintf(ferr, "The kinetic energy blew up. Stopping simulation\n");
      fp = fopen (logFile, "a");
      fprintf(fp, "The kinetic energy blew up. Stopping simulation\n");
      fclose(fp);
      dump(file=dumpFile);
      return 1;
    }
  }
  if (ke < 1e-6 && i > 1e1){
    if (pid() == 0){
      fprintf(ferr, "kinetic energy too small now! Stopping!\n");
      dump(file=dumpFile);
      fp = fopen (logFile, "a");
      fprintf(fp, "kinetic energy too small now! Stopping!\n");
      fclose(fp);
      return 1;
    }
  }
  }
}
