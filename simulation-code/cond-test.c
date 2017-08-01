/* cond-test.c: quantitative test of the anisotropic conduction routine


   Commentary:

   sets up a 2D problem with circular field lines using code from the
   field_loop test problem.

   adds an circularly symmetric Gaussian temperature perturbation;
   this should not evolve if conduction is perfectly anisotropic, but
   will evolve due to numeric errors...  we want to look at this
   evolution and quantify it as an effective perpendicular
   conductivity as a function of resolution.

   using an analytic solution for the evolution of this perturbation,
   we find that the central temperature reaches a value of 1.01839
   after one half of a conduction time, or 1.03570 after a full
   conduction time, where the conduction time is defined as
   (sigma^2)/kappa_perp.

   we kill the simulation when the central temperature reaches this
   threshold, dump a temperature profile to compare to the analytic
   solution, and compute an effective kappa_perp from the time to
   reach this central temperature.


   Usage:

   note that this problem should be run WITHOUT THE MHD INTEGRATOR.
   comment out calls to the integrator in main.c; we only want to test
   the conduction routine.

   compile using either

   ./configure --with-problem=cond-test --enable-conduction

   or

   ./configure --with-problem=cond-test --enable-conduction --enable-sts

   code also works with MPI

   note that the code should be run with tlim = nlim = infinity; the
   code runs until the central temperature reaches some threshold and
   then quits on its own.


   Output:

   dumps initial and final temperature profiles, and reports the
   effective conduction time via an ath_error() call.

   see the included run.rb for an example of how to use this problem.

*/
#include "copyright.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"



extern Real kappa_aniso;
extern Real kappa_iso;

static Real r0    = 0.5;
static Real sigma = 0.15;

static void bc_ix1(GridS *pGrid);
static void bc_ox1(GridS *pGrid);

static void bc_ix2(GridS *pGrid);
static void bc_ox2(GridS *pGrid);

void dump_temp_profile(DomainS *pD, char *fname);



void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;

  int i, j, k;
  int is, ie, js, je, ks, ke;
  int nx1, nx2;

  Real x1c, x2c, x3c;
  Real x1f, x2f, x3f;

  Real **az;
#ifdef MHD
  int ku;
#endif

  Real r, press;


  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;

  if ( (je-js) == 0 || (ke-ks) > 0 )
    ath_error("[cond_test]: This problem can only be run in 2D.\n");


  az = (Real**) calloc_2d_array(nx2, nx1, sizeof(Real));


#ifdef THERMAL_CONDUCTION
  kappa_iso   = 0.0;
  kappa_aniso = 1.0;
#else
  ath_error("[cond_test]: This problem requires conduction.\n");
#endif


  /* vector potential for circular field loops */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {

        fc_pos(pGrid, i, j, k, &x1f, &x2f, &x3f);

        r = sqrt(SQR(x1f) + SQR(x2f));

        if (r < 1.0)
          az[j][i] = 1.0 - r;
        else
          az[j][i] = 0.0;
      }
    }
  }


  /* initialize density and momentum */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->U[k][j][i].d  = 1.0;
        pGrid->U[k][j][i].M1 = 0.0;
        pGrid->U[k][j][i].M2 = 0.0;
        pGrid->U[k][j][i].M3 = 0.0;
      }
    }
  }


  /* curl of vector potential to get interface B field */
#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie+1; i++) {
        pGrid->B1i[k][j][i] = (az[j+1][i] - az[j][i])/pGrid->dx2;
      }
    }
  }
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->B2i[k][j][i] = -1.0*(az[j][i+1] - az[j][i])/pGrid->dx1;
      }
    }
  }
  ku = ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->B3i[k][j][i] = 0.0;
      }
    }
  }
#endif
  free_3d_array((void***)az);


  /* cell-centered B and magnetic energy */
#if defined MHD || !defined ISOTHERMAL
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid, i, j, k, &x1c, &x2c, &x3c);

        r = sqrt(SQR(x1c) + SQR(x2c));

#ifdef MHD
        pGrid->U[k][j][i].B1c = 0.5*(pGrid->B1i[k][j][i  ] +
                                     pGrid->B1i[k][j][i+1]);

        pGrid->U[k][j][i].B2c = 0.5*(pGrid->B2i[k][j  ][i] +
                                     pGrid->B2i[k][j+1][i]);

        pGrid->U[k][j][i].B3c = pGrid->B3i[k][j][i];
#endif

#ifndef ISOTHERMAL
        /* Gaussian temperature perturbation */
        press = 1.0 + 0.1 * exp(-SQR(r-r0) / (2*SQR(sigma)));
        pGrid->U[k][j][i].E = press/Gamma_1
#ifdef MHD
          + 0.5*(SQR(pGrid->U[k][j][i].B1c) + SQR(pGrid->U[k][j][i].B2c)
                 + SQR(pGrid->U[k][j][i].B3c))
#endif
          + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2)
                 + SQR(pGrid->U[k][j][i].M3))/pGrid->U[k][j][i].d;
#endif /* ISOTHERMAL */
      }
    }
  }
#endif


  /* set custom BCs which hold the temperature fixed along the boundary */
  bvals_mhd_fun(pDomain, left_x1,  bc_ix1);
  bvals_mhd_fun(pDomain, right_x1, bc_ox1);

  bvals_mhd_fun(pDomain, left_x2,  bc_ix2);
  bvals_mhd_fun(pDomain, right_x2, bc_ox2);


  /* dump an initial temperature profile */
  dump_temp_profile(pDomain, "ic.dat");

  return;
}


void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  DomainS *pDomain = &(pM->Domain[0][0]);

#ifdef THERMAL_CONDUCTION
  kappa_iso   = 0.0;
  kappa_aniso = 1.0;
#else
  ath_error("[cond_test]: This problem requires conduction.\n");
#endif

  bvals_mhd_fun(pDomain, left_x1,  bc_ix1);
  bvals_mhd_fun(pDomain, right_x1, bc_ox1);

  bvals_mhd_fun(pDomain, left_x2,  bc_ix2);
  bvals_mhd_fun(pDomain, right_x2, bc_ox2);

  return;
}


#ifndef BAROTROPIC
static Real temp(const GridS *pG, const int i, const int j, const int k)
{
  PrimS W;

  W = Cons_to_Prim(&(pG->U[k][j][i]));

  return (W.P/W.d);
}
#endif

ConsFun_t get_usr_expr(const char *expr)
{
#ifndef BAROTROPIC
  if(strcmp(expr, "temp")==0) return temp;
#endif
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}


/* monitor the central temperature and kill the simulation when it
   reaches a threshold corresponding to an effective conduction
   time */
void Userwork_in_loop(MeshS *pM)
{
  GridS *pG = pM->Domain[0][0].Grid;

  int i, j, k;
  int is, ie, js, je, ks, ke;

  Real x1, x2, x3, r, dx;

  Real mean_t[2];
#ifdef MPI_PARALLEL
  Real global_mean_t[2];
  int ierr;
#endif  /* MPI_PARALLEL */

  is = pG->is; ie = pG->ie;
  js = pG->js; je = pG->je;
  ks = pG->ks; ke = pG->ke;

  dx = MAX(pG->dx1, pG->dx2);


  /* calculate the mean temperature in the center of the domain */
  mean_t[0] = mean_t[1] = 0.0;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pG, i, j, k, &x1, &x2, &x3);
        r = sqrt(SQR(x1) + SQR(x2));

        if (r < 2*dx){
          mean_t[0] += temp(pG, i, j, k);
          mean_t[1] += 1.0;
        }
      }
    }
  }

  /* sum over processors, if necessary */
#ifdef MPI_PARALLEL
  ierr = MPI_Allreduce(&mean_t, &global_mean_t, 2,
                       MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(ierr)
    ath_error("[Userwork_in_loop]: MPI_Allreduce returned error %d\n", ierr);

  mean_t[0] = global_mean_t[0];
  mean_t[1] = global_mean_t[1];
#endif  /* MPI_PARALLEL */

  mean_t[0] = mean_t[0] / mean_t[1];

  ath_pout(0, "time=%f\tcentral-temp=%f\n",
           pM->time, mean_t[0]);



  /* using the analytic solution, after x% of a conduction time, the
     central temperature reaches y: */

  /*   X       Y    */
  /* -------------- */
  /*   5%   1.00159 */
  /*  10%   1.00286 */
  /*  50%   1.01839 */
  /* 100%   1.03570 */


  /* kill the simulation when the central temperature reaches
     1.01839... report this as half of an effective conduction time */
  if(mean_t[0] > 1.01839){
    dump_temp_profile(&(pM->Domain[0][0]), "end.dat");
    ath_error("simulation finished at time t = %f\n", pG->time);
  }

}

void Userwork_after_loop(MeshS *pM)
{
  return;
}


/* write out an azimuthally-averaged radial temperature profile */
void dump_temp_profile(DomainS *pD, char *fname)
{
  GridS *pG=pD->Grid;

  int i, j, k, s;
  int is, ie, js, je, ks, ke;

  int num;

  Real x1, x2, x3, r, dx1;

  Real **tdata;
#ifdef MPI_PARALLEL
  Real **tdata_global;
  int ierr;
#endif

  char fmt[80];
  FILE *outfile;

  sprintf(fmt," %%12.8e"); /* Use a default format */

  /* assume cubic cells! */
  dx1 = pD->dx[1];


  num = MIN(pD->Nx[0], pD->Nx[1]);
  num = ceil(num / sqrt(2));

  tdata = (Real**) calloc_2d_array(num, 2, sizeof(Real));
  for(i=0; i<num; i++){
    tdata[i][0] = tdata[i][1] = 0.0;
  }

#ifdef MPI_PARALLEL
  tdata_global = (Real**) calloc_2d_array(num, 2, sizeof(Real));
  for(i=0; i<num; i++){
    tdata_global[i][0] = tdata_global[i][1] = 0.0;
  }
#endif

  is = pG->is; ie = pG->ie;
  js = pG->js; je = pG->je;
  ks = pG->ks; ke = pG->ke;

  /* profile along 45 degrees */
  for(k=ks; k<=ke; k++){
    for(j=js; j<=je; j++){
      for(i=is; i<=ie; i++){
        cc_pos(pG, i, j, k, &x1, &x2, &x3);

        r = sqrt(x1*x1 + x2*x2);
        s = (int) floor(r/dx1);

        tdata[s][0] += temp(pG, i, j, k);
        tdata[s][1] += 1.0;
      }
    }
  }

#ifdef MPI_PARALLEL
  ierr = MPI_Allreduce(&tdata[0][0], &tdata_global[0][0], 2*num,
                       MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(ierr)
    ath_error("[dump_temp_profile]: MPI_Allreduce returned error %d\n", ierr);

  for(i=0; i<num; i++){
    tdata[i][0] = tdata_global[i][0];
    tdata[i][1] = tdata_global[i][1];
  }
#endif  /* MPI_PARALLEL */

  for(i=0; i<num; i++){
    if (tdata[i][1] > 0.0)
      tdata[i][0] = tdata[i][0] / tdata[i][1];
  }


  outfile = fopen(fname, "w");
  if (outfile == NULL)
    ath_error("[dump_temp_profile]: could not open file %s for writing.\n", fname);

  fprintf(outfile, "# N = %d\n", num);
  fprintf(outfile, "# temperature profile at Time = %g\n", pG->time);

  for(i=0; i<num; i++){
    fprintf(outfile, "%12d", i);
    fprintf(outfile, fmt, i*dx1);
    fprintf(outfile, fmt, tdata[i][0]);
    fprintf(outfile, "\n");
  }

  fclose(outfile);

  return;
}



/* fixed boundary conditions.  since we don't call the MHD integrator
   (right??), we only need to fix the temperature. */

static void bc_ix1(GridS *pGrid)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  Real x1, x2, x3, r, press;

  cc_pos(pGrid, i, j, k, &x1, &x2, &x3);

  r = sqrt(SQR(x1) + SQR(x2));

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
#ifndef ISOTHERMAL
        press = 1.0 + 0.1 * exp(-SQR(r-r0) / (2*SQR(sigma)));
        pGrid->U[k][j][i].E = press/Gamma_1
#ifdef MHD
          + 0.5*(SQR(pGrid->U[k][j][i].B1c) + SQR(pGrid->U[k][j][i].B2c)
                 + SQR(pGrid->U[k][j][i].B3c))
#endif
          + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2)
                 + SQR(pGrid->U[k][j][i].M3))/pGrid->U[k][j][i].d;
#endif /* ISOTHERMAL */
      }
    }
  }

  return;
}

static void bc_ox1(GridS *pGrid)
{
  int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  Real x1, x2, x3, r, press;

  cc_pos(pGrid, i, j, k, &x1, &x2, &x3);

  r = sqrt(SQR(x1) + SQR(x2));

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=ie; i<=ie+nghost; i++) {
#ifndef ISOTHERMAL
        press = 1.0 + 0.1 * exp(-SQR(r-r0) / (2*SQR(sigma)));
        pGrid->U[k][j][i].E = press/Gamma_1
#ifdef MHD
          + 0.5*(SQR(pGrid->U[k][j][i].B1c) + SQR(pGrid->U[k][j][i].B2c)
                 + SQR(pGrid->U[k][j][i].B3c))
#endif
          + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2)
                 + SQR(pGrid->U[k][j][i].M3))/pGrid->U[k][j][i].d;
#endif /* ISOTHERMAL */
      }
    }
  }

  return;
}

static void bc_ix2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  Real x1, x2, x3, r, press;

  cc_pos(pGrid, i, j, k, &x1, &x2, &x3);

  r = sqrt(SQR(x1) + SQR(x2));

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
#ifndef ISOTHERMAL
        press = 1.0 + 0.1 * exp(-SQR(r-r0) / (2*SQR(sigma)));
        pGrid->U[k][j][i].E = press/Gamma_1
#ifdef MHD
          + 0.5*(SQR(pGrid->U[k][j][i].B1c) + SQR(pGrid->U[k][j][i].B2c)
                 + SQR(pGrid->U[k][j][i].B3c))
#endif
          + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2)
                 + SQR(pGrid->U[k][j][i].M3))/pGrid->U[k][j][i].d;
#endif /* ISOTHERMAL */
      }
    }
  }

  return;
}

static void bc_ox2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  Real x1, x2, x3, r, press;

  cc_pos(pGrid, i, j, k, &x1, &x2, &x3);

  r = sqrt(SQR(x1) + SQR(x2));

  for (k=ks; k<=ke; k++) {
    for (j=je; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
#ifndef ISOTHERMAL
        press = 1.0 + 0.1 * exp(-SQR(r-r0) / (2*SQR(sigma)));
        pGrid->U[k][j][i].E = press/Gamma_1
#ifdef MHD
          + 0.5*(SQR(pGrid->U[k][j][i].B1c) + SQR(pGrid->U[k][j][i].B2c)
                 + SQR(pGrid->U[k][j][i].B3c))
#endif
          + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2)
                 + SQR(pGrid->U[k][j][i].M3))/pGrid->U[k][j][i].d;
#endif /* ISOTHERMAL */
      }
    }
  }

  return;
}

/* end cond-test.c */
