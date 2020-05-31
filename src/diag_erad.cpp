/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "diag_erad.h"
#include "app.h"
#include "app_erad.h"
#include "comm_lattice.h"
#include "timer.h"
#include "error.h"
#include "memory.h"

using namespace SPPARKS_NS;

enum{inter,floater};                              // data type 
enum{ZERO,FE,VACANCY,I1,I2,I3,I4,I5,I6};        // diagnosis terms   
enum{hFE=11,hVAC,hI1,hI2,hI3,hI4,hI5,hI6};               // hop steps for each element
enum{react=20};                                 //reactions 
enum{sink=30};                                    // number of sink absorption   
enum{RcFe=41,RcVAC,RcI1,RcI2,RcI3,RcI4,RcI5,RcI6};                               // number of recombination    
enum{FPair=51,nclst};                             // Frenkle Pairs created by ballistic  
enum{dFE=61,dVAC,dI1,dI2,dI3,dI4,dI5,dI6}; // MSD for each element 
enum{cFE=71,cVAC,cI1,cI2,cI3,cI4,cI5,cI6}; // time everaged concentration  
enum{energy=81,rclst,csia,percolation};                        // energy and realistic time  
/* ---------------------------------------------------------------------- */

DiagErad::DiagErad(SPPARKS *spk, int narg, char **arg) : Diag(spk,narg,arg)
{
  if (strcmp(app->style,"erad") != 0)
    error->all(FLERR,"Diag_style erad requires app_style erad");

  nlist = 0;
 
  int iarg = iarg_child;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"list") == 0) {
      nlist = narg - iarg - 1;
      list = new char*[nlist];
      int j = 0;
      for (int i = iarg+1; i < narg; i++) {
	int n = strlen(arg[i]) + 1;
	list[j] = new char[n];
	strcpy(list[j],arg[i]);
	j++;
      }
      iarg = narg;
    } else error->all(FLERR,"Illegal diag_style erad command");
  }

  if (nlist == 0) error->all(FLERR,"Illegal diag_style erad command");
  which = new int[nlist];
  index = new int[nlist];
  itype = new int[nlist];
  ivector = new int[nlist];
  dvector = new double[nlist+1]; //nelement + total energy 
}

/* ---------------------------------------------------------------------- */

DiagErad::~DiagErad()
{
  for (int i = 0; i < nlist; i++) delete [] list[i];
  delete [] list;
  delete [] which;
  delete [] index;
  delete [] itype;
  delete [] ivector;
  delete [] dvector;
}

/* ---------------------------------------------------------------------- */

void DiagErad::init()
{
  apperad = (AppErad *) app;

  siteflag = 0;
  hopflag = 0;
  clstflag = 0;
  csiteflag = 0;

  for (int i = 0; i < nlist; i++) {
    if (strcmp(list[i],"fe") == 0) which[i] = FE; //total sites
    else if (strcmp(list[i],"vac") == 0) which[i] = VACANCY;
    else if (strcmp(list[i],"i1") == 0) which[i] = I1;
    else if (strcmp(list[i],"i2") == 0) which[i] = I2;
    else if (strcmp(list[i],"i3") == 0) which[i] = I3;
    else if (strcmp(list[i],"i4") == 0) which[i] = I4;
    else if (strcmp(list[i],"i5") == 0) which[i] = I5;
    else if (strcmp(list[i],"i6") == 0) which[i] = I6;
   
    else if (strcmp(list[i],"hfe") == 0) which[i] = hFE;//total hop events
    else if (strcmp(list[i],"hvac") == 0) which[i] = hVAC;
    else if (strcmp(list[i],"hi1") == 0) which[i] = hI1;
    else if (strcmp(list[i],"hi2") == 0) which[i] = hI2;
    else if (strcmp(list[i],"hi3") == 0) which[i] = hI3;
    else if (strcmp(list[i],"hi4") == 0) which[i] = hI4;
    else if (strcmp(list[i],"hi5") == 0) which[i] = hI5;
    else if (strcmp(list[i],"hi6") == 0) which[i] = hI6;
 
    else if (strcmp(list[i],"dfe") == 0) which[i] = dFE;//MSD 
    else if (strcmp(list[i],"dvac") == 0) which[i] = dVAC; 
    else if (strcmp(list[i],"di1") == 0) which[i] = dI1;
    else if (strcmp(list[i],"di2") == 0) which[i] = dI2;
    else if (strcmp(list[i],"di3") == 0) which[i] = dI3;
    else if (strcmp(list[i],"di4") == 0) which[i] = dI4;
    else if (strcmp(list[i],"di5") == 0) which[i] = dI5;
    else if (strcmp(list[i],"di6") == 0) which[i] = dI6; 
 
    else if (strcmp(list[i],"cfe") == 0) which[i] = cFE;//concentration  
    else if (strcmp(list[i],"cvac") == 0) which[i] = cVAC; 
    else if (strcmp(list[i],"ci1") == 0) which[i] = cI1;
    else if (strcmp(list[i],"ci2") == 0) which[i] = cI2;
    else if (strcmp(list[i],"ci3") == 0) which[i] = cI3;
    else if (strcmp(list[i],"ci4") == 0) which[i] = cI4;
    else if (strcmp(list[i],"ci5") == 0) which[i] = cI5;
    else if (strcmp(list[i],"ci6") == 0) which[i] = cI6; 

    else if (strcmp(list[i],"rcvac") == 0) which[i] = RcVAC; 
    else if (strcmp(list[i],"rci1") == 0) which[i] = RcI1;
    else if (strcmp(list[i],"rci2") == 0) which[i] = RcI2;
    else if (strcmp(list[i],"rci3") == 0) which[i] = RcI3;
    else if (strcmp(list[i],"rci4") == 0) which[i] = RcI4;
    else if (strcmp(list[i],"rci5") == 0) which[i] = RcI5;
    else if (strcmp(list[i],"rci6") == 0) which[i] = RcI6; 
    else if (strcmp(list[i],"energy") == 0) which[i] = energy;
    else if (strcmp(list[i],"nfp") == 0) which[i] = FPair;
    else if (strcmp(list[i],"nclst") == 0) which[i] = nclst; 
    else if (strcmp(list[i],"rclst") == 0) which[i] = rclst; 
    else if (strcmp(list[i],"csia") == 0) which[i] = csia; 
    else if (strcmp(list[i],"percolation") == 0) which[i] = percolation; 

    else if (list[i][0] == 'r' && list[i][1] == 'e' && list[i][2] == 'c' && list[i][3] == 't') {
      int id = list[i][4] - '0';
      which[i] = react + id; 
    }

    else if (list[i][0] == 's' && list[i][1] == 'i' && list[i][2] == 'n' && list[i][3] == 'k') {
      int id = list[i][4] - '0';
      which[i] = sink + id; 
    }

    else error->all(FLERR,"Invalid value setting in diag_style erad");
  }

  for (int i = 0; i < nlist; i++) {
    if (which[i] == FE || which[i] == VACANCY || which[i] == I1 || which[i] == I2 || which[i] == I3 || which[i] == I4 || which[i] == I5 || which[i] == I6)
      siteflag = 1;
    if (which[i] == hFE || which[i] == hVAC || which[i] == hI1 || which[i] == hI2 || which[i] == hI3 || which[i] == hI4 || which[i] == hI5 || which[i] == hI6)
      hopflag = 1;
    if (which[i] == cFE || which[i] == cVAC || which[i] == cI1 || which[i] == cI2 || which[i] == cI3 || which[i] == cI4 || which[i] == cI5 || which[i] == cI6)
      csiteflag = 1;
    if (which[i] == nclst || which[i] == rclst) clstflag = 1; // do cluster analysis 
  }

  for (int i = 0; i < nlist; i++) { ivector[i] = 0;
    dvector[i] = 0.0;
  }

  for (int i=0; i < nlist; i++) { itype[i] = inter;
    if(which[i] >= dFE) itype[i] = floater;
  }
}

/* ---------------------------------------------------------------------- */

void DiagErad::compute()
{
  int ninter,nfloater;
  int sites[10],nhop[10],ivalue; // int data
  int nlocal = apperad->nlocal;
  int nelement = apperad->nelement;
  double *csites; // concentration 
  double dvalue ; // double data 

  ninter = nfloater = 0;
  dvalue = 0.0;

  if(csiteflag) {csites = apperad->ct;}
  /*// print the recombination vectors, test only  
  for (int ii=0; ii<361; ii++) {
  for (int ij=0; ij<181; ij++) {
      if(apperad->reccount[ii][ij] == 0) continue;  
      fprintf(screen," %d %d %d \n", ii, ij, apperad->reccount[ii][ij]); 
  } 
  }*/ 
 
  if (siteflag) {
    sites[FE] = sites[VACANCY] = sites[I1] = sites[I2] = 0;
    sites[I3] = sites[I4] = sites[I5] = sites[I6] = 0;
    int *element = apperad->element;
    for (int i = 0; i < nlocal; i++) sites[element[i]]++;
  }

  if(hopflag) { // count absoprtion by external sinks 
    for(int i = 1; i < nelement+1; i++) {nhop[i] = 0; nhop[i] = apperad->nhmfp[i];}   
  }

  //if(clstflag) apperad->cluster(); // cluster analysis 

  for (int i = 0; i < nlist; i++) {
    if (which[i] == FE) ivalue = sites[FE]; //total sites
    else if (which[i] == VACANCY) ivalue = sites[VACANCY];
    else if (which[i] == I1) ivalue = sites[I1];
    else if (which[i] == I2) ivalue = sites[I2];
    else if (which[i] == I3) ivalue = sites[I3]; 
    else if (which[i] == I4) ivalue = sites[I4]; 
    else if (which[i] == I5) ivalue = sites[I5];
    else if (which[i] == I6) ivalue = sites[I6];
    /*
    else if (which[i] == mVACANCY) ivalue = monomer_local[VACANCY];
    else if (which[i] == mI1) ivalue = monomer_local[I1];
    else if (which[i] == mI2) ivalue = monomer_local[I2];
    else if (which[i] == mI3) ivalue = monomer_local[I3]; 
    else if (which[i] == mP) ivalue = monomer_local[P]; 
    else if (which[i] == mC) ivalue = monomer_local[C];
    else if (which[i] == mI4A) ivalue = monomer_local[I4A];
    */
    else if (which[i] == hFE) ivalue = nhop[FE]; //total hop events
    else if (which[i] == hVAC) ivalue = nhop[VACANCY];
    else if (which[i] == hI1) ivalue = nhop[I1];
    else if (which[i] == hI2) ivalue = nhop[I2];
    else if (which[i] == hI3) ivalue = nhop[I3]; 
    else if (which[i] == hI4) ivalue = nhop[I4]; 
    else if (which[i] == hI5) ivalue = nhop[I5];
    else if (which[i] == hI6) ivalue = nhop[I6];

    else if (which[i] == cFE) dvalue = csites[FE]; //time_averaged concentration 
    else if (which[i] == cVAC) dvalue = csites[VACANCY];
    else if (which[i] == cI1) dvalue = csites[I1];
    else if (which[i] == cI2) dvalue = csites[I2];
    else if (which[i] == cI3) dvalue = csites[I3];
    else if (which[i] == cI4) dvalue = csites[I4];
    else if (which[i] == cI5) dvalue = csites[I5];
    else if (which[i] == cI6) dvalue = csites[I6];

    else if (which[i] == RcVAC) ivalue = apperad->nrecombine[VACANCY-1]; //number of reocmbination 
    else if (which[i] == RcI1) ivalue = apperad->nrecombine[I1-1]; //number of reocmbination 
    else if (which[i] == RcI2) ivalue = apperad->nrecombine[I2-1]; //number of reocmbination 
    else if (which[i] == RcI3) ivalue = apperad->nrecombine[I3-1]; //number of reocmbination 
    else if (which[i] == RcI4) ivalue = apperad->nrecombine[I4-1]; //number of reocmbination 
    else if (which[i] == RcI5) ivalue = apperad->nrecombine[I5-1]; //number of reocmbination 
    else if (which[i] == RcI6) ivalue = apperad->nrecombine[I6-1]; //number of reocmbination 

    else if (which[i] == FPair) ivalue = apperad->nFPair; //total Frenkel pairs generater  

    else if (which[i] == energy) dvalue = apperad->total_energy(); //system energy   
    else if (which[i] == nclst) ivalue = apperad->ncluster; //# of cluster
    else if (which[i] == rclst) dvalue = apperad->rcluster; //size of cluster 
    else if (which[i] == csia) dvalue = apperad->csia; //sia concentration
    else if (which[i] == percolation) dvalue = apperad->percolation(); //percolation rate  
  
    else if (which[i] > react && which[i] < sink) {
      int id = which[i] - react; 
      ivalue = apperad->rcount[id-1];
    }
 
    else if (which[i] > sink && which[i] < RcFe) {
      int id = which[i] - sink; 
      ivalue = apperad->nabsorption[id-1];
    }
 
    if(which[i] >= dFE) nfloater++;
    else ninter++;

    MPI_Allreduce(&ivalue,&ivector[ninter],1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&dvalue,&dvector[nfloater],1,MPI_DOUBLE,MPI_SUM,world);
  }
}

/* ---------------------------------------------------------------------- */

void DiagErad::stats(char *str)
{
  int ninter, nfloater;

  ninter = nfloater = 0;
  for (int i = 0; i < nlist; i++) {
    if(itype[i] == inter) { ninter++;
      sprintf(str," %d",ivector[ninter]);
    }
    else if(itype[i] == floater) { nfloater++;
      sprintf(str," %14.8g",dvector[nfloater]);
    }
    str += strlen(str);
  }
}

/* ---------------------------------------------------------------------- */

void DiagErad::stats_header(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str," %s",list[i]);
    str += strlen(str);
  }
}
