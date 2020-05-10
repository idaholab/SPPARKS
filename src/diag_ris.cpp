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
#include "stdlib.h"
#include "string.h"
#include "diag_ris.h"
#include "app.h"
#include "app_ris.h"
#include "comm_lattice.h"
#include "timer.h"
#include "error.h"
#include "memory.h"

using namespace SPPARKS_NS;

enum{inter,floater};                              // data type 
enum{FE=0,CU,NI,VACANCY,I1,I2,I3,I4,I5,I6};       // diagnosis terms   
enum{hFE=11,hCU,hNI,hMN,hSi,hP,hC};               // hop steps for each element
enum{sink=20};                                    // number of sink absorption   
enum{recombine=31,FPair};                        // number of recombination    
enum{cFE=40,cCU,cNI,cVACANCY,cI1,cI2,cI3,cI4,cI5,cI6};        // time averaged concentration    
enum{dFE=51,dCU,dNI,dVACANCY,dI1,dI2,dI3,dI4,dI5,dI6}; // MSD for each element 
enum{energy=61,treal,fvt};                        // energy and realistic time  
enum{ris=70};                                     // number of ris   
enum{lij=80};                                     // onsager coefficient   
/* ---------------------------------------------------------------------- */

DiagRis::DiagRis(SPPARKS *spk, int narg, char **arg) : Diag(spk,narg,arg)
{
  if (strcmp(app->style,"ris") != 0)
    error->all(FLERR,"Diag_style ris requires app_style ris");

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
    } else error->all(FLERR,"Illegal diag_style ris command");
  }

  if (nlist == 0) error->all(FLERR,"Illegal diag_style ris command");
  which = new int[nlist];
  index = new int[nlist];
  itype = new int[nlist];
  ivector = new int[nlist];
  dvector = new double[nlist+1]; //nelement + total energy 
}

/* ---------------------------------------------------------------------- */

DiagRis::~DiagRis()
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

void DiagRis::init()
{
  appris = (AppRis *) app;

  siteflag = 0;
  csiteflag = 0;
  hopflag = 0;
  msdflag = 0;
  risflag = appris->ris_flag;; 

  for (int i = 0; i < nlist; i++) {
    if (strcmp(list[i],"fe") == 0) which[i] = FE; //total sites
    else if (strcmp(list[i],"cu") == 0) which[i] = CU;
    else if (strcmp(list[i],"ni") == 0) which[i] = NI;
    else if (strcmp(list[i],"vac") == 0) which[i] = VACANCY;
    else if (strcmp(list[i],"i1") == 0) which[i] = I1;
    else if (strcmp(list[i],"i2") == 0) which[i] = I2;
    else if (strcmp(list[i],"i3") == 0) which[i] = I3;
    else if (strcmp(list[i],"i4") == 0) which[i] = I4;
    else if (strcmp(list[i],"i5") == 0) which[i] = I5;
    else if (strcmp(list[i],"i6") == 0) which[i] = I6;

    else if (strcmp(list[i],"cfe") == 0) which[i] = cFE; //time averaged concentration 
    else if (strcmp(list[i],"ccu") == 0) which[i] = cCU;
    else if (strcmp(list[i],"cni") == 0) which[i] = cNI;
    else if (strcmp(list[i],"cvac") == 0) which[i] = cVACANCY;
    else if (strcmp(list[i],"ci1") == 0) which[i] = cI1;
    else if (strcmp(list[i],"ci2") == 0) which[i] = cI2;
    else if (strcmp(list[i],"ci3") == 0) which[i] = cI3;
    else if (strcmp(list[i],"ci4") == 0) which[i] = cI4;
    else if (strcmp(list[i],"ci5") == 0) which[i] = cI5;
    else if (strcmp(list[i],"ci6") == 0) which[i] = cI6;

   /*
    else if (strcmp(list[i],"hfe") == 0) which[i] = hFE;//total hop events
    else if (strcmp(list[i],"hcu") == 0) which[i] = hCU;
    else if (strcmp(list[i],"hni") == 0) which[i] = hNI;
    else if (strcmp(list[i],"hmn") == 0) which[i] = hMN;
    else if (strcmp(list[i],"hp") == 0) which[i] = hP;
    else if (strcmp(list[i],"hc") == 0) which[i] = hC;
   */ 
    else if (strcmp(list[i],"dfe") == 0) which[i] = dFE;//MSD 
    else if (strcmp(list[i],"dcu") == 0) which[i] = dCU;
    else if (strcmp(list[i],"dni") == 0) which[i] = dNI;
    else if (strcmp(list[i],"dvac") == 0) which[i] = dVACANCY; 
    else if (strcmp(list[i],"di1") == 0) which[i] = dI1;
    else if (strcmp(list[i],"di2") == 0) which[i] = dI2;
    else if (strcmp(list[i],"di3") == 0) which[i] = dI3;
    else if (strcmp(list[i],"di4") == 0) which[i] = dI4;
    else if (strcmp(list[i],"di5") == 0) which[i] = dI5;
    else if (strcmp(list[i],"di6") == 0) which[i] = dI6;
 /*
    else if (strcmp(list[i],"mfe") == 0) which[i] = mFE;//MSD 
    else if (strcmp(list[i],"mvac") == 0) which[i] = mVACANCY; 
    else if (strcmp(list[i],"mcu") == 0) which[i] = mCU;
    else if (strcmp(list[i],"mni") == 0) which[i] = mNI;
    else if (strcmp(list[i],"mmn") == 0) which[i] = mMN;
    else if (strcmp(list[i],"mp") == 0) which[i] = mP;
    else if (strcmp(list[i],"mc") == 0) which[i] = mC;
    else if (strcmp(list[i],"msia") == 0) which[i] = mSIA;
 
    else if (strcmp(list[i],"treal") == 0) which[i] = treal;
    else if (strcmp(list[i],"fvt") == 0) which[i] = fvt;
*/
    else if (strcmp(list[i],"recombine") == 0) which[i] = recombine;
    else if (strcmp(list[i],"nfp") == 0) which[i] = FPair;
    else if (strcmp(list[i],"energy") == 0) which[i] = energy;
    else if (list[i][0] == 's' && list[i][1] == 'i' && list[i][2] == 'n' && list[i][3] == 'k') {
      int id = list[i][4] - '0';
      which[i] = sink + id; 
    }

    // ris of element i
    else if (list[i][0] == 'r' && list[i][1] == 'i' && list[i][2] == 's') {
      int id = list[i][3] - '0';
      which[i] = ris + id;
    }

    // onsager coefficients  
    else if (list[i][0] == 'l' && list[i][1] == 'i' && list[i][2] == 'j') {
      int id1 = list[i][3] - '0';
      int id2 = list[i][4] - '0';
      which[i] = lij + id1*10 + id2;
    }

    else error->all(FLERR,"Invalid value setting in diag_style ris");
  }

  for (int i = 0; i < nlist; i++) {
    if (which[i] == FE || which[i] == CU || which[i] == NI || which[i] == VACANCY || which[i] == I1 || which[i] == I2 || which[i] == I3 || which[i] == I4 || which[i] == I5 || which[i] == I6)
      siteflag = 1;
    if (which[i] == cFE || which[i] == cCU || which[i] == cNI || which[i] == cVACANCY || which[i] == cI1 || which[i] == cI2 || which[i] == cI3 || which[i] == cI4 || which[i] == cI5 || which[i] == cI6)
      csiteflag = 1;
  //  if (which[i] == hFE || which[i] == hCU || which[i] == hNI || which[i] == hMN || which[i] == hP || which[i] == hC)
  //    hopflag = 1;
    if(which[i] >= dFE && appris->diffusionflag == 1) msdflag = 1;
  //  if(which[i] >= dFE && which[i] <= dI6 && msdflag == 0) error->all(FLERR,"MSD calculation need displacement calculated in appris"); 
  }

  for (int i = 0; i < nlist; i++) {ivector[i] = 0;
    dvector[i] = 0.0;
  }

  for (int i=0; i < nlist; i++) { itype[i] = inter;
    if(which[i] >= cFE) itype[i] = floater;
  }
}

/* ---------------------------------------------------------------------- */

void DiagRis::compute()
{
  int ninter,nfloater;
  int sites[10],nhop[10],ivalue; // int data
  int nlocal = appris->nlocal;
  int nelement = appris->nelement;
  double dvalue; // double data 
  double *csites; 
  double msd[10];

  ninter = nfloater = 0;
  dvalue = 0.0;
 
  if (risflag) {
     appris->ris_time(); 
  } 

  if (siteflag) {
    sites[FE] = sites[VACANCY] = sites[CU] = sites[NI] = 0;
    sites[I1] = sites[I2] = sites[I3] = sites[I4] = sites[I5] = sites[I6] = 0;
    int *element = appris->element;
    for (int i = 0; i < nlocal; i++) sites[element[i]]++;
  }

  if(csiteflag) {csites = appris->ct;} 
/*
  if(hopflag) {// hop event of each element 
    for(int i = 1; i < nelement+1; i++) {nhop[i] = 0; nhop[i] = appris->hcount[i];}   
  }
*/
  if(msdflag) {// MSD calculation 
    int *element = appris->element;
    msd[FE] = msd[VACANCY] = msd[CU] = msd[NI] = 0.0;
    msd[I1] = msd[I2] = msd[I3] = msd[I4] =msd[I5] =msd[I6]  =0.0;
    if(siteflag == 0) {//need to count total sites if not calculated above 
      sites[FE] = sites[VACANCY] = sites[CU] = sites[NI] = 0;
      sites[I1] = sites[I2] = sites[I3] = sites[I4] = sites[I5] = sites[I6] = 0;
      for (int i = 0; i < nlocal; i++) sites[element[i]]++;
    }

    for(int i = 0; i < nlocal; i++) {
       msd[element[i]] += appris->disp[3][i];
    }   
  }

    //fprintf(screen, "%f %f \n", csites[1],csites[2]);
  for (int i = 0; i < nlist; i++) {
    if (which[i] == FE) ivalue = sites[FE]; //total sites
    else if (which[i] == CU) ivalue = sites[CU];
    else if (which[i] == NI) ivalue = sites[NI];
    else if (which[i] == VACANCY) ivalue = sites[VACANCY];
    else if (which[i] == I1) ivalue = sites[I1]; 
    else if (which[i] == I2) ivalue = sites[I2]; 
    else if (which[i] == I3) ivalue = sites[I3];
    else if (which[i] == I4) ivalue = sites[I4];
    else if (which[i] == I5) ivalue = sites[I5];
    else if (which[i] == I6) ivalue = sites[I6];

    else if (which[i] == cFE) dvalue = csites[FE]; //time_averaged concentration 
    else if (which[i] == cCU) dvalue = csites[CU];
    else if (which[i] == cNI) dvalue = csites[NI];
    else if (which[i] == cVACANCY) dvalue = csites[VACANCY];
    else if (which[i] == cI1) dvalue = csites[I1]; 
    else if (which[i] == cI2) dvalue = csites[I2]; 
    else if (which[i] == cI3) dvalue = csites[I3];
    else if (which[i] == cI4) dvalue = csites[I4];
    else if (which[i] == cI5) dvalue = csites[I5];
    else if (which[i] == cI6) dvalue = csites[I6];
    /*
    else if (which[i] == mVACANCY) ivalue = monomer_local[VACANCY];
    else if (which[i] == mCU) ivalue = monomer_local[CU];
    else if (which[i] == mNI) ivalue = monomer_local[NI];
    else if (which[i] == mMN) ivalue = monomer_local[MN]; 
    else if (which[i] == mP) ivalue = monomer_local[P]; 
    else if (which[i] == mC) ivalue = monomer_local[C];
    else if (which[i] == mSIA) ivalue = monomer_local[SIA];
    
    else if (which[i] == hFE) ivalue = nhop[FE]; //total hop events
    else if (which[i] == hCU) ivalue = nhop[CU];
    else if (which[i] == hNI) ivalue = nhop[NI];
    else if (which[i] == hMN) ivalue = nhop[MN]; 
    else if (which[i] == hP) ivalue = nhop[P]; 
    else if (which[i] == hC) ivalue = nhop[C];
    */
    else if (which[i] == dFE && sites[FE] > 0) dvalue = msd[FE]/sites[FE];//MSD 
    else if (which[i] == dVACANCY && sites[VACANCY] > 0) dvalue = msd[VACANCY]/sites[VACANCY]; 
    else if (which[i] == dCU && sites[CU] > 0) dvalue = msd[CU]/sites[CU]; 
    else if (which[i] == dNI && sites[NI] > 0) dvalue = msd[NI]/sites[NI]; 
    else if (which[i] == dI1 && sites[I1] > 0) dvalue = msd[I1]/sites[I1]; 
    else if (which[i] == dI2 && sites[I2] > 0) dvalue = msd[I2]/sites[I2]; 
    else if (which[i] == dI3 && sites[I3] > 0) dvalue = msd[I3]/sites[I3]; 
    else if (which[i] == dI4 && sites[I4] > 0) dvalue = msd[I4]/sites[I4]; 
    else if (which[i] == dI5 && sites[I5] > 0) dvalue = msd[I5]/sites[I5]; 
    else if (which[i] == dI6 && sites[I6] > 0) dvalue = msd[I6]/sites[I6]; 
   /*
    else if (which[i] == treal) dvalue = appris->realtime; //realistic time  
    else if (which[i] == fvt) dvalue = appris->fvt; //realistic time  
   */
    else if (which[i] == energy) dvalue = appris->total_energy(); //system energy 
    else if (which[i] > sink && which[i] < recombine) { // to be updateed lated 
     // int id = which[i] - sink; 
     // ivalue = appris->nabsorption[id-1];
     ivalue = 0; 
    }
    else if (which[i] == recombine) ivalue = appris->nrecombine[VACANCY]; //number of reocmbination 
    else if (which[i] == FPair) ivalue = appris->nFPair; //number of reocmbination 
    else if (which[i] >= ris && which[i] < lij) { // ris  
      int id = which[i] - ris; 
      dvalue = appris->ris_total[id];
    }
    else if (which[i] >= lij) { // ris  
      int id2 = (which[i] - lij)%10; 
      int id1 = (which[i] - lij - id2)/10; 
      dvalue = appris->Lij[id1][id2]; // calcualted in appris->onsager()
    }

    if(which[i] >= cFE) nfloater++;
    else ninter++;
    MPI_Allreduce(&ivalue,&ivector[ninter],1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&dvalue,&dvector[nfloater],1,MPI_DOUBLE,MPI_SUM,world);
  }
}

/* ---------------------------------------------------------------------- */

void DiagRis::stats(char *str)
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

void DiagRis::stats_header(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str," %s",list[i]);
    str += strlen(str);
  }
}
