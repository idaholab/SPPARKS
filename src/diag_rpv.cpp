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
#include "diag_rpv.h"
#include "app.h"
#include "app_rpv.h"
#include "comm_lattice.h"
#include "timer.h"
#include "error.h"
#include "memory.h"

using namespace SPPARKS_NS;

enum{inter,floater};                          // data type 
enum{ZERO,FE,VACANCY,CU,NI,MN,Si,P,C,SIA};      // diagnosis terms   
enum{hFE=11,hCU,hNI,hMN,hSi,hP,hC}; // hop steps for each element
//enum{mFE=21,mVACANCY,mCU,mNI,mMN,mP,mC,mSIA}; // number of solute (monomers) in the matrix  
enum{sink=30}; // number of sink absorption   
enum{recombine=41}; // number of sink absorption   
enum{dFE=61,dVACANCY,dCU,dNI,dMN,dSi,dP,dC,dSIA}; // MSD for each element 
enum{energy=71,treal,fvt}; //energy and realistic time  
/* ---------------------------------------------------------------------- */

DiagRpv::DiagRpv(SPPARKS *spk, int narg, char **arg) : Diag(spk,narg,arg)
{
  if (strcmp(app->style,"rpv") != 0)
    error->all(FLERR,"Diag_style rpv requires app_style rpv");

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
    } else error->all(FLERR,"Illegal diag_style rpv command");
  }

  if (nlist == 0) error->all(FLERR,"Illegal diag_style rpv command");
  which = new int[nlist];
  index = new int[nlist];
  itype = new int[nlist];
  ivector = new int[nlist];
  dvector = new double[nlist+1]; //nelement + total energy 
}

/* ---------------------------------------------------------------------- */

DiagRpv::~DiagRpv()
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

void DiagRpv::init()
{
  apprpv = (AppRpv *) app;

  siteflag = 0;
  hopflag = 0;
  msdflag = 0;

  for (int i = 0; i < nlist; i++) {
    if (strcmp(list[i],"fe") == 0) which[i] = FE; //total sites
    else if (strcmp(list[i],"vac") == 0) which[i] = VACANCY;
    else if (strcmp(list[i],"cu") == 0) which[i] = CU;
    else if (strcmp(list[i],"ni") == 0) which[i] = NI;
    else if (strcmp(list[i],"mn") == 0) which[i] = MN;
    else if (strcmp(list[i],"p") == 0) which[i] = P;
    else if (strcmp(list[i],"c") == 0) which[i] = C;
    else if (strcmp(list[i],"sia") == 0) which[i] = SIA;
   
    else if (strcmp(list[i],"hfe") == 0) which[i] = hFE;//total hop events
    else if (strcmp(list[i],"hcu") == 0) which[i] = hCU;
    else if (strcmp(list[i],"hni") == 0) which[i] = hNI;
    else if (strcmp(list[i],"hmn") == 0) which[i] = hMN;
    else if (strcmp(list[i],"hp") == 0) which[i] = hP;
    else if (strcmp(list[i],"hc") == 0) which[i] = hC;
 
    else if (strcmp(list[i],"dfe") == 0) which[i] = dFE;//MSD 
    else if (strcmp(list[i],"dvac") == 0) which[i] = dVACANCY; 
    else if (strcmp(list[i],"dcu") == 0) which[i] = dCU;
    else if (strcmp(list[i],"dni") == 0) which[i] = dNI;
    else if (strcmp(list[i],"dmn") == 0) which[i] = dMN;
    else if (strcmp(list[i],"dp") == 0) which[i] = dP;
    else if (strcmp(list[i],"dc") == 0) which[i] = dC;
    else if (strcmp(list[i],"dsia") == 0) which[i] = dSIA; 
/* 
    else if (strcmp(list[i],"mfe") == 0) which[i] = mFE;//MSD 
    else if (strcmp(list[i],"mvac") == 0) which[i] = mVACANCY; 
    else if (strcmp(list[i],"mcu") == 0) which[i] = mCU;
    else if (strcmp(list[i],"mni") == 0) which[i] = mNI;
    else if (strcmp(list[i],"mmn") == 0) which[i] = mMN;
    else if (strcmp(list[i],"mp") == 0) which[i] = mP;
    else if (strcmp(list[i],"mc") == 0) which[i] = mC;
    else if (strcmp(list[i],"msia") == 0) which[i] = mSIA;
*/ 
    else if (strcmp(list[i],"recombine") == 0) which[i] = recombine;
    else if (strcmp(list[i],"energy") == 0) which[i] = energy;
    else if (strcmp(list[i],"treal") == 0) which[i] = treal;
    else if (strcmp(list[i],"fvt") == 0) which[i] = fvt;

    else if (list[i][0] == 's' && list[i][1] == 'i' && list[i][2] == 'n' && list[i][3] == 'k') {
      int id = list[i][4] - '0';
      which[i] = sink + id; 
    }

    else error->all(FLERR,"Invalid value setting in diag_style rpv");
  }

  for (int i = 0; i < nlist; i++) {
    if (which[i] == FE || which[i] == VACANCY || which[i] == SIA || which[i] == CU || which[i] == NI || which[i] == MN || which[i] == P || which[i] == C)
      siteflag = 1;
    if (which[i] == hFE || which[i] == hCU || which[i] == hNI || which[i] == hMN || which[i] == hP || which[i] == hC)
      hopflag = 1;
    if(which[i] >= dFE && apprpv->diffusionflag == 1) msdflag = 1;
    if(which[i] >= dFE && which[i] <= dSIA && msdflag == 0) error->all(FLERR,"MSD calculation need displacement calculated in apprpv"); 
  }

  for (int i = 0; i < nlist; i++) { ivector[i] = 0;
    dvector[i] = 0.0;
  }

  for (int i=0; i < nlist; i++) { itype[i] = inter;
    if(which[i] >= dFE) itype[i] = floater;
  }
}

/* ---------------------------------------------------------------------- */

void DiagRpv::compute()
{
  int ninter,nfloater;
  int sites[10],nhop[10],ivalue; // int data
  int estyle = apprpv->engstyle; // energy style 
  int nlocal = apprpv->nlocal;
//  int *monomer_local = apprpv->monomer_local; // monomers
  double dvalue ; // double data 
  double msd[10] ;

  ninter = nfloater = 0;
  dvalue = 0.0;
 
  if (siteflag) {
    sites[FE] = sites[VACANCY] = sites[CU] = sites[NI] = 0;
    sites[MN] = sites[P] = sites[C] = sites[SIA] = 0;
    int *element = apprpv->element;
    for (int i = 0; i < nlocal; i++) sites[element[i]]++;
  }

  if(hopflag) {
    for(int i = 0; i < nlist; i++) {nhop[i] = 0; nhop[i] = apprpv->hcount[i];}   
  }

  if(msdflag) {// MSD calculation 
    int *element = apprpv->element;
    msd[FE] = msd[VACANCY] = msd[CU] = msd[NI] = 0.0;
    msd[MN] = msd[P] = msd[C] = msd[SIA] =0.0;
    if(siteflag == 0) {//need to count total sites if not calculated above 
      sites[FE] = sites[VACANCY] = sites[CU] = sites[NI] = 0;
      sites[MN] = sites[P] = sites[C] = sites[SIA] = 0;
      for (int i = 0; i < nlocal; i++) sites[element[i]]++;
    }

    for(int i = 0; i < nlocal; i++) {
       msd[element[i]] += apprpv->disp[3][i];
    }   
  }

  for (int i = 0; i < nlist; i++) {
    if (which[i] == FE) ivalue = sites[FE]; //total sites
    else if (which[i] == VACANCY) ivalue = sites[VACANCY];
    else if (which[i] == CU) ivalue = sites[CU];
    else if (which[i] == NI) ivalue = sites[NI];
    else if (which[i] == MN) ivalue = sites[MN]; 
    else if (which[i] == P) ivalue = sites[P]; 
    else if (which[i] == C) ivalue = sites[C];
    else if (which[i] == SIA) ivalue = sites[SIA];
    /*
    else if (which[i] == mVACANCY) ivalue = monomer_local[VACANCY];
    else if (which[i] == mCU) ivalue = monomer_local[CU];
    else if (which[i] == mNI) ivalue = monomer_local[NI];
    else if (which[i] == mMN) ivalue = monomer_local[MN]; 
    else if (which[i] == mP) ivalue = monomer_local[P]; 
    else if (which[i] == mC) ivalue = monomer_local[C];
    else if (which[i] == mSIA) ivalue = monomer_local[SIA];
    */
    else if (which[i] == hFE) ivalue = nhop[FE]; //total hop events
    else if (which[i] == hCU) ivalue = nhop[CU];
    else if (which[i] == hNI) ivalue = nhop[NI];
    else if (which[i] == hMN) ivalue = nhop[MN]; 
    else if (which[i] == hP) ivalue = nhop[P]; 
    else if (which[i] == hC) ivalue = nhop[C];

    else if (which[i] == dFE && sites[FE] > 0) dvalue = msd[FE]/sites[FE];//MSD 
    else if (which[i] == dVACANCY && sites[VACANCY] > 0) dvalue = msd[VACANCY]/sites[VACANCY]; 
    else if (which[i] == dCU && sites[CU] > 0) dvalue = msd[CU]/sites[CU]; 
    else if (which[i] == dNI && sites[NI] > 0) dvalue = msd[NI]/sites[NI]; 
    else if (which[i] == dMN && sites[MN] > 0) dvalue = msd[MN]/sites[MN]; 
    else if (which[i] == dP && sites[P] > 0) dvalue = msd[P]/sites[P]; 
    else if (which[i] == dC && sites[C] > 0) dvalue = msd[C]/sites[C]; 
    else if (which[i] == dSIA && sites[SIA] > 0) dvalue = msd[SIA]/sites[SIA]; 
   
    else if (which[i] == recombine) dvalue = apprpv->nrecombine; //number of reocmbination 
    else if (which[i] == energy) dvalue = apprpv->total_energy(); //system energy 
    else if (which[i] == treal) dvalue = apprpv->realtime; //realistic time  
    else if (which[i] == fvt) dvalue = apprpv->fvt; //realistic time  
   
    else if (which[i] > sink && which[i] < recombine) {
      int id = which[i] - sink; 
      ivalue = apprpv->nabsorption[id-1];
    }
 
    if(which[i] >= dFE) nfloater++;
    else ninter++;

    MPI_Allreduce(&ivalue,&ivector[ninter],1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&dvalue,&dvector[nfloater],1,MPI_DOUBLE,MPI_SUM,world);
  }
}

/* ---------------------------------------------------------------------- */

void DiagRpv::stats(char *str)
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

void DiagRpv::stats_header(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str," %s",list[i]);
    str += strlen(str);
  }
}
