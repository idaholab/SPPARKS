/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Class AppPottsPhaseField - added by Eric Homer, ehomer@sandia.gov
   Mar 31, 2011 - Most recent version.  Most of this was copied from
   AppPotts and AppPottsNeighOnly.

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "app_potts_phasefield_marmot.h"
#include "solve.h"
#include "random_park.h"
#include "error.h"
#include "memory.h"
#include "domain.h"
#include "lattice.h"
#include "comm_lattice.h"
#include "timer.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"

using namespace SPPARKS_NS;

// same as in create_sites.cpp and diag_cluster.cpp and lattice.cpp
enum{NONE,LINE_2N,SQ_4N,SQ_8N,TRI,SC_6N,SC_26N,FCC,BCC,DIAMOND,
  FCC_OCTA_TETRA,RANDOM_1D,RANDOM_2D,RANDOM_3D};

/* ---------------------------------------------------------------------- */

AppPottsPhaseFieldMarmot::AppPottsPhaseFieldMarmot(SPPARKS *spk, int narg, char **arg) :
  AppPottsNeighOnly(spk,narg,arg)
{
  ninteger = 2;
  ndouble = 1;
  allow_kmc = 0;
  allow_rejection = 1;
  allow_app_update = 1;
  allow_masking = 0;
  numrandom = 2;
  delpropensity = 2;//need full neighbor lists of the 1st layer ghost sites

  //add the double array
  recreate_arrays();

  //check the number of spins
  if (nspins % 2)
    error->all(FLERR,"AppPottsPhaseFieldMarmot must have even # of spins");

  phaseChangeInt = nspins/2;

  if (narg != 11)
    error->all(FLERR,"Illegal app_style command - AppPottsPhaseFieldMarmot");

  //dt_phasefield is set as a multiple of dt_rkmc in setup_end_app();
  //set the multiple here. ( dt_phasefield = dt_rkmc / dt_phasefield_mult )
  dt_phasefield_mult = atof(arg[2]);

  //set the variables for the energetics and evolution
  gamma=atof(arg[3]);
  M_c=atof(arg[4]);
  kappaC=atof(arg[5]);
  a_1 = a_2 = atof(arg[6]);

  c_1 = atof(arg[7]);
  c_2 = atof(arg[8]);
  c_3 = atof(arg[9]);
  c_4 = atof(arg[10]);

  //set other default values
  nlocal_app=0;
  cnew=NULL;
  pf_resetfield=false;
  warn_concentration_deviation=0;
  warn_concentration_deviation_all=0;
  cmap_ready=false;
  print_cmap=false;
  enforceConcentrationLimits=false;
  initialize_values=false;
  dimension=0;
  cmap=NULL;

}

/* ---------------------------------------------------------------------- */

AppPottsPhaseFieldMarmot::~AppPottsPhaseFieldMarmot()
{
  if (pf_resetfield) {
    memory->sfree(pf_resetlist);
    memory->sfree(pf_resetlistvals);
  }
  if (nlocal_app) {
    memory->sfree(cnew);
    nlocal_app=0;
  }
  if (cmap)
    memory->sfree(cmap);
}

/* ---------------------------------------------------------------------- */

void AppPottsPhaseFieldMarmot::input_app(char *command, int narg, char **arg)
{

  if (strcmp(command,"pottspfmarmot/command") == 0) {

    if (narg != 2) error->all(FLERR,"Invalid command for app_style");

    if (strcmp(arg[0],"reset_phasefield") == 0)
      set_phasefield_resetfield(narg,arg);

    else if (strcmp(arg[0],"print_connectivity") == 0) {
      if (strcmp(arg[1],"yes") == 0) print_cmap=true;
      else print_cmap=false;
    }

    else if (strcmp(arg[0],"initialize_values") == 0) {
      if (strcmp(arg[1],"yes") == 0) initialize_values=true;
      else initialize_values=false;
    }

    else if (strcmp(arg[0],"enforce_concentration_limits") == 0) {
      if (strcmp(arg[1],"yes") == 0) {
        enforceConcentrationLimits=true;
        //since I'm forcing the concentration limits
        //set the warn flags to 1 so that it doesn't check
        warn_concentration_deviation=1;
        warn_concentration_deviation_all=1;
      }
      else {
        enforceConcentrationLimits=false;
        //if someone turns this off reset the flags to check
        // for concentration deviations
        warn_concentration_deviation=0;
        warn_concentration_deviation_all=0;
      }
    }
    else
      error->all(FLERR,"Invalid command for app_style");
  }
  else
    error->all(FLERR,"Invalid command for app_style");
}

/* ---------------------------------------------------------------------- */

void AppPottsPhaseFieldMarmot::init_app()
{
  //init the parent class first
  AppPottsNeighOnly::init_app();

  //setup the connectivity map
  if (!cmap_ready) setup_connectivity_map();

  //initialize values if command is called
  if (initialize_values) init_values();

}

/* ---------------------------------------------------------------------- */

void AppPottsPhaseFieldMarmot::setup_end_app()
{
  //set the time step for phase field
  dt_phasefield = dt_rkmc / dt_phasefield_mult;

  //add code here at a later date to check the stability of the phase field
}

/* ----------------------------------------------------------------------
 set site value ptrs each time iarray/darray are reallocated
 ------------------------------------------------------------------------- */

void AppPottsPhaseFieldMarmot::grow_app()
{
  //set integer pointers
  spin = iarray[0];
  phase = iarray[1];

  //setup the initial phase values
  for (int i=0; i<nlocal+nghost; i++)
    set_site_phase(i);

  //setup double pointers
  c = darray[0];

  //grow cnew locally so it doesn't have to be communicated
  if (nlocal_app < nlocal) {
    nlocal_app = nlocal;

    cnew = (double *)
    memory->srealloc(cnew,nlocal_app*sizeof(double),"app_potts_pf:cnew");
  }
}

/* ----------------------------------------------------------------------
 compute chemical potential of site, dG/dc + kappaC*del^2(c)
 ------------------------------------------------------------------------- */

double AppPottsPhaseFieldMarmot::site_chem_pot(int i)
{
  int j;
  double qi_alpha,qi_beta;
  double del_sq_c = 0;

  if (phase[i]==0) {
    qi_alpha=1.0;
    qi_beta=0.0;
  }
  else {
    qi_alpha=0.0;
    qi_beta=1.0;
  }

  //Calculate dGdc.
  double C=c[i];
  double dGdc = 2.0*gamma*(C-c_1+C-c_2) + 2.0*a_1*qi_alpha*(C-c_3) + 2.0*a_2*qi_beta*(C-c_4);

  // Calculate del^2(c) = c(i-1) - 2*c(i) + c(i+1)
  for (j=0; j<2*dimension; j++) {

    //neighbor
    int s=neighbor[i][cmap[j]];
    del_sq_c += c[s];
  }
  del_sq_c -= 2*dimension*c[i];
  return dGdc - kappaC*del_sq_c;
}

/* ----------------------------------------------------------------------
 compute energy of site
 ------------------------------------------------------------------------- */

double AppPottsPhaseFieldMarmot::site_energy(int i)
{
  //get the energy without the gradient term
  double energy = site_energy_no_gradient(i);

  //Add the gradient energy to obtain total energy, (del c)^2
   double val=0.0;
   //perform finite difference on all the cells in the list
   for (int j=0; j<2*dimension; j++) {

     //site
     int s=neighbor[i][cmap[j]];

     if (j%2==0)
       val -= c[s];
     else
       val += c[s];
   }

   energy += 2.0 * kappaC * pow(val / 2.0,2.0);

  return energy;
}

/* ----------------------------------------------------------------------
 compute energy of site without the gradient term for efficient site event rejection
 ------------------------------------------------------------------------- */

double AppPottsPhaseFieldMarmot::site_energy_no_gradient(int i)
{
  int isite = spin[i];
  int eng = 0;
  for (int j = 0; j < numneigh[i]; j++)
    if (isite != spin[neighbor[i][j]]) eng++;


  // Now add the energy from deltaG
  double energy = (double) eng;
  double C=c[i];

  energy+= (gamma)*( pow(C-c_1,2.0) + pow(c_2-C,2.0) );
  //now add the energy from the alpha or beta phase
  if (phase[i]==0) {
    energy+= (a_1)*pow(C-c_3,2.0);
  } else {
    energy+= (a_2)*pow(c_4-C,2.0);
  }

  /*--------------------------------------------------
   The gradient energy term is not calculated here for
   computational efficiency.  This function is called
   by site_event rejection where the gradient term
   won't change and therefore this is more
   computationally efficient.  The real site energy
   term is calculated in site_energy.
  --------------------------------------------------*/

  return energy;
}

/* ----------------------------------------------------------------------
 rKMC method
 perform a site event with no null bin rejection
 flip to random neighbor spin without null bin
 technically this is an incorrect rejection-KMC algorithm
 ------------------------------------------------------------------------- */

void AppPottsPhaseFieldMarmot::site_event_rejection(int i, RandomPark *random)
{
  int oldstate = spin[i];
  double einitial = site_energy_no_gradient(i);

  // events = spin flips to neighboring site different than self

  int j,m,value;
  int nevent = 0;

  for (j = 0; j < numneigh[i]; j++) {
    value = spin[neighbor[i][j]];
    if (value == spin[i]) continue;
    for (m = 0; m < nevent; m++)
      if (value == unique[m]) break;
    if (m < nevent) continue;
    unique[nevent++] = value;
  }

  if (nevent == 0) return;
  int iran = (int) (nevent*random->uniform());
  if (iran >= nevent) iran = nevent-1;
  spin[i] = unique[iran];
  set_site_phase(i);
  double efinal = site_energy_no_gradient(i);

  // accept or reject via Boltzmann criterion

  if (efinal <= einitial) {
  } else if (temperature == 0.0) {
    spin[i] = oldstate;
    set_site_phase(i);
  } else if (random->uniform() > exp((einitial-efinal)*t_inverse)) {
    spin[i] = oldstate;
    set_site_phase(i);
  }

  if (spin[i] != oldstate) naccept++;
}

/* ----------------------------------------------------------------------
 iterate through the phase field solution
 ------------------------------------------------------------------------- */

void AppPottsPhaseFieldMarmot::user_update(double stoptime)
{
  double localtime=0.0;

  //communicate all sites to make sure it's up-to-date when it starts
  timer->stamp();
  comm->all();
  timer->stamp(TIME_COMM);
}

/* ---------------------------------------------------------------------- */

void AppPottsPhaseFieldMarmot::setup_connectivity_map()
{
  int i,j;

  //this check is redundant but I'm leaving it anyway
  if (domain->lattice->nbasis > 1)
    error->all(FLERR,
      "only single basis units are allowed for app_style potts/phasefield_marmot");

  //set the dimension variable
  dimension = domain->dimension;

  if (!cmap)
    cmap = (int *)
      memory->smalloc(2*dimension*sizeof(int),"app_potts_pf:cmap");


  //setup the connectivity map for the appropriate style
  int style=domain->lattice->style;

  if (style == LINE_2N) {
    cmap[0]=0;// x-1
    cmap[1]=1;// x+1
  }
  else if (style == SQ_4N) {
    cmap[0]=0;// x-1
    cmap[1]=3;// x+1
    cmap[2]=1;// y-1
    cmap[3]=2;// y+1
  }
  else if (style == SQ_8N) {
    cmap[0]=1;// x-1
    cmap[1]=6;// x+1
    cmap[2]=3;// y-1
    cmap[3]=4;// y+1
  }

  else if (style == SC_6N){
    cmap[0]=0;// x-1
    cmap[1]=5;// x+1
    cmap[2]=1;// y-1
    cmap[3]=4;// y+1
    cmap[4]=2;// z-1
    cmap[5]=3;// z+1
  }
  else if (style == SC_26N) {
    cmap[0]=4; // x-1
    cmap[1]=21;// x+1
    cmap[2]=10;// y-1
    cmap[3]=15;// y+1
    cmap[4]=12;// z-1
    cmap[5]=13;// z+1
  }
  else
    error->all(FLERR,
      "Lattice style not compatible with app_style potts/phasefield_marmot");

  //The connectivity map defined above should not change unless
  //create_sites changes. However, the following is an error check
  //to ensure that the connectivity map is correct.

  double lim[3][2];
  lim[0][0]=domain->boxxlo;
  lim[0][1]=domain->boxxhi;
  lim[1][0]=domain->boxylo;
  lim[1][1]=domain->boxyhi;
  lim[2][0]=domain->boxzlo;
  lim[2][1]=domain->boxzhi;

  int done,id=-1;
  for (i=0; i<nlocal; i++) {
    done=0;
    for (j=0; j<dimension; j++)
      if (xyz[i][j]!=lim[j][0] && xyz[i][j]!=lim[j][1])
        done++;

    if (done==dimension) {
      id=i;
      break;
    }
  }

  if (id==-1)
    error->all(FLERR,
      "Error checking connectivity map for app_style potts/phasefield_marmot");

  //now check the ith site because it's not on any of the boundaries
  double pos;
  for (i=0; i<dimension; i++) {
    for (j=0; j<2; j++) {
      if (j==0)
        pos=xyz[id][i]-1;
      else
        pos=xyz[id][i]+1;
      int s=neighbor[id][cmap[2*i+j]];
      if (pos != xyz[s][i]) {
      //  print_connectivity_map();
      //  fprintf(screen,"pos %f xyz[s][i] %f, i = %d j = %d \n", pos, xyz[s][i],i,j);
      //  error->all(FLERR,
      //    "Invalid connectivity map for app_style potts/phasefield_marmot");
      }
    }
  }

  cmap_ready=true;

  //print the connectivity map if it has been asked for
  if (cmap_ready && print_cmap && me==0) print_connectivity_map();
}

/* ---------------------------------------------------------------------- */

void AppPottsPhaseFieldMarmot::print_connectivity_map()
{
  char str[256];
  char *strptr;
  int i,j;

  strptr=str;

  if (dimension==1)
    sprintf(strptr,"  +/-     x\n");
  else if (dimension==2)
    sprintf(strptr,"  +/-     x     y\n");
  else
    sprintf(strptr,"  +/-     x     y     z\n");
  strptr += strlen(strptr);

  //cycle through the appropriate dimensions
  for (j=0; j<2; j++) {
    if (j==0)
      sprintf(strptr,"    -");
    else
      sprintf(strptr,"    +");
    strptr += strlen(strptr);

    for (i=0; i<dimension; i++) {
      sprintf(strptr," %5d",cmap[2*i+j]);
      strptr += strlen(strptr);
    }
    sprintf(strptr,"\n");
    strptr += strlen(strptr);
  }

  if (screen)
    fprintf(screen,
            "Connectivity map for app_style: potts/phasefield_marmot\n%s\n",str);
  if (logfile)
    fprintf(logfile,
            "Connectivity map for app_style: potts/phasefield_marmot\n%s\n",str);
}

/* ---------------------------------------------------------------------- */

void AppPottsPhaseFieldMarmot::set_site_phase(int i)
{
  if ( (spin[i] - phaseChangeInt) > 0 )
    phase[i]=0;//Alpha
  else
    phase[i]=1;//Beta
}

/* ---------------------------------------------------------------------- */

void AppPottsPhaseFieldMarmot::set_phasefield_resetfield(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal pottspfmarmot/command reset_phasefield");
  if (strcmp(arg[1],"yes") == 0)
    pf_resetfield = true;
  else if (strcmp(arg[1],"no") == 0) //redundant but allows error catching
    pf_resetfield = false;
  else
    error->all(FLERR,"Illegal pottspfmarmot/command reset_phasefield - yes/no");

  //setup array for easy reset
  if (pf_resetfield) {
    //get the count
    pf_nresetlist=0;
    for (int i=0;i<nlocal;i++)
      if ( ( xyz[i][0] <= domain->boxxlo ) ||
           ( xyz[i][0] >= (domain->boxxhi - 1) ) )
        pf_nresetlist++;

    //setup lists
    pf_resetlist = (int *) memory->smalloc(
        pf_nresetlist*sizeof(int),"app_potts_pf:pf_resetlist");
    pf_resetlistvals = (double *) memory->smalloc(
        pf_nresetlist*sizeof(double),"app_potts_pf:pf_resetlistvals");

    //set values in the list
    int counter=0;
    for (int i=0;i<nlocal;i++) {
      if ( xyz[i][0] <= domain->boxxlo ) {
        pf_resetlist[counter]=i;
        pf_resetlistvals[counter++]=0.0;
      } else if ( xyz[i][0] >= (domain->boxxhi - 1) ) {
        pf_resetlist[counter]=i;
        pf_resetlistvals[counter++]=1.0;
      }
    }
    if (counter != pf_nresetlist)
      error->all(FLERR,"Problem setting up pfsolver reset field lists");
  }
}

/* ---------------------------------------------------------------------- */

void AppPottsPhaseFieldMarmot::init_values()
{
  double q_alpha,q_beta,cAlpha,cBeta;

  //phase=0 - alpha (a_1,c_3)
  q_alpha = 1.0; q_beta = 0.0;
  cAlpha=(gamma * (c_1 + c_2) + a_1*c_3*q_alpha + a_2*c_4*q_beta )
    / (2.0*gamma + a_1*q_alpha + a_2*q_beta);

  //phase=1 - beta  (a_2,c_4)
  q_alpha = 0.0; q_beta = 1.0;
  cBeta =(gamma * (c_1 + c_2) + a_1*c_3*q_alpha + a_2*c_4*q_beta )
    / (2.0*gamma + a_1*q_alpha + a_2*q_beta);

  for (int i=0;i<nlocal;i++) {
    set_site_phase(i);
    if (phase[i]==0)
      c[i] = cAlpha;
    else
      c[i] = cBeta;
  }
  initialize_values=false;//ensure this doesn't get called again
}
