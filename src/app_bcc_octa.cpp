/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
*************************************************************************************
   This application does the vacancy diffusion in RPV steel, via hop in a BCC lattice
   and recombination with SIAs. Barriers for diffusion was calcualted by saddle point
   bond energies. Contributer: Yongfeng Zhang, yongfeng.zhang@inl.gov
------------------------------------------------------------------------- */

#include "math.h"
#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "app_bcc_octa.h"
#include "domain.h"
#include "solve.h"
#include "random_park.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

enum{NOOP,BCC,OCTA};                                // all sites are BCC or OCTA
enum{ZERO,FE,VACANCY,SB,HE,VO,I1,I2,I3,I4,I5,I6};       // same as DiagBCCOCTA; element

#define DELTAEVENT 100000
#define MAX2NN 18  // max 2NN of BCC lattice
#define BIGNUMBER 1e18 // define a big number 

/* ---------------------------------------------------------------------- */

AppBccOcta::AppBccOcta(SPPARKS *spk, int narg, char **arg) :
  AppLattice(spk,narg,arg)
{
  ninteger = 3; // first for lattice type; second for element; third for cluster analysis  
  ndouble = 0;
  delpropensity = 2;
  delevent = 1;
  allow_kmc = 1;
  allow_rejection = 0;

  engstyle = 1; //1 for 1NN interaction, 2 for 2NN interaction; default 1
  rhop = 1; // default only 1NN hop 
  rrecombine = 1; // recombination radius, default 1
  cflag = clst_flag = 0; // if activate concentration dependent VV bond, default 0 (no); 
  if (narg < 1) error->all(FLERR,"Illegal app_style command");
  if (narg >= 2) engstyle = atoi(arg[1]);
  if (narg >= 3) rhop = atoi(arg[2]);
  if (narg >= 4) rrecombine = atoi(arg[3]);
  if (rrecombine == 4) { 
     if (narg < 5) error->all(FLERR,"Recombination radius is required but not supplied");
     rrec = atof(arg[4]);
     if (narg >= 6) cflag = atoi(arg[5]);
  } else if (narg >= 5) cflag = atoi(arg[4]);
  if (engstyle == 2) delpropensity += 1;// increase delpropensity for 2NN interaction
  if (rhop == 3) delpropensity +=1; 

  //ninteger ++; // for recombination vector or cluster analysis, test only  
  create_arrays();

  dmigration = 3;
  dratio = 1.0; 
  firsttime = 1;
  esites = NULL;
  echeck = NULL;
  events = NULL;
  maxevent = 0;
  firstevent = NULL;

  // set random number generator
  ranbccocta = new RandomPark(ranmaster->uniform());
  double seed = ranmaster->uniform();
  ranbccocta->reset(seed,me,100);

  // flags for bond interactions
  ebond1 = NULL;
  ebond2 = NULL;
  mbarrier = NULL;
  hcount = NULL; //numner of hop events for mobile elements
  nn1flag = nn2flag = mfpflag = barrierflag = 0; //flags for bond energy and migration barriers

  // flags and parameters for sinks, dislocations, reactions and ballistic mixing
  sink_flag = reaction_flag = 0; //flags for sink dislocation and vacancy
  nsink = nreaction = nballistic = 0;

  // arrays for sinks
  sink_shape = sink_type = sink_normal = sink_segment = nabsorption = NULL;
  sink_strength = sink_radius = sink_mfp = NULL;
  xsink = NULL;

  // arrays for reactions
  rsite = rinput = routput = rcount = renable = rtarget = NULL;
  nsites_local = target_local = target_global =NULL; //number of reactions
  rbarrier = rrate = NULL;

  // arrays for ballistic mixing
  bfreq = NULL; 
  time_old = time_new = NULL;
  min_bfreq = BIGNUMBER;

  // 2NN neighbor information
  ibox = NULL;
  numneigh2 = numneigh3 = numneigh4 = NULL;
  neighbor2 = neighbor3 = neighbor4 = NULL;

  // cluster analysis 
  ncelem = cid = csize = NULL;

  // He occupation of a vacancy 
  iHe=NULL;
}

/* ---------------------------------------------------------------------- */

AppBccOcta::~AppBccOcta()
{
  delete [] esites;
  delete [] echeck;
  delete ranbccocta;

  memory->sfree(events);
  memory->destroy(firstevent);
  memory->destroy(ebond1);
  memory->destroy(ebond2);
  memory->destroy(mbarrier);
  memory->destroy(nsites_local);
  memory->destroy(hcount);
  memory->destroy(iHe);

  if (engstyle == 2) {// memory use for 2NNs
    memory->destroy(numneigh2);
    memory->destroy(neighbor2);
    if(rhop == 3 || rrecombine == 3) memory->destroy(numneigh3);
    if(rhop == 3 || rrecombine == 3) memory->destroy(neighbor3);
  }

  if (rrecombine == 4) {
    memory->destroy(numneigh4);
    memory->destroy(neighbor4);
  }

  if (sink_flag) { // memory use related to sink
    memory->destroy(sink_shape);
    memory->destroy(sink_type);
    memory->destroy(sink_strength);
    memory->destroy(xsink);
    memory->destroy(isink);
    memory->destroy(sink_normal);
    memory->destroy(sink_segment);
    memory->destroy(sink_radius);
    memory->destroy(sink_mfp);
    memory->destroy(nabsorption);
  }

  if (ballistic_flag) { // memory use related to ballistic mixing
    memory->destroy(bfreq);
    memory->destroy(time_old);
    memory->destroy(time_new);
  }

  if (reaction_flag) {// memory use related to reaction
    memory->destroy(rsite);
    memory->destroy(rinput);
    memory->destroy(routput);
    memory->destroy(rcount);
    memory->destroy(renable);
    memory->destroy(rtarget);
    memory->destroy(rbarrier);
    memory->destroy(rrate);
    memory->destroy(target_local);
    memory->destroy(target_global);
  }

  if (mfpflag) {// mean free path for absorption  
    memory->destroy(rhmfp);
    memory->destroy(nhmfp);
    memory->destroy(mfp);
  }

}

/* ----------------------------------------------------------------------
  setup bond energies for 1NN, 2NN, and SPs
------------------------------------------------------------------------- */

void AppBccOcta::input_app(char *command, int narg, char **arg)
{
  int i,j,ibond;
  int nlattice = nlocal + nghost;

  // initial coordiation of each atom, for recombination vector analysis, test only 
  //memory->create(xyz0,nlocal,3,"app/bccocta:xyz0");
  memory->create(iHe,nlocal,"app/bccocta:iHe");

  // 1NN bond energy taken in the order of 11 12 ... 1N; 22 ... 2N; ...; NN
  if (strcmp(command,"ebond1") == 0) {

    if (narg < 3) error->all(FLERR,"Illegal ebond1 command");
    nelement = atoi(arg[0]);   // num of elements

    memory->create(nsites_local,nelement,"app/bccocta:nsites_local");
    memory->create(ebond1,nelement+1,nelement+1,"app/bccocta:ebond1");
    memory->create(hcount,nlattice,"app/bccocta:hcount");

    if(narg != nelement*(nelement+1)/2+1) error->all(FLERR,"Illegal ebond command");
    nn1flag = 1;

    for (i = 1; i <= nelement; i++ ) {
      for (j = i; j <= nelement; j++ ) {
        ibond = ibonde(i,j,nelement);
        ebond1[i][j] = atof(arg[ibond]);
        if(j > i) ebond1[j][i] = ebond1[i][j];
      }
    }
  }

  // 2NN bond energy taken in the order of 11 12 ... 1N; 22 ... 2N; ...; NN
  else if (strcmp(command,"ebond2") == 0) {

    nelement = atoi(arg[0]);   // num of elements
    memory->create(ebond2,nelement+1,nelement+1,"app/bccocta:ebond2");
    if(narg != nelement*(nelement+1)/2+1) error->all(FLERR,"Illegal ebond command");

    nn2flag = 1;
    for (i = 1; i <= nelement; i++ ) {
      for (j = i; j <= nelement; j++ ) {
        ibond = ibonde(i,j,nelement);
        ebond2[i][j] = atof(arg[ibond]);
        if (j > i) ebond2[j][i] = ebond2[i][j];
      }
    }
  }

  // migration barriers for each element
  else if (strcmp(command, "migbarrier") ==0) {

    if (narg < 2 || narg % 2 != 0) error->all(FLERR,"Illegal migbarrier command");
    barrierflag = 1;
    memory->create(mbarrier,nelement+1,"app/bccocta:mbarrier");

    for (i=0; i<narg-1; i++) {
      if(i % 2 == 0){ j = atoi(arg[i]);
        mbarrier[j] = atof(arg[i+1]);
      }
    }
  }

  // define sinks to defects, could be dislocations or interfaces or 3D regions
  else if (strcmp(command, "sink") ==0) {

    if(narg != 10) error->all(FLERR,"illegal sink command");
    if(sink_flag == 0) {
      memory->create(isink, nlattice, nelement,"app/bccocta:isink");
      for(i = 0; i < nlattice; i++) {
        for(j = 0; j < nelement; j++)  isink[i][j] = 0;
      }
    }

    sink_flag = 1;
    grow_sinks();

    sink_type[nsink] = atoi(arg[0]); // sink to certain element
    sink_strength[nsink] = atof(arg[1]); // thickness of sink
    sink_shape[nsink] = atoi(arg[2]); // 1 dislocation, 2 interface, 3 3D region
    xsink[nsink][0] = atof(arg[3]); // coordinaiton of sink center
    xsink[nsink][1] = atof(arg[4]);
    xsink[nsink][2] = atof(arg[5]);
    sink_normal[nsink] = atoi(arg[6]); // normal of planar sinks
    sink_radius[nsink] = atof(arg[7]); // radius for circular or polygonal sinks
    sink_segment[nsink] = atoi(arg[8]); // # of segment for polygon sinks
    sink_mfp[nsink] = atof(arg[9]); // mean free path in this sink
    nabsorption[nsink] = 0; // initialize number of absorptions

    sink_creation(nsink); //create the nth sink, can overlap with other sinks
    nsink ++;
  }

  // reactions for absorption and emission
  else if (strcmp(command, "reaction") ==0) {

    if(narg != 6) error->all(FLERR,"illegal reaction command");
    reaction_flag = 1;
    grow_reactions(); // grow reation list

    rsite[nreaction] = atoi(arg[0]); // reaciton site: type of lattice site
    rinput[nreaction] = atoi(arg[1]); // input element
    routput[nreaction] = atoi(arg[2]); // output element
    rbarrier[nreaction] = atof(arg[3]); // reaction barrier
    rrate[nreaction] = atof(arg[4]); // reaction rate
    rtarget[nreaction] = atoi(arg[5]); // target number of output element

    nreaction ++;
  }

  else if (strcmp(command, "ballistic") ==0) {
    if(narg < 2) error->all(FLERR,"illegal ballistic command");
    ballistic_flag = 1;
    grow_ballistic();

    double dose_rate=atof(arg[0]);// dose rate 
    bfreq[nballistic] = 1e12/nlocal/dose_rate; // time interval to introduce an FP 
    if(min_bfreq > bfreq[nballistic]) min_bfreq = bfreq[nballistic];
    gas_rate = atof(arg[1]); // gas/V ratio 
    nballistic ++; // number of mixing events
  }

  else if (strcmp(command, "migration_vector") ==0) {
    if(nelement < 3 || narg != (nelement-I1+1)*3 + 2) error->all(FLERR,"illegal migration vector command");
    dmigration = atoi(arg[0]); // dimension of migration, default 3D  
    dratio = atof(arg[1]); // ratio of anisotropic diffusion   
    j = I1-1;
    for (i=0; i<narg-2; i++) {
      int k = i % 3;
      if(k == 0) j++;
      Vmig[j][k] = atof(arg[i+2]);
    }

    for (i=I1; i<= nelement; i++) {
      Vmig[i][3] = Vmig[i][0]*Vmig[i][0] + Vmig[i][1]*Vmig[i][1] + Vmig[i][2]*Vmig[i][2];
    }
  }

  // meam hop steps before absorption by external sink 
  else if (strcmp(command, "mfp") ==0) {

    if (narg < 2 || narg % 2 == 0) error->all(FLERR,"Illegal mfp command");

    mfpflag = 1;
    memory->create(mfp,nelement+1,"app/bccocta:mfp");
    memory->create(nhmfp,nelement+1,"app/bccocta:nhmfp");
    memory->create(rhmfp,nelement+1,"app/bccocta:rhmfp");

    for (i=0; i<=nelement; i++) mfp[i] = -1.0; 

    sigmamfp = atof(arg[0]); 
    for (i=1; i<narg-1; i++) {
      if(i % 2 == 1) { j = atoi(arg[i]);
        mfp[j] = atof(arg[i+1]);
      }
    }
  }
  // meam hop steps before absorption by external sink 
  else if (strcmp(command, "cluster") ==0) {

    if (narg < 2) error->all(FLERR,"Illegal mfp command");
    
    clst_flag = 1;
    nctype = atoi(arg[0]);
    ncsize = atoi(arg[1]);;
    memory->create(ncelem,nctype,"app/bccocta:ncelem");
    for (i=0; i<nctype; i++) ncelem[i] = atoi(arg[i+2]);
  }

  else error->all(FLERR,"Unrecognized command");
}


/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
   This is before set commands are executed.  
------------------------------------------------------------------------- */

void AppBccOcta::grow_app()
{
  // i1 needs to be redefined to separete bcc and octa sites
  type = iarray[0];   // lattice type; i1 in input
  element = iarray[1];  // element type; i2 in input
  for (int i = 0; i < nlocal; i++) {
      if(numneigh[i] == 32) type[i] = BCC;
      if(numneigh[i] == 18) type[i] = OCTA; // include 2NN OCTA and BCC 
  } 
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppBccOcta::init_app()
{
  int i,j;

  // separate second nearest neighbors from the neighbor list, required here
  define_2NN();

  if (firsttime) {
    firsttime = 0;

    echeck = new int[nlocal];
    memory->create(firstevent,nlocal,"app:firstevent");

    // esites must be large enough for 3 sites and their 1st neighbors
    esites = new int[3 + 3*maxneigh];
  }

  int flag = 0;
  for (i = 0; i < nelement; i++) nsites_local[i] = 0;
  // setup neighbor list for large recombination redius 
  // recombine V and I if created on 1NN sites 
  if(rrecombine == 4) recombine_list(rrec);
  for (i = 0; i < nelement; i++) nrecombine[i] = 0;
  for (i = 0; i < nlocal; i++) {iHe[i] = 0; recombine(i);} 

  // site validity
  for (i = 0; i < nlocal; i++) {
    if (type[i] < BCC || type[i] > OCTA) flag = 1;
    if (element[i] < FE || element[i] > I6) flag = 1;
    if (element[i] == FE && type[i] == OCTA) error->all(FLERR,"site value incorrect");
    nsites_local[element[i]-1]++;
  }

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");

  // check if reactions need to be enabled or disabled
  if(reaction_flag) {
    for(i = 0; i< nreaction; i++) {
      rcount[i] = 0;
      renable[i] = 0;
      target_global[i] = 0;
      target_local[i] = 0;
    }
  }

  // initialize the time_list for ballistic mixing
  if(ballistic_flag) {
    nFPair = 0;
    csia = 0.0;
    for(i = 0; i < nballistic; i ++) {
       time_old[i] = 0;
       time_new[i] = 0;
    }
  }

  // initialize mfp calculations 
  if(mfpflag) {
    for (i = 0; i < nlocal+nghost; i++) hcount[i] = 0;
    for (i = 0; i <= nelement; i++) {
        rhmfp[i] = 0.0;
        nhmfp[i] = 0;
    }

    varmfp = 2.0*sigmamfp; //sqrt(2*pi)*sigma 
    sigmamfp *= sigmamfp;
    sigmamfp *= 2.0; //2*sigma^2   
  }

}

/* ---------------------------------------------------------------------- */
/* ---------------------------initiation of looing-up list -------------- */

void AppBccOcta::setup_app()
{
      
  for (int i = 0; i < nlocal; i++) echeck[i] = 0;

  // clear event list
  nevents = 0;
  for (int i = 0; i < nlocal; i++) firstevent[i] = -1;
  for (int i = 0; i < maxevent; i++) events[i].next = i+1;
  freeevent = 0;

  // set propensities from rates
  if(temperature == 0.0) error->all(FLERR,"Temperature cannot be 0.0 for app bccocta");
  if(nn1flag == 0) error->all(FLERR, "First neareast neighbor bond not defined: AppBccOcta");
  if(nn2flag == 0 && engstyle == 2) error->all(FLERR, "Second nearest neighbor bond not defined: AppBccOcta");
  if(barrierflag == 0) error->warning(FLERR, "Diffusion barrier not defined: AppBccOcta");

  double KB = 0.00008617;
  KBT = temperature * KB;

}

/* ----------------------------------------------------------------------
   separate the total neighbors into 1NN and 2NN. For bcc sites, there are 
   8 bcc and 6 octa 1NNs, and 6 bcc and 12 octa 2NNs. For octa sites, there 
   are 2 bcc and 4 octa 1NNs, and 4 bcc and 8 octa 2NNs. All are defined as 
   1NN in create_site.cpp     
------------------------------------------------------------------------- */

void AppBccOcta::define_2NN()
{ 
  int i,j,jd,n1nn,n2nn;
  double rij[4],cutoff,cutoff1nn[2][2];

  memory->create(numneigh2,nlocal,"app/bccocta, numneigh2");
  memory->create(neighbor2,nlocal,MAX2NN,"app/bccocta, neighbor2");

  double epsln = 0.001;
  cutoff1nn[0][0] = 3.0/4.0; // squared, rij*rij; needs to be multiplied lattice parameter square if other than 1.0  
  cutoff1nn[0][1] = 0.25;   
  cutoff1nn[1][0] = 0.25;   
  cutoff1nn[1][1] = 0.25;
   
  for (i = 0; i < nlocal; i++) {
      numneigh2[i] = 0;
      n1nn = n2nn = 0;

      for (j = 0; j < numneigh[i]; j++) {
          jd = neighbor[i][j];
          distanceIJ(i,jd, rij);
          
          if (type[i] == 1 && type[jd] == 1) cutoff = cutoff1nn[0][0];
          else if (type[i] == 1 && type[jd] == 2) cutoff = cutoff1nn[0][1];
          else if (type[i] == 2 && type[jd] == 1) cutoff = cutoff1nn[1][0];
          else if (type[i] == 2 && type[jd] == 2) cutoff = cutoff1nn[1][1];
          else error->all(FLERR,"Lattice type defined incorrectly!");

          rij[3] = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
          if (rij[3] < cutoff + epsln) { 
             neighbor[i][n1nn] = jd; 
             n1nn ++;
          } else {
             neighbor2[i][n2nn] = jd; 
             n2nn ++;
          } 
      } 

      numneigh[i] = n1nn;
      numneigh2[i] = n2nn;
      if(n2nn > MAX2NN) error->all(FLERR,"Total 2NN exceeds MAX2NN!");
  }
}

/* ----------------------------------------------------------------------
   setup neighbor list for large recombination radius
   currently for serial runs only  
------------------------------------------------------------------------- */

void AppBccOcta::recombine_list(double rrec)
{
  int i,j,id; 
  int ****ibox; 
  int ncell[3];
  int index[3],cell[3];
  double rij[4],lprd[3],xlo[3],rcell[3];

  int max4nn = 0;
  rrec *= rrec;
  // calculate # of neighbors within rrec using site 0  
  for(i = 0; i < nlocal; i++) {
     distanceIJ(i,0,rij);
     rij[3] = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
     if(rij[3] > rrec) continue;
     max4nn ++;
  }   
 
  memory->create(numneigh4,nlocal,"app/bccocta, numneigh4");
  memory->create(neighbor4,nlocal,max4nn,"app/bccocta, neighbor4");

  for(i = 0; i < nlocal; i++) {
     numneigh4[i] = 0;
     for(j = 0; j < max4nn; j++) neighbor4[i][j] =0;
  }

  // setup list for each lattice site
  lprd[0] = domain->xprd;
  lprd[1] = domain->yprd;
  lprd[2] = domain->zprd;
  for(j = 0; j < 3; j++) {
     ncell[j] = static_cast<int> (lprd[j] / sqrt(rrec));
     rcell[j] = lprd[j]/ncell[j];
  }  
  
  xlo[0] = domain->boxxlo;
  xlo[1] = domain->boxylo;
  xlo[2] = domain->boxzlo;

  memory->create(ibox,ncell[0],ncell[1],ncell[2],2*max4nn,"appbccocta: ibox");
  for(i = 0; i < nlocal; i++) {
     for(j = 0; j < 3; j++) index[j] = static_cast<int>((xyz[i][j] - xlo[j]) / rcell[j]);
     ibox[index[0]][index[1]][index[2]][0] ++;
     int k = ibox[index[0]][index[1]][index[2]][0];
     ibox[index[0]][index[1]][index[2]][k] = i;
  } 

  for(i = 0; i < nlocal; i++) {
     for(j = 0; j < 3; j++) index[j] = static_cast<int>(xyz[i][j] - xlo[j]) / rcell[j];
     
     for(int m = -1; m < 2; m++){  
     for(int n = -1; n < 2; n++){  
     for(int l = -1; l < 2; l++){ 
        cell[0] = index[0] + m;
        cell[1] = index[1] + n;
        cell[2] = index[2] + l;
  
        for(j = 0; j < 3; j++) {
           if(cell[j] < 0) cell[j] = ncell[j]-1; 
           if(cell[j] > ncell[j]-1) cell[j] = 0; 
        }

        for(j = 1; j <= ibox[cell[0]][cell[1]][cell[2]][0]; j ++){
           id = ibox[cell[0]][cell[1]][cell[2]][j];
           distanceIJ(i,id,rij);
           rij[3] = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
           if(rij[3] > rrec) continue; 
           neighbor4[i][numneigh4[i]] = id; 
           numneigh4[i] ++;
        } 
     }
     }
     }
  }

  memory->destroy(ibox); 
  fprintf(screen,"End setup recombination list \n");
}  

/* ----------------------------------------------------------------------
   define virtual function site_energy(int)
------------------------------------------------------------------------- */

double AppBccOcta::site_energy(int)
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   compute energy of site i
------------------------------------------------------------------------- */

double AppBccOcta::sites_energy(int i, int estyle)
{
  int j,jd,n1nn;
  double eng = 0.0;
  double ci, cj;

  //energy from 1NN bonds
  n1nn = numneigh[i];  //num of 1NN
  ci = cj = 0.0;
  for (j = 0; j < n1nn; j++) {
    jd = neighbor[i][j];
    if(cflag && element[i] == element[jd] && element[i] == VACANCY) {
      if(ci == 0.0) ci = site_concentration(i,1); 
      cj = site_concentration(jd,1);

      eng += ci*cj*ebond1[element[i]][element[jd]]; 
    } else 
    eng += ebond1[element[i]][element[jd]];
  }

  //energy from 2NN bonds
  if(estyle == 2) {
    ci = cj = 0.0;
    int n2nn = numneigh2[i];
    for (j = 0; j < n2nn; j++) {
      jd = neighbor2[i][j];
      if(cflag && element[i] == element[jd] && element[i] == VACANCY) {
        if(ci == 0.0) ci = site_concentration(i,2);
        cj = site_concentration(jd,2);

        eng += ci*cj*ebond2[element[i]][element[jd]];
      } else
      eng += ebond2[element[i]][element[jd]];
    }
  }

  return eng/2.0;
}

/* ----------------------------------------------------------------------
  compute local concentration to adjust V-V bond energy  
------------------------------------------------------------------------- */

double AppBccOcta::site_concentration(int i, int estyle)
{ 
  double ci = 0.0;
 
  if(estyle == 1) {
    int n1nn = numneigh[i];  //num of 1NN
    for(int j = 0; j < n1nn; j++) {
       int jd = neighbor[i][j];
       if(element[jd]!=VACANCY) ci += 1.0/(n1nn-1); 
    }
    return ci;
  }

  if(estyle == 2) {
    int n1nn = numneigh2[i];  //num of 1NN
    for(int j = 0; j < n1nn; j++) {
       int jd = neighbor2[i][j];
       if(element[jd]!=VACANCY) ci += 1.0/(n1nn-1);
    }
    return ci;
  }

  return 0.0;
}
/* ----------------------------------------------------------------------
  compute barriers for an exchange event between i & j
------------------------------------------------------------------------- */

double AppBccOcta::site_SP_energy(int i, int j, int estyle)
{
  double eng = 0.0;
  double ci,cj,eng0i, eng0j, eng1i, eng1j; //energy before and after jump
  int m, k, jd;

  int iele = element[i];
  int jele = element[j];

  // no excahange of same type of elements   
  if(iele == jele) return -1.0;
 
  //SIA diffusion on bcc lattice only  
  if(iele >= I1) {
    if(type[j] == OCTA) return -1.0; // SIA does not move to octahedral sublattice 
    if(jele >= I1) return -1.0; // No SIA-SIA exchange 
    double ran1 = ranbccocta->uniform();
    if(dmigration < 3 && ran1 < dratio) { // ratio of 1D/3D diffusion
      double dij[4];
      dij[0] = dij[1] = dij[2] = 0.0;
      distanceIJ(i,j,dij);
      dij[3] = dij[1]*dij[1] + dij[2]*dij[2] + dij[0] * dij[0];
      double prdt = dij[1]*Vmig[iele][1] + dij[2]*Vmig[iele][2] + dij[0]*Vmig[iele][0];
      if (dmigration == 1 && prdt*prdt/dij[3]/Vmig[iele][3] < 0.99) return -1.0; // 1D diffusion
      if (dmigration == 2 && prdt*prdt/dij[3]/Vmig[iele][3] > 0.01) return -1.0; // 2D diffusion
    }
    return mbarrier[iele];
  } 

  //vacancy and SB diffusion on bcc lattice only 
  if(iele == VACANCY || iele == SB) { 
    if(type[j] == OCTA) return -1.0; // no SB dissociation currently (except when a recombination occurs), to be improved later  
    if(iele == SB and jele == VACANCY) return -1.0; // avoid double counting of SB and vacancy exchange    

    // Chenge in bonding energy due to vavancy move  
    eng0i = sites_energy(i,estyle); //broken bond with i initially,
    eng0j = sites_energy(j,estyle); //broken bond with j initially

    // switch the element and recalculate the site energy 
    element[i] = jele;
    element[j] = iele; 
    eng1i = sites_energy(i,estyle); //broken bond with i initially,
    eng1j = sites_energy(j,estyle); //broken bond with j initially

    // switch back 
    element[j] = jele; 
    element[i] = iele; 

    eng = mbarrier[iele] + (eng1i + eng1j - eng0i -eng0j);
    if(eng < 0.0) eng = 0.0; // zero barrier event
    return eng;
  } 

  //He diffuses on both lattices 
  if(iele == HE) {
    // switch with an VO 
    eng0i = sites_energy(i,estyle); //broken bond with i initially,
    eng0j = sites_energy(j,estyle); //broken bond with j initially

    if(jele == VO) {// He move on the octahedral sublattice  
      element[i] = jele;
      element[j] = iele;
    } else if (jele == VACANCY) {// He going to a vacancy creates a VO and an SB  
      element[i] = VO;
      element[j] = SB; 
    } else {// No other He moves allowed 
      return -1.0;
    }
 
    eng1i = sites_energy(i,estyle); //broken bond with i initially,
    eng1j = sites_energy(j,estyle); //broken bond with j initially

    // switch back 
    element[j] = jele; 
    element[i] = iele; 

    eng = mbarrier[iele] + (eng1i + eng1j - eng0i -eng0j);
    if(eng < 0.0) eng = 0.0; // zero barrier event
    return eng;
  } 

  return -1.0; 
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppBccOcta::site_propensity(int i)
{
  int j, iid, jid;

  // valid hop and recombination tabulated lists
  // propensity for each event is input by user
  clear_events(i);
  double prob_reaction = 0.0;
  double prob_hop = 0.0;
  double ebarrier = 0.0;
  double hpropensity = 0.0;

  // propensity for reactions, only when flagged and enabled
  // zero barrier event are dealted with separately in site_events()
  // barriers are from input
  if(reaction_flag) { //reaction flag
    for(j = 0; j < nreaction; j++) {
      if(renable[j] == 0) continue;
      if(element[i] == rinput[j] && type[i] == rsite[j]) {
        iid = rinput[j];
        jid = routput[j];

        if(!(sink_flag && isink[i][jid -1] == 1)) {//production at sinks allowed
          ebarrier = rbarrier[j];
          hpropensity = rrate[j] * exp(-ebarrier/KBT);
          add_event(i,jid,2,j,hpropensity);
          prob_reaction += hpropensity;
        }
      }
    }
  }

  // for hop events, vacancy and SIA only currently
  // propensity calculated in site_SP_energy();
  if (element[i] != VACANCY && element[i] != SB && element[i] != HE  && element[i] < I1) return prob_reaction;

  // for hop events, vacancy and SIA only 
  // propensity calculated in site_SP_energy();
  for (j = 0; j < numneigh[i]; j++) {
      jid = neighbor[i][j];
      ebarrier = site_SP_energy(i,jid,engstyle); // diffusion barrier
      if(ebarrier < 0) continue;
      hpropensity = exp(-ebarrier/KBT);
      add_event(i,jid,1,-1,hpropensity);
      prob_hop += hpropensity;
  }

  if (rhop < 2) return prob_hop + prob_reaction;

  // 2NN hop for <100> 1D SIA diffusion in bcc 
  for (j = 0; j < numneigh2[i]; j++) {
      jid = neighbor2[i][j];
      ebarrier = site_SP_energy(i,jid,engstyle); // diffusion barrier
      if(ebarrier < 0) continue;
      hpropensity = exp(-ebarrier/KBT);
      add_event(i,jid,1,-1,hpropensity);
      prob_hop += hpropensity;
  }

  if (rhop < 3) return prob_hop + prob_reaction;

  // 3NN hop for <110> 1D SIA diffusion in bcc, not used here, to be updated later  
  for (j = 0; j < numneigh3[i]; j++) {
      jid = neighbor3[i][j];
      ebarrier = site_SP_energy(i,jid,engstyle); // diffusion barrier
      if(ebarrier < 0) continue;
      hpropensity = exp(-ebarrier/KBT);
      add_event(i,jid,1,-1,hpropensity);
      prob_hop += hpropensity;
  }

  return prob_hop + prob_reaction;
}

/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppBccOcta::site_event(int i, class RandomPark *random)
{
  int j,jd,k,l,m,n,ii;

  // perform events with non_zero barriers
  // pick one event from total propensity by accumulating its probability
  // compare prob to threshhold, break when reach it to select event
  double threshhold = random->uniform() * propensity[i2site[i]];
  double proball = 0.0;

  // check if trapped by a solute
  // find the event to perform
  int ievent = firstevent[i];
  while (1) {
    proball += events[ievent].propensity;
    if (proball >= threshhold) break;
    ievent = events[ievent].next;
  }

  // perform hop or reaction event
  int rstyle = events[ievent].style;
  int which = events[ievent].which; // type of reactions or neighbor id for acceleration
  j = events[ievent].jpartner;

  // switch element between site i and jpartner for hop diffusion
  if(rstyle == 1) { 
    k = element[i];
    jd = element[j];

    if(k == HE && jd == VACANCY) {// He becomes substitution 
      hcount[i] = 0; // reset due to the reaction of He + V = SB 
      hcount[j] = 0;
      element[i] = VO;
      element[j] = SB; 
      nsites_local[VO-1] ++;  
      nsites_local[SB-1] ++; 
 
      nsites_local[HE-1] --;  
      nsites_local[VACANCY-1] --;  
    } else {
      hcount[i] ++; // will be switch for i & j later
      element[i] = element[j];
      element[j] = k;
    }

    // update He occupancy 
    if(element[i] == SB) {
      iHe[i] = 1;
      for(int nn = 0; nn < numneigh[i]; nn ++) {
         if(element[neighbor[i][nn]] == HE) iHe[i] ++;
      }
    } else {iHe[i] = 0;}

    // update He occupancy 
    if(element[j] == SB) {
      iHe[j] = 1;
      for(int nn = 0; nn < numneigh[j]; nn ++) {
         if(element[neighbor[j][nn]] == HE) iHe[j] ++;
      }
    } else {iHe[j] = 0;} 
 
    /*// switch global atomic id for recombination vector analysis, test only 
    k = aid[i];
    aid[i] = aid[j];
    aid[j] = aid[i];

    double xi = xyz0[i][0];
    double yi = xyz0[i][1];
    double zi = xyz0[i][2];

    xyz0[i][0] = xyz0[j][0];
    xyz0[i][1] = xyz0[j][1];
    xyz0[i][2] = xyz0[j][2];

    xyz0[j][0] = xi;
    xyz0[j][1] = yi;
    xyz0[j][2] = zi;
    // end test */

    if(mfpflag && mfp[element[j]] > 0.0) { // This mean-field absorption works only for vacancy and SIA 
      k = hcount[i]; 
      hcount[i] = hcount[j]; 
      hcount[j] = k; 

      //double r1 = (hcount[j] - mfp[element[j]]);
      //double fx = exp(-r1*r1/sigmamfp)/varmfp; 
      //if(ranbccocta->uniform() < fx) {
      if(hcount[j] > mfp[element[j]]) { 
        nhmfp[element[j]] ++; 
        rhmfp[element[j]] += hcount[j]; 
 
        nsites_local[element[j]-1] --;  
        nsites_local[FE-1] ++;  
        element[j] = FE; //absorption 
        hcount[j] = 0; 
      }
    } 

  } else { // reaction, element[i]<-j; rstyle 2  
    k = element[i];
    element[i] = j;
    rcount[which] ++;
    nsites_local[k-1] --;
    nsites_local[j-1] ++;

    if(mfpflag) hcount[i] = 0; 
    // update reaction target number
    for(ii = 0; ii < nreaction; ii++) {
      if(routput[ii] != j) continue;
      target_local[ii] --;
    }
  }

  // perform zero_barrier events: absorption and recombination
  int rid = recombine(i); // recombine site i with its neighbor rid 
  if(rid >= 0) update_propensity(rid);
 
  // calculate the vector between the two recombined site before diffusing into recombination radius 
  // if(rid >= 0 && rstyle == 1) vec_recombine(i,rid); // i was at j site before diffusion 
 
  if(rstyle == 1) {
    //recombine j with its neighbors too after hopping; 
    int rid = recombine(j); // recombine site i with its neighbor rid 
    if(rid >= 0) update_propensity(rid);
  
    // calculate the vector between the two recombined site before diffusing into recombination radius 
    //if(rid >= 0) vec_recombine(j,rid); // j was at i site before diffusion 

    //sink absorption of element[j] after hopping,for hop only since no element produced at its sinks
    if(nsink > 0) {
      absorption(i);
      absorption(j); 
    }
  } 

  // compute propensity changes for participating sites i & j and their neighbors
  update_propensity(i);
  if(rstyle == 1) update_propensity(j);

  // check if any active reactions needs to be disabled
  if(reaction_flag == 1) {
    for(ii = 0; ii < nreaction; ii ++) {
      if(target_local[ii] <= 0 && renable[ii] == 1) {
        renable[ii] = 0;
        reset_propensity();
        return;
      }
    }
  }

  // check site validity in case of errors  
    //if (element[i] == FE && type[i] == OCTA) fprintf(screen,"%d %d %d %d \n", i,j,which,rstyle);
    //if (element[j] == FE && type[j] == OCTA) fprintf(screen,"%d %d %d %d \n", i,j,which,rstyle);
    //if (element[i] != HE && element[i] != VO && type[i] == OCTA) error->all(FLERR,"site i value invalid");
    //if ((element[i] == HE || element[i] == VO) && type[i] == BCC) error->all(FLERR,"site i value invalid");
    //if (element[j] != HE && element[j] != VO && type[j] == OCTA) error->all(FLERR,"site j value invalid");
    //if ((element[j] == HE || element[j] == VO) && type[j] == BCC) error->all(FLERR,"site j value invalid");
}

/* ----------------------------------------------------------------------
   absorption of an element at site i  
------------------------------------------------------------------------- */
void AppBccOcta::absorption(int i)
{
  int k = element[i];
  for (int n = 0; n < nsink; n++) {
      if(isink[i][k-1] == 1) {
        double rand_me = ranbccocta->uniform();
        if(rand_me <= 1.0/sink_mfp[n]) {
          nabsorption[n] ++;
          nsites_local[k-1] --;
          nsites_local[FE-1] ++;
          element[i] = FE;

          if(mfpflag) hcount[i] = 0; 
          // update reaction target number
          if(reaction_flag == 1) {
            for(int ii = 0; ii < nreaction; ii++) {
               if(routput[ii] == k) target_local[ii] ++;
            }
          }
        }
      }
   }
}

/* ----------------------------------------------------------------------
   calculate the vector between two recombined atoms i and j in vector form. 
------------------------------------------------------------------------- */
/*void AppBccOcta::vec_recombine (int i, int j)
{  
   double dij[3], lprd[3], phi, eta; 
   int m, n, periodicity[3];

   periodicity[0] = domain->xperiodic;
   periodicity[1] = domain->yperiodic;
   periodicity[2] = domain->zperiodic;
   lprd[0] = domain->xprd;
   lprd[1] = domain->yprd;
   lprd[2] = domain->zprd;

   for (int k = 0; k < 3; k++) { 
       dij[k] = xyz0[j][k] - xyz0[i][k];
       if (periodicity[k] && dij[k] >= lprd[k]/2.0) dij[k] -= lprd[k];
       if (periodicity[k] && dij[k] <= -lprd[k]/2.0) dij[k] += lprd[k];
   } 
  
   double x = dij[0];
   double y = dij[1];
   double z = dij[2];
   double xy = sqrt(x*x + y*y);
   double xyz = sqrt(xy*xy+z*z); 
 
   if(xyz > 0.0) { // do not count recombination in short distance to reduce initial effect 
     eta = phi = 0.0; 
     double pi = acos(-1.0); 
     if (xy == 0) { 
        phi = 0;
        eta = pi/2; 
        if (z <= 0) eta = -pi/2;
     } else {
        eta = atan(z/xy);
        phi = acos(x/xy);
        if(y < 0) phi += pi; 
     }  

     m = static_cast<int>(phi*180/pi);
     n = static_cast<int>((eta+pi/2.0)*180/pi);
 
     //fprintf(screen, "%d %d %d %d %f %f vec_rec \n",i, j, m, n, x, y); 
     if(m<0 || n<0) error->all(FLERR,"recombination vector calculated wrongly!");    
     reccount[m][n] ++;
   }

}*/

/* ----------------------------------------------------------------------
   calculate the distance between i and j in vector form. 
------------------------------------------------------------------------- */
void AppBccOcta::distanceIJ(int i, int j, double dij[3])
{
   //update and switch displacement
   int periodicity[3];
   double lprd[3];

   periodicity[0] = domain->xperiodic;
   periodicity[1] = domain->yperiodic;
   periodicity[2] = domain->zperiodic;
   lprd[0] = domain->xprd;
   lprd[1] = domain->yprd;
   lprd[2] = domain->zprd;

   for (int k = 0; k < 3; k++) { //update
       dij[k] = xyz[j][k] - xyz[i][k];
       if (periodicity[k] && dij[k] >= lprd[k]/2.0) dij[k] -= lprd[k];
       if (periodicity[k] && dij[k] <= -lprd[k]/2.0) dij[k] += lprd[k];
   } 
}
/* ----------------------------------------------------------------------
   recombine site i with one of its neighbor if needed.
   return site id which just recombined with i to update propensity.  
   reset target numbers of reactions if recombined 
------------------------------------------------------------------------- */
int AppBccOcta::recombine(int i)
{ 
  if(type[i] == OCTA || element[i] == FE) return -1;

  // recombine by radius 
  int jd = -1; // id of the neighbor to be recombined

  if(rrecombine == 4) { 
    for(int n = 0; n < numneigh4[i]; n++) {
       int m = neighbor4[i][n];

       if(type[m] == OCTA || iHe[m] > 2) continue;
       if(element[i] != VACANCY && element[m] != VACANCY && element[i] != SB && element[m] != SB) continue; // None vacancy or substititute 
       if(element[i] < I1 && element[m] < I1) continue; // None SIA 
 
       jd = m; 
       if(jd >= 0) break;
    } 
  } else if(rrecombine == 2) {  
    for(int n = 0; n < numneigh[i]; n++) {
       int m = neighbor[i][n];

       if(type[m] != BCC || iHe[m] > 2) continue;
       //if(element[i] != VACANCY && element[m] != VACANCY) continue; // None vacancy 
       if(element[i] != VACANCY && element[m] != VACANCY && element[i] != SB && element[m] != SB) continue; // None vacancy 
       if(element[i] < I1 && element[m] < I1) continue; // None SIA 
   
       jd = m; 
       if(jd >= 0) break;
    } 

    if(jd < 0) {
      for(int n = 0; n < numneigh2[i]; n++) {
         int m = neighbor2[i][n];
 
         if(type[m] != BCC || iHe[m] > 2) continue;
         if(element[i] != VACANCY && element[m] != VACANCY && element[i] != SB && element[m] != SB) continue; // None vacancy 
         if(element[i] < I1 && element[m] < I1) continue; // None SIA 

         jd = m; 
         if(jd >= 0) break;
      }
    } 
  } else {
    for(int n = 0; n < numneigh[i]; n++) {
       int m = neighbor[i][n];

       if(type[m] != BCC || iHe[m] > 2) continue;
       if(element[i] != VACANCY && element[m] != VACANCY && element[i] != SB && element[m] != SB) continue; // None vacancy 
       if(element[i] < I1 && element[m] < I1) continue; // None SIA 
   
       jd = m; 
       if(jd >= 0) break;
    } 
  }
  
  if(jd < 0) return jd; 
 
  nrecombine[element[i]-1] ++;
  nrecombine[element[jd]-1] ++;
  nsites_local[element[i]-1] --;
  nsites_local[element[jd]-1] --;
  nsites_local[FE-1] += 2;

  int sd  = -1; 
  if(element[i] == SB) sd = i;
  if(element[jd] == SB) sd = jd;

  element[i] = FE;
  element[jd] = FE;

  if(mfpflag) {hcount[i] = 0; hcount[jd-1] = 0;}
  
  // update reaction target number
  if(reaction_flag == 1) {
    for(int k = 0; k < nreaction; k++) {
      if(routput[k] == element[i] || routput[k] == element[jd])  target_local[k] ++;
    }
  }

  if(sd < 0) return jd; 

  // add an He of an SB is recombined 
  for(int n = 0; n < numneigh[sd]; n ++) {
     int sn = neighbor[sd][n];
     if(element[sn] != VO) continue;

     element[sn] = HE;
     break; 
  }

  return jd;

}

/* ----------------------------------------------------------------------
   update propensity for site i and its neighbors after exchange
   ignore update of sites with i2site < 0; updated during each sweep
   use echeck[] to avoid resetting propensity of same site
------------------------------------------------------------------------- */
void AppBccOcta::update_propensity(int i)
{
  int m,n;
  int nsites = 0;
  int isite = i2site[i];

  if (isite < 0) return;

  propensity[isite] = site_propensity(i); // update site i
  esites[nsites++] = isite;
  echeck[isite] = 1;

  for (n = 0; n < numneigh[i]; n++) { // update propensity for 1NN of site i
    m = neighbor[i][n];
    isite = i2site[m];
    if (isite >= 0 && echeck[isite] == 0) {
      propensity[isite] = site_propensity(m);
      esites[nsites++] = isite;
      echeck[isite] = 1;
    }
  }

  if(engstyle == 2 ) {
    for (n = 0; n < numneigh2[i]; n++) { // update propensity for 2NN of site i
      m = neighbor2[i][n];
      isite = i2site[m];
      if(isite >= 0 && echeck[isite] == 0) {
        propensity[isite] = site_propensity(m);
        esites[nsites++] = isite;
        echeck[isite] = 1;
      }
    }
  }
  solve->update(nsites,esites,propensity);

  // clear echeck array
  for (m = 0; m < nsites; m++)
    echeck[esites[m]] = 0;
}

/* ----------------------------------------------------------------------
   check if any reaction is need to be enabled or disabled
------------------------------------------------------------------------- */
void AppBccOcta::check_reaction()
{
  int i,m,n,n_me,n_total,flag_local,nsites_global;
  double nprcs = (double)(domain->nprocs);

  //update the statistics of each elements locally
  for (i = 0; i < nelement; i++) nsites_local[i] = 0;
  for (i = 0; i < nlocal; i++) nsites_local[element[i]-1]++;

  //compute how many atoms need to be generated globally
  for (i = 0; i < nelement; i++) {
    nsites_global = 0;
    MPI_Allreduce(&nsites_local[i],&nsites_global,1,MPI_INT,MPI_SUM,world);

    for(n = 0; n < nreaction; n++){
      if(i+1 == routput[n]) target_global[n] = rtarget[n] - nsites_global;
    }
  }

  //assign the expected atoms evenly to local processors, turn on reaction locally if needed
  flag_local = 0;
  for(n = 0; n < nreaction; n++){
    target_local[n] = 0;
    m = renable[n];
    renable[n] = 0;

    if(target_global[n] > 0) {
      for( i = 0; i < target_global[n]; i++){
        n_total = 0;

        while(n_total < 1) {
          n_me = 0;
          double rand_me = ranbccocta->uniform();
          if(rand_me < 1.0/nprcs) n_me += 1;
          target_local[n] += n_me;
          MPI_Allreduce(&n_me,&n_total,1,MPI_INT,MPI_SUM,world);

          if(n_total > 1) { //only increase by one allowed in each round
            n_total = 0;
            target_local[n] -= n_me;
          }
        }
      }
    }

    if(target_local[n] > 0) renable[n] = 1;
    flag_local += abs(renable[n] - m);
  }

  // reset site propensity if any enabling or disabling occurs
  if(flag_local > 0) reset_propensity();

}

/* ----------------------------------------------------------------------
   if any reaction is enabled or disabled, reset the propensity
   for all local lattice sites
------------------------------------------------------------------------- */
void AppBccOcta::reset_propensity()
{
  for (int i = 0; i < nset; i++) { //need update all sectors
    for (int m = 0; m < set[i].nlocal; m++){
      set[i].propensity[m] = site_propensity(set[i].site2i[m]);
      set[i].solve->update(m,set[i].propensity);
    }
  }
}

/* ----------------------------------------------------------------------
   clear all events out of list for site I
   add cleared events to free list
------------------------------------------------------------------------- */

void AppBccOcta::clear_events(int i)
{
  int next;
  int index = firstevent[i];
  while (index >= 0) {
    next = events[index].next;
    events[index].next = freeevent;
    freeevent = index;
    nevents--;
    index = next;
  }
  firstevent[i] = -1;
}

/* ----------------------------------------------------------------------
  add an event to list for site I
  event = exchange with site J with probability = propensity
------------------------------------------------------------------------- */

void AppBccOcta::add_event(int i, int j, int rstyle, int which, double propensity)
{
  // grow event list and setup free list
  if (nevents == maxevent) {
    maxevent += DELTAEVENT;
    events =
      (Event *) memory->srealloc(events,maxevent*sizeof(Event),"app:events");
    for (int m = nevents; m < maxevent; m++) events[m].next = m+1;
    freeevent = nevents;
  }

  int next = events[freeevent].next;

  //for hop, rstyle = 1, and switch element at sites I and which
  //for reaction, rstyle = 2, and switch element at site I to which
  //events[freeevent].kpartner = kpartner;
  events[freeevent].style = rstyle;
  events[freeevent].jpartner = j;
  events[freeevent].which = which;
  events[freeevent].propensity = propensity;

  events[freeevent].next = firstevent[i];
  firstevent[i] = freeevent;
  freeevent = next;
  nevents++;
}

/* ----------------------------------------------------------------------
   track sia concentration on this procs: (c(t)*t))
------------------------------------------------------------------------- */
void AppBccOcta::sia_concentration(double t)
{
  double ci = 0.0;
  for (int i = 0; i < nlocal; i ++ ) {if(element[i] >= I1) ci ++; }  
  csia += ci*t/nlocal;
  
}

/* ----------------------------------------------------------------------
   check if perform ballistic mixing
------------------------------------------------------------------------- */
void AppBccOcta::check_ballistic(double t)
{
  int nmix = 0;
  for(int i = 0; i < nballistic; i ++) {
     time_new[i] = static_cast<int>(t/bfreq[i]);
     nmix = time_new[i] - time_old[i];

     if(nmix > 100) fprintf(screen,"Too many Frenkle Pairs generated one time, %d \n", nmix);
     while (nmix) {  //perform mixing nmix times
       nmix --;
       ballistic(i);
       if(nmix == 0) time_old[i] = time_new[i];  //update time
    }
  }
  return;
}

/* ----------------------------------------------------------------------
  Create Frenkel pairs randomly. may need to scale the dose rate 
  by the # of processors when  work in parallel   
------------------------------------------------------------------------- */

void AppBccOcta::ballistic(int n)
{
  // creat an vacancy
  if(nsites_local[FE-1] == 0) error->all(FLERR, "No bcc sites available for FP generation!");      
  if(nsites_local[VO-1] == 0) error->all(FLERR, "No octahedron sites available for He generation!");      

  nFPair ++;   
  int findv = 1;   
  while (findv) { 
    int id = static_cast<int> (nlocal*ranbccocta->uniform()); 
    if(id < nlocal && element[id] == FE) { 
      element[id] = VACANCY; 
      nsites_local[VACANCY-1] ++; 
      nsites_local[FE-1] --;     
      
      // recalculate the propensity if defects are generated 
      update_propensity(id);
      findv = 0; 
    }
  } 
  
  //return; // spk_v no interstitial creation for test 
  //create an interstitial  
  int findi = 1;   
  while (findi) { 
    int id = static_cast<int> (nlocal*ranbccocta->uniform()); 
    if(id < nlocal && element[id] == FE) { 
      if(nelement == VACANCY) error->all(FLERR, "simulation contains no SIAs"); 
      double dx = 1.0/(nelement - I1 + 1);  // equal probability for each type of SIA 
      double r1 = ranbccocta->uniform(); 

      for (int i = 1; i <= nelement - I1 + 1; i++) {
          if(findi == 0) continue; // SIA already generated  
          if(r1 < i*dx) { 
            element[id] = I1 -1  + i; 
            nsites_local[I1-1+i-1] ++; 
            nsites_local[FE-1] --;  
   
            // recalculate the propensity if defects are generated 
            update_propensity(id);      
            findi = 0;
          } 
      } 
    }
  }

  // gas generation
  if(gas_rate == 0) return; // no gas generation 
  if(ranbccocta->uniform() > gas_rate) return; // gas/V_octa ratio no big than one 
  int findg = 1;   
  while (findg) { 
    int id = static_cast<int> (nlocal*ranbccocta->uniform()); 
    if (id < nlocal && element[id] == VO) { 
       element[id] = HE; 
       nsites_local[HE-1] ++; 
       nsites_local[VO-1] --;  
   
       // recalculate the propensity if defects are generated 
       update_propensity(id);      
       findg = 0;
       return;  
    } 
  }
 
}

/* ----------------------------------------------------------------------
  calculating total energy
------------------------------------------------------------------------- */

double AppBccOcta::total_energy( )
{
  double penergy = 0.0;
  for(int i = 0; i < nlocal; i++) penergy += sites_energy(i,engstyle);

  return penergy;
}

/* ----------------------------------------------------------------------
  map n by n matrix to a vector
------------------------------------------------------------------------- */

int AppBccOcta::ibonde(int a, int b, int c)
{
  return ((a-1)*c + b - a*(a-1)/2);
}

 
/* ----------------------------------------------------------------------
  match element. Return 0 if match, 1 otherwise 
------------------------------------------------------------------------- */ 
int AppBccOcta::match(int i)
{ 
  int j,iele,imatch,jtype;
  iele = element[i];
  imatch = 1;
  for (j=0; j < nctype; j++) imatch *=(iele - ncelem[j]);
  return imatch;
 
}

/* ----------------------------------------------------------------------
  perform cluster analysis to obtain precipitation kinetics, currently 
  works in serial only  
------------------------------------------------------------------------- */

void AppBccOcta::cluster()
{
  int i,j,n1nn,jd;

  csize = new int[nlocal];
  cid = new int[nlocal];
  for (i=0; i<nlocal; i++) cid[i] = -1; 
  for (i=0; i<nlocal; i++) csize[i] = 0; 
  
  // check 1NN atoms 
  for (i=0; i<nlocal; i++) {
      if(cid[i] >=0) continue; // already assigned to a cluster
      if(match(i)) continue; // type not match 
      cid[i] = i;
      csize[i] ++;
      
      n1nn = numneigh[i]; // 1NN linkage
      for (j=0; j<n1nn; j++) {
          jd = neighbor[i][j];
          if(cid[jd] >= 0) continue; // already assigned to a cluster
          if(match(jd)) continue; // type not match 
          cid[jd] = i;
          csize[i] ++;
      } 
  }    

  // iterate to link all clusters until done  
  int ccflag = 1;
  int cmin = 0;
  int n = 0;
  while(ccflag) {
      n++;
      ccflag = 0;
      for (i=0; i<nlocal; i++) {
          if(cid[i] < 0) continue;
          cmin = cid[i];

          n1nn = numneigh[i];
          for (j=0; j<n1nn; j++) {
              jd = neighbor[i][j];
              if(cid[jd] < 0) continue;
              if(cid[jd] < cmin) {
                ccflag = 1; 
                cmin = cid[jd];
              } 
          }

          if(ccflag) {
            csize[cid[i]] --;
            csize[cmin] ++;
            cid[i] = cmin;
            for (j=0; j<n1nn; j++) {
                jd = neighbor[i][j];
                if(cid[jd] < 0 || cid[jd] == cmin) continue;
                csize[cid[jd]] --;
                csize[cmin] ++;
                cid[jd] = cmin;
            } 
          }
      }
  }          
       
  // calculate average size and total numbers
  ncluster = 0;
  rcluster = 0.0; 
  for (i=0; i<nlocal; i++) {
      iarray[2][i] = -1;
      if(cid[i] < 0) continue;
      if(csize[cid[i]] < ncsize) continue;
      iarray[2][i] = cid[i];
      if(cid[i] != i) continue;
      rcluster +=csize[i]; // total atoms in all clusters 
      ncluster ++; // each cluster is represented by its member with smallest atom id 
  }

  if(ncluster > 0) rcluster /= ncluster;
  delete [] cid;  
  delete [] csize;  
}

/* ----------------------------------------------------------------------
  grow memory for sink
------------------------------------------------------------------------- */

void AppBccOcta::grow_sinks()
{
  int n = nsink + 1;
  memory->grow(sink_shape,n,"app/bccocta:sink_shape");
  memory->grow(sink_strength,n,"app/bccocta:sink_strength");
  memory->grow(xsink,n,3,"app/bccocta:xsink");
  memory->grow(sink_type,n,"app/bccocta:sink_type");
  memory->grow(sink_normal,n,"app/bccocta:sink_normal");
  memory->grow(sink_segment,n,"app/bccocta:sink_segmant");
  memory->grow(sink_radius,n,"app/bccocta:sink_radius");
  memory->grow(sink_mfp,n,"app/bccocta:sink_mfp");
  memory->grow(nabsorption,n,"app/bccocta:nabsorption");
}

/* ----------------------------------------------------------------------
  grow memory for reaction
------------------------------------------------------------------------- */

void AppBccOcta::grow_reactions()
{
  int n = nreaction + 1;
  memory->grow(rsite,n,"app/bccocta:rsite");
  memory->grow(rinput,n,"app/bccocta:rinput");
  memory->grow(routput,n,"app/bccocta:routput");
  memory->grow(rcount,n,"app/bccocta:rcount");
  memory->grow(renable,n,"app/bccocta:renable");
  memory->grow(rtarget,n,"app/bccocta:rtarget");
  memory->grow(rbarrier,n,"app/bccocta:rbarrier");
  memory->grow(rrate,n,"app/bccocta:rrate");
  memory->grow(target_local,n,"app/bccocta:target_local");
  memory->grow(target_global,n,"app/bccocta:target_global");
}

/* ----------------------------------------------------------------------
  grow memory for ballistic mixing
------------------------------------------------------------------------- */

void AppBccOcta::grow_ballistic()
{
  int n = nballistic + 1;
  int m = 3;
  memory->grow(bfreq,n,"app/bccocta:bfreq");
  memory->grow(time_old,n,"app/bccocta:time_old");
  memory->grow(time_new,n,"app/bccocta:time_new");
}

/* ----------------------------------------------------------------------
  create sinks
------------------------------------------------------------------------- */

void AppBccOcta::sink_creation(int n)
{
  int i,j,periodicity[3];
  int shape = sink_shape[n];
  int normal = sink_normal[n];
  int ntype = sink_type[n];
  int segment = sink_segment[n];
  int nlattice = nlocal + nghost;
  double dx,dij[3],rik,rjk,lprd[3];
  double radius = sink_radius[n];
  double strength = sink_strength[n]*sink_strength[n];

  // get periodicity and box length
  periodicity[0] = domain->xperiodic;
  periodicity[1] = domain->yperiodic;
  periodicity[2] = domain->zperiodic;

  lprd[0] = domain->xprd;
  lprd[1] = domain->yprd;
  lprd[2] = domain->zprd;

  // shape =1: straight line sink, e.g., dislocations
  if(shape == 1) {
    for(i = 0; i < nlattice; i++){
      for(j = 0; j < 3; j ++) {
        if(j != normal) dij[j] = xyz[i][j]-xsink[n][j];
        else dij[j] = 0.0;

        if(periodicity[j] && dij[j] >= lprd[j]/2.0) dij[j] -= lprd[j];
        if(periodicity[j] && dij[j] <= -lprd[j]/2.0) dij[j] += lprd[j];
      }

      dx = dij[0]*dij[0]+dij[1]*dij[1]+dij[2]*dij[2];
      if(dx < strength) isink[i][ntype-1] = 1;
    }
  }

  // circular or polygon sinks, e.g., dislocation loops
  else if (shape == 2) {

    if( segment == 0) {
      for(i = 0; i < nlattice; i++){
        dx = 0.0;
        for( j = 0; j < 3; j ++) {
          dij[j] = xyz[i][j]-xsink[n][j];
          if(periodicity[j] && dij[j] >= lprd[j]/2.0) dij[j] -= lprd[j];
          if(periodicity[j] && dij[j] <= -lprd[j]/2.0) dij[j] += lprd[j];
          if(j != normal) dx += dij[j]*dij[j];
        }

        rik = sqrt(dx);
        rjk = (rik-radius)*(rik-radius) + dij[normal]*dij[normal];
        if( rjk < strength) isink[i][ntype-1] = 1;
      }
    }

    else {

      if( segment < 4) error->all(FLERR,"wrong number of segment for polygon sinks");
      double pi = acos(-1.0);
      double theta = 2*pi/segment;
      double theta_xyz;

      int b = normal + 1;
      int c = b + 1;
      if( b > 2) b -= 3;
      if( c > 2) c -= 3;

      for(i = 0; i < nlattice; i++){
        dx = 0.0;
        for( j = 0; j < 3; j ++) {
          dij[j] = xyz[i][j]-xsink[n][j];
          if(periodicity[j] && dij[j] >= lprd[j]/2.0) dij[j] -= lprd[j];
          if(periodicity[j] && dij[j] <= -lprd[j]/2.0) dij[j] += lprd[j];
          if(j != normal) dx += dij[j]*dij[j];
        }

        if(xyz[i][c] == 0 && xyz[i][b] >= 0) theta_xyz = 0.0;
        else if(xyz[i][c] == 0 && xyz[i][b] < 0) theta_xyz = pi;
        else if(xyz[i][b] == 0) theta_xyz = pi/2;
        else theta_xyz = atan(xyz[i][c] / xyz[i][b]);
        if ( theta_xyz < 0 ) theta_xyz += pi;

        int d = static_cast<int>(theta_xyz/theta);
        double phi = fabs(theta_xyz - (d+0.5)*theta);
        rik = sqrt(dx);
        rjk = (rik-radius/cos(phi))*(rik-radius/cos(phi)) + dij[normal]*dij[normal];
        if( rjk < strength) isink[i][ntype-1] = 1;
      }
    }
  }

  // planar sinks, e.g., grain boundaries
  else if (shape == 3) {
    for(i = 0; i < nlattice; i++){
      dx = xyz[i][normal]-xsink[n][normal];
      if(periodicity[normal] && dx >= lprd[normal]/2.0) dx -= lprd[normal];
      if(periodicity[normal] && dx <= -lprd[normal]/2.0) dx += lprd[normal];
      if(dx*dx < strength) isink[i][ntype-1] = 1;
    }
  }

  // 3D spherical sinks, e.g., precipitates
  else {
    for(i = 0; i < nlattice; i++){
      for( j = 0; j < 3; j ++) {
        dij[j] = xyz[i][j]-xsink[n][j];
        if(periodicity[j] && dij[j] >= lprd[j]/2.0) dij[j] -= lprd[j];
        if(periodicity[j] && dij[j] <= -lprd[j]/2.0) dij[j] += lprd[j];
      }

      dx = dij[0]*dij[0]+dij[1]*dij[1]+dij[2]*dij[2];
      if(dx < strength) isink[i][ntype-1] = 1;
    }
  }
}

