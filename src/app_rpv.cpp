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
#include "app_rpv.h"
#include "domain.h"
#include "solve.h"
#include "random_park.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS; 

enum{NOOP,BCC,NBCC};                          // all sites are BCC except for SIAs; type 
enum{ZERO,FE,VACANCY,CU,NI,MN,Si,P,C,SIA};       // same as DiagRpv; element  

#define DELTAEVENT 100000
#define MAX2NN 6 // max 2NN of BCC lattice 

/* ---------------------------------------------------------------------- */

AppRpv::AppRpv(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg)
{
  ninteger = 2; // first for lattice type,second for element 
  ndouble = 0;   
  delpropensity = 2;
  delevent = 1;
  allow_kmc = 1;
  allow_rejection = 0;

  engstyle = 1; //1 for 1NN interaction, 2 for 2NN interaction; default 1  
  diffusionflag = 0; //flag for MSD calculations, 1 for yes, 0 for no; default 0  
  if (narg < 1) error->all(FLERR,"Illegal app_style command");
  if (narg >= 2) engstyle = atoi(arg[1]);  
  if (narg >= 3) diffusionflag = atoi(arg[2]);  
  if (narg >= 4) concentrationflag = atoi(arg[3]);  
  if (engstyle == 2) delpropensity += 1;// increase delpropensity for 2NN interaction

  // darray 1-4 for msd if activated, followed by concentrations  
  if (diffusionflag == 1) {ninteger++; ndouble += 4;} 
  // calculate concentration fiels for certain elements 
  if (concentrationflag) {ndouble += concentrationflag + 1;}    
  ndiffusion = diffusionflag*4; 
 
  create_arrays();

  firsttime = 1;
  esites = NULL;
  echeck = NULL;
  events = NULL;
  maxevent = 0;
  firstevent = NULL;

  // set random number generator 
  ranrpv = new RandomPark(ranmaster->uniform());
  double seed = ranmaster->uniform();
  ranrpv->reset(seed,me,100);

  // flags for bond interactions 
  ebond1 = NULL;
  ebond2 = NULL;
  mbarrier = NULL;
  hcount = NULL; //numner of vacancy switch events 
  nn1flag = nn2flag = barrierflag = time_flag = 0; //flags for bond energy and migration barriers
  
  // flags and parameters for sinks, dislocations, reactions and ballistic mixing    
  sink_flag = elastic_flag = moduli_flag = dislocation_flag = reaction_flag = 0; //flags for sink dislocation and vacancy
  nsink = ndislocation = nreaction = nballistic = 0; 

  // arrays for dislocations 
  dislocation_type = line_vector = nsegment = NULL;
  burgers = xdislocation = NULL;
  dislocation_radius = NULL; 
  
  // arrays for sinks 
  sink_shape = sink_type = sink_normal = sink_segment = nabsorption = NULL;
  sink_strength = sink_radius = sink_mfp = NULL;
  xsink = NULL;   

  // arrays for reactions 
  rsite = rinput = routput = rcount = renable = rtarget = NULL;
  nsites_local = target_local = target_global =NULL; //number of reactions 
  rbarrier = rrate = NULL; 

  // arrays for ballistic mixing 
  rdamp = pn_local = pn_global = NULL;
  bfreq = time_old = time_new = NULL;
  xmix = pmix = NULL; 
 
  // 2NN neigbor information 
  numneigh2 = NULL; 
  neighbor2 = NULL;

  // number of recombinations
  nrecombine = 0;  
}

/* ---------------------------------------------------------------------- */

AppRpv::~AppRpv()
{
  delete [] esites;
  delete [] echeck;
  delete [] hcount; // number of vacancy switch   
  delete ranrpv; 

  memory->sfree(events);
  memory->destroy(firstevent);
  memory->destroy(ebond1);
  memory->destroy(ebond2);
  memory->destroy(mbarrier);

  if (engstyle == 2) {// memory use for 2NNs
    memory->destroy(numneigh2);
    memory->destroy(neighbor2);
  } 

  if (dislocation_flag) {// memory use related to dislocation 
    memory->destroy(stress); 
    memory->destroy(dislocation_type); 
    memory->destroy(burgers); 
    memory->destroy(dislocation_radius); 
    memory->destroy(line_vector); 
    memory->destroy(nsegment); 
    memory->destroy(xdislocation);
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
    memory->destroy(rdamp);
    memory->destroy(pn_local);
    memory->destroy(pn_global);
    memory->destroy(xmix);
    memory->destroy(pmix);
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
    memory->destroy(nsites_local);
  } 
}

/* ---------------------------------------------------------------------- 
  setup bond energies for 1NN, 2NN, and SPs 
------------------------------------------------------------------------- */ 

void AppRpv::input_app(char *command, int narg, char **arg)
{
  int i,j,ibond;
  int nlattice = nlocal + nghost; 

  // 1NN bond energy taken in the order of 11 12 ... 1N; 22 ... 2N; ...; NN 
  if (strcmp(command,"ebond1") == 0) {

    if (narg < 3) error->all(FLERR,"Illegal ebond1 command");
    nelement = atoi(arg[0]);   // num of elements

    memory->create(nsites_local,nelement,"app/rpv:nsites_local"); 
    memory->create(ebond1,nelement+1,nelement+1,"app/rpv:ebond1");
   
    hcount = new int [nelement]; // total numner of switching with a vacancy;
 
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
    memory->create(ebond2,nelement+1,nelement+1,"app/rpv:ebond2");
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
    memory->create(mbarrier,nelement+1,"app/rpv:mbarrier");
   
    for (i=0; i<narg-1; i++) {
      if(i % 2 == 0){ j = atoi(arg[i]);
        mbarrier[j] = atof(arg[i+1]);  
      } 
    }
  }
 
  // time intervals used to estimate solute trapping   
  else if (strcmp(command, "time_tracer") ==0) {
  
    if (narg < 1 ) error->all(FLERR,"Illegal time_tracer command");
    time_flag = 1;
    dt_interval = atof(arg[0]); 
  }
  
  // elastic moduli for stress calculation  
  else if (strcmp(command, "moduli") ==0) {

    moduli_flag = 1; 
    if(narg != 5) error->all(FLERR,"illegal moduli command"); 
    c11 = atof(arg[0]);
    c12 = atof(arg[1]);
    c44 = atof(arg[2]);
    ninteg = atoi(arg[3]);
    dcore = atof(arg[4]);
  }
  
  // calculating stress field of dislocations or loops  
  else if (strcmp(command, "dislocation") ==0) {

    if(narg < 8 || moduli_flag == 0) error->all(FLERR,"illegal dislocation command"); 
    if(dislocation_flag == 0) {
      memory->create(stress,nlattice,6,"app/rpv:stress"); 
      for(i = 0; i < nlattice; i++) {
        for(j = 0; j < 6; j++) stress[i][j] = 0; 
      }
    }

    dislocation_flag = 1; 
    grow_dislocations(); 

    // dislocation type, 1 straight, others loop  
    dislocation_type[ndislocation] = atoi(arg[0]);
    burgers[ndislocation][0] = atof(arg[1]); //burgers vector 
    burgers[ndislocation][1] = atof(arg[2]); 
    burgers[ndislocation][2] = atof(arg[3]);
    xdislocation[ndislocation][0] = atof(arg[4]); //center of location or loop 
    xdislocation[ndislocation][1] = atof(arg[5]); 
    xdislocation[ndislocation][2] = atof(arg[6]); 
    line_vector[ndislocation] = atoi(arg[7]); //line vector or plane normal for loop, X(0) or Y(1) or Z(2) only 

    // extra parameters for loops 
    if(narg > 8) {
      dislocation_radius[ndislocation] = atof(arg[8]); //loop radius, for loop only  
      nsegment[ndislocation] = atoi(arg[9]); //loop segments(shape), for loop only 
    } 

    // compute the stress field 
    stress_field(ndislocation);
    ndislocation ++; 
  }  

  // elastic interaction with stress field  
  else if (strcmp(command, "elastic_interaction") ==0) {

    elastic_flag = 1; 
    if(narg < 2 || dislocation_flag == 0) error->all(FLERR,"illegal elastic interaction command"); 
    
    int iarg = 0;
    int itype = 0; 

    for (i = 1; i <= nelement; i++ ) evol[i] = 0.0; 
 
    while(iarg < narg) {
      itype = atoi(arg[iarg]);
      evol[itype] = atof(arg[iarg+1]); 
      iarg += 2;  
    } 
  }

  // define sinks to defects, could be dislocations or interfaces or 3D regions 
  else if (strcmp(command, "sink") ==0) {

    if(narg != 10) error->all(FLERR,"illegal sink command");
    if(sink_flag == 0) {
      memory->create(isink, nlattice, nelement,"app/rpv:isink");   
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
    if(diffusionflag) error->warning(FLERR,"MSD calculated with reactions");  
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
    
    bfreq[nballistic] = atoi(arg[0]); // dose rate 
    rdamp[nballistic] = atof(arg[1]); // damage range 

    nballistic ++; // number of mixing events 

  }

  else error->all(FLERR,"Unrecognized command");
}  

/* ----------------------------------------------------------------------
   define 2NN list based on 1NN list. In BCC lattice, each 2NN of site i is shared
   by 4 1NNs of site i as their 1NNs. 
------------------------------------------------------------------------- */

void AppRpv::define_2NN()
{ int candidate[64],frequency[64];
  int i,j,k,jd,kd,n1nn,n2nn,njnn,ncandidate;

  memory->create(numneigh2,nlocal+nghost,"app, numneigh2");
  memory->create(neighbor2,nlocal+nghost,MAX2NN,"app, neighbor2");
  for (i = 0; i < nlocal+nghost; i++) {
    for (j = 0; j < 64; j++) candidate[j] = 0;  
    for (j = 0; j < 64; j++) frequency[j] = 0;  
  
    ncandidate = 0;
    n1nn = numneigh[i];
    n2nn = 0;
    for (j =0; j < n1nn; j++) 
    { 
      jd = neighbor[i][j]; 
      njnn = numneigh[jd];
      for (k = 0; k < njnn; k++) 
      { 
        kd = neighbor[jd][k];
        
        //skip i itself
        if (kd != i) { 
          candidate[ncandidate] = kd; 
          ncandidate++;
        }
      }
    }

    for (j = 0; j < ncandidate; j++) {
      jd = candidate[j];
      frequency[j]++;  
      for (k = j+1; k < ncandidate; k++) {
        kd = candidate[k];
        if(kd == jd) frequency[j]++;
      }  
      if(frequency[j] == 4) {
        if (n2nn == MAX2NN) error->all(FLERR, "Two many 2nd nearest neighbors defined, please expand the simulation cell dimensions or MAX2NN");  
        neighbor2[i][n2nn] = jd;
        n2nn++; // selected if shared by 4 NNs
      }
    }
    numneigh2[i] = n2nn;
  }
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppRpv::grow_app()
{
  type = iarray[0];   // lattice type; i1 in input 
  element = iarray[1];  // element type; i2 in input
  
  if(diffusionflag) {
    aid = iarray[2]; // initially set as global ID, must use set i3 unique in command line 
    disp = darray; // msd; zero initially 
  }
}

/* ----------------------------------------------------------------------
   define virtual function site_energy(int) 
------------------------------------------------------------------------- */

double AppRpv::site_energy(int)
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppRpv::init_app()
{ 
  int i,j;
 
  // define second nearest neighbor 
  if (engstyle ==2) define_2NN(); 
  
  if (firsttime) {
    firsttime = 0;

    echeck = new int[nlocal];
    memory->create(firstevent,nlocal,"app:firstevent");

    // esites must be large enough for 3 sites and their 1st neighbors
    esites = new int[3 + 3*maxneigh];
  }

  // site validity
  int flag = 0;
  for ( i = 0; i < nelement; i++) nsites_local[i] = 0;
 
  for ( i = 0; i < nlocal; i++) {
    if (type[i] < BCC || type[i] > NBCC) flag = 1;
    if (element[i] < FE || element[i] > SIA) flag = 1;
    nsites_local[element[i]-1]++;
  }

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");
 
  // check if reactions need to be enabled or disabled 
  if(reaction_flag) {
    for( i = 0; i< nreaction; i++) {
      rcount[i] = 0; 
      renable[i] = 0;
      target_global[i] = 0; 
      target_local[i] = 0; 
    }
  }

  // initialize the time_list for ballistic mixing 
  if(ballistic_flag) {
    for(i = 0; i < nballistic; i ++) {
       time_old[i] = 0;
       time_new[i] = 0; 
    }
  }

  // initialize the time_list for ballistic mixing 
  if(diffusionflag) {
    for(i = 0; i < ndiffusion; i ++) {
      for(j = 0; j < nlocal; j++){
         disp[i][j] = 0.0; 
      }
    }
  }

  // initialize the time_list for ballistic mixing 
  if(concentrationflag) {
    for(i = ndiffusion; i < ndiffusion + concentrationflag; i ++) {
      for(j = 0; j < nlocal; j++){
         disp[i][j] = 0.0; 
      }
    }
  }

  // initiate parameters for vacancy trapping 
  if(time_flag) {
    double nprcs = (double)(domain->nprocs);  
    fvt = 1.0/nprcs; 
    itrap = 0; 
    itime_current = itime_old = 0;
    treal_me = takmc_me = 0.0;
    dt_real = dt_akmc = 0.0; 
  }
}

/* ---------------------------------------------------------------------- */
/* ---------------------------initiation of looing-up list -------------- */

void AppRpv::setup_app()
{
  for (int i = 0; i < nlocal; i++) echeck[i] = 0;

  // clear event list
  nevents = 0;
  for (int i = 0; i < nlocal; i++) firstevent[i] = -1;
  for (int i = 0; i < maxevent; i++) events[i].next = i+1;
  freeevent = 0;

  // set propensities from rates
  if(temperature == 0.0) error->all(FLERR,"Temperature cannot be 0.0 for app rpv");
  if(nn1flag == 0) error->all(FLERR, "First neareast neighbor bond not defined: AppRpv");
  if(nn2flag == 0 && engstyle == 2) error->all(FLERR, "Second nearest neighbor bond not defined: AppRpv");
  if(barrierflag == 0) error->warning(FLERR, "Diffusion barrier not defined: AppRpv");

  double KB = 0.00008617;
  KBT = temperature * KB;  
}

/* ----------------------------------------------------------------------
   compute energy of site i  
------------------------------------------------------------------------- */

double AppRpv::sites_energy(int i, int estyle)
{
  int j,jd,n1nn;
  double eng = 0.0; 

  //energy from 1NN bonds
  n1nn = numneigh[i];  //num of 1NN
  for (j = 0; j < n1nn; j++) {
    jd = neighbor[i][j];
    eng += ebond1[element[i]][element[jd]] / 2.0;      
  }

  //energy from 2NN bonds 
  if (estyle == 2) {     
    int n2nn = numneigh2[i];
    for (j = 0; j < n2nn; j++) {
      jd = neighbor2[i][j];
      eng += ebond2[element[i]][element[jd]] / 2.0;
    }     
  }

  return eng; 
}

/* ----------------------------------------------------------------------
   compute elastic energy   
------------------------------------------------------------------------- */

double AppRpv::elastic_energy(int i, int itype)
{
  double eng = 0.0;
  double pressure = 0.0; 
  
  pressure = -(stress[i][0] + stress[i][1] + stress[i][2])/3.0; 
  eng = pressure * evol[itype];

  return eng; 
} 
/* ----------------------------------------------------------------------
  compute barriers for an exchange event between i & j  
------------------------------------------------------------------------- */

double AppRpv::site_SP_energy(int i, int j, int estyle)
{
  double eng = 0.0; 
  double eng0i, eng0j, eng1i, eng1j; //energy before and after jump 
  int m, jd;

  eng0i = 2 * sites_energy(i,estyle); //broken bond with i initially, 
  eng0j = 2 * sites_energy(j,estyle); //broken bond with j initially 
  eng1i = eng1j = ebond1[element[i]][element[j]]; //bond between i&j, 
  
  //bond formed with j after switch 
  for (m = 0; m < numneigh[i]; m++) 
  {
    jd = neighbor[i][m];
    if (jd != j) 
       eng1i += ebond1[element[j]][element[jd]]; 
  }  

  if (estyle == 2) 
  {
    for(m = 0; m < numneigh2[i]; m++) 
    {
      jd = neighbor2[i][m]; 
      eng1i += ebond2[element[j]][element[jd]]; 
    }
  }    

  //bond formed with i after switch 
  for (m = 0; m < numneigh[j]; m++) 
  {
    jd = neighbor[j][m];
    if (jd != i) 
       eng1j += ebond1[element[i]][element[jd]];
  }  

  if (estyle == 2) 
  {
    for (m = 0; m < numneigh2[j]; m++) 
    {
      jd = neighbor2[j][m]; 
      eng1j += ebond2[element[i]][element[jd]]; 
    }
  }
    
  //barrier = migbarrier + (eng_after - eng_before)/2.0; 
  eng = mbarrier[element[j]] + (eng1i + eng1j - eng0i -eng0j) / 2.0;

  //add elastic contribution if applicable 
  if(elastic_flag) {
    int itype = element[i];
    int jtype = element[j];
    double eng_ei = elastic_energy(j,itype) - elastic_energy(i,itype);
    double eng_ej = elastic_energy(i,jtype) - elastic_energy(j,jtype);
 
    eng += (eng_ei + eng_ej)/2.0;
  } 

  return eng; 
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppRpv::site_propensity(int i)
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
          if(elastic_flag) ebarrier += elastic_energy(i,jid) - elastic_energy(i,iid); 
          hpropensity = rrate[j] * exp(-ebarrier/KBT);
          add_event(i,jid,2,j,hpropensity);
          prob_reaction += hpropensity;
        }
      } 
    }
  }

  // for hop events, vacancy only currently  
  // propensity calculated in site_SP_energy(); 
  if (element[i] != VACANCY) return prob_reaction;

  for (j = 0; j < numneigh[i]; j++) {
    jid = neighbor[i][j];
    if(element[jid] != VACANCY) { // no vacancy-vacancy switch 
      ebarrier = site_SP_energy(i,jid,engstyle); // diffusion barrier
      hpropensity = exp(-ebarrier/KBT);
      add_event(i,jid,1,-1,hpropensity);
      prob_hop += hpropensity;
    } 
  }
 
  return prob_hop + prob_reaction;
}

/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppRpv::site_event(int i, class RandomPark *random)
{
  int j,k,m,n,ii;
  
  // perform events with non_zero barriers 
  // pick one event from total propensity by accumulating its probability
  // compare prob to threshhold, break when reach it to select event
  double threshhold = random->uniform() * propensity[i2site[i]];
  double proball = 0.0;
 
  // check if trapped by a solute 
  if(time_flag) itrap = vacancy_trap(i);
  // find the event to perform 
  int ievent = firstevent[i];
  while (1) {
    proball += events[ievent].propensity;
    if (proball >= threshhold) break;
    ievent = events[ievent].next;
  }

  // perform hop or reaction event
  int rstyle = events[ievent].style;
  int which = events[ievent].which;
  j = events[ievent].jpartner;

  // switch element between site i and jpartner for hop diffusion 
  if(rstyle == 1) {  
    k = element[i];
    element[i] = element[j];
    element[j] = k;
    hcount[element[i]] ++;

    // calculate MSD for each atom if activated 
    if(diffusionflag) {  
      // switch global atomic id 
      k = aid[i];
      aid[i] = aid[j];
      aid[j] = aid[i];

      //update and switch displacement
      int periodicity[3]; 
      double dij[3],lprd[3];
      periodicity[0] = domain->xperiodic;
      periodicity[1] = domain->yperiodic;
      periodicity[2] = domain->zperiodic;
      lprd[0] = domain->xprd;
      lprd[1] = domain->yprd;
      lprd[2] = domain->zprd;

      for (k = 0; k < 3; k++) { //update 
        dij[k] = xyz[j][k] - xyz[i][k];   
        if (periodicity[k] && dij[k] >= lprd[k]/2.0) dij[k] -= lprd[k];
        if (periodicity[k] && dij[k] <= -lprd[k]/2.0) dij[k] += lprd[k];
        disp[k][i] += dij[k]; 
        disp[k][j] -= dij[k]; 
      }
 
      for (k = 0; k < 3; k++) { //switch  
         dij[k] = disp[k][i];   
         disp[k][i] = disp[k][j]; 
         disp[k][j] = dij[k]; 
      }
      disp[3][i] = disp[0][i]*disp[0][i] + disp[1][i]*disp[1][i] + disp[2][i]*disp[2][i];  
      disp[3][j] = disp[0][j]*disp[0][j] + disp[1][j]*disp[1][j] + disp[2][j]*disp[2][j];  
    }

  } else {
  
 // reaction events with non-zero barriers
    k = element[i];
    element[i] = j;
    rcount[which-1] ++;
    nsites_local[k-1] --; 
    nsites_local[j-1] ++;

    // update reaction target number 
    for(ii = 0; ii < nreaction; ii++) {
      if(routput[ii] != j) continue; 
      target_local[ii] --; 
    }
  }

  // perform zero_barrier events: absorption and recombination 
  if(rstyle == 1) {
  
    //sink absorption,for hop only since no element produced at its sinks  
    if(nsink > 0) {
      k = element[j]; 
      for (n = 0; n < nsink; n++) {  
        if(isink[j][k-1] == 1) { 
          double rand_me = ranrpv->uniform();
          if(rand_me <= 1.0/sink_mfp[n]) {
            nabsorption[n] ++;  
            nsites_local[k-1] --; 
            element[j] = FE;
         
            // update reaction target number
            if(reaction_flag == 1) { 
              for(ii = 0; ii < nreaction; ii++) {
                if(routput[ii] != k) continue; 
                target_local[ii] ++;
              } 
            }
          }
        }
      }
    }

  //check recombinations after hopping; recombine with first SIA found   
    for (n = 0; n < numneigh[j]; n++) { 
      m = neighbor[j][n];
      if(element[j] == VACANCY && element[m] == SIA) {
        nrecombine ++; 
        nsites_local[VACANCY-1] --; 
        nsites_local[SIA-1] --; 

        element[j] = FE;
        element[m] = FE;
         
        // update reaction target number 
        if(reaction_flag == 1) { 
          for(ii = 0; ii < nreaction; ii++) {
            if(routput[ii] != VACANCY && routput[ii] != SIA) continue; 
            target_local[ii] ++;
          }
        }
      }
    }

  } else {
  
    //for a vavancy produced (rstyle 2), check if recombine with neighbor SIAs  
    for (n = 0; n < numneigh[i]; n++) { 
      m = neighbor[i][n];
      if(element[i] == VACANCY && element[m] == SIA) {
        nrecombine ++; 
        nsites_local[VACANCY-1] --; 
        nsites_local[SIA-1] --; 

        element[i] = FE;
        element[m] = FE;

        // update reaction target number 
        if(reaction_flag == 1) { 
          for(ii = 0; ii < nreaction; ii++) {
            if(routput[ii] != VACANCY && routput[ii] != SIA) continue; 
            target_local[ii] ++;
          }
        }
      }
    }
  }
  
  // compute propensity changes for participating sites i & j and their neighbors
  update_propensity(i); 
  update_propensity(j); 

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
} 

/* ----------------------------------------------------------------------
   update propensity for site i,j and their neighbors after exchange 
   ignore update of sites with i2site < 0; updated during each sweep  
   use echeck[] to avoid resetting propensity of same site
------------------------------------------------------------------------- */
void AppRpv::update_propensity(int i)
{ 
  int m,n;  
  int nsites = 0;
  int isite = i2site[i];

  if (isite < 0) return;
  
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
void AppRpv::check_reaction()
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
          double rand_me = ranrpv->uniform();
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
void AppRpv::reset_propensity()
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

void AppRpv::clear_events(int i)
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

void AppRpv::add_event(int i, int j, int rstyle, int which, double propensity)
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
   check if perform ballistic mixing  
------------------------------------------------------------------------- */
void AppRpv::check_ballistic(double t)
{
  int nmix = 0;  
  for(int i = 0; i < nballistic; i ++) {  
     time_new[i] = static_cast<int>(t/bfreq[i]); 
     nmix = time_new[i] - time_old[i]; 

     while (nmix) {  //perform mixing nmix times 
       nmix --; 
       ballistic(i); 
       if(nmix == 0) time_old[i] = time_new[i];  //update time 
    } 
  }
} 

/* ----------------------------------------------------------------------
  calculating total energy  
------------------------------------------------------------------------- */

void AppRpv::ballistic(int n)
{
  int i,iwhich,iid,jid,ibtype,jbtype,ibtype_all,jbtype_all;
  int iproc,jproc,found_me,found_all;
  double rand_me,xiid[3]; 

  // find atom j for the exchange 
  iproc = jproc = -1; 
  ibtype = jbtype = ibtype_all = jbtype_all = 0;  
  xmix[n][0] = xmix[n][1] = xmix[n][2] = 0.0;  
  xiid[0] = xiid[1] = xiid[2] = 0.0; 
  
  found_me = found_all = 0; 
  while (found_all == 0) {
    rand_me = ranrpv->uniform();
    if(rand_me < static_cast<double> (1.0/nprocs)) found_me ++; 
    MPI_Allreduce(&found_me, &found_all,1,MPI_INT,MPI_SUM,world); 

    // only one processor can be selected each time 
    if(found_all > 1) {
      found_me = 0; 
      found_all = 0; 
    } 
  } 

  iid = -1; 
  if(found_me) {
    iproc = me; 
    iwhich = -1; 
    while(iwhich < 0 || iwhich >= nlocal) {
      rand_me = ranrpv->uniform();
      iwhich = static_cast<int> (nlocal*rand_me);}
       
    iid = iwhich;
    ibtype = element[iid];  
    xiid[0] = xyz[iid][0];     
    xiid[1] = xyz[iid][1];     
    xiid[2] = xyz[iid][2];    
  }
 
  MPI_Allreduce(&xiid[0],&xmix[n][0],1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&xiid[1],&xmix[n][1],1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&xiid[2],&xmix[n][2],1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&ibtype,&ibtype_all,1,MPI_INT,MPI_SUM,world);

  // calculate the total exchange propensity
  ballistic_probability(n); 
   
  // find proc j and then jid based on the location of iid  
  found_me = found_all = 0; 
  while (found_all == 0) {
    rand_me = ranrpv->uniform();
    if(rand_me < pn_local[n]/pn_global[n]) found_me ++; 
    MPI_Allreduce(&found_me, &found_all,1,MPI_INT,MPI_SUM,world); 

    // only one processor can be selected each time 
    if(found_all > 1) {
      found_me = 0; 
      found_all = 0; 
    } 
  } 

  // find atom j on selected processor 
  jid = -1;
  if(found_me) {
    i = 0; 

    rand_me = ranrpv->uniform()*pn_local[n]; 
    while(jid < 0 && i < nlocal) {
      if(rand_me <= pmix[n][i]) jid = i;
      i ++;  
    }

    if(jid < 0) error->all(FLERR,"could not find atom j for mixing !");  
    
    jproc = me; // proc owns jid 
    jbtype = element[jid];
  } 
 
  MPI_Allreduce(&jbtype,&jbtype_all,1,MPI_INT,MPI_SUM,world);


  // exchange type of iid and jid and update propensity 
  if(me == iproc) {
    element[iid] = jbtype_all; 
    update_propensity(iid);
  } 

  if(me == jproc) {  
    element[jid] = ibtype_all; 
    update_propensity(jid);
  } 

}

/* ----------------------------------------------------------------------
  calculating total exchange propensity   
------------------------------------------------------------------------- */

void AppRpv::ballistic_probability(int n) 
{
  int i,k,periodicity[3]; 
  double dij[4],lprd[3];
 
  periodicity[0] = domain->xperiodic;
  periodicity[1] = domain->yperiodic;
  periodicity[2] = domain->zperiodic;
  lprd[0] = domain->xprd;
  lprd[1] = domain->yprd;
  lprd[2] = domain->zprd;

  pn_local[n] = pn_global[n] = 0.0; 
  for(i = 0; i < nlocal; i ++) {pmix[n][i] = 0.0;} 
  for(i = 0; i < nlocal; i ++) {
 
     for(k = 0; k < 3; k++) {  
        dij[k] = xyz[i][k] - xmix[n][k];   
        if(periodicity[k] && dij[k] >= lprd[k]/2.0) dij[k] -= lprd[k];
        if(periodicity[k] && dij[k] <= -lprd[k]/2.0) dij[k] += lprd[k];
     }

     dij[3] = sqrt(dij[0]*dij[0]+dij[1]*dij[1]+dij[2]*dij[2]); 
     double pb = 0.0; 
     if(dij[3] > 0.0) pb = exp(-fabs(dij[3])/rdamp[n]); 

     pn_local[n] += pb;
     pmix[n][i] = pn_local[n]; 
  }           

  MPI_Allreduce(&pn_local[n],&pn_global[n],1,MPI_DOUBLE,MPI_SUM,world); 
  if(pn_global[n] == 0.0) {
    error->all(FLERR,"total exchange probability returns 0!"); 
  }
}

/* ----------------------------------------------------------------------
  calculating total energy  
------------------------------------------------------------------------- */

double AppRpv::total_energy( )
{
  double penergy = 0.0;
  for(int i = 0; i < nlocal; i++) penergy += sites_energy(i,engstyle);
  if(elastic_flag) {
    for(int j = 0; j < nlocal; j++) {
      int jtype = element[j]; 
      penergy += elastic_energy(j,jtype);
    }
  }
  
  return penergy;
}

/* ----------------------------------------------------------------------
  calculating "average" concentration at each lattice site 
  for 1NN interaction (engstyle==1): Ci(a)=0.5*ni(a)+0.5*sum(nj(a))/sum(nj)
  for 2NN interaction (engstyle==2): Ci(a)=0.25*ni(a)+0.5*sum(nj(a))/sum(nj)+0.25*sum(nk(a))/sum(nk))
  subscript i for corners atoms and j for 1NN and k for 2NN of i 
  The contributions of corner and center atoms (of bcc cell) are equal 
------------------------------------------------------------------------- */

void AppRpv::concentration_field( )
{ 
  int i,j,jd;

  double *ci = new double[nelement];  
  for(i = 0; i < nlocal; i++) { 
    for(j = 0; j < nelement; j++) ci[j] = 0.0; 

    ci[element[i]-1] += 0.25; 
    if(engstyle == 1) ci[element[i]-1] += 0.25; 
   
    for(j = 0; j < numneigh[i]; j ++) {
      jd = neighbor[i][j]; 
      ci[element[jd]-1] += 0.5*1.0/numneigh[i];
    }

    if(engstyle == 2) { 
      for(j = 0; j < numneigh2[i]; j ++) {
        jd = neighbor2[i][j]; 
        ci[element[jd]-1] += 0.25*1.0/numneigh2[i];
      }
    }

    for(j = ndiffusion; j < ndiffusion + concentrationflag; j++) { // count for concentrationflag elements 
      disp[j][i] = ci[j-ndiffusion];
    } 
  } 

  delete [] ci; 
}

/* ----------------------------------------------------------------------
  check if the vacancy is trapped by solute by counting its  
  first and second NNs, return 0 if any of those is non-Fe 
------------------------------------------------------------------------- */

int AppRpv::vacancy_trap(int i)
{
  int j,jd; 

  //energy from 1NN bonds
  for (j = 0; j < numneigh[i]; j++) {
    jd = neighbor[i][j];
    if(element[jd] > 2) return 0;      
  }

  //energy from 2NN bonds 
  if (engstyle == 2) {     
    for (j = 0; j < numneigh2[i]; j++) {
    jd = neighbor2[i][j];
    if(element[jd] > 2) return 0;      
    }     
  }

  return 1; 
}

/* ----------------------------------------------------------------------
  record KMC time and real time  
------------------------------------------------------------------------- */
void AppRpv::time_tracer(double dt)
{
  dt_akmc += dt; 
  dt_real += dt*itrap; // not contribute to real time if trapped by solute 
} 

/* ----------------------------------------------------------------------
  check if the vacancy is trapped by solute, contribute 0 to physical 
  time if so, return zero if itrap == 1; 1 otherwise 
------------------------------------------------------------------------- */

double AppRpv::real_time(double t)
{
  double treal_all,takmc_all;  
  double dtreal_all,dtakmc_all; 
  double nprcs = (double)(domain->nprocs);  

  dtreal_all = dtakmc_all = 0.0; 
  MPI_Allreduce(&dt_real,&dtreal_all,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&dt_akmc,&dtakmc_all,1,MPI_DOUBLE,MPI_SUM,world);
  
  //compute fvt every dt_interval 
  itime_current = static_cast<int> (t/dt_interval); 
  treal_me += dt_real; 
  takmc_me += dt_akmc; 

  if(itime_current > itime_old)  {
    treal_all = takmc_all = 0.0; 
    MPI_Allreduce(&treal_me,&treal_all,1,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&takmc_me,&takmc_all,1,MPI_DOUBLE,MPI_SUM,world);
    if(takmc_all > 0.0) fvt = treal_all/takmc_all/nprcs; 
    itime_old = itime_current;
    treal_me = 0.0; 
    takmc_me = 0.0;
    dt_real = 0.0; 
    dt_akmc = 0.0;
  }

  if(dtakmc_all == 0.0) return 1.0;  
  return dtreal_all/dtakmc_all; //averaged real time over all processors 
}

/* ----------------------------------------------------------------------
  map n by n matrix to a vector 
------------------------------------------------------------------------- */

int AppRpv::ibonde(int a, int b, int c)
{
  return ((a-1)*c + b - a*(a-1)/2);
}

/* ----------------------------------------------------------------------
  grow memory for dislocation 
------------------------------------------------------------------------- */

void AppRpv::grow_dislocations()
{
  int n = ndislocation + 1;
  memory->grow(dislocation_type,n,"app/rpv:dislocation_type");
  memory->grow(burgers,n,3,"app/rpv:burgers");
  memory->grow(xdislocation,n,3,"app/rpv:xdislocation");
  memory->grow(line_vector,n,"app/rpv:line_vector");
  memory->grow(dislocation_radius,n,"app/rpv:dislocation_radius");
  memory->grow(nsegment,n,"app/rpv:nsegment");
}

/* ----------------------------------------------------------------------
  grow memory for sink 
------------------------------------------------------------------------- */

void AppRpv::grow_sinks()
{
  int n = nsink + 1;
  memory->grow(sink_shape,n,"app/rpv:sink_shape");
  memory->grow(sink_strength,n,"app/rpv:sink_strength");
  memory->grow(xsink,n,3,"app/rpv:xsink");
  memory->grow(sink_type,n,"app/rpv:sink_type");
  memory->grow(sink_normal,n,"app/rpv:sink_normal");
  memory->grow(sink_segment,n,"app/rpv:sink_segmant");
  memory->grow(sink_radius,n,"app/rpv:sink_radius");
  memory->grow(sink_mfp,n,"app/rpv:sink_mfp");
  memory->grow(nabsorption,n,"app/rpv:nabsorption");
}

/* ----------------------------------------------------------------------
  grow memory for reaction
------------------------------------------------------------------------- */

void AppRpv::grow_reactions()
{
  int n = nreaction + 1;
  memory->grow(rsite,n,"app/rpv:rsite");
  memory->grow(rinput,n,"app/rpv:rinput");
  memory->grow(routput,n,"app/rpv:routput");
  memory->grow(rcount,n,"app/rpv:rcount");
  memory->grow(renable,n,"app/rpv:renable");
  memory->grow(rtarget,n,"app/rpv:rtarget");
  memory->grow(rbarrier,n,"app/rpv:rbarrier");
  memory->grow(rrate,n,"app/rpv:rrate");
  memory->grow(target_local,n,"app/rpv:target_local");
  memory->grow(target_global,n,"app/rpv:target_global");
}

/* ----------------------------------------------------------------------
  grow memory for ballistic mixing 
------------------------------------------------------------------------- */

void AppRpv::grow_ballistic()
{
  int n = nballistic + 1;
  int m = 3; 
  memory->grow(bfreq,n,"app/rpv:bfreq");
  memory->grow(time_old,n,"app/rpv:time_old");
  memory->grow(time_new,n,"app/rpv:time_new");
  memory->grow(rdamp,n,"app/rpv:rdamp");
  memory->grow(pn_local,n,"app/rpv:pn_local");
  memory->grow(pn_global,n,"app/rpv:pn_global");
  memory->grow(xmix,n,m,"app/rpv:xmix");
  memory->grow(pmix,n,nlocal,"app/rpv:pmix");
  
}

/* ----------------------------------------------------------------------
  create sinks  
------------------------------------------------------------------------- */

void AppRpv::sink_creation(int n)
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

/* ----------------------------------------------------------------------
  calculate dislocation stess field  
------------------------------------------------------------------------- */

void AppRpv::stress_field(int n)
{  
  if(dislocation_type[n] == 1) stress_dislocation(n); //straight dislocation
  else stress_loop(n); //loop 
}

/* ----------------------------------------------------------------------
  calculate stess field for straight dislocations  
------------------------------------------------------------------------- */

void AppRpv::stress_dislocation(int n)
{ int i,j,k,l,ii;
  int nlattice = nlocal + nghost;
  int periodicity[3];

  double pix; 
  double ni[3],m[3],nn[3][3],nm[3][3],nni[3][3];
  double temp1[3][3],temp2[3][3];
  double bmatx[3][3],qmatx[3][3],smatx[3][3];
  double sigma[3][3],stran[3][3],stres[3][3]; 
  double tvect[3],dij[3];
  double lprd[3]; 

  // get periodicity and box length 
  periodicity[0] = domain->xperiodic;
  periodicity[1] = domain->yperiodic;
  periodicity[2] = domain->zperiodic;

  lprd[0] = domain->xprd;
  lprd[1] = domain->yprd;
  lprd[2] = domain->zprd;

  // calculate the elastic tensor
  elastic_tensor();
 
  for (i = 0; i < 3; i ++) {
    if(i == line_vector[n]) tvect[i] = 1.0;
    else tvect[i] = 0.0;   
  } 

  pix = acos(-1.0); 

  vector_normalize(tvect);
  stroh(tvect,qmatx,bmatx,smatx); 
  
  // calculate stress at each lattice site
  for( int natom = 0; natom < nlattice; natom ++) {   
    for( j = 0; j < 3; j ++) { 
      dij[j] = xyz[natom][j] - xdislocation[n][j]; 
      if(periodicity[j] && dij[j] > lprd[j]/2.0) dij[j] -= lprd[j];
      if(periodicity[j] && dij[j] <= -lprd[j]/2.0) dij[j] += lprd[j];
    }
 
    double norm = dij[0]*tvect[0] + dij[1]*tvect[1] + dij[2]*tvect[2]; 
    for(i = 0; i < 3; i ++) m[i] = dij[i] - tvect[i]*norm; 

    double d = sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]); 
    
    if(d == 0.0) {
      for(i = 0; i < 6; i ++) stress[natom][i] = 0.0; 
    } else {
 
      vector_normalize(m); 

      cross_product(tvect,m,ni);  
   
      for(j = 0; j < 3; j ++) { 
        for(k = 0; k < 3; k ++) { 
          nn[j][k] = 0.0;
          nm[j][k] = 0.0;
        }
      }   

      for(i = 0; i < 3; i ++) { 
        for(j = 0; j < 3; j ++) { 
          for(k = 0; k < 3; k ++) { 
            for(l = 0; l < 3; l ++) { 
              nn[i][j] += ni[k]*cijkl[i][j][k][l]*ni[l];
              nm[i][j] += ni[k]*cijkl[i][j][k][l]*m[l];
            }
          }
        }
      }

      matrix_inversion(nn,nni);

      // inteimediate matrices
      for(i = 0; i < 3; i ++) {
        for(j = 0; j < 3; j ++) {
          temp1[i][j] = bmatx[i][j]; 
          for(k = 0; k < 3; k ++) {
            temp1[i][j] += nm[i][k]*smatx[k][j]; 
          }
       }
      }   
 
      for(i = 0; i < 3; i ++) {
        for(j = 0; j < 3; j ++) {
          temp2[i][j] = 0.0; 
          for(k = 0; k < 3; k ++) {
            temp2[i][j] += nni[i][k]*temp1[k][j]; 
          }
        }
      }   

      // form sigma angular 
      for(i = 0; i < 3; i ++) {
        for(j = 0; j < 3; j ++) {
          sigma[i][j] = 0.0; 
          stran[i][j] = 0.0; 
          stres[i][j] = 0.0; 
        } 
      }

      for(i = 0; i < 3; i ++) { 
        for(j = 0; j < 3; j ++) { 
          for(k = 0; k < 3; k ++) { 
            for(l = 0; l < 3; l ++) { 
              for(ii = 0; ii < 3; ii ++) { 
                sigma[i][j] += cijkl[i][j][k][l]*burgers[n][ii]*
                   (-m[l]*smatx[k][ii] + ni[l]*temp2[k][ii])/2/pix; 
              }
            }
          }
        }
      }

      for(i = 0; i < 3; i ++) {
        for(j = 0; j < 3; j ++) {
          // stres[i][j] = sigma[i][j]/d; 
          stres[i][j] = sigma[i][j]*d/(d*d + dcore*dcore); //stress converge at core distance 
        } 
      }

    // calculate strain and stress 
      for(i = 0; i < 3; i ++) {
        for(j = 0; j < 3; j ++) {
          for(ii = 0; ii < 3; ii ++) { 
            // stran[i][j] += burgers[n][ii]*(-m[j]*smatx[i][ii] + ni[j]*temp2[i][ii])/d; 
            stran[i][j] += burgers[n][ii]*(-m[j]*smatx[i][ii] + ni[j]*temp2[i][ii])*d/(d*d + dcore*dcore); 
          }
        } 
      }

      for(i = 0; i < 3; i ++) { 
        for(j = 0; j < 3; j ++) { 
          stres[i][j] = 0.0; 

          for(k = 0; k < 3; k ++) { 
            for(l = 0; l < 3; l ++) { 
              stres[i][j] += cijkl[i][j][k][l]*stran[k][l]; 
            }
          }
        }
      }

      // record stress at lattice point 
      stress[natom][0] = stres[0][0]; 
      stress[natom][1] = stres[1][1]; 
      stress[natom][2] = stres[2][2]; 

      stress[natom][3] = (stres[0][1]+stres[1][0])/2.0; 
      stress[natom][4] = (stres[0][2]+stres[2][0])/2.0; 
      stress[natom][5] = (stres[1][2]+stres[2][1])/2.0; 
    
    } 
  }

}

/* ----------------------------------------------------------------------
  calculate stess field for dislocation loops  
------------------------------------------------------------------------- */

void AppRpv::stress_loop(int n)
{   
  int i,j;
  int normal,nseg,iseg,jseg;
  int nlattice = nlocal + nghost;
  int periodicity[3];

  double pix,theta; 
  double dij[3],xseg[40][3]; //maximum 40 segment to represent a circle 
  double A[3],B[3],P[3];
  double bv[3],rloop;  // burgers vector and loop normal & radius
  double stres[3][3],sstres[3][3]; // stress tensor at lattice site 
  double lprd[3]; 

  // get periodicity and box length 
  periodicity[0] = domain->xperiodic;
  periodicity[1] = domain->yperiodic;
  periodicity[2] = domain->zperiodic;

  lprd[0] = domain->xprd;
  lprd[1] = domain->yprd;
  lprd[2] = domain->zprd;

  elastic_tensor(); 

  for(i = 0; i < 3; i ++) {
    bv[i] = burgers[n][i];
  } 
 
  rloop = dislocation_radius[n]; 
  normal = line_vector[n]; 
  nseg = nsegment[n];
  pix = acos(-1.0);
  theta = 2*pix/nseg;  

  int b = normal + 1;
  int c = b + 1;  
  if( b > 2) b -= 3;                   
  if( c > 2) c -= 3;                   

  // get vertice at intersections of segments 
  for (iseg = 0; iseg < nseg+1; iseg ++) {
    xseg[iseg][b] = rloop*cos(theta*iseg);
    xseg[iseg][c] = rloop*sin(theta*iseg);
    xseg[iseg][normal] = 0.0; 
  }
  
  // calculate stress at each lattice site
  for( int natom = 0; natom < nlattice; natom ++) {   
    for( j = 0; j < 3; j ++) { 
      dij[j] = xyz[natom][j] - xdislocation[n][j]; 
      if(periodicity[j] && dij[j] >= lprd[j]/2.0) dij[j] -= lprd[j];
      if(periodicity[j] && dij[j] <= -lprd[j]/2.0) dij[j] += lprd[j];
    }

    for( i = 0; i < 3; i++) {
      P[i] = dij[i]; 
      for( j = 0; j < 3; j ++) stres[i][j] = 0.0; 
    } 

    for( iseg = 0; iseg < nseg; iseg ++) {
      jseg = iseg+1; 
      if(jseg == nseg) jseg = 0; 

      for(i = 0; i < 3; i ++) {
        A[i] = xseg[iseg][i]; 
        B[i] = xseg[jseg][i]; 
        
        for( j = 0; j < 3; j ++) sstres[i][j] = 0.0; 
      }

      seg_stress(A,B,P,bv,sstres); 

      for(i = 0; i < 3; i ++) {
        for(j = 0; j < 3; j ++) stres[i][j] += sstres[i][j]; 
      }
    }
 
    // record stress at lattice point 
    stress[natom][0] = stres[0][0]; 
    stress[natom][1] = stres[1][1]; 
    stress[natom][2] = stres[2][2]; 

    stress[natom][3] = (stres[0][1]+stres[1][0])/2.0; 
    stress[natom][4] = (stres[0][2]+stres[2][0])/2.0; 
    stress[natom][5] = (stres[1][2]+stres[2][1])/2.0;
  }
 
}
/* ----------------------------------------------------------------------
  Calculate stress field due to segment AB     
------------------------------------------------------------------------- */

void AppRpv::seg_stress( double A[3], double B[3], double P[3], double bv[3], double sstres[3][3])
{ int i,j; 
  
  double xi[3],PA[3],PB[3],x[3],ni[3]; 
  double norm,eps,dotx,beta1,beta2; 
  double tau1[3],tau2[3],m1[3],m2[3];
  double sigma1[3][3],sigma2[3][3],sigma_p1[3][3],sigma_p2[3][3];  
  
  eps = 1.0e-6; 

  for(i = 0; i < 3; i ++) {
    x[i] = 0.0; 
    xi[i] = B[i] - A[i]; 
  }
  vector_normalize(xi);   

  for(i = 0; i < 3; i ++) {
    PA[i] = P[i] - A[i]; 
    PB[i] = P[i] - B[i]; 
  } 

  dotx = 0.0; 
  for(i = 0; i < 3; i ++) dotx += PA[i]*xi[i]; 
  for(i = 0; i < 3; i ++) x[i] += PA[i] - dotx*xi[i]; 

  double normx = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  // return zero if co-linear 
  if(normx < eps) {
    for(i = 0; i < 3; i ++) {
      for(j = 0; j < 3; j ++) sstres[i][j] = 0.0; 
    } 

    return;

    // if normx < dcore, scale normx to dcore 
  } else if (normx < dcore) {

    for(i = 0; i < 3; i ++) {
      PA[i] = dotx*xi[i] + x[i]*dcore/normx;
      PB[i] = PA[i] + A[i] - B[i];
    }
  } 

  // ABP plane normal to ni 
  vector_normalize(x); 
  
  norm = sqrt(PA[0]*PA[0] + PA[1]*PA[1] + PA[2]*PA[2]);
  for(i = 0; i < 3; i ++) tau1[i] = PA[i]/norm; 

  norm = sqrt(PB[0]*PB[0] + PB[1]*PB[1] + PB[2]*PB[2]);
  for(i = 0; i < 3; i ++) tau2[i] = PB[i]/norm; 

  // ni & m1 & m2 vactor 
  cross_product(xi,x,ni); 
  cross_product(ni,tau1,m1);
  cross_product(ni,tau2,m2);

  // angles
  dotx = 0.0; 
  for(i = 0; i < 3; i ++)  dotx += tau1[i]*xi[i]; 
  beta1 = acos(dotx); 
  
  dotx = 0.0; 
  for(i = 0; i < 3; i ++)  dotx += tau2[i]*xi[i]; 
  beta2 = acos(dotx); 
  
  // Augular factors & stress factors 

  sigma_A(tau1, m1, bv, sigma1); 
  sigma_A(tau2, m2, bv, sigma2); 
  sigma_P(tau1, ni, bv, sigma_p1); 
  sigma_P(tau2, ni, bv, sigma_p2); 

  // calculate the stress field due to AB segment 

  for(i = 0; i < 3; i ++) { 
    for(j = 0; j < 3; j ++) { 
      sstres[i][j] = (cos(beta1)*sigma1[i][j] - cos(beta2)*sigma2[i][j])*normx/2
                     /(normx*normx + dcore*dcore) + (sin(beta2)*sigma_p2[i][j] - sin(beta1)* 
                     sigma_p1[i][j])*normx/2/(normx*normx + dcore*dcore); 
    }
  } 
}

/* ----------------------------------------------------------------------
  Calculate sigma_A 
------------------------------------------------------------------------- */

void AppRpv::sigma_A(double t[3], double m[3], double bv[3], double sigma[3][3])
{ 
  int i,j,k,l,ii; 
  double pix;
  double ni[3],nx[3],tx[3]; 
  double bmatx[3][3],qmatx[3][3],smatx[3][3];
  double bpmat[3][3],qpmat[3][3],spmat[3][3];
  double nn[3][3],nm[3][3],nni[3][3];
  double temp[3][3],ma[3][3];

  pix = acos(-1.0); 
  vector_normalize(t); 
  
  cross_product(t,m,ni); 
  vector_normalize(ni);

  // qpmatx,bpmatx,spmatx
  for(i = 0; i < 3; i ++) {
    nx[i] = ni[i];
    tx[i] = t[i];
  }

  stroh_p(tx,nx,qpmat,bpmat,spmat); 

  // qmatx,bmatx,smatx
  for(i = 0; i < 3; i ++) tx[i] = t[i]; 
  stroh(tx,qmatx,bmatx,smatx); 
 
  // nn & nm 
  for(i = 0; i < 3; i ++) {
    for(j = 0; j < 3; j ++) {
      nn[i][j] = 0.0; 
      nm[i][j] = 0.0; 
    }
  } 

  for(j = 0; j < 3; j ++) { 
    for(k = 0; k < 3; k ++) { 
      for(i = 0; i < 3; i ++) { 
        for(l = 0; l < 3; l ++) { 
          nn[j][k] += ni[i]*cijkl[i][j][k][l]*ni[l];
          nm[j][k] += ni[i]*cijkl[i][j][k][l]*m[l];
        }
      }
    }
  }

  matrix_inversion(nn,nni);
  
  //temporaty matrix ma 
  for(i = 0; i < 3; i ++) {
    for(j = 0; j < 3; j ++) {
      temp[i][j] = 0.0 ; 
      for(k = 0; k < 3; k ++)  temp[i][j] += nm[i][k]*smatx[k][j]; 
    }
  }
   
  for(i = 0; i < 3; i ++) {
    for(j = 0; j < 3; j ++)  temp[i][j] += bmatx[i][j];
  }  

  for(i = 0; i < 3; i ++) {
    for(j = 0; j < 3; j ++) {
      ma[i][j] = 0.0 ; 
      for(k = 0; k < 3; k ++)  ma[i][j] += nni[i][k]*temp[k][j]; 
    }
  }

  // calculate the angular stress factors 
  for(i = 0; i < 3; i ++) {
    for(j = 0; j < 3; j ++)  sigma[i][j] = 0.0;
  }  

  for(i = 0; i < 3; i ++) { 
    for(j = 0; j < 3; j ++) { 
      for(k = 0; k < 3; k ++) { 
        for(l = 0; l < 3; l ++) { 
          for(ii = 0; ii < 3; ii ++) { 
            sigma[i][j] += cijkl[i][j][k][l]*bv[ii]*
                 (-m[l]*smatx[k][ii] + ni[l]*ma[k][ii])/2/pix; 
          }
        }
      }
    }
  }

}

/* ----------------------------------------------------------------------
  Calculate sigma_P 
------------------------------------------------------------------------- */

void AppRpv::sigma_P(double t[3], double ni[3], double bv[3], double sigmap[3][3])
{ int i, j, k, l, ii; 
 
  double pix;
  double nn[3][3],nm[3][3],nt[3][3],nni[3][3],mi[3];
  double bmatx[3][3],qmatx[3][3],smatx[3][3];
  double bpmat[3][3],qpmat[3][3],spmat[3][3];
  double temp1[3][3],temp2[3][3],temp3[3][3];
  double ma[3][3]; 

  pix = acos(-1.0); 
  vector_normalize(t); 
  vector_normalize(ni);

  cross_product(ni,t,mi); 
  stroh_p(t,ni,qpmat,bpmat,spmat); 
  stroh(t,qmatx,bmatx,smatx); 
 
  // nn & nm 
  for(i = 0; i < 3; i ++) {
    for(j = 0; j < 3; j ++) {
      nn[i][j] = 0.0; 
      nm[i][j] = 0.0; 
      nt[i][j] = 0.0; 
    }
  } 

  for(i = 0; i < 3; i ++) { 
    for(j = 0; j < 3; j ++) { 
      for(k = 0; k < 3; k ++) { 
        for(l = 0; l < 3; l ++) { 
          nn[j][k] += ni[i]*cijkl[i][j][k][l]*ni[l];
          nm[j][k] += ni[i]*cijkl[i][j][k][l]*mi[l];
          nt[j][k] += ni[i]*cijkl[i][j][k][l]*t[l];
        }
      }
    }
  }

  matrix_inversion(nn,nni);
  
  // temporaty matrix ma 
  for(i = 0; i < 3; i ++) {
    for(j = 0; j < 3; j ++) {
      temp1[i][j] = 0.0 ; 
      for(k = 0; k < 3; k ++)  temp1[i][j] += nt[i][k]*smatx[k][j]; 
    }
  }
   
  for(i = 0; i < 3; i ++) {
    for(j = 0; j < 3; j ++) {
      temp2[i][j] = 0.0 ; 
      for(k = 0; k < 3; k ++)  temp2[i][j] += nm[i][k]*spmat[k][j]; 
    }
  }  

  for(i = 0; i < 3; i ++) {
    for(j = 0; j < 3; j ++)    temp3[i][j] = bpmat[i][j] + temp2[i][j] - temp1[i][j]; 
  } 
 
  for(i = 0; i < 3; i ++) {
    for(j = 0; j < 3; j ++) {
      ma[i][j] = 0.0 ; 
      for(k = 0; k < 3; k ++)  ma[i][j] += nni[i][k]*temp3[k][j]; 
    }
  }

  // calculate the angular stress factors 

  for(i = 0; i < 3; i ++) {
    for(j = 0; j < 3; j ++)    sigmap[i][j] = 0.0;
  }  

  for(i = 0; i < 3; i ++) { 
    for(j = 0; j < 3; j ++) { 
      for(k = 0; k < 3; k ++) { 
        for(l = 0; l < 3; l ++) { 
          for(ii = 0; ii < 3; ii ++) { 
            sigmap[i][j] += cijkl[i][j][k][ii]*bv[l]*
                 (t[ii]*smatx[k][l] - mi[ii]*spmat[k][l] + ni[ii]*ma[k][l])/2/pix; 
          }
        }
      }
    }
  }

}

/* ----------------------------------------------------------------------
  stress calculation based on stroh formulism        
------------------------------------------------------------------------- */

void AppRpv::stroh( double tvect[3], double qmatx[3][3], double bmatx[3][3], double smatx[3][3])
{
  int i,j,k,l; 
  double pix,omega,domega; 
  double mvect[3],nvect[3],mmvect[3],nnvect[3];
  double nn[3][3],mm[3][3],nm[3][3],mn[3][3];
  double nni[3][3],nn2[3][3],nn3[3][3];
  
  pix = acos(-1.0); 
  domega = 2*pix/ninteg; 

  right_hand_coord(mvect,nvect,tvect); 
  
  for(i = 0; i < 3; i ++) {
    for(j = 0; j < 3; j ++) {
      qmatx[i][j] = 0.0;
      smatx[i][j] = 0.0;
      bmatx[i][j] = 0.0;
    } 
  }
  
  for(int integ = 0; integ < ninteg; integ ++) {

    for(i = 0; i < 3; i ++) {
      for(j = 0; j < 3; j ++) {
        nn[i][j] = 0.0;
        nm[i][j] = 0.0;
        mn[i][j] = 0.0;
        mm[i][j] = 0.0;
      } 
    }

    omega = integ*domega; 
   
    for(i = 0; i < 3; i ++) {
      mmvect[i] = mvect[i]*cos(omega)+nvect[i]*sin(omega);
      nnvect[i] = -mvect[i]*sin(omega)+nvect[i]*cos(omega);
    }

    // different with Bulent's version, double check 
    for(i = 0; i < 3; i++) {
      for(j = 0; j < 3; j++) {
        for(k = 0; k < 3; k++) { 
          for(l = 0; l < 3; l++) {
             nn[i][j] += nnvect[k]*cijkl[i][j][k][l]*nnvect[l]; 
             nm[i][j] += nnvect[k]*cijkl[i][j][k][l]*mmvect[l]; 
             mm[i][j] += mmvect[k]*cijkl[i][j][k][l]*mmvect[l]; 
          }
        }
      }
    }


    for(i = 0; i < 3; i++) {
      for(j = 0; j < 3; j++) mn[j][i] = nm[i][j]; 
    } 

    matrix_inversion(nn,nni); 
 
    // get nn2 
    for(i = 0; i < 3; i++) {
      for(j = 0; j < 3; j++) {
        nn2[i][j] = 0.0; 
        for(k = 0; k < 3; k++)  nn2[i][j] += nni[i][k]*nm[k][j]; 
      } 
    }
 
    // get nn3 
    for(i = 0; i < 3; i++) {
      for(j = 0; j < 3; j++) {
        nn3[i][j] = mm[i][j]; 
        for(k = 0; k < 3; k++)  nn3[i][j] -= mn[i][k]*nn2[k][j]; 
      } 
    } 

    // calculate qmatx, smatx and bmatx 
    for(i = 0; i < 3; i++) {
      for(j = 0; j < 3; j++) {
        qmatx[i][j] -= nni[i][j]/ninteg;
        smatx[i][j] -= nn2[i][j]/ninteg;
        bmatx[i][j] += nn3[i][j]/ninteg;
      }
    }
  } // loop for integ = 0 - ninteg 
}

/* ----------------------------------------------------------------------
  angular factors in stress calculation based on stroh formulism        
------------------------------------------------------------------------- */

void AppRpv::stroh_p( double t[3], double n0[3], double qpmat[3][3], double bpmat[3][3], double spmat[3][3])
{ int i,j,k,l,ii;

  double pix,value;
  double N1[3],M1[3],ni[3],mi[3]; 
  double domega,omega,omegas[100]; 
  double qmatx[3][3],bmatx[3][3],smatx[3][3];
  double nn[3][3],mm[3][3],mn[3][3],nm[3][3],nni[3][3]; 
  double nt[3][3],tn[3][3],mt[3][3],tm[3][3]; 
  double F[3][3],Bi[3][3],Si[3][3],Qi[3][3];
  double temp1[3][3],temp2[3][3],temp3[3][3],temp4[3][3]; 
  double BP[9][100],SP[9][100],QP[9][100]; 
   
  pix = acos(-1.0); 
  
  vector_normalize(t);
  vector_normalize(n0);

  for(i = 0; i < 3; i ++) N1[i] = n0[i];
  cross_product(N1,t,M1); 
  stroh(t,qmatx,bmatx,smatx); 

  // angles 
  domega = 2.0*pix/ninteg; 
  for(i = 0; i < ninteg+1; i ++)  omegas[i] = i*domega; 
  
  for(i = 0; i < 9; i ++) {
    for(j = 0; j < ninteg+1; j ++) {
      BP[i][j] = 0.0; 
      SP[i][j] = 0.0; 
      QP[i][j] = 0.0; 
    }
  } 

  for(i = 0; i < ninteg+1; i ++) {
    omega = omegas[i];

    for(j = 0; j < 3; j ++) {
      mi[j] = M1[j]*cos(omega) + N1[j]*sin(omega);
      ni[j] = -M1[j]*sin(omega) + N1[j]*cos(omega);
    }

    for(ii = 0; ii < 3; ii ++) {
      for(j = 0; j < 3; j ++) {
        nn[ii][j] = 0.0; 
        nm[ii][j] = 0.0; 
        mm[ii][j] = 0.0; 
        mn[ii][j] = 0.0; 
        nt[ii][j] = 0.0; 
        tn[ii][j] = 0.0; 
        mt[ii][j] = 0.0; 
        tm[ii][j] = 0.0; 
      }
    } 

    for(ii = 0; ii < 3; ii ++){
      for(j = 0; j < 3; j ++){
        for(k = 0; k < 3; k ++){
          for(l = 0; l < 3; l ++){
            nn[j][k] += ni[ii]*cijkl[ii][j][k][l]*ni[l];  
            nm[j][k] += ni[ii]*cijkl[ii][j][k][l]*mi[l];  
            mm[j][k] += mi[ii]*cijkl[ii][j][k][l]*mi[l];  
            mn[j][k] += mi[ii]*cijkl[ii][j][k][l]*ni[l];  
            nt[j][k] += ni[ii]*cijkl[ii][j][k][l]*t[l];  
            tn[j][k] += t[ii]*cijkl[ii][j][k][l]*ni[l];  
            mt[j][k] += mi[ii]*cijkl[ii][j][k][l]*t[l];  
            tm[j][k] += t[ii]*cijkl[ii][j][k][l]*mi[l];  
          }
        }
      }
    }

    matrix_inversion(nn,nni);

    // define F matrix 
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        temp1[j][k] = nt[j][k] + tn[j][k]; 
      }
    }

    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        temp2[j][k] = 0.0; 
        for(l = 0; l < 3; l ++) temp2[j][k] += temp1[j][l]*nni[l][k]; 
      }
    }

    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        F[j][k] = 0.0; 
        for(l = 0; l < 3; l ++) F[j][k] += nni[j][l]*temp2[l][k]; 
      }
    }

    // define Si matrix 
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        temp1[j][k] = tm[j][k]*sin(omega) - nt[j][k]*cos(omega); 
      }
    }

    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        temp2[j][k] = 0.0; 
        for(l = 0; l < 3; l ++) temp2[j][k] += nni[j][l]*temp1[l][k]; 
      }
    }

    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        temp1[j][k] = 0.0; 
        for(l = 0; l < 3; l ++) temp1[j][k] += F[j][l]*nm[l][k]; 
      }
    }

    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        Si[j][k] = -temp1[j][k]*sin(omega) + temp2[j][k]; 
      }
    }

    // define Qi 
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        Qi[j][k] = -F[j][k]*sin(omega); 
      }
    }

    // define Bi 
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        temp1[j][k] = tm[j][k]*sin(omega) - nt[j][k]*cos(omega); 
      }
    }

    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        temp2[j][k] = 0.0; 
        for(l = 0; l < 3; l ++) temp2[j][k] += nni[j][l]*temp1[l][k]; 
      }
    }

    // keep this temp1
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        temp1[j][k] = 0.0; 
        for(l = 0; l < 3; l ++) temp1[j][k] += mn[j][l]*temp2[l][k]; 
      }
    }

    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        temp2[j][k] = 0.0; 
        for(l = 0; l < 3; l ++) temp2[j][k] += nni[j][l]*nm[l][k]; 
      }
    }

    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        temp3[j][k] = tn[j][k]*cos(omega) - mt[j][k]*sin(omega); 
      }
    }

    // keep this temp4
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        temp4[j][k] = 0.0; 
        for(l = 0; l < 3; l ++) temp4[j][k] += temp3[j][l]*temp2[l][k]; 
      }
    }

    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        temp2[j][k] = 0.0; 
        for(l = 0; l < 3; l ++) temp2[j][k] += F[j][l]*nm[l][k]*sin(omega); 
      }
    }

    // keep this temp3
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        temp3[j][k] = 0.0; 
        for(l = 0; l < 3; l ++) temp3[j][k] += mn[j][l]*temp2[l][k]; 
      }
    }

    // keep this temp2
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        temp2[j][k] = mt[j][k] + tm[j][k]; 
      }
    }

    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        Bi[j][k] = -temp2[j][k]*cos(omega) + temp3[j][k] + temp4[j][k] - temp1[j][k]; 
      }
    }

    // define BP, SP, QP 
    
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        BP[3*j+k][i] = Bi[j][k]; 
        SP[3*j+k][i] = Si[j][k]; 
        QP[3*j+k][i] = Qi[j][k]; 
      }
    }
  } // end loop for i = 0 - ninteg  

  value = 0.0; 
  for(j = 0; j < 3; j ++) {
    for(k = 0; k < 3; k ++) {
      trapozidal(omegas,BP,j,k,value); 
      bpmat[j][k] = value/2/pix;  

      trapozidal(omegas,SP,j,k,value); 
      spmat[j][k] = value/2/pix;  

      trapozidal(omegas,QP,j,k,value); 
      qpmat[j][k] = value/2/pix;  
    }
  }
}

/* ----------------------------------------------------------------------
  establish elastic tensor matrix from c11 c12 & c44    
------------------------------------------------------------------------- */

void AppRpv::elastic_tensor()
{   
  int i,j,k,l; 
  double anisotropy,delta[3][3]; 

  // shift c44 
  anisotropy = fabs(2*c44 + c12 - c11);
  if(anisotropy < 1.0e-6) c44 += 1.0e-6; 

  // identity matrix 
  for(i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
      if(i == j) delta[i][j] = 1; 
      else delta[i][j] = 0; 
    }
  }

  // cijkl tensor 
  anisotropy = 2*c44 + c12 - c11;
  for(i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
      for(k = 0; k < 3; k++) { 
        for(l = 0; l < 3; l++) {
          cijkl[i][j][k][l] = c44*(delta[i][k]*delta[j][l] + delta[i][l]*delta[j][k])
                              + c12*delta[i][j]*delta[k][l]
                              -anisotropy*delta[i][j]*delta[k][l]*delta[i][k];
        }
      }
    }
  }
   
}

/* ----------------------------------------------------------------------
  3x3 matrix inversion     
------------------------------------------------------------------------- */

void AppRpv::matrix_inversion(double mm[3][3], double nn[3][3])
{
   double det = mm[0][0]*(mm[2][2]*mm[1][1] - mm[2][1]*mm[1][2]) 
                - mm[1][0]*(mm[2][2]*mm[0][1] - mm[2][1]*mm[0][2])
                + mm[2][0]*(mm[1][2]*mm[0][1] - mm[1][1]*mm[0][2]);

   if(det <= 0) error->all(FLERR,"matrix not invertible!"); 
  
   // get the inverse matrix 

   nn[0][0] = (mm[2][2]*mm[1][1] - mm[2][1]*mm[1][2])/det; 
   nn[0][1] = -(mm[2][2]*mm[0][1] - mm[2][1]*mm[0][2])/det; 
   nn[0][2] = (mm[1][2]*mm[0][1] - mm[1][1]*mm[0][2])/det; 
   nn[1][0] = -(mm[2][2]*mm[1][0] - mm[2][0]*mm[1][2])/det; 
   nn[1][1] = (mm[2][2]*mm[0][0] - mm[2][0]*mm[0][2])/det; 
   nn[1][2] = -(mm[1][2]*mm[0][0] - mm[1][0]*mm[0][2])/det; 
   nn[2][0] = (mm[2][1]*mm[1][0] - mm[2][0]*mm[1][1])/det; 
   nn[2][1] = -(mm[2][1]*mm[0][0] - mm[2][0]*mm[0][1])/det; 
   nn[2][2] = (mm[1][1]*mm[0][0] - mm[1][0]*mm[0][1])/det; 

}

/* ----------------------------------------------------------------------
  cross product of double vectors     
------------------------------------------------------------------------- */

void AppRpv::cross_product( double m[3], double n[3], double l[3])
{
  l[0] = m[1]*n[2] - m[2]*n[1];
  l[1] = m[2]*n[0] - m[0]*n[2];
  l[2] = m[0]*n[1] - m[1]*n[0];

}

/* ----------------------------------------------------------------------
  cross product of int vectors     
------------------------------------------------------------------------- */

void AppRpv::cross_product( int m[3], int n[3], int l[3])
{
  l[0] = m[1]*n[2] - m[2]*n[1];
  l[1] = m[2]*n[0] - m[0]*n[2];
  l[2] = m[0]*n[1] - m[1]*n[0];

}

/* ----------------------------------------------------------------------
  normalize double vector      
------------------------------------------------------------------------- */

void AppRpv::vector_normalize( double m[3])
{
  double norm = sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]); 
  if(norm == 0) error->all(FLERR,"can not normalize zero vector"); 
  for(int i = 0; i < 3; i ++) m[i] /= norm; 

}

/* ----------------------------------------------------------------------
  right-hand coordination system, m x n x t        
------------------------------------------------------------------------- */

void AppRpv::right_hand_coord( double m[3], double n[3], double t[3])
{
  double taux[3];
  // random vector taux to find m  
  taux[0] = ranrpv->uniform();
  taux[1] = ranrpv->uniform();
  taux[2] = ranrpv->uniform();
  
  vector_normalize(t);
  vector_normalize(taux); 
  
  // m = t x taux
  cross_product(t,taux,m); 
  vector_normalize(m); 

  // n = t x m
  cross_product(t,m,n); 
  vector_normalize(n); 

}

/* ----------------------------------------------------------------------
  trapozidal integration     
------------------------------------------------------------------------- */

void AppRpv::trapozidal(double omegas[100], double XP[9][100], int j, int k, double value)
{ int i, ii;
  double dy,domega;  

  for(i = 0; i < ninteg; i ++) { 
    ii = i + 1; 
    domega = omegas[ii] - omegas[i]; 
    dy = XP[3*j+k][ii] + XP[3*j+k][i]; 
    value += dy*domega/2.0; 
  }
}
 
