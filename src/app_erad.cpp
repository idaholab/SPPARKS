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
#include "app_erad.h"
#include "domain.h"
#include "solve.h"
#include "random_park.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

enum{NOOP,BCC,NBCC};                                // all sites are BCC or FCC; type
enum{ZERO,FE,VACANCY,I1, I2, I3, I4, I5, I6};       // same as DiagErad; element

#define DELTAEVENT 100000
#define MAX2NN 6  // max 2NN of BCC lattice
#define MAX3NN 12 // max 3NN of BCC lattice
#define BIGNUMBER 1e18 // define a big number 

/* ---------------------------------------------------------------------- */

AppErad::AppErad(SPPARKS *spk, int narg, char **arg) :
  AppLattice(spk,narg,arg)
{
  ninteger = 2; // first for lattice type; second for element
  ndouble = 0;
  delpropensity = 2;
  delevent = 1;
  allow_kmc = 1;
  allow_rejection = 0;

  engstyle = 1; //1 for 1NN interaction, 2 for 2NN interaction; default 1
  concentrationflag = 0; //flag for concentration calculation  
  percolationflag = 0; //flag for concentration calculation  
  rhop = 1; // default only 1NN hop 
  rrecombine = 1; // recombination radius, default 1
  cflag = clst_flag = 0; // if activate concentration dependent VV bond, default 0 (no); 
  if (narg < 1) error->all(FLERR,"Illegal app_style command");
  if (narg >= 2) engstyle = atoi(arg[1]);
  if (narg >= 3) rhop = atoi(arg[2]);
  if (narg >= 4) concentrationflag = atoi(arg[3]);
  if (narg >= 5) rrecombine = atoi(arg[4]);
  if (rrecombine == 4) { 
     rrec = atof(arg[5]);
     if(narg >= 7) cflag = atoi(arg[6]);
  } else if (narg >= 6) cflag = atoi(arg[5]);
  //if (concentrationflag) ndouble += concentrationflag;
  if (engstyle == 2) delpropensity += 1;// increase delpropensity for 2NN interaction
  if (rhop == 3) delpropensity +=1; 

  ninteger ++; // for percolation  roc cluster analysis, test only  
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
  ranerad = new RandomPark(ranmaster->uniform());
  double seed = ranmaster->uniform();
  ranerad->reset(seed,me,100);

  // flags for bond interactions
  ebond1 = NULL;
  ebond2 = NULL;
  mbarrier = NULL;
  hcount = NULL; //numner of hop events for mobile elements
  nn1flag = nn2flag = mfpflag = barrierflag = 0; //flags for bond energy and migration barriers

  // flags and parameters for sinks, dislocations, reactions and ballistic mixing
  sink_flag = elastic_flag = moduli_flag = dislocation_flag = reaction_flag = acceleration_flag = 0; //flags for sink dislocation and vacancy
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
  bfreq = NULL; 
  time_old = time_new = NULL;
  min_bfreq = BIGNUMBER;

  // 2NN neighbor information
  ibox = NULL;
  numneigh2 = numneigh3 = numneigh4 = NULL;
  neighbor2 = neighbor3 = neighbor4 = NULL;

  // cluster analysis 
  ncelem = cid = csize = NULL;
}

/* ---------------------------------------------------------------------- */

AppErad::~AppErad()
{
  delete [] esites;
  delete [] echeck;
  delete ranerad;

  memory->sfree(events);
  memory->destroy(firstevent);
  memory->destroy(ebond1);
  memory->destroy(ebond2);
  memory->destroy(mbarrier);
  memory->destroy(nsites_local);
  memory->destroy(hcount);

  if (engstyle == 2) {// memory use for 2NNs
    memory->destroy(numneigh2);
    if(rhop == 3 || rrecombine ==3) memory->destroy(numneigh3);
    memory->destroy(neighbor2);
    if(rhop == 3 || rrecombine ==3) memory->destroy(neighbor3);
  }

  if (rrecombine == 4) {
    memory->destroy(numneigh4);
    memory->destroy(neighbor4);
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

  if (mfpflag) {// memory use for acceleraton
    memory->destroy(rhmfp);
    memory->destroy(nhmfp);
    memory->destroy(mfp);
  }

  if (acceleration_flag) {// memory use for acceleraton
    memory->destroy(trap_type);
  }
}

/* ----------------------------------------------------------------------
  setup bond energies for 1NN, 2NN, and SPs
------------------------------------------------------------------------- */

void AppErad::input_app(char *command, int narg, char **arg)
{
  int i,j,ibond;
  int nlattice = nlocal + nghost;

  // initial coordiation of each atom, for recombination vector analysis, test only 
  //memory->create(xyz0,nlocal,3,"app/erad:xyz0");

  // 1NN bond energy taken in the order of 11 12 ... 1N; 22 ... 2N; ...; NN
  if (strcmp(command,"ebond1") == 0) {

    if (narg < 3) error->all(FLERR,"Illegal ebond1 command");
    nelement = atoi(arg[0]);   // num of elements

    memory->create(nsites_local,nelement,"app/erad:nsites_local");
    memory->create(ebond1,nelement+1,nelement+1,"app/erad:ebond1");
    memory->create(hcount,nlattice,"app/erad:hcount");
    if(concentrationflag) {
      memory->create(ct,nelement+1,"app/erad:ct"); //time averaged concentration 
      memory->create(ct_new,nelement+1,"app/erad:ct_new"); //time averaged concentration 
    } 

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
    memory->create(ebond2,nelement+1,nelement+1,"app/erad:ebond2");
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
    memory->create(mbarrier,nelement+1,"app/erad:mbarrier");

    for (i=0; i<narg-1; i++) {
      if(i % 2 == 0){ j = atoi(arg[i]);
        mbarrier[j] = atof(arg[i+1]);
      }
    }
  }

  // elastic moduli for stress calculation
  else if (strcmp(command, "percolation") ==0) {

    percolationflag = 1;
    if(narg != 1) error->all(FLERR,"illegal percolation command");
    dpercolation = atof(arg[0]); // 2D or 3D 
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
      memory->create(stress,nlattice,6,"app/erad:stress");
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
      memory->create(isink, nlattice, nelement+1,"app/erad:isink");
      for(i = 0; i < nlattice; i++) {
        for(j = 0; j < nelement+1; j++)  isink[i][j] = 0;
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
    if(narg < 1) error->all(FLERR,"illegal ballistic command");
    ballistic_flag = 1;
    grow_ballistic();

    double dose_rate=atof(arg[0]);// dose rate 
    bfreq[nballistic] = 1e12/nlocal/dose_rate; // time interval to introduce an FP 
    if(min_bfreq > bfreq[nballistic]) min_bfreq = bfreq[nballistic];
    nballistic ++; // number of mixing events
  }

  else if (strcmp(command, "migration_vector") ==0) {
    if(nelement < 3 || narg != (nelement-2)*3 + 2) error->all(FLERR,"illegal migration vector command");
    dmigration = atoi(arg[0]); // dimension of migration, default 3D  
    dratio = atof(arg[1]); // ratio of anisotropic diffusion   
    j = 2;
    for (i=0; i<narg-2; i++) {
      int k = i % 3;
      if(k == 0) j++;
      Vmig[j][k] = atof(arg[i+2]);
    }

    for (i=3; i<= nelement; i++) {
      Vmig[i][3] = Vmig[i][0]*Vmig[i][0] + Vmig[i][1]*Vmig[i][1] + Vmig[i][2]*Vmig[i][2];
    }
  }

  // meam hop steps before absorption by external sink 
  else if (strcmp(command, "mfp") ==0) {

    if (narg < 2 || narg % 2 == 0) error->all(FLERR,"Illegal mfp command");

    mfpflag = 1;
    memory->create(mfp,nelement+1,"app/erad:mfp");
    memory->create(nhmfp,nelement+1,"app/erad:nhmfp");
    memory->create(rhmfp,nelement+1,"app/erad:rhmfp");

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
    memory->create(ncelem,nctype,"app/erad:ncelem");
    for (i=0; i<nctype; i++) ncelem[i] = atoi(arg[i+2]);
  }

  else error->all(FLERR,"Unrecognized command");
}


/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppErad::grow_app()
{
  type = iarray[0];   // lattice type; i1 in input
  element = iarray[1];  // element type; i2 in input
  ipercolation = iarray[2]; // percolation analysis 
  //aid = iarray[2]; // initially set as global ID, must use set i3 unique in command line
  //if(concentrationflag) disp = darray; // msd; zero initially
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppErad::init_app()
{
  int i,j;

  // define second nearest neighbor
  if (engstyle == 2) define_2NN();

  if (firsttime) {
    firsttime = 0;

    echeck = new int[nlocal];
    memory->create(firstevent,nlocal,"app:firstevent");

    // esites must be large enough for 3 sites and their 1st neighbors
    esites = new int[3 + 3*maxneigh];
  }

  int flag = 0;
  for (i = 0; i < nelement; i++) nsites_local[i] = 0;
  // ecombine V and I if created on 1NN sites 
  // setup neighbor list for large recombination redius 
  if(rrecombine == 4) recombine_list(rrec);
  for (i = 0; i < nelement; i++) nrecombine[i] = 0;
  for (i = 0; i < nlocal; i++) recombine(i); 

  // site validity
  for (i = 0; i < nlocal; i++) {
    if (type[i] < BCC || type[i] > NBCC) flag = 1;
    if (element[i] < FE || element[i] > I6) flag = 1;
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

   //initialize the time_list for ballistic mixing
  if(concentrationflag) {
    //concentration_field();
    dt_new = 0.0; 
    for(i = 0; i <= nelement; i ++) {
       ct[i] = 0.0; 
       ct_new[i] = 0.0; 
    } 
    //for(i = 0; i < concentrationflag; i ++) {
    //   for(j = 0; j < nlocal; j++){
    //      disp[i][j] = 0.0;
    //   }
    //}
  }

  // initialize percolation list 
  if(percolationflag) {
    for (i = 0; i < nlocal; i++) ipercolation[i] = 0;
  }  

  /*/ initialize recombination vector analysis arrays, test only 
  for (i = 0; i < nlocal; i ++) {
  for (j = 0; j < 3; j ++) {
      xyz0[i][j] = xyz[i][j]; 
  }
  }
  for (i = 0; i < 361; i ++) {
  for (j = 0; j < 181; j ++) {
      reccount[i][j] = 0;
  }
  }
  // end test */
}

/* ---------------------------------------------------------------------- */
/* ---------------------------initiation of looing-up list -------------- */

void AppErad::setup_app()
{
  for (int i = 0; i < nlocal; i++) echeck[i] = 0;

  // clear event list
  nevents = 0;
  for (int i = 0; i < nlocal; i++) firstevent[i] = -1;
  for (int i = 0; i < maxevent; i++) events[i].next = i+1;
  freeevent = 0;

  // set propensities from rates
  if(temperature == 0.0) error->all(FLERR,"Temperature cannot be 0.0 for app erad");
  if(nn1flag == 0) error->all(FLERR, "First neareast neighbor bond not defined: AppErad");
  if(nn2flag == 0 && engstyle == 2) error->all(FLERR, "Second nearest neighbor bond not defined: AppErad");
  if(barrierflag == 0) error->warning(FLERR, "Diffusion barrier not defined: AppErad");

  double KB = 0.00008617;
  KBT = temperature * KB;
}

/* ----------------------------------------------------------------------
   define 2NN and 3NN list based on 1NN list. In BCC lattice, each 2NN of 
   site i is shared by 4 1NNs of site i as their 1NNs, and 3NN by 2. 
------------------------------------------------------------------------- */

void AppErad::define_2NN()
{ int candidate[144],frequency[144];
  int i,j,k,jd,kd,n1nn,n2nn,n3nn,njnn,ncandidate;
  int ndimension,n2freq,n3freq;

  ndimension = domain->dimension;
  n2freq = 4;
  if(ndimension == 2) n2freq = 2; // for both square and hex  
  n3freq = 2;
  if(ndimension == 2) n3freq = 1;
 
  memory->create(numneigh2,nlocal+nghost,"app, numneigh2");
  memory->create(neighbor2,nlocal+nghost,MAX2NN,"app, neighbor2");
  if(rhop == 3 || rrecombine == 3) memory->create(numneigh3,nlocal+nghost,"app, numneigh3");
  if(rhop == 3 || rrecombine == 3) memory->create(neighbor3,nlocal+nghost,MAX3NN,"app, neighbor3");
  for (i = 0; i < nlocal+nghost; i++) {
    for (j = 0; j < 144; j++) candidate[j] = 0;
    for (j = 0; j < 144; j++) frequency[j] = 0;

    ncandidate = 0;
    n1nn = numneigh[i];
    n2nn = 0;
    n3nn = 0;
    for (j =0; j < n1nn; j++)
    {
      jd = neighbor[i][j];
      njnn = numneigh[jd];
      for (k = 0; k < njnn; k++)
      {
        kd = neighbor[jd][k];

        if (kd != i) {
          candidate[ncandidate] = kd;
          ncandidate++;
        }
      }
    }

    for (j = 0; j < ncandidate; j++) {
      if(frequency[j] < 0) continue; // jd already observed earlier  
      jd = candidate[j];
      frequency[j]++;
      for (k = j+1; k < ncandidate; k++) {
        kd = candidate[k];
        if(kd == jd) {frequency[j]++; frequency[k] = -1;}
      }
    } 

    for (j = 0; j < ncandidate; j++) {
      jd = candidate[j];
      if(frequency[j] == n2freq) {
        int ifbreak = 0;
  	for (int m = 0; m < n1nn; m ++) { // To screen 1NNs for fcc crystals
	    if(jd == neighbor[i][m]) ifbreak = 1;
        }
        if (ifbreak) continue;

        if (n2nn == MAX2NN) error->all(FLERR, "Two many 2nd nearest neighbors defined, please expand the simulation cell dimensions or MAX2NN");
        neighbor2[i][n2nn] = jd;
        n2nn++; // selected if shared by 4 NNs
      } 
      else if ((rhop == 3 || rrecombine == 3) && frequency[j] == n3freq) { // 3NN 
        if (n3nn == MAX3NN) error->all(FLERR, "Two many 3rd nearest neighbors defined, please expand the simulation cell dimensions or MAX3NN");
        neighbor3[i][n3nn] = jd;
        n3nn++; // selected if shared by 4 NNs
      }
    }
    numneigh2[i] = n2nn;
    if(rhop == 3 || rrecombine == 3) numneigh3[i] = n3nn;
  }
}

/* ----------------------------------------------------------------------
   setup neighbor list for large recombination radius
   currently for serial runs only  
------------------------------------------------------------------------- */

void AppErad::recombine_list(double rrec)
{
  int i,j,id; 
  int ****ibox; 
  int ncell[3];
  int index[3],cell[3];
  double rij[4],lprd[3],xlo[3],rcell[3];

  int max4nn = 0;
  rrec *= rrec;
  // calculate # of neighbors within rrec 
  for(i = 0; i < nlocal; i++) {
     distanceIJ(i,0,rij);
     rij[3] = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
     if(rij[3] > rrec) continue;
     max4nn ++;
  }   
 
  memory->create(numneigh4,nlocal,"app/erad, numneigh4");
  memory->create(neighbor4,nlocal,max4nn,"app/erad, neighbor4");

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

  memory->create(ibox,ncell[0],ncell[1],ncell[2],2*max4nn,"apperad: ibox");
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

double AppErad::site_energy(int)
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   compute energy of site i
------------------------------------------------------------------------- */

double AppErad::sites_energy(int i, int estyle)
{
  int j,k,jd,n1nn;
  double eng = 0.0;
  double ci, cj;

  if(estyle == 0) return eng; // no thermodynamic interaction 
  //energy from 1NN bonds
  n1nn = numneigh[i];  //num of 1NN
  ci = cj = 0.0;
  int iele = element[i];
  for (j = 0; j < n1nn; j++) {
    jd = neighbor[i][j];
    if(cflag && iele == element[jd] && iele == VACANCY) {
      if(ci == 0.0) ci = site_concentration(i,1); 
      cj = site_concentration(jd,1);

      eng += ci*cj*ebond1[iele][element[jd]]; 
    } else 
    eng += ebond1[iele][element[jd]];
  }

  //energy from 2NN bonds
  if(estyle == 2) {
    ci = cj = 0.0;
    int n2nn = numneigh2[i];
    for (j = 0; j < n2nn; j++) {
      jd = neighbor2[i][j];
      if(cflag && iele == element[jd] && iele == VACANCY) {
        if(ci == 0.0) ci = site_concentration(i,2);
        cj = site_concentration(jd,2);

        eng += ci*cj*ebond2[iele][element[jd]];
      } else
      eng += ebond2[iele][element[jd]];
    }
  }

  return eng/2.0;
}

/* ----------------------------------------------------------------------
   compute elastic energy
------------------------------------------------------------------------- */

double AppErad::elastic_energy(int i, int itype)
{
  double eng = 0.0;
  double pressure = 0.0;

  pressure = -(stress[i][0] + stress[i][1] + stress[i][2])/3.0;
  eng = pressure * evol[itype];

  return eng;
}

/* ----------------------------------------------------------------------
  compute local concentration to adjust V-V bond energy  
------------------------------------------------------------------------- */

double AppErad::site_concentration(int i, int estyle)
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

double AppErad::site_SP_energy(int i, int j, int estyle)
{
  double eng = 0.0;
  double eng0i, eng0j, eng1i, eng1j; //energy before and after jump

  int iele = element[i];
  int jele = element[j];
  // no Fe hop event  
  if(iele < VACANCY) return -1.0;
  // no V-V or V-I or I-I switch
  if(jele >= VACANCY) return -1.0;

  // 1D % 2D SIA hop; no SIA bonding currently 
  if(iele > VACANCY) {
    double ran1 = ranerad->uniform(); 
    if(dmigration < 3 && ran1 < dratio) {
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
    
/*
  int jdele; 
  //bond formed between j and i's neighbors after switch, j can not be a vacancy 
  for (m = 0; m < numneigh[i]; m++) {
    jd = neighbor[i][m];
    jdele = element[jd];   
    if (jd == j) jdele = element[i]; // i & j switched  
    eng1i += ebond1[jele][jdele];
  }

  if (estyle == 2) {
    for(m = 0; m < numneigh2[i]; m++) {
      jd = neighbor2[i][m];
      jdele = element[jd];   
      if (jd == j) jdele = element[i]; // i & j switched  
      eng1i += ebond2[jele][jdele];
    }
  }

  //bond formed between i and j's neighbors after switch, i is a vavancy 
  if(cflag) ci = site_concentration(i,1); 
  for (m = 0; m < numneigh[j]; m++) {
    jd = neighbor[j][m];
    jdele = element[jd];   
    if (jd == i) jdele = element[j]; // i & j switched  

    if(cflag && jdele == VACANCY) { // only CC bond are concentration dependendent 
      cj = site_concentration(jd,1);
      eng1j += ci*cj*ebond1[iele][jeele];
    } else {
      eng1j += ebond1[iele][jdele];
    }
  }

  if (estyle == 2){
    if(cflag) ci = site_concentration(i,2); 
    for (m = 0; m < numneigh2[j]; m++) {
      jd = neighbor2[j][m];
      jdele = element[jd];   
      if (jd == i) jdele = element[j]; // i & j switched  

      if(cflag && jdele == VACANCY) {
        cj = site_concentration(jd,2);
        eng1j += ci*cj*ebond2[iele][jdele];
      } else
      eng1j += ebond2[iele][jdele];
    }
  }
*/

  eng = mbarrier[iele] + (eng1i + eng1j - eng0i -eng0j);

  //add elastic contribution if applicable
  if(elastic_flag) {
    double eng_ei = elastic_energy(j,iele) - elastic_energy(i,iele);
    double eng_ej = elastic_energy(i,jele) - elastic_energy(j,jele);

    eng += (eng_ei + eng_ej)/2.0;
  }

  return eng;
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppErad::site_propensity(int i)
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

  // for hop events, vacancy and SIA only currently
  // propensity calculated in site_SP_energy();
  if (element[i] < VACANCY) return prob_reaction;

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

  // 3NN hop for <110> 1D SIA diffusion in bcc 
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

void AppErad::site_event(int i, class RandomPark *random)
{
  int j,k,l,m,n,ii;

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
    hcount[i] ++;

    element[i] = element[j];
    element[j] = k;

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

    if(mfpflag && mfp[element[j]] > 0.0) {
      k = hcount[i]; 
      hcount[i] = hcount[j]; 
      hcount[j] = k; 

      //double r1 = (hcount[j] - mfp[element[j]]);
      //double fx = exp(-r1*r1/sigmamfp)/varmfp; 
      //if(ranerad->uniform() < fx) {
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
    if(sink_flag && isink[i][element[i]-1] == 1)   {absorption(i); }
    if(sink_flag && isink[j][element[j]-1] == 1)   {absorption(j); }
    //if(fabs(xyz[j][2]) < 2.0) fprintf(screen,"i %d %d %d, %f \n",j, element[j], isink[16791][2],xyz[j][2]); 
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
}

/* ----------------------------------------------------------------------
   absorption of an element at site i  
------------------------------------------------------------------------- */
void AppErad::absorption(int i)
{
  int k = element[i];

  double rand_me = ranerad->uniform();
  //if(rand_me <= 1.0/sink_mfp[n]) {
     //nabsorption[n] ++;
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
  //}
//	fprintf(screen,"%d %d %d %d %d %d %f,\n",k,isink[i][1],isink[i][2],isink[i][4],isink[i][5],isink[i][k-1],xyz[i][2]);
}

/* ----------------------------------------------------------------------
   calculate the vector between two recombined atoms i and j in vector form. 
------------------------------------------------------------------------- */
/*void AppErad::vec_recombine (int i, int j)
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
void AppErad::distanceIJ(int i, int j, double dij[3])
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
int AppErad::recombine(int i)
{ 
  // recombine by radius 
  if(rrecombine == 4) { 
    for (int n = 0; n < numneigh4[i]; n++) {
        int m = neighbor4[i][n];

        if(element[i] != VACANCY && element[m] != VACANCY) continue; // None vacancy 
        if(element[i] <= VACANCY && element[m] <= VACANCY) continue; // None SIA 
        //if(element[m] == VACANCY && element[i] <= VACANCY) return; // None SIA 
  
        nrecombine[element[i]-1] ++;
        nrecombine[element[m]-1] ++;
        nsites_local[element[i]-1] --;
        nsites_local[element[m]-1] --;
        nsites_local[FE-1] += 2;

        element[i] = FE;
        element[m] = FE;

        if(mfpflag) {hcount[i] = 0; hcount[m] = 0;}
  
       // update reaction target number
       if(reaction_flag == 1) {
         for(int k = 0; k < nreaction; k++) {
         if(routput[k] == element[i] || routput[k] == element[m])  target_local[k] ++;
         }
       }

       return m;
    }

    return -1;
  }

  for (int n = 0; n < numneigh[i]; n++) {
      int m = neighbor[i][n];

      if(element[i] != VACANCY && element[m] != VACANCY) continue; // None vacancy 
      if(element[i] <= VACANCY && element[m] <= VACANCY) continue; // None SIA 
      //if(element[m] == VACANCY && element[i] <= VACANCY) return; // None SIA 
         
      nrecombine[element[i]-1] ++;
      nrecombine[element[m]-1] ++;
      nsites_local[element[i]-1] --;
      nsites_local[element[m]-1] --;
      nsites_local[FE-1] += 2;

      element[i] = FE;
      element[m] = FE;

      if(mfpflag) {hcount[i] = 0; hcount[m] = 0;} 
      
     // update reaction target number
     if(reaction_flag == 1) {
       for(int k = 0; k < nreaction; k++) {
       if(routput[k] == element[i] || routput[k] == element[m])  target_local[k] ++;
       }
     }
     //fprintf(screen, "%d %d rec \n",i, m); 
     return m;  
  }

  // 2NN recombination 
  //if (rrecombine < 2) return -1; // no recombination takes place
  for (int n = 0; n < numneigh2[i]; n++) {
      int m = neighbor2[i][n];

      if(element[i] != VACANCY && element[m] != VACANCY) continue; // None vacancy 
      if(element[i] <= VACANCY && element[m] <= VACANCY) continue; // None SIA 
         
      nrecombine[element[i]-1] ++;
      nrecombine[element[m]-1] ++;
      nsites_local[element[i]-1] --;
      nsites_local[element[m]-1] --;
      nsites_local[FE-1] += 2;

      element[i] = FE;
      element[m] = FE;

      if(mfpflag) {hcount[i] = 0; hcount[m] = 0;} 
      
     // update reaction target number
     if(reaction_flag == 1) {
       for(int k = 0; k < nreaction; k++) {
       if(routput[k] == element[i] || routput[k] == element[m])  target_local[k] ++;
       }
     }

     return m;  
  }
  
  // 3NN recombination 
  if (rrecombine < 3) return -1; 
  for (int n = 0; n < numneigh3[i]; n++) {
      int m = neighbor3[i][n];

      if(element[i] != VACANCY && element[m] != VACANCY) continue; // None vacancy 
      if(element[i] <= VACANCY && element[m] <= VACANCY) continue; // None SIA 
         
      nrecombine[element[i]-1] ++;
      nrecombine[element[m]-1] ++;
      nsites_local[element[i]-1] --;
      nsites_local[element[m]-1] --;
      nsites_local[FE-1] += 2;

      element[i] = FE;
      element[m] = FE;

      if(mfpflag) {hcount[i] = 0; hcount[m] = 0;} 
      
     // update reaction target number
     if(reaction_flag == 1) {
       for(int k = 0; k < nreaction; k++) {
       if(routput[k] == element[i] || routput[k] == element[m])  target_local[k] ++;
       }
     }

     return m;  
  }

  return -1; 
}

/* ----------------------------------------------------------------------
   update propensity for site i and its neighbors after exchange
   ignore update of sites with i2site < 0; updated during each sweep
   use echeck[] to avoid resetting propensity of same site
------------------------------------------------------------------------- */
void AppErad::update_propensity(int i)
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
void AppErad::check_reaction()
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
          double rand_me = ranerad->uniform();
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
void AppErad::reset_propensity()
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

void AppErad::clear_events(int i)
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

void AppErad::add_event(int i, int j, int rstyle, int which, double propensity)
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
void AppErad::sia_concentration(double t)
{
  double ci = 0.0;
  for (int i = 0; i < nlocal; i ++ ) {if(element[i] >= I1) ci ++; }  
  csia += ci*t/nlocal;
  
}

/* ----------------------------------------------------------------------
   check if perform ballistic mixing
------------------------------------------------------------------------- */
void AppErad::check_ballistic(double t)
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
}

/* ----------------------------------------------------------------------
  Create Frenkel pairs randomly. may need to scale the dose rate 
  by the # of processors when  work in parallel   
------------------------------------------------------------------------- */

void AppErad::ballistic(int n)
{
  // creat an vacancy
  if(nsites_local[FE-1] == 0) error->all(FLERR, "No matrix sites available for FP generation!");      

  nFPair ++;   
  int findv = 1;   
  while (findv) { 
    int id = static_cast<int> (nlocal*ranerad->uniform()); 
    if(id < nlocal && element[id] == 1) { 
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
    int id = static_cast<int> (nlocal*ranerad->uniform()); 
    if(id < nlocal && element[id] == 1) { 
      if(nelement == VACANCY) error->all(FLERR, "simulation contains no SIAs"); 
      double dx = 1.0/(nelement - VACANCY);  // equal probability for each type of SIA 
      double r1 = ranerad->uniform(); 

      for (int i = 1; i <= nelement - VACANCY; i++) { 
          if(r1 < i*dx) { 
            element[id] = VACANCY + i; 
            nsites_local[VACANCY+i-1] ++; 
            nsites_local[FE-1] --;  
   
            // recalculate the propensity if defects are generated 
            update_propensity(id);      
            findi = 0;
            return;  
          } 
      } 
    }
  }  
}

/* ----------------------------------------------------------------------
  calculating total energy
------------------------------------------------------------------------- */
double AppErad::total_energy( )
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
  calculating percolation rate  
------------------------------------------------------------------------- */
double AppErad::percolation( )
{
  int i,j,jd,nvac,nperco,nnew; 
  double prate = 0.0;
  double xlo[3],xhi[3];

  xlo[0] = domain->boxxlo; 
  xhi[0] = domain->boxxhi; 
  xlo[1] = domain->boxylo; 
  xhi[1] = domain->boxyhi; 
  xlo[2] = domain->boxzlo; 
  xhi[2] = domain->boxzhi; 

  //reinitializa each time of calculations 
  nvac = nperco = 0; 
  for(i = 0; i < nlocal; i++) ipercolation[i] = 0;
 
  // find vacancies on boundaries 
  for(i = 0; i < nlocal; i++) { 
     if(element[i] != VACANCY) continue;  
     nvac ++; 

     if(xyz[i][0]-xlo[0] == 0 || xyz[i][0]-xhi[0] == 0 || xyz[i][1]-xlo[1] == 0 || xyz[i][1]-xhi[1] == 0) {       
       ipercolation[i] = 1;
       nperco ++;
     }

     if(dpercolation == 3 && ipercolation[i] == 0 && (xyz[i][2]-xlo[2] == 0 || xyz[i][2]-xhi[2] == 0)) {      
       ipercolation[i] = 1;
       nperco ++;
     }
  }

  //calculate transient percolation rate until no percolated vacancies can be found  
  nnew = 1; 
  ivisit = new int[nlocal];
  for(i = 0; i < nlocal; i++) ivisit[i] = 0;
 
  while(nnew > 0) {
     nnew = 0; 
     for(i = 0; i < nlocal; i++) { 
        if(ipercolation[i] != 1) continue;  
        if(ivisit[i] == 1) continue;
        for(j = 0; j < numneigh[i]; j++) { 
           jd = neighbor[i][j];
           if(element[jd] != VACANCY || ipercolation[jd] == 1) continue;
           nnew ++; 
           ipercolation[jd] = 1;
           nperco ++; 
        }
        ivisit[i] = 1;
     } 
  }

  delete [] ivisit;
  if(nvac > 0) prate = nperco*1.0/nvac; 
  return prate;
}

/* ----------------------------------------------------------------------
  Integrate c*t at each site for fractional occupancy over time 
------------------------------------------------------------------------- */
void AppErad::concentration_field(double dt)
{
  dt_new += dt; // update time interval 
  for(int i = 0; i < nlocal; i++) {
     //disp[element[i]][i] += dt;
     ct_new[element[i]] += dt;      
  } // c[element[i]][i] = 1; disp += c[element[i][i] * dt  }
}

/* ----------------------------------------------------------------------
  update time averaged concentrations 
------------------------------------------------------------------------- */
void AppErad::time_averaged_concentration()
{
  for(int i = 1; i <= nelement; i++) { // element starts from FE=1  
     if(dt_new > 0.0) ct[i] = ct_new[i]/dt_new/nlocal; 
     ct_new[i] = 0.0;
  }
  dt_new = 0.0;   
}

/* ----------------------------------------------------------------------
  map n by n matrix to a vector
------------------------------------------------------------------------- */
int AppErad::ibonde(int a, int b, int c)
{
  return ((a-1)*c + b - a*(a-1)/2);
}

 
/* ----------------------------------------------------------------------
  match element. Return 0 if match, 1 otherwise 
------------------------------------------------------------------------- */ 
int AppErad::match(int i)
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

void AppErad::cluster()
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
  grow memory for dislocation
------------------------------------------------------------------------- */

void AppErad::grow_dislocations()
{
  int n = ndislocation + 1;
  memory->grow(dislocation_type,n,"app/erad:dislocation_type");
  memory->grow(burgers,n,3,"app/erad:burgers");
  memory->grow(xdislocation,n,3,"app/erad:xdislocation");
  memory->grow(line_vector,n,"app/erad:line_vector");
  memory->grow(dislocation_radius,n,"app/erad:dislocation_radius");
  memory->grow(nsegment,n,"app/erad:nsegment");
}

/* ----------------------------------------------------------------------
  grow memory for sink
------------------------------------------------------------------------- */

void AppErad::grow_sinks()
{
  int n = nsink + 1;
  memory->grow(sink_shape,n,"app/erad:sink_shape");
  memory->grow(sink_strength,n,"app/erad:sink_strength");
  memory->grow(xsink,n,3,"app/erad:xsink");
  memory->grow(sink_type,n,"app/erad:sink_type");
  memory->grow(sink_normal,n,"app/erad:sink_normal");
  memory->grow(sink_segment,n,"app/erad:sink_segmant");
  memory->grow(sink_radius,n,"app/erad:sink_radius");
  memory->grow(sink_mfp,n,"app/erad:sink_mfp");
  memory->grow(nabsorption,n,"app/erad:nabsorption");
}

/* ----------------------------------------------------------------------
  grow memory for reaction
------------------------------------------------------------------------- */

void AppErad::grow_reactions()
{
  int n = nreaction + 1;
  memory->grow(rsite,n,"app/erad:rsite");
  memory->grow(rinput,n,"app/erad:rinput");
  memory->grow(routput,n,"app/erad:routput");
  memory->grow(rcount,n,"app/erad:rcount");
  memory->grow(renable,n,"app/erad:renable");
  memory->grow(rtarget,n,"app/erad:rtarget");
  memory->grow(rbarrier,n,"app/erad:rbarrier");
  memory->grow(rrate,n,"app/erad:rrate");
  memory->grow(target_local,n,"app/erad:target_local");
  memory->grow(target_global,n,"app/erad:target_global");
}

/* ----------------------------------------------------------------------
  grow memory for ballistic mixing
------------------------------------------------------------------------- */

void AppErad::grow_ballistic()
{
  int n = nballistic + 1;
  int m = 3;
  memory->grow(bfreq,n,"app/erad:bfreq");
  memory->grow(time_old,n,"app/erad:time_old");
  memory->grow(time_new,n,"app/erad:time_new");
}

/* ----------------------------------------------------------------------
  create sinks
------------------------------------------------------------------------- */

void AppErad::sink_creation(int n)
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

void AppErad::stress_field(int n)
{
  if(dislocation_type[n] == 1) stress_dislocation(n); //straight dislocation
  else stress_loop(n); //loop
}

/* ----------------------------------------------------------------------
  calculate stess field for straight dislocations
------------------------------------------------------------------------- */

void AppErad::stress_dislocation(int n)
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

void AppErad::stress_loop(int n)
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

void AppErad::seg_stress( double A[3], double B[3], double P[3], double bv[3], double sstres[3][3])
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

void AppErad::sigma_A(double t[3], double m[3], double bv[3], double sigma[3][3])
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

void AppErad::sigma_P(double t[3], double ni[3], double bv[3], double sigmap[3][3])
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

void AppErad::stroh( double tvect[3], double qmatx[3][3], double bmatx[3][3], double smatx[3][3])
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

void AppErad::stroh_p( double t[3], double n0[3], double qpmat[3][3], double bpmat[3][3], double spmat[3][3])
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

void AppErad::elastic_tensor()
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

void AppErad::matrix_inversion(double mm[3][3], double nn[3][3])
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

void AppErad::cross_product( double m[3], double n[3], double l[3])
{
  l[0] = m[1]*n[2] - m[2]*n[1];
  l[1] = m[2]*n[0] - m[0]*n[2];
  l[2] = m[0]*n[1] - m[1]*n[0];

}

/* ----------------------------------------------------------------------
  cross product of int vectors
------------------------------------------------------------------------- */

void AppErad::cross_product( int m[3], int n[3], int l[3])
{
  l[0] = m[1]*n[2] - m[2]*n[1];
  l[1] = m[2]*n[0] - m[0]*n[2];
  l[2] = m[0]*n[1] - m[1]*n[0];

}

/* ----------------------------------------------------------------------
  normalize double vector
------------------------------------------------------------------------- */

void AppErad::vector_normalize( double m[3])
{
  double norm = sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
  if(norm == 0) error->all(FLERR,"can not normalize zero vector");
  for(int i = 0; i < 3; i ++) m[i] /= norm;

}

/* ----------------------------------------------------------------------
  right-hand coordination system, m x n x t
------------------------------------------------------------------------- */

void AppErad::right_hand_coord( double m[3], double n[3], double t[3])
{
  double taux[3];
  // random vector taux to find m
  taux[0] = ranerad->uniform();
  taux[1] = ranerad->uniform();
  taux[2] = ranerad->uniform();

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

void AppErad::trapozidal(double omegas[100], double XP[9][100], int j, int k, double value)
{ int i, ii;
  double dy,domega;

  for(i = 0; i < ninteg; i ++) {
    ii = i + 1;
    domega = omegas[ii] - omegas[i];
    dy = XP[3*j+k][ii] + XP[3*j+k][i];
    value += dy*domega/2.0;
  }
}
