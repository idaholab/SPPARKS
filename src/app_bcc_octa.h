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

#ifdef APP_CLASS
AppStyle(bccocta,AppBccOcta)
#else

#ifndef SPK_APP_BCC_OCTA_H
#define SPK_APP_BCC_OCTA_H

#include "app_lattice.h"

namespace SPPARKS_NS {

class AppBccOcta : public AppLattice {
  friend class DiagBccOcta;

 public:
  AppBccOcta(class SPPARKS *, int, char **);
  ~AppBccOcta();
  void input_app(char *, int, char **);
  void grow_app();
  void init_app();
  void setup_app();

  double site_energy(int);
  void site_event_rejection(int, class RandomPark *) {}
  double site_propensity(int);
  void site_event(int, class RandomPark *);

 private:
  int engstyle,rhop,rrecombine,nn1flag,nn2flag,barrierflag,mfpflag,cflag; // 1NN or 2NN bonds
  int *type,*element,*iHe; // variables on each lattice site
  int firsttime;

  int *esites;
  int *echeck;
 
  int *numneigh2,*numneigh3,*numneigh4; // number of 2NN and 3NN
  int **neighbor2,**neighbor3,**neighbor4; // index of 2NN and 3NN
  int ****ibox;//used for defining recombination radius
  int nelement,dmigration,nFPair; //total element, dimension of migraiton
  int *nsites_local; //statics of local element 
  double KBT,dratio,sigmamfp,varmfp,rrec; //rate and propensity for hop diffusion 
  double **ebond1,**ebond2; //bond energy 
  double *mbarrier,*mfp,*rhmfp; //migration barriers
  double Vmig[12][4]; //migration vector
  int *hcount,*nhmfp; 

  int nrecombine[10]; //number of recombination event for each element 

  class RandomPark *ranbccocta; //random number generator 

//parameter for sinks 
  int nsink,sink_flag;
  int *sink_type,*sink_shape,*sink_segment,*sink_normal,*nabsorption;
  int **isink;
  double *sink_strength,*sink_radius,*sink_mfp;
  double **xsink;  

//parameter for reaction
  int nreaction;  
  int *rsite,*rinput,*routput,*rcount,*renable,*rtarget;
  int *target_local,*target_global; 
  double *rbarrier,*rrate; // rrate is scaled by the attempt rate of atom hopping 

//paramter for vacancy trapping  
//  int itrap,itime_current,itime_old;
//  double dt_interval,fvt,dt_real,dt_akmc,treal_me,takmc_me;  

//parameter for ballistic mixing 
  int nballistic;
  int *time_old,*time_new;
  double *bfreq,gas_rate;

//parameter for cluster analysis  
  int nctype,ncsize,ncluster;
  int *ncelem,*cid,*csize;
  double rcluster,csia;
 
//parameter for acceleration 
  int ntrap; 
  int *trap_type;

  struct Event {           // one event for an owned site
    int style;             // reaction style = HOP,RECOMBINE 
    int which;             // which reaction of this type
    int jpartner;          // which J neighbors of I are part of event
    int next;              // index of next event for this site
    double propensity;     // propensity of this event
  };

  Event *events;           // list of events for all owned sites
  int nevents;             // # of events for all owned sites
  int maxevent;            // max # of events list can hold
  int *firstevent;         // index of 1st event for each owned site
  int freeevent;           // index of 1st unused event in list

  int ibonde(int, int, int);  //list to matrix bond energy
  void define_2NN();
  void clear_events(int);
  void add_event(int, int, int, int, double);
  void update_propensity(int);
  double add_acceleration_event(int, int); // add acceleration event and return a probability
  double total_energy();
  double sites_energy(int, int);
  double site_SP_energy(int, int, int);
  double site_concentration(int, int);  
  void distanceIJ(int, int, double []); 
  void absorption(int);  
  void recombine_list(double);
  void vec_recombine(int,int);
  int recombine(int);  
  
  void grow_reactions(); //reactions
  void check_reaction(); 
  void reset_propensity(); 

  void grow_ballistic();// ballistic mixing 
  void ballistic(int); 
  void check_ballistic(double); 
  void sia_concentration(double); // track sia concentration 

  void grow_sinks(); //sink
  void sink_creation(int);

  int match(int); //cluster
  void cluster();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Unrecognized command

The command is assumed to be application specific, but is not
known to SPPARKS.  Check the input script.

E: One or more sites have invalid values

The application only allows sites to be initialized with specific
values.

E: Temperature cannot be 0.0 for appbccocta 

UNDOCUMENTED

*/
