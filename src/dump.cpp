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
#include "string.h"
#include "stdio.h"
#include "stdlib.h"
#include "dump.h"
#include "app.h"
#include "app_lattice.h"
#include "app_off_lattice.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

/* ---------------------------------------------------------------------- */

Dump::Dump(SPPARKS *spk, int narg, char **arg) : Pointers(spk)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // parse dump args

  if (narg < 4) error->all(FLERR,"Illegal dump command");

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  n = strlen(arg[1]) + 1;
  style = new char[n];
  strcpy(style,arg[1]);

  delta = atof(arg[2]);
  if (delta <= 0.0) error->all(FLERR,"Illegal dump command");

  n = strlen(arg[3]) + 1;
  filename = new char[n];
  strcpy(filename,arg[3]);

  // parse filename for special syntax
  // if contains '%', write one file per proc and replace % with proc-ID
  // if contains '*', write one file per timestep and replace * with timestep
  // check file suffixes
  //   if ends in .bin = binary file
  //   else if ends in .gz = gzipped text file
  //   else ASCII text file

  compressed = 0;
  binary = 0;
  multifile = 0;
  multiproc = 0;

  char *ptr;
  if ((ptr = strchr(filename,'%'))) {
    multiproc = 1;
    char *extend = new char[strlen(filename) + 16];
    *ptr = '\0';
    sprintf(extend,"%s%d%s",filename,me,ptr+1);
    delete [] filename;
    n = strlen(extend) + 1;
    filename = new char[n];
    strcpy(filename,extend);
    delete [] extend;
  }

  if (strchr(filename,'*')) multifile = 1;

  char *suffix = filename + strlen(filename) - strlen(".bin");
  if (suffix > filename && strcmp(suffix,".bin") == 0) binary = 1;
  suffix = filename + strlen(filename) - strlen(".gz");
  if (suffix > filename && strcmp(suffix,".gz") == 0) compressed = 1;

  // on-lattice or off-lattice app

  if (app->appclass == App::LATTICE) {
    applattice = (AppLattice *) app;
    latticeflag = 1;
  } else if (app->appclass == App::OFF_LATTICE) {
    appoff = (AppOffLattice *) app;
    latticeflag = 0;
  } else
    error->all(FLERR,"Dump command can only be used for spatial applications");

  // dump params

  boxxlo = domain->boxxlo;
  boxxhi = domain->boxxhi;
  boxylo = domain->boxylo;
  boxyhi = domain->boxyhi;
  boxzlo = domain->boxzlo;
  boxzhi = domain->boxzhi;

  logfreq = 0;
  delay = 0.0;
  flush_flag = 1;
  padflag = 0;

  maxbuf = 0;
  buf = NULL;

  idump = 0;
}

/* ---------------------------------------------------------------------- */

Dump::~Dump()
{
  delete [] id;
  delete [] style;
  delete [] filename;

  memory->sfree(buf);

  if (multifile == 0 && fp != NULL) {
    if (compressed) {
      if (multiproc) pclose(fp);
      else if (me == 0) pclose(fp);
    } else {
      if (multiproc) fclose(fp);
      else if (me == 0) fclose(fp);
    }
  }
}

/* ---------------------------------------------------------------------- */

void Dump::init()
{
  init_style();
}

/* ----------------------------------------------------------------------
   dump a snapshot of site values as atom coords
------------------------------------------------------------------------- */

void Dump::write(double time)
{
  // if file per timestep, open new file

  if (multifile) openfile();

  // nmine = # of dump lines this proc will contribute to dump
  // ntotal = total # of dump lines
  // nmax = max # of dump lines on any proc

  int nme = count();
  bigint bnme = nme;

  bigint ntotal;
  int nmax;
  if (multiproc) nmax = nme;
  else {
    MPI_Allreduce(&bnme,&ntotal,1,MPI_SPK_BIGINT,MPI_SUM,world);
    MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);
  }

  // write timestep header

  if (multiproc) write_header(bnme,time);
  else write_header(ntotal,time);

  idump++;

  // grow communication buffer if necessary

  if (nmax*size_one > maxbuf) {
    maxbuf = nmax*size_one;
    memory->sfree(buf);
    buf = (double *) memory->smalloc(maxbuf*sizeof(double),"dump:buf");
  }

  // pack my data into buf

  pack();

  // multiproc = 1 = each proc writes own data to own file 
  // multiproc = 0 = all procs write to one file thru proc 0
  //   proc 0 pings each proc, receives it's data, writes to file
  //   all other procs wait for ping, send their data to proc 0

  if (multiproc) write_data(nme,buf);
  else {
    int tmp,nlines;
    MPI_Status status;
    MPI_Request request;

    if (me == 0) {
      for (int iproc = 0; iproc < nprocs; iproc++) {
	if (iproc) {
	  MPI_Irecv(buf,maxbuf,MPI_DOUBLE,iproc,0,world,&request);
	  MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	  MPI_Wait(&request,&status);
	  MPI_Get_count(&status,MPI_DOUBLE,&nlines);
	  nlines /= size_one;
	} else nlines = nme;

	write_data(nlines,buf);
      }
      if (flush_flag) fflush(fp);
      
    } else {
      MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
      MPI_Rsend(buf,nme*size_one,MPI_DOUBLE,0,0,world);
    }
  }

  // if file per timestep, close file

  if (multifile) {
    if (compressed) {
      if (multiproc) pclose(fp);
      else if (me == 0) pclose(fp);
    } else {
      if (multiproc) fclose(fp);
      else if (me == 0) fclose(fp);
    }
  }
}

/* ---------------------------------------------------------------------- */

void Dump::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all(FLERR,"Illegal dump_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"flush") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      if (strcmp(arg[iarg+1],"yes") == 0) flush_flag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) flush_flag = 0;
      else error->all(FLERR,"Illegal dump_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"delta") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      delta = atof(arg[iarg+1]);
      if (delta <= 0.0) error->all(FLERR,"Illegal dump_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"logfreq") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal dump_modify command");
      nrepeat = atoi(arg[iarg+1]);
      scale = atof(arg[iarg+2]);
      if (nrepeat < 0) error->all(FLERR,"Illegal dump_modify command");
      if (nrepeat == 0) logfreq = 0;
      else logfreq = 1;
      iarg += 3;
    } else if (strcmp(arg[iarg],"loglinfreq") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal dump_modify command");
      nrepeat = atoi(arg[iarg+1]);
      scale = atof(arg[iarg+2]);
      if (nrepeat < 0) error->all(FLERR,"Illegal dump_modify command");
      if (nrepeat == 0) logfreq = 0;
      else logfreq = 2;
      iarg += 3;
    } else if (strcmp(arg[iarg],"delay") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      delay = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"pad") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump_modify command");
      padflag = atoi(arg[iarg+1]);
      if (padflag < 0) error->all(FLERR,"Illegal dump_modify command");
      iarg += 2;

    } else {
      int n = modify_param(narg-iarg,&arg[iarg]);
      if (n == 0) error->all(FLERR,"Illegal dump_modify command");
      iarg += n;
    }
  }
}

/* ----------------------------------------------------------------------
   generic opening of a dump file
   ASCII or binary or gzipped
------------------------------------------------------------------------- */

void Dump::openfile()
{
  // if one file per timestep, replace '*' with current timestep

  char *filecurrent;
  if (multifile == 0) filecurrent = filename;
  else {
    filecurrent = new char[strlen(filename) + 16];
    char *ptr = strchr(filename,'*');
    *ptr = '\0';
    if (padflag == 0) 
      sprintf(filecurrent,"%s%d%s",filename,idump,ptr+1);
    else {
      char bif[8],pad[16];
      strcpy(bif,"%d");
      sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
      sprintf(filecurrent,pad,filename,idump,ptr+1);
    }
    *ptr = '*';
  }

  // open one file on proc 0 or file on every proc

  if (me == 0 || multiproc) {
    if (compressed) {
#ifdef SPPARKS_GZIP
      char gzip[128];
      sprintf(gzip,"gzip -6 > %s",filecurrent);
      fp = popen(gzip,"w");
#else
      error->one(FLERR,"Cannot open gzipped file");
#endif
    } else if (binary) {
      fp = fopen(filecurrent,"wb");
    } else {
      fp = fopen(filecurrent,"w");
    }

    if (fp == NULL) error->one(FLERR,"Cannot open dump file");
  } else fp = NULL;

  // delete string with timestep replaced

  if (multifile) delete [] filecurrent;
}
