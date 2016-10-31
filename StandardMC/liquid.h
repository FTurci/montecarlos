#ifndef __LIQUID_H__
#define __LIQUID_H__

#define maxNrpart 50000    // maximum number of particles read from file

#include "particle.h"
#include "box.h"
#include "cell.h"

#include <random>
#include <vector>
#include <fstream>
// #include  <gsl/gsl_errno.h>
// #include <gsl/gsl_math.h>
// #include <gsl/gsl_vector.h>     
// gsl vector stuff
// #include <gsl/gsl_multimin.h>
// gsl multidimensional minimization 

class liquid {

 public:


  const char *ParameterFile, *ConfigFile;  // Input files

  std::mt19937 gen;
  static std::uniform_real_distribution<double> u;

  void usage();
  void run();
  double uniform();

  Particle *part;
  Particle **inWell;
  Cell *cells;
  Box *box;


  int N;             // number of particles
  int NA,NB;
  double Ninv;       // 1/N
  double beta;       // 1/k_BT
  double LJcutoff;   // cutoff for Lennard Jones potential
  double shift;      // shifted to zero at cutoff
  double interactionRange2; // LJcutoff^2

  bool withCells;    // If systems is smaller than 3 times cutoff radius in 
                     // any dimension, cell system is not used

  // Free energy calculation is done for linear potential well
  // and tanh blending function

  double rcut;     // cutoff for reference potential well
  double rcutinv, rcut2, rcut3, Dcut;
  double rcutAA, rcutBB, rcutAB;
  double rcutAA2, rcutBB2, rcutAB2;

  // staged insertion
  bool stageOnOff;
  int numStages=4;
  int stage;
  int  ghostIndex;
  long int stageAccept;
  double ene_min[4]={ -1.5, -1.5, -1.5, -1.5};
  double ene_max[4]={  10,  3  , 11, 6};

  double cumulatedEnergy;
  double cumulatedCount;
  std::vector<double> Energies;
  std::vector<long int> timestamps;
  long int tictac;

  int transAccept, transMoves;
  int relocAccept, relocMoves;
  int swapAccept, swapMoves;
  int eqSteps;          
  int eqEvalSteps;      // # blockaverages taken during equilibration.
  int eqSnapSteps;      // # configuration snapshots taken during equil.
  int mcSteps;
  int mcEvalSteps;      // # measurements for averages
  int mcSnapSteps;      // # snapshots during mc

  double maxDisplace;

  long seed;          //seed for ran3
  long *iran;         //pointer to seed

  int measureCount;
  double *meanHLJ, *meanHref; // <H_LJ>, <H_ref>
  double *meandHdl;   //   < dblend1/dlambda H_ref + dblend2/dlambda H_LJ>
  double mHLJ, mHref, mdHdl;
  double HLJ, Href, dHdl, Htotal;
  double inW;

  void setUpFromFile();
  void setUpFromXYZ();
  void initialize();
  void setupList();
  void setUpNeighbours();


  void equil();
  void MC();

  void particleMove();
  void stageMove(bool sampling);


  void computeAllEnergies(bool countGhost);
  double PotentialEnergy(bool countGhost);    

  void displaceParticle(double newP[3], double oldP[3]);
  void relocateParticle(Particle &p);

  double interactionAll(Particle &p);
  double pairInteraction(Particle &p1, Particle &p2,bool countGhost);

  // io
  std::ofstream trjfile;
  void writeOutConfig(char name[16]);
  void appendOutConfig();
  void writeOutSimData(char name[16]);
  void writeOutAllSim(char name[16]);

  double minimise();
  // static double my_f ( int N3,double X[]);

  // Kob-Andersen constants (cutoff 2.5 sigmaij)
  double shiftAA;
  double shiftBB;
  double shiftAB;

};


template< typename F >
  class gsl_function_pp : public gsl_function {
  public:
  gsl_function_pp(const F& func) : _func(func) {
    function = &gsl_function_pp::invoke;
    params=this;
  }
  private:
  const F& _func;
  static double invoke(double x, void *params) {
    return static_cast<gsl_function_pp*>(params)->_func(x);
  }
};
#endif
