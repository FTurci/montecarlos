#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#include "cell.h"
class Particle {

 public:
  int index; //particle index
  int type;
  
  double R[3];     // center of mass
  double Rmin[3];     // center of mass
  Particle *next, *prev; // next, previous in cell list
  Cell *cell;
  
  bool ghost=false;

  Particle(){
      next = 0;
      prev = 0;
      cell = 0;
  }

  void insertToCell(Cell &c){

    prev = 0;               // initialize pointer
    next = c.firstParticle; // -> firstParticle must be initialized to 0!
    if (next) next->prev=this;
    c.firstParticle=this;
    cell=&c;
    
  };

  void moveBetweenCells(Cell &from, Cell &to){

    // first remove the particle from the old cell
    if (prev) prev->next=next;
    else from.firstParticle=next;
    if (next) next->prev=prev;

    // now add particle
    prev=0;
    next=to.firstParticle;
    if (next) next->prev=this;
    to.firstParticle=this;
    cell=&to;
  };


};
#endif