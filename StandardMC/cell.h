#ifndef __CELL_H__
#define __CELL_H__
#include "particle.h"

class Particle;

class Cell {

 public:

  Particle *firstParticle;
  Cell *neighbours[26];      // neighbouring cells 
  Cell(void){
      firstParticle = 0;
      for (int i=0; i<26; i++){
      neighbours[i] = 0;
      }
  }

};
#endif