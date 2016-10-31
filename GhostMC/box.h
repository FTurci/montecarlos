#ifndef __BOX_H__
#define __BOX_H__
class Box {

 public:

  double x[3], halfx[3], halfxInv[3];  // dimensions
  double V;      // Volume
  double xCell[3], xCellInv[3]; //Cell dimensions
  int nxCell[3]; // number of cells
  int nCells;

  void update(){
      for (int i=0; i<3; i++){
      halfx[i] = 0.5*x[i];
      halfxInv[i] = 1.0/halfx[i];
      }
      V = x[0]*x[1]*x[2];
  }

};
#endif