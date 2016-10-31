
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "liquid.h"
#include "particle.h"
#include "box.h"
#include "cell.h"
# include "compass_search.hpp"


using namespace std;

const double Pi = M_PI;
std::uniform_real_distribution<double> liquid::u(0,1);


double ljpow(double r, double epsilon, double sigma){
  return 4.0*epsilon*(pow(sigma/r,12) - pow(sigma/r,6));
}
double liquid::uniform(){
  return ((double) rand() / (RAND_MAX));
  // return u(gen);
}

void liquid::run(){

  setUpFromXYZ();
  setupList();
  cout << "# Eqsteps = " << eqSteps << " MCsteps = "<<mcSteps<<endl;
  computeAllEnergies(false);
  cout<<"# Starting energy: "<<Htotal<<endl;
  cout<<Htotal/N<<endl;

  double minimum=minimise();
  cout<<"The minimum is "<<minimum<<endl;
  // writeOutConfig("Start");
  // equil();
  // MC();
  cout<<"Summarising..."<<endl;
  // writeOutConfig("Config");
  // writeOutSimData("Stats");
}

void liquid::equil(){

    char outstring[16];
    char temp[16];
    int eqInterval;
    if (eqEvalSteps > 0 )
      eqInterval = eqSteps/eqEvalSteps;
    else eqInterval = 0;
    int measuring = eqInterval;
    int snapInterval;
    if (eqSnapSteps > 0)
  snapInterval = eqSteps/eqSnapSteps;
    else snapInterval = 0;
    int snapshot = snapInterval;
    

    transAccept = 0; transMoves = 0;
    relocAccept = 0; relocMoves = 0;
    swapAccept = 0; swapMoves = 0;

    // int oldeqSteps=eqSteps;
    // if(iter==0)eqSteps*=1;
    // else eqSteps=oldeqSteps;
    for (int steps=1; steps<=eqSteps; steps++){
  
  measuring--;
  snapshot--;

  for (int i=0; i<N; i++){  
    // cout<<"moving"<<endl;    
    particleMove();
    
  }
  if(stageOnOff==true)stageMove(false);
  tictac++;

  if (measuring == 0){
    computeAllEnergies(false);
    cout << "# Eq  Energy StageAccept Stage" << endl;
    cout << steps << " "
         << Htotal/N << " "
         << stageAccept/(double)(steps)<<" "
         << stage 
         <<  endl;
      measuring = eqInterval;
  }
  
  if(snapshot == 0){
      snapshot = snapInterval;
      strcpy(outstring, "CoEQ");
      sprintf(temp,"%.*d", 0, int(steps/snapInterval));
      strcat(outstring, temp);
      writeOutConfig(outstring);
  }
    }
    cout << "# Equilibration ended" << endl;
    
}



void liquid::MC(){
  
  char outstring[16];
  char temp[16];
  int mcInterval;
  if (mcEvalSteps > 0 )
      mcInterval = mcSteps/mcEvalSteps;
  else mcInterval = 0;
  int measuring = mcInterval;
  int snapInterval;
  if (mcSnapSteps > 0)
      snapInterval = mcSteps/mcSnapSteps;
  else snapInterval = 0;
  int snapshot = snapInterval;
  
  measureCount = 0;
  inW = 0.0;

  for (int steps=1; steps<=mcSteps; steps++){
    
    measuring--;
    snapshot--;
    for (int i=0; i<N; i++){ 
      particleMove();
     
    }

    if(stageOnOff==true) stageMove(true);
    tictac++;

    if (measuring == 0) {
      computeAllEnergies(false);
      measureCount++;
      if(stage==numStages) 
        appendOutConfig();
      
      if (measuring == 0){
  cout << "# MC  Energy StageAccept Stage" 
       << endl;
  cout << steps + eqSteps << " " 
       << Htotal /N << " "
       << stageAccept/(double)(steps + eqSteps)/N<<" "
       << stage
       << endl;
  measuring = mcInterval;
      }
      
      
      if(snapshot == 0) {
        snapshot = snapInterval;
        strcpy(outstring, "CoMC");
        sprintf(temp,"%.*d", 0, int(steps/snapInterval));
        strcat(outstring, temp);
        writeOutConfig(outstring);
      }
    }
  }
  printf("measure count %i\n", measureCount);
}


void liquid::particleMove(){
  
  int pNr;
  double Rold[3];
  int newCell;
  Cell *oldCell;
  double oldE, newE;
  Particle *p;


  pNr = int((N)* uniform());
  p = &part[pNr];
  
  for (int i=0; i<3; i++){
    Rold[i] = p->R[i];
  }
  

  if(withCells){
    oldCell = p->cell;
    oldE = interactionAll(*p);
    displaceParticle(p->R, Rold); 
  
    newCell = int((p->R[0]+box->halfx[0])*box->xCellInv[0]) + 
      int((part[pNr].R[1]+box->halfx[1])*box->xCellInv[1])*box->nxCell[0] +
      int((part[pNr].R[2]+box->halfx[2])*box->xCellInv[2])*
      box->nxCell[0]*box->nxCell[1];
    
    if (&cells[newCell] != oldCell){
      p->moveBetweenCells(*oldCell, cells[newCell]);
    }
    
    newE = interactionAll(*p);
    if ( uniform() < exp(beta*(oldE- newE))) transAccept ++;
  
    else {
      for (int i=0; i<3; i++){
       p->R[i] = Rold[i];
      }
      if (&cells[newCell] != oldCell){
       p->moveBetweenCells(cells[newCell], *oldCell);
      }
    }
    
    transMoves ++;
  }

  else{ //if not using the cell system
    oldE = 0;
    for (int i=0; i<N; i++){
     if (&part[i] != p){
       oldE += pairInteraction(part[i],*p, true);

      }
    }
    
    displaceParticle(p->R, Rold);
    
    newE = 0;
    for (int i=0; i<N; i++){
      if (&part[i] != p){
        newE += pairInteraction(part[i],*p, true);
      }
    }

    if (log( uniform()) < beta*(oldE- newE)) 
    {
      transAccept ++;
    }
   
    else { // restore positions
      for (int i=0; i<3; i++){
        p->R[i] = Rold[i];
      }
    }
    transMoves ++;
  }

}

void liquid::stageMove(bool sampling){
  double eghostOld=0, eghostNew=0;

 
  int stageOld=this->stage;
  // cout<<ghostIndex<<endl;

  for (int j=0; j<N; j++) {
    if(j !=ghostIndex)eghostOld+=pairInteraction(part[ghostIndex], part[j],true);
  }
        
  // propose stage move
  if (stage>0){
    if(uniform()<0.5) stage++;
    else stage--;
  }
  else stage++;

  if (stage==numStages){ // the ghost particle becomes real
    part[ghostIndex].ghost=false;
    if(sampling){
    computeAllEnergies(false);
    // cmulate the energy of a system of only real particles
    cumulatedEnergy+=Htotal/N;
    cumulatedCount++;
    Energies.push_back(Htotal/N);
    timestamps.push_back(tictac);
  }

  }
  else if(stage>numStages){
        // select another particle randomly and make it ghost
    int selected=(int)(uniform()*NA+NB);
    part[selected].ghost=true;
    ghostIndex=selected;
    // resetting the stage

    stage=0;
  }
  else{

    for (int j=0; j<N; j++) {
      if(j !=ghostIndex) eghostNew+=pairInteraction(part[ghostIndex], part[j],true);
    }
    //  accept:
    if(log(uniform())<beta*(eghostOld-eghostNew)){
      stageAccept++;
    }
    else{//reject
      stage=stageOld;
    }
  }

}

void liquid::usage(){

  cerr << "Usage is monte -f<ConfigFile> -p<ParameterFile> \n";
  exit(8);

}




void liquid::displaceParticle(double newP[3], double oldP[3]) {

 double displace;

  for (int i=0; i<3; i++){
 
    displace = 2*( uniform() - 0.5)*maxDisplace;

    newP[i] =oldP[i]+displace;

    if(newP[i]>box->x[i]) newP[i]-=box->x[i];
    else if (newP[i]<0) newP[i]+=box->x[i];

  }

}

// computes LJ interaction energy of particle p with all particles from cell list

double liquid::interactionAll(Particle &p) {

  Cell *curCell;
  Particle *curPart;
  double ene = 0;
  
  curCell = p.cell;
  curPart = curCell->firstParticle;

  // cell of p
  while (curPart){ 
      if (curPart != &p) 
  ene += pairInteraction(p, *curPart, true);
      curPart = curPart->next;
  }

  // neighbouring cells
  for(int j=0; j<26; j++){
      curCell = p.cell->neighbours[j];
      curPart = curCell->firstParticle;
      while (curPart){
    ene += pairInteraction(p, *curPart,true);
    curPart = curPart->next; 
      }
  }

  return ene;
} 


// computes interaction energy between two particles
double liquid::pairInteraction(Particle &p1, Particle &p2, bool countGhost) {

  double d[3], d2,id2, d6inv, d12inv;
  double ene = 0.0;
  for (int j=0; j<3; j++){
      d[j] = fabs(p1.R[j] - p2.R[j]);
      if(d[j]>box->halfx[j]) d[j]-=box->x[j];
  }

  d2 = d[0]*d[0]+d[1]*d[1]+d[2]*d[2];

  id2 = 1.0/d2;
  d6inv = id2*id2*id2;
  d12inv = d6inv*d6inv;



  if(p1.type==2 && p2.type==2){

    // cout<<p1.index<<" is B and "<<p2.index<<" is B"<<endl;
    if(d2<rcutBB2){
       ene= 4.0*0.5*(d12inv*.215671156 - d6inv*.464404087)-shiftBB;
    }

  }
  else if((p1.type==2 && p2.type>=1) || (p1.type==1 && p2.type==2) )
  {

    if(d2<rcutAB2){
       ene = 4.0*1.5*(d12inv*.068719477 - d6inv*.262144)-shiftAB;
      if(countGhost)  {
          if((p1.ghost || p2.ghost) )
              { 
                // staged potential
                if (d2<0.8*0.8)
                  ene=min(ene, ene_max[stage]);
                else
                  ene=max(ene, ene_min[stage]);
              }

               //if neither is a ghost, ene is unperturbed
     }
    }
  }
  else{
     if(d2<rcutAA2){
       ene = 4.0*(d12inv - d6inv) - shiftAA;
    }
  }

      // Lennardjonesium
      // ene = 4.0*(d12inv - d6inv) - shift;
  

  return ene;
}




double liquid::PotentialEnergy(bool countGhost) {

  double ene=0.0;

  if (withCells){
    for(int i=0; i<N; i++){
      ene+=interactionAll(part[i]);
    }
    return 0.5*ene;
  }
  
  else{

    for(int i=0; i<N; i++){
      for (int j=i+1; j<N; j++){
        ene+=pairInteraction(part[i], part[j], countGhost);
      }
    }
    return ene; 
  } 

}


void liquid::computeAllEnergies(bool countGhost){

  HLJ = PotentialEnergy( countGhost);
  Htotal = HLJ;

  return;
}





void liquid::initialize(){
    char filename[256];
    sprintf(filename, "trajectory-T%g-stage%d.xyz",1./beta, stageOnOff);
    cout<<"Opening file "<<filename<<endl;

    trjfile.open(filename);
    Ninv = 1.0/float(N);
    NA=ceil(N*0.8); //number of A particles in the Kob-Andersen Mixture
    NB=N-NA;
    //truncation and shift
    // interactionRange2 = LJcutoff*LJcutoff;

    // shift = 4*(pow(LJcutoff,-12) - pow(LJcutoff,-6));
    // cout << "# cutoff " << LJcutoff << " shift " << shift << endl; 

    // rcutinv = 1.0/rcut;
    // rcut2 = rcut*rcut;
    // rcut3 = rcut2*rcut;
    // Dcut = 2.0*rcut;

    if(withCells) setupList();

    rcutAA=2.5;
    rcutBB=2.2;//2.5*0.88;
    rcutAB=2.;//2.5*0.8;

    rcutAA2=rcutAA*rcutAA;
    rcutBB2=rcutBB*rcutBB;
    rcutAB2=rcutAB*rcutAB;
    shiftAA=ljpow(rcutAA,1,1);
    shiftAB=ljpow(rcutAB,1.5,0.8);
    shiftBB=ljpow(rcutBB,0.5, 0.88);


    
    int selected=(int)(uniform()*NA+NB);
    if(stageOnOff==true){
      part[selected].ghost=true;
      // set an intermediate stage
      stage=ceil(numStages/2.);
      ghostIndex=selected;
    }
    else{
      ghostIndex=-1;
      stage=numStages;
    }
    cumulatedCount=0;
    cumulatedEnergy=0;
}




void liquid::setupList(){

  int icell;
  double longestRange = LJcutoff;
  if(rcut > LJcutoff) longestRange = rcut;

  box->update();

  // dimensions
  for (int i=0;i<3;i++){
    box->nxCell[i] = int(box->x[i]/longestRange);
    box->xCell[i] = box->x[i]/box->nxCell[i];
    box->xCellInv[i] = 1.0/box->xCell[i];
  }
  box->nCells = box->nxCell[0]*box->nxCell[1]*box->nxCell[2];

   cout<<"The box dimensions are "<<box->x[0]<<" "<<box->x[1]<<" "<<box->x[2]<<endl;

  if ((box->nxCell[0] < 3)||(box->nxCell[1] < 3)||(box->nxCell[2] < 3)){
    cout << "# Box dimension is less than 3*largest cutoff." << endl;
    cout << "# Not using cell system." << endl;
    withCells = false;
  }
  else {

    cout << "# Box dimension is larger than 3*largest cutoff." << endl;
    cout << "# Using cell system." << endl;
    withCells = true;
  }


  if (withCells==true){
  cells = new Cell[box->nCells];

  // put particles into cells

  for (int i=0; i<N; i++){
      icell = int((part[i].R[0]+box->halfx[0])*box->xCellInv[0]) + 
    int((part[i].R[1]+box->halfx[1])*box->xCellInv[1])*box->nxCell[0] +
    int((part[i].R[2]+box->halfx[2])*box->xCellInv[2])*
    box->nxCell[0]*box->nxCell[1];
      if (icell < box->nCells) part[i].insertToCell(cells[icell]);
      else {
    cerr << "Error: problem with cell assignment!" << endl;
    cerr << "Positions: " << part[i].R[0] << " " << part[i].R[1] << " " 
         << part[i].R[2] << endl;
    cerr << "Box: " << box->x[0] << " " << box->x[1] << " "
         << box->x[2] << endl;
    cerr << "icell: " << icell << " ncells " << box->nCells << endl;
    exit(8);
      }
    }
  
  // find neighbouring cells
  setUpNeighbours();
  }
}



void liquid::setUpNeighbours(){

  int x,y,z,count,cl,cln;
   
  for (int k=0; k<box->nxCell[2]; k++){
    for (int j=0; j<box->nxCell[1]; j++){
      for (int i=0; i<box->nxCell[0]; i++){
  count = 0;
  for (int zs = -1; zs <= 1; zs++){
    // Modulokonstr. period. Randbed.
    z = (k+zs + 10*box->nxCell[2])%box->nxCell[2]; 
    for (int ys = -1; ys <= 1; ys++){
      y = (j+ys + 10*box->nxCell[1])%box->nxCell[1];
      for (int xs = -1; xs <= 1; xs++){
        // don't count the middle cell
        if (!((xs==0)&&(ys==0)&&(zs==0))){
    x = (i+xs + 10*box->nxCell[0])%box->nxCell[0];
    cl = i + box->nxCell[0]*j + box->nxCell[1]*box->nxCell[0]*k;
    cln = z*box->nxCell[1]*box->nxCell[0] + y*box->nxCell[0] + x;
    cells[cl].neighbours[count] = &cells[cln];
    count ++;
        }
      }
    }
  }
      }
    }
  }
} 




// reads simulation parameters in case -f<configFile> ("set up from 
// configuration file") is chosen and reads particle positions from configFile
  
void liquid::setUpFromFile() {

   ifstream infile, infile2;
   char line[200];
   char quantity[200];
   float number;
   int check;
   // read new simulation parameters from ParameterFile
   infile.open(ParameterFile);
   if (infile.bad()) {
       cerr << "Error " << ParameterFile << "not found\n";
       exit(8);
   }
   while (infile.peek() != EOF) {
       infile.getline(line,sizeof(line));
       check = sscanf(line, "%s%f", quantity, &number);
       cout<<line<<endl;
       if (check != 2) {
     cerr << "Error: wrong line format in " << ParameterFile << endl;
     exit(8);
       }

       if (strcmp(quantity,"Temperature")==0) beta=1.0/number; 
       else if (strcmp(quantity,"McSteps")==0) mcSteps=int(number);
       else if (strcmp(quantity,"McEvalSteps")==0) mcEvalSteps=int(number);
       else if (strcmp(quantity,"McSnapSteps")==0) mcSnapSteps=int(number);
       else if (strcmp(quantity,"EqSteps")==0) eqSteps=int(number); 
       else if (strcmp(quantity,"EqEvalSteps")==0) eqEvalSteps=int(number);
       else if (strcmp(quantity,"EqSnapSteps")==0) eqSnapSteps=int(number);
       else if (strcmp(quantity,"MaxDisplacement")==0) maxDisplace=number;
       else if (strcmp(quantity,"Seed")==0) seed=int(number);
       else if (strcmp(quantity,"LJCutoff")==0) LJcutoff=number;
       else if (strcmp(quantity,"StageOnOff")==0) stageOnOff=number>0;
   }
   
  infile.close();
  
  // read configuration from ConfigFile

  int allSet = 0;
  int count = 0;
  double Rin[3];
  double dump[maxNrpart][3];

  infile2.open(ConfigFile);

  if (infile2.bad()) {
      cerr << "Error " << ConfigFile << "not found\n";
      exit(8);
  }

  box = new Box;
 
  while (infile2.peek() != EOF) {
      infile2.getline(line,sizeof(line));
      // cout<<"ALL "<<allSet<<endl;
      if(allSet == 5){

      check = sscanf(line, "%le %le %le", &Rin[0], &Rin[1], &Rin[2]);
      
      // cerr<<"check "<<check<<" "<<line<<endl;
      
      if (check == 3) { 
          for (int j=0; j<3; j++){
              dump[count][j] = Rin[j];
          }
          count++; 
      }
      else { 
          cerr << "Error: wrong line format in body of " 
         << ConfigFile << endl;
          exit(8);
        }
    } 
      
      else {
    check = sscanf(line, "%s%f", quantity, &number);
    if (check != 2) {
        cerr << "Error: wrong line format in header of " 
       << ConfigFile << endl;
        exit(8);
    }
    if (strcmp(quantity,"#Number")==0){
      N=int(number); 
      if (N > maxNrpart){
        cerr << "N particles is " << N << endl;
        cerr << "Too many particles in " << ConfigFile << "!" << endl;
        cerr << "Please change maxNrpart in headerfile." << endl;
        exit(8);
      }
      allSet++;
    }
    else if (strcmp(quantity,"#Pressure")==0) {allSet++;}
    else if (strcmp(quantity,"#Density")==0) {allSet++;}
    else if (strcmp(quantity,"#Boxx")==0) {box->x[0]=number; allSet++;}
    else if (strcmp(quantity,"#Boxy")==0) {box->x[1]=number; allSet++;}
    else if (strcmp(quantity,"#Boxz")==0) {box->x[2]=number; allSet++;}
      }
  }
  infile2.close();

  part = new Particle[N];

  for (int i=0; i<N; i++){
    //assign index
    part[i].index=i;
    //get coordinates
    for (int j=0; j<3; j++){
      part[i].R[j] = dump[i][j];
    }
  }

  if (count==0) {
    cerr << "Error: parameter missing in header of " << ConfigFile << endl;
    exit(8);
  }

  else if (count!=N) {
    cerr << "Error: wrong number of particle coordinates in " 
   << ConfigFile << endl;
    exit(8);
  }

  initialize();

  cout << "# read parameters and configuration" << endl;
  cout << "# N " << N <<" Temperature "<<1./beta<< " beta " << beta << endl;

}
void liquid::setUpFromXYZ() {

   ifstream infile, infile2;
   char line[200];
   char quantity[200];
   float number;
   int check;
   // read new simulation parameters from ParameterFile
   infile.open(ParameterFile);
   if (infile.bad()) {
       cerr << "Error " << ParameterFile << "not found\n";
       exit(8);
   }
   while (infile.peek() != EOF) {
       infile.getline(line,sizeof(line));
       check = sscanf(line, "%s%f", quantity, &number);
       cout<<line<<endl;
       if (check != 2) {
     cerr << "Error: wrong line format in " << ParameterFile << endl;
     exit(8);
       }

       if (strcmp(quantity,"Temperature")==0) beta=1.0/number; 
       else if (strcmp(quantity,"McSteps")==0) mcSteps=int(number);
       else if (strcmp(quantity,"McEvalSteps")==0) mcEvalSteps=int(number);
       else if (strcmp(quantity,"McSnapSteps")==0) mcSnapSteps=int(number);
       else if (strcmp(quantity,"EqSteps")==0) eqSteps=int(number); 
       else if (strcmp(quantity,"EqEvalSteps")==0) eqEvalSteps=int(number);
       else if (strcmp(quantity,"EqSnapSteps")==0) eqSnapSteps=int(number);
       else if (strcmp(quantity,"MaxDisplacement")==0) maxDisplace=number;
       else if (strcmp(quantity,"Seed")==0) seed=int(number);
       else if (strcmp(quantity,"LJCutoff")==0) LJcutoff=number;
       else if (strcmp(quantity,"StageOnOff")==0) stageOnOff=number>0;
   }
   
  infile.close();
  
  // read configuration from ConfigFile
  double Rin[3];
  double dump[maxNrpart][3];

  infile2.open(ConfigFile);

  if (infile2.bad()) {
      cerr << "Error " << ConfigFile << "not found\n";
      exit(8);
  }

  box = new Box;
 char t;

 infile2.getline(line,sizeof(line));
 sscanf(line, "%d", &N);
 cout<<"The number of particles is "<<N<<endl;
 infile2.getline(line,sizeof(line));
 sscanf(line, "%s %le %le %le", quantity, &box->x[0],&box->x[1],&box->x[2]);

 box->update();

 for (int i = 0; i < N; ++i)
 {
      infile2.getline(line,sizeof(line));
      check = sscanf(line, "%s %le %le %le", &t, &Rin[0], &Rin[1], &Rin[2]);
      
          for (int j=0; j<3; j++){
              dump[i][j] = Rin[j];
              
          }
  }
  infile2.close();

  part = new Particle[N];

  for (int i=0; i<N; i++){
    //assign index
    part[i].index=i;
    //get coordinates
    for (int j=0; j<3; j++){
      part[i].R[j] = dump[i][j];
         
    }
 
  }


  initialize();

  cout << "# read parameters and configuration" << endl;
  cout << "# N " << N <<" Temperature "<<1./beta<< " beta " << beta << endl;

}



void liquid::writeOutConfig(char name[16]){

  char temp[16];
  char allname[50];  strcpy(allname, name);
  strcat(allname, "_");
  strcat(allname, temp);
  strcat(allname, "_");
  sprintf(temp,"%.*lf", 2, 1.0/beta);//Temperature
  strcat(allname, temp);
  strcat(allname, "_");
  sprintf(temp,"%.*lf", 2, N/box->V);
  strcat(allname, temp);



  ofstream outfile(allname, ios::out);
  outfile.precision(10);

  outfile<<N<<"\nBox "<< box->x[0] <<" "<< box->x[1] <<" "<< box->x[2] <<endl;

  for (int i=0; i<N; i++){
    if (i<NB)
    outfile <<"B "<< part[i].R[0] << " " << part[i].R[1]  << " "<< part[i].R[2] << endl;
    else
    outfile <<"A "<< part[i].R[0] << " " << part[i].R[1]  << " " << part[i].R[2] << endl;
  }

  outfile.close();
}


void liquid::appendOutConfig(){


  trjfile<<N<<"\nBox "<< box->x[0] <<" "<< box->x[1] <<" "<< box->x[2] <<endl;

  for (int i=0; i<N; i++){
    if (ghostIndex==i)
      trjfile <<"G "<< part[i].R[0] << " " << part[i].R[1]  << " "<< part[i].R[2] << endl;
    else if (i<NB)
    trjfile <<"B "<< part[i].R[0] << " " << part[i].R[1]  << " "<< part[i].R[2] << endl;
    else
    trjfile <<"A "<< part[i].R[0] << " " << part[i].R[1]  << " " << part[i].R[2] << endl;
  }

  
}



void liquid::writeOutSimData(char name[16]){

  char temp[16];
  char allname[50];
  
  strcpy(allname, name);
  strcat(allname, "_");
  strcat(allname, temp);
  strcat(allname, "_");
  sprintf(temp,"%.*lf", 2, 1.0/beta);
  strcat(allname, temp);
  strcat(allname, "_");
  sprintf(temp,"%.*lf", 2, N/box->V);
  strcat(allname, temp);
  
  
  ofstream outfile(allname, ios::out);
  
  outfile << "# Trans Rate: " 
    << float(transAccept)/float(transMoves) << endl;

  outfile.close();
}

void liquid::writeOutAllSim(char name[16]){

  char temp[16];
  char allname[50];
  
  strcpy(allname, name);
  strcat(allname, "_");
  sprintf(temp,"%.*lf", 2, 1.0/beta);
  strcat(allname, temp);
  strcat(allname, "_");
  sprintf(temp,"%.*lf", 2, N/box->V);
  strcat(allname, temp);

  
  ofstream outfile(allname, ios::out);
  
  outfile << "# Free energy of ideal gas (per particle) = " << 1/beta*log(N/box->V) << endl;
  outfile << "# lambda <HLJ>/N <Href>/N <dHdlambda>/N" << endl;
  
  outfile.close();
}


 double  my_f ( int N3,double X[])
{  
  int NN=N3/3;
  Particle *parts=new Particle[NN];
  for (int i = 0; i < NN; ++i)
  {
    for (int k = 0; k < 3; ++k)
    {
      parts[i].R[k] = X[i*3+k];
    }
    
  }

  double ene=0;
    for(int i=0; i<NN; i++){
      for (int j=i+1; j<NN; j++){
        ene+=pairInteraction(parts[i], parts[j], false);
        cout<<ene<<endl;
      }
    }

    return ene; 

}
double liquid::minimise(){
    int dim=3*N;
    double *X=new double[dim];
    double *Xmin=new double [dim];
    for (int i = 0; i < N; ++i)
    {
      for (int k = 0; k < 3; ++k)
      {
        X[3*i+k]=part[i].R[k];
      }
    }
    int steps_taken;
    double value;

    

    Xmin=compass_search(my_f, dim, X, 0.0000000001, 0.01,10000,value, steps_taken, this);
    cout<<"Value "<<value<<" steps taken "<<steps_taken<<endl;
    return value;

  // return status;
}





