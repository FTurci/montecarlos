#include "liquid.h"
#include <iostream>
#include <random>
#include <fstream>

using namespace std;

int main(int argc, char *argv[]){
  
  std::mt19937 gen (123);
 
  liquid Model;

  if (argc < 3) {
    cerr << "Too few arguments \n";
    Model.usage();
  }

  Model.gen=gen;

  while ((argc > 1) && (argv[1][0] == '-')) {

    switch (argv[1][1]) {
    
    case 'f':
      Model.ConfigFile=&argv[1][2];
      break;
    case 'p' :
      Model.ParameterFile=&argv[1][2];
      break;
    default:
      cerr << "Bad option " << argv[1] << '\n';
      Model.usage();
    }

    ++argv;
    --argc;
  }

  cout << "# Parameterfile: " << Model.ParameterFile << endl;


  Model.run();

  char fname[256];
  sprintf(fname,"energies-T%g-stage%d.txt", 1./Model.beta, Model.stageOnOff);
  ofstream fout(fname,ios::out);
  for (unsigned int i = 0; i < Model.Energies.size(); ++i)
  {
    fout<<Model.timestamps.at(i)<<" "<<Model.Energies.at(i)<<endl;
  }
  fout.close();
  Model.trjfile.close();

}