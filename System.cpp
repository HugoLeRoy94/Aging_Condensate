#include "Header.h"

using namespace std;

System::System(){}
System::System(double ell_tot,double distance_anchor,double rho0,double temperature){
  // ---------------------------------------------------------------------------
  // -----------------------constant of the simulation--------------------------
  ell = ell_tot;
  D = distance_anchor;
  rho0 = rho;
  kBT = temperature;
  // ---------------------------------------------------------------------------
  //-----------------------------initialize loops-------------------------------
  // *reserve* makes the vector to never reallocate memory, and fasten the push_backs
  loops.reserve(ell_tot); // maximum number of loop that can be inserted
  array<double,3> Rf={D,0.,0.},R0={0.,0.,0.};
  loops.push_back(new Loop(R0,Rf,ell)); // create the first loop.
}
System::~System(){
  for(auto& it : loops){delete it;}
  loops.clear();
}
// -----------------------------------------------------------------------------
// ----------------------------Main function------------------------------------
// -----------------------------------------------------------------------------
double System::evolve(){
  double time(0);
  return time;
}
// -----------------------------------------------------------------------------
// -----------------------------accessor----------------------------------------
// -----------------------------------------------------------------------------
int System::get_N() const{return loops.size();}

void System::get_R(double* R, int size) const{
  if(size!=3*loops.size()){
    throw invalid_argument("invalid size in System::get_R");}
  // fill the vector R with the value of R of each loops
  for(int n=0;n<loops.size();n++){
    R[3*n] = loops[n]->get_Rright()[0];
    R[3*n+1] = loops[n]->get_Rright()[1];
    R[3*n+2] = loops[n]->get_Rright()[2];
  }
}

void System::get_ell(double* ells,int size) const{
  if(size!=loops.size()){throw invalid_argument("invalid size in System::get_ell");}
  for(int n=0;n<size;n++){
    ells[n] = loops[n]->get_ell();
  }
}
