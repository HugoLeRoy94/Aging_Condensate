#include "Header.h"
using namespace std;
int main(int argc, char* argv[]){
  double ell_tot(100.);
  double distance_anchor(50.);
  double rho0(0.1);
  double temperature(0.1);
  System* S = new System(ell_tot,
                         distance_anchor,
                         rho0,
                         temperature);
  int size(S->get_N());
  double ell[size];
  double R[3*size];
  try{
    S->get_R(R,3*size);
  }
  catch(invalid_argument& message){
    cout<<message.what()<<endl;
  }
  try{
    S->get_ell(ell,size);
  }
  catch(invalid_argument& message){
    cout<<message.what()<<endl;
  }
  for(int n=0;n<size;n++)
  {
    cout<<ell[n]<<" ";
  }
  cout<<"\n";
  for(int n=0;n<size*3;n+=3)
  {
    cout<<R[n]<<" "<<R[n+1]<<" "<<R[n+2]<<"\n";
  }
  return 0;
}
