#include "Header.h"
using namespace std;
int main(int argc, char* argv[]){
  double ell_tot(100.);
  double distance_anchor(50.);
  double rho0(0.01);
  double temperature(0.1);
  System* S = new System(ell_tot,distance_anchor,rho0,temperature,198730);
  //S->Print_Loop_positions();
  for(int i=0;i<20;i++){
    cout<<i<<endl;
    //try{
      S->evolve();
    }
    //catch(const exception& e){cout<<e.what()<<endl;}
  //}
  //S->Print_Loop_positions();
  return 0;
}
