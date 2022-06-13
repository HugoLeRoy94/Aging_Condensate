#include "Header.h"
using namespace std;
int main(int argc, char* argv[]){
  double ell_tot(100.);
  double distance_anchor(50.);
  double rho0(0.01);
  double temperature(0.00001);
  System* S = new System(ell_tot,distance_anchor,rho0,temperature,198730);
  //S->Print_Loop_positions();
  for(int i=0;i<10000;i++){
    cout<<i<<endl;
    //try{
      S->evolve();
      cout<<"number of loop in the system = "<<S->get_N()<<endl;
      cout<<"number of crosslinkers in the system = "<<S->get_r_size()<<endl;
      cout<<endl;
    }
    //catch(const exception& e){cout<<e.what()<<endl;}
  //}
  //S->Print_Loop_positions();
  return 0;
}
