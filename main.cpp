#include "Header.h"
using namespace std;
int main(int argc, char* argv[]){
  double ell_tot(100.);
  double distance_anchor(50.);
  double rho0(0.1);
  double temperature(0.1);
  System* S = new System(ell_tot,distance_anchor,rho0,temperature);
  S->Print_Loop_positions();
  S->evolve();
  S->Print_Loop_positions();
  return 0;
}
