#include "Header.h"
using namespace std;
int main(int argc, char* argv[]){
  int t_tot(10000);
  double ell_tot(100.);
  double distance_anchor(50.);
  double rho0(0.1);
  double temperature(0.1);
  bool bind(true);
  System* S = new System(ell_tot,distance_anchor,rho0,temperature,198730);
  for(int n = 0; n<t_tot;n++){
    cout<<n<<endl;
  cout<<S->evolve(&bind)<<endl;
  cout<<bind<<endl;}
  return 0;
}
