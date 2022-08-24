#include "Header.h"
using namespace std;
int main(int argc, char* argv[]){
  int t_tot(2000);
  double ell_tot(2000.);
  //double distance_anchor(1000.);
  double rho0(0.00001);
  double temperature(0.05);
  bool bind(true);
  System* S = new System(ell_tot,rho0,temperature,198730,false);
//for(auto& it : S->get_r()){for(auto r : it){cout<<r<<" ";}cout<<endl;}
  //cout<<S->evolve(&bind)<<endl;
for(int n = 0; n<t_tot;n++){
    cout<<"step : "<<n<<endl;
    cout<<S->evolve(&bind)<<endl;
    cout<<bind<<endl;}
}
