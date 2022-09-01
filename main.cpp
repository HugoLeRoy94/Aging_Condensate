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
  
  /* 
  map3d<double, double ,double, double> M3D;
  M3D(0,0,0) = 1.5;
  M3D(1.1,0,0) = 9100.;
  M3D(0,1.2,0) = 9010.;
  M3D(0,0,1.3) = 9001.;
  M3D(2,0,0) = 9200.;
  M3D(0,2,0) = 9020.;
  M3D(0,0,2) = 9002.;
  M3D.cut_slice(0,1.5,0,1.,0,1.);
  cout<<endl;
  for(auto& it : M3D.underlying_array()){for(auto& it2 : it.second){for(auto& it3 : it2.second){cout<<it3.second<<endl;}}}
  */
}
