#include "Header.h"
using namespace std;
int main(int argc, char* argv[]){
  int t_tot(2000);
  double ell_tot(2000.);
  //double distance_anchor(1000.);
  double rho0(0.0001);
  double temperature(0.1);
  bool bind(true);
  System* S = new System(ell_tot,rho0,temperature,19830,false);
  S->evolve(&bind);
//for(auto& it : S->get_r()){for(auto r : it){cout<<r<<" ";}cout<<endl;}
  //cout<<S->evolve(&bind)<<endl;
  //S->print_random_stuff();
/*for(int n = 0; n<t_tot;n++){
    cout<<"step : "<<n<<endl;
    cout<<S->evolve(&bind)<<endl;
    cout<<bind<<endl;}*/
  
   
   /*
  map3d<double, double ,double, double> M3D;
  M3D(0,0,0) = 1.5;
  double to_insert(9100.);
  //double* address(&to_insert);
  M3D(1.1,0,0) = to_insert;
  //double* address(&M3D(1.1,0,0));
  M3D(0,1.2,0) = 9010.;
  M3D(0,0,1.3) = 9001.;
  M3D(2,0,0) = 9200.;
  M3D(0,2,0) = 9020.;
  M3D(0,0,2) = 9002.;
  M3D.cut_slice(0,1.,0,1.,0,1.);
  double* address(M3D.add_return_address(3,4,5,9345));
  for(auto& it : M3D.underlying_array()){for(auto& it2 : it.second){for(auto& it3 : it2.second){cout<<it3.second<<endl;}}}
  cout<<M3D(3,4,5)<<endl;
  cout<<&M3D(3,4,5)<< " "<<address<<endl;*/
  }
