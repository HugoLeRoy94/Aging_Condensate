#include "Header.h"
using namespace std;
int main(int argc, char* argv[]){
  int t_tot(100);
  double ell_tot(200.);
  //double distance_anchor(1000.);
  double rho0(0.1/ell_tot);
  double BindingEnergy(-10);
  bool bind(true);
  double* R;
  System* S = new System(ell_tot,rho0,BindingEnergy,19874,false);
  for(int n(0);n<t_tot;n++){
  S->evolve(&bind);
  cout<<S->get_F()<<endl;
  }
  //cout<<bind<<endl;

  /*
  R = new double[3*S->get_N_strand()];
  S->get_R(R,3*S->get_N_strand());
  cout<<bind<<endl;
  cout<<"binding points :"<<endl;
  for(int n =0;n<3*S->get_N_strand();n++)
  {
    cout<<R[n]<<" ";
    if(n%3==2){cout<<endl;}
  }
  double* r = new double[S->get_r_system_size()];
  S->get_r_system(r,S->get_r_system_size());
  cout<<"linkers :"<<endl;
  for(int n =0;n<S->get_r_system_size();n++)
  {
    cout<<r[n]<<" ";
    if(n%3==2){cout<<endl;}
  }
  */
   
  //for(int n = 0;n<100;n++){cout<<n<<endl;S->evolve(&bind);}
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
