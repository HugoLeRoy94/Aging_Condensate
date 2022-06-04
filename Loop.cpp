#include "Header.h"
using namespace std;

Loop::Loop(array<double,3> R0,array<double,3> R1, double ell_loop,double rho){
  // the right anchoring point of the loop
  rho0 = rho;
  Rleft=R0;
  Rright=R1;
  ell = ell_loop;
  // Compute the volume covered by the polymer
  if(diff(Rleft,Rright)<0.1*ell){
    V =4/3*Pi*pow(ell,1.5);
    a=sqrt(ell)*0.5;
    b=sqrt(ell)*0.5;
  }
  else{
    V = Pi/6*diff(Rleft,Rright)*ell;
    a=diff(Rleft,Rright)*0.5;
    b=sqrt(ell)*0.5;
  }
  //cout<<"a and b = "<<a<<" "<<b<<endl;
  // create r by randmoly generating positions
  generate_binding_sites();
  // compute the array of binding rate to any crosslinker at any length
  for(auto& rlink : r){
    vector<double> rates_ell;
    rates_ell.reserve((int)ell);
    for(int ELL=0;ELL<(int)ell;ELL++){
      rates_ell.push_back(Omega(Rleft,rlink,(double)ELL)*Omega(rlink,Rright,ell-(double)ELL));
    }
    rates.push_back(rates_ell);
  }
  // compute the total binding rate:
  for(auto& ell_r : rates){
    for(auto& rate: ell_r){
      total_rates+=rate;
    }
  }
}
Loop::~Loop(){}
// ---------------------------------------------------------------------------
//-----------------------------------accessor---------------------------------
// ---------------------------------------------------------------------------
std::array<double,3> Loop::get_Rright() const{return Rright;}
std::array<double,3> Loop::get_Rleft() const{return Rleft;}
std::vector<array<double,3>> Loop::get_r() const{return r;}
double Loop::get_theta()const{return atan2(0.5*(Rright[1]-Rleft[1]),0.5*(Rright[0]-Rleft[0]));}
double Loop::get_phi()const{return atan2(0.5*(Rright[0]-Rleft[0]),0.5*(Rright[2]-Rleft[2]))-Pi/2.;}
array<double,3> Loop::get_Rg()const{return {0.5*(Rright[0]+Rleft[0]),0.5*(Rright[1]+Rleft[1]),0.5*(Rright[2]+Rleft[2])};}
double Loop::get_ell() const{return ell;}
double Loop::get_V()const{return V;}
double Loop::get_total_binding_rates() const{return total_rates;}
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
double Loop::compute_binding_rate(int r_index, double ell){

  return 0;
}
void Loop::select_link_length(double& length, array<double,3>& r_selected,double& time){
  // let select a link and a length uniformly for now
  uniform_int_distribution<int> r_distrib(0,r.size());
  r_selected =r[r_distrib(generator)];

  //cout<<"find a length between"<<diff(Rleft,r_selected)<<"and "<<ell-diff(r_selected,Rright)<<endl;
  uniform_real_distribution<double> ell_distrib(diff(Rleft,r_selected),ell-diff(r_selected,Rright));
  length =ell_distrib(generator);

  time = 1/(Omega(Rleft,r_selected,length)*Omega(r_selected,Rright,ell-length));

}
void Loop::generate_binding_sites(){
  default_random_engine generator;
  poisson_distribution<int> distribution(rho0*V);
  int N_crosslinker=distribution(generator);
  for(int n=0;n<N_crosslinker;n++){
    r.push_back(random_in_ellipse(a,b,b,0.5*(Rright[0]+Rleft[0]),0.5*(Rright[1]+Rleft[1]),0.5*(Rright[2]+Rleft[2])));
  }
  //cout<<"size of r "<<r.size()<<endl;
  //for(auto& it : r){cout<<it[0]<<" "<<it[1]<<" "<<it[2]<<endl;}
}
double Loop::Omega(array<double,3> r1,array<double,3> r2,double ell){
  //return pow(4*Pi,ell)*pow(3/(2*Pi*ell),1.5)*exp(-3/2*(get_square_diff(Rleft,Rright))/ell);
  return pow(3/(2*Pi*ell),1.5)*exp(-3/2*(get_square_diff(Rleft,Rright))/ell);
}
