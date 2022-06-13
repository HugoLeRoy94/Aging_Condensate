#include "Header.h"
using namespace std;


Loop::Loop(double ell_loop,double rho){
  // a looop with a single crosslink, will be defined as a normal loop
  // with ell loop = 2x a normal loop and both anchorring point at the same
  // place
  IF(true){cout<<"Loop : creator with a single crosslink"<<endl;}
  //generator.seed(seed);

  rho0 = rho;
  Rright = {0.,0.,0.};
  Rleft={0.,0.,0.};
  // Rright is not initialized
  ell = ell_loop;
  IF(true){cout<<"Loop : compute the volume covered by the polymer"<<endl;}
  V =4/3*Pi*pow(ell,1.5);
  a=sqrt(ell)*0.5;
  b=sqrt(ell)*0.5;
  // create r by randmoly generating positions
  IF(true){cout<<"Loop : generate the binding sites"<<endl;}
  generate_binding_sites();
  // compute the array of binding rate to any crosslinker at any length
  IF(true){cout<<"loop : compute the unbinding rates for every ell"<<endl;}
  for(auto& rlink : r){
    vector<double> rates_ell;
    rates_ell.reserve((int)ell);
    for(int ELL=1;ELL<(int)ell;ELL++){
      rates_ell.push_back(Omega(Rleft,rlink,(double)ELL)*Omega(rlink,Rright,ell-(double)ELL));
    }
    rates.push_back(rates_ell);
  }
  // compute the total binding rate:
  IF(true){cout<<"compute the total unbinding rates"<<endl;}
  total_rates=0; // make sure it's 0 if there is no bond
  for(auto& ell_r : rates){
    for(auto& rate: ell_r){
      //cout<<rate<<" ";
      total_rates+=rate;
    }
    //cout<<endl;
  }
}
Loop::Loop(array<double,3> R0,array<double,3> R1, double ell_loop,double rho){
  IF(true){cout<<"Loop : creator"<<endl;}
  //generator.seed(seed);
  // the right anchoring point of the loop
  rho0 = rho;
  Rleft=R0;
  Rright=R1;
  ell = ell_loop;
  // Compute the volume covered by the polymer
  IF(true){cout<<"Loop : compute the volume covered by the polymer"<<endl;}
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
  IF(true){cout<<"Loop : generate the binding sites"<<endl;}
  generate_binding_sites();
  // compute the array of binding rate to any crosslinker at any length
  IF(true){cout<<"loop : compute the unbinding rates for every ell"<<endl;}

}
Loop::Loop(const Loop& loop){
  Rleft = loop.Rleft;
  Rright = loop.Rright;
  rho0 = loop.rho0;
  ell = loop.ell;
  V = loop.V;
  a = loop.a;
  b = loop.b;
  generate_binding_sites();
  compute_all_rates();
}
Loop::~Loop(){}
bool Loop::operator<(const Loop& loop_right) const{
  return Rright[0]<loop_right.Rright[0];
}
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
void Loop::compute_all_rates(){
  for(auto& rlink : r){
    vector<double> rates_ell;
    rates_ell.reserve((int)ell);
    for(int ELL=1;ELL<(int)ell;ELL++){
      rates_ell.push_back(Omega(Rleft,rlink,(double)ELL)*Omega(rlink,Rright,ell-(double)ELL));
    }
    rates.push_back(rates_ell);
  }
  // compute the total binding rate:
  IF(true){cout<<"compute the total unbinding rates"<<endl;}
  total_rates=0; // make sure it's 0 if there is no bond
  for(auto& ell_r : rates){
    for(auto& rate: ell_r){
      //cout<<rate<<" ";
      total_rates+=rate;
    }
    //cout<<endl;
  }
};
double Loop::compute_binding_rate(int r_index, double ell){

  return 0;
}
void Loop::select_link_length(double& length, array<double,3>& r_selected,double& time) const{
  // check if there is a neighboring crosslinker, it shouldn't be possible
  if(r.size() == 0){throw out_of_range("No crosslinker in the vicinity");}
  IF(true){cout<<"Loop : start selecting a link and a length"<<endl;}
  // let select a link and a length uniformly for now
  IF(true){cout<<"Loop : select a r"<<endl;}
  uniform_int_distribution<int> r_distrib(0,r.size()-1);
  int index(r_distrib(generator));
  //cout<<r.size()<<" "<<index<<endl;
  r_selected =r[index];

  IF(true){cout<<"Loop : select a length"<<endl;}
  //cout<<"find a length between"<<diff(Rleft,r_selected)<<"and "<<ell-diff(r_selected,Rright)<<endl;
  uniform_real_distribution<double> ell_distrib(diff(Rleft,r_selected),ell-diff(r_selected,Rright));
  length =ell_distrib(generator);
  IF(true){cout<<"Loop : compute the associated time"<<endl;}
  time = 1/(Omega(Rleft,r_selected,length)*Omega(r_selected,Rright,ell-length));
  time=1.;

}
void Loop::generate_binding_sites(){
  poisson_distribution<int> distribution(rho0*V);
  int N_crosslinker=distribution(generator);
  for(int n=0;n<N_crosslinker;n++){
    r.push_back(random_in_ellipse(a,b,b,0.5*(Rright[0]+Rleft[0]),0.5*(Rright[1]+Rleft[1]),0.5*(Rright[2]+Rleft[2])));
  }
  //cout<<"size of r "<<r.size()<<endl;
  //for(auto& it : r){cout<<it[0]<<" "<<it[1]<<" "<<it[2]<<endl;}
}
double Loop::Omega(array<double,3> r1,array<double,3> r2,double ell) const{
  //return pow(4*Pi,ell)*pow(3/(2*Pi*ell),1.5)*exp(-3/2*(get_square_diff(Rleft,Rright))/ell);
  return pow(3/(2*Pi*ell),1.5)*exp(-3/2*(get_square_diff(Rleft,Rright))/ell);
}

array<double,3> Loop::random_in_ellipse(double a,double b,double c,double xg,double yg,double zg){
  /*
    draw three random number that  must lie within an ellipse of
    revolution of big axe a and two small axes b.
  */
  // draw a x,y,z x in [-a,a] and y,z in [-b,b]
  bool OUT(true);
  double x(0),y(0),z(0);
  int counter(0);
  while(OUT){
    counter++;
    if(counter>100000){throw runtime_error("cannot find a place to put a crosslinker");}
    uniform_real_distribution<double> distribution(-1,1); //doubles from -1 to 1

    x = distribution(generator)*a;
    y = distribution(generator)*b;
    z = distribution(generator)*c;

    /*x = (double)xint/(double)randmax*a;
    y = (double)yint/randmax*b;
    z = (double)zint/randmax*c;*/
    //cout<<"x,y,z"<<x<<","<<y<<","<<z<<endl;
    //cout<<pow(x/a,2)+pow(y/b,2)+pow(z/c,2)<<endl;
    //cout<<a<<" "<<b<<" "<<c<<endl;
    if(pow(x/a,2)+pow(y/b,2)+pow(z/c,2)<=1){OUT=false;}
  }
  double theta(atan2(yg,xg)),phi(atan2(xg,zg)-Pi/2.);
  //cout<<"xg yg zg ="<<xg<<" "<<yg<<" "<<zg<<endl;
  //cout<< "theta phi ="<<theta<<" "<<phi<<endl;
  //cout<<x<<" "<<y<<" "<<z<<endl;
  //cout<<x<<" "<<y<<" "<<z<<endl;
  array<double,3> res{cos(phi)*(cos(theta)*x-sin(theta)*y+sin(phi)*z),
                      sin(theta)*x+cos(theta)*y,
                      -sin(phi)*(cos(theta)*x-sin(theta)*y)+cos(phi)*z};
  res[0]+=xg;
  res[1]+=yg;
  res[2]+=zg;
  //cout<<res[0]<<" "<<res[1]<<" "<<res[2]<<endl;
  return res;
}
int compare(const pair<int,double> & t1, const pair<int,double>& t2){
  return (t1.second<t2.second);
}
