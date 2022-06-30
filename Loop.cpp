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
  compute_all_rates();
}
Loop::Loop(array<double,3> R0,array<double,3> R1, double ell_loop,double rho,double dsl,bool rho_adjust){
  IF(true){cout<<"Loop : creator"<<endl;}
  //generator.seed(seed);
  // the right anchoring point of the loop
  Rleft=R0;
  Rright=R1;
  ell = ell_loop;
  if(rho_adjust){rho0 = rho+(0.1-rho)*(diff(Rright,Rleft)/ell-dsl);}
  else{rho0 = rho;}
  // Compute the volume covered by the polymer
  IF(true){cout<<"Loop : compute the volume covered by the polymer"<<endl;}
  if(ell<2.){V=0.;}
  else if(diff(Rleft,Rright)<0.1*ell){
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
  compute_all_rates();
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
std::vector<std::vector<double>> Loop::get_rates() const{return rates;}
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
    vector<double> rates_ell,cum_rates_ell; // rates is just a temporary vector to append the double vector of rates, cum_rates is just the corresponding cumulative array.
    rates_ell.reserve((int)ell);
    // fill the rates vector in which all rates associated to the binding of any crosslinker
    // at any length ell
    double unbound_term(1.5*get_square_diff(Rleft,Rright)/ell); // intermediate fastener computation
    double rate_sum_l(0); // initialize the double that just gives the binding rate transition to whatever length
    for(int ELL=1;ELL<(int)ell;ELL++){
      double li(ELL);
      rates_ell.push_back(exp(1.5*log(3*ell/(2*Pi*li*(ell-li)))-log(4*Pi)-1.5*(get_square_diff(Rleft,rlink)/li+get_square_diff(rlink,Rright)/(ell-li))+unbound_term)); // push back the rate
      // ad it to the cumulative vector
      if(cum_rates_ell.size()== 0){cum_rates_ell.push_back(rates_ell.back());}
      else{cum_rates_ell.push_back(cum_rates_ell.back()+rates_ell.back());}
      rate_sum_l+=rates_ell.back();
    }
    rates.push_back(rates_ell);
    cum_rates.push_back(cum_rates_ell);
    if(sum_l_cum_rates.size() == 0){sum_l_cum_rates.push_back(rate_sum_l);}
    else{sum_l_cum_rates.push_back(sum_l_cum_rates.back()+rate_sum_l);}
  }
  // compute the total binding rate, which is the sum of the rates vector:
  IF(true){cout<<"compute the total unbinding rates"<<endl;}
  total_rates=0; // make sure it's 0 if there is no bond
  for(auto& ell_r : rates){for(auto& rate: ell_r){total_rates+=rate;}}
}
void Loop::select_link_length(double& length, array<double,3>& r_selected) const{
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////select a crosslinker////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  // check if there is a neighboring crosslinker, it shouldn't be possible
  if(r.size() == 0){throw out_of_range("No crosslinker in the vicinity");}
  IF(true){cout<<"Loop : start selecting a link and a length"<<endl;}
  IF(true){cout<<"Loop : select a r"<<endl;}
  uniform_real_distribution<double> distribution(0,sum_l_cum_rates.back());
  double pick_rate = distribution(generator);
  vector<double> copy_sum_l_cum_rates(sum_l_cum_rates);
  vector<double>::iterator rate_selec = lower_bound(copy_sum_l_cum_rates.begin(),copy_sum_l_cum_rates.end(),pick_rate);
  int rindex(distance(copy_sum_l_cum_rates.begin(),rate_selec));
  /*cout<<"np.array([";
  for(auto& it : copy_sum_l_cum_rates){cout<<it<<",";}cout<<"])";
  cout<<endl<<rindex<<endl;*/
  r_selected = r[rindex];
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////select a length/////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  IF(true){cout<<"Loop : select a length"<<endl;}
  uniform_real_distribution<double> distribution_ell(0,cum_rates[rindex].back());
  double pick_rate_ell = distribution_ell(generator);
  vector<double> copy_cum_rates_rindex(cum_rates[rindex]);
  vector<double>::iterator rate_ell_selec = lower_bound(copy_cum_rates_rindex.begin(),copy_cum_rates_rindex.end(),pick_rate_ell);
  int ell_index(distance(copy_cum_rates_rindex.begin(),rate_ell_selec)+1);
  length=(double)ell_index;
  /*cout<<"np.array([";
  for(auto& it :copy_cum_rates_rindex){cout<<it<<",";}cout<<"])";
  cout<<endl<<ell_index<<endl;*/
  IF(true){cout<<"Loop : access the associated rate"<<endl;}

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
  return pow(3/(2*Pi*ell),1.5)*exp(-3/2*get_square_diff(Rleft,Rright)/ell);
}

array<double,3> Loop::random_in_ellipse(double a,double b,double c,double xg,double yg,double zg){
  /*
    draw three random number that  must lie within an ellipse of
    revolution of big axe a and two small axes b.
  */
  // draw a x,y,z x in [-a,a] and y,z in [-b,b]
  draw:
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
  if(res[0]>anchor[0] | res[0]<0.){goto draw;}
  return res;
}
int compare(const pair<int,double> & t1, const pair<int,double>& t2){
  return (t1.second<t2.second);
}
