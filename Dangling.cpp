#include "Header.h"
using namespace std;

// dsl : distance sur ell : ratio of the two usefull to know how to adjust the concentration of linkers
Dangling::Dangling(std::array<double,3> R0, double ell, double rho,double dsl,bool rho_adjust){
  IF(true){cout<<"Dangling : creator"<<endl;}
  Rleft = R0;
  radius = sqrt(ell);
  if(rho_adjust){rho0 = rho+(0.1-rho)*(radius/ell-dsl);}
  else{rho0 = rho;}
  V = 4/3*Pi*pow(ell,1.5);
  IF(true){cout<<"Dangling : generate the binding sites"<<endl;}
  generate_binding_sites();
  IF(true){cout<<"Dangling : compute the binding rates" <<endl;}
  compute_all_rates();
}
Dangling::Dangling(const Dangling& dangling){
  Rleft = dangling.Rleft;
  rho0 = dangling.rho0;
  ell = dangling.ell;
  V=dangling.V;
  radius = dangling.radius;
  generate_binding_sites();
  compute_all_rates();
}
Dangling::~Dangling(){}
// ---------------------------------------------------------------------------
//-----------------------------------accessor---------------------------------
// ---------------------------------------------------------------------------
std::array<double,3> Dangling::get_Rleft() const{return Rleft;}
std::vector<array<double,3>> Dangling::get_r() const{return r;}
std::vector<std::vector<double>> Dangling::get_rates() const{return rates;}
double Dangling::get_ell() const{return ell;}
double Dangling::get_V()const{return V;}
double Dangling::get_total_binding_rates() const{return total_rates;}
double Dangling::get_S()const{return pow((4*Pi),ell);}
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
void Dangling::compute_all_rates(){
  for(auto& rlink : r){
    vector<double> rates_ell,cum_rates_ell; // rates is just a temporary vector to append the double vector of rates, cum_rates is just the corresponding cumulative array.
    rates_ell.reserve((int)ell);
    // fill the rates vector in which all rates associated to the binding of any crosslinker
    // at any length ell
    //double unbound_term(1.5*get_square_diff(Rleft,Rright)/ell); // intermediate fastener computation
    double rate_sum_l(0); // initialize the double that just gives the binding rate transition to whatever length
    for(int ELL=1;ELL<(int)ell;ELL++){
      double li(ELL);
      //rates_ell.push_back(exp(
      //  1.5*log(3*ell/(2*Pi*li*(ell-li)))-log(4*Pi)-1.5*(get_square_diff(Rleft,rlink)/li+get_square_diff(rlink,Rright)/(ell-li))+unbound_term)); // push back the rate
      rates_ell.push_back(exp(1.5*log(1.5/(Pi*li))-1.5*get_square_diff(Rleft,rlink)/li));
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
  IF(true){cout<<"Dangling : compute the total unbinding rates"<<endl;}
  total_rates=0; // make sure it's 0 if there is no bond
  for(auto& ell_r : rates){for(auto& rate: ell_r){total_rates+=rate;}}
}
void Dangling::select_link_length(double& length, array<double,3>& r_selected) const{
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////select a crosslinker////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  // check if there is a neighboring crosslinker, it shouldn't be possible
  if(r.size() == 0){throw out_of_range("No crosslinker in the vicinity");}
  IF(true){cout<<"Dangling : start selecting a link and a length"<<endl;}
  IF(true){cout<<"Dangling : select a r"<<endl;}
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
  IF(true){cout<<"Dangling : select a length"<<endl;}
  uniform_real_distribution<double> distribution_ell(0,cum_rates[rindex].back());
  double pick_rate_ell = distribution_ell(generator);
  vector<double> copy_cum_rates_rindex(cum_rates[rindex]);
  vector<double>::iterator rate_ell_selec = lower_bound(copy_cum_rates_rindex.begin(),copy_cum_rates_rindex.end(),pick_rate_ell);
  int ell_index(distance(copy_cum_rates_rindex.begin(),rate_ell_selec)+1);
  length=(double)ell_index;
}
void Dangling::generate_binding_sites(){
  poisson_distribution<int> distribution(rho0*V);
  int N_crosslinker=distribution(generator);
  for(int n=0;n<N_crosslinker;n++){
    r.push_back(random_in_sphere(Rleft[0], Rleft[1],Rleft[2]));
  }
}
double Dangling::Omega(double ell) const{
  return pow(4*Pi,ell);
}
array<double,3> Dangling::random_in_sphere(double xg,double yg,double zg){
  bool OUT(true);
  double x(0),y(0),z(0);
  int counter(0);
  while(OUT){
    counter++;
    if(counter>100000){throw runtime_error("cannot find a place to put a crosslinker");}
    uniform_real_distribution<double> distribution(-1,1); //doubles from -1 to 1

    x = distribution(generator)*radius;
    y = distribution(generator)*radius;
    z = distribution(generator)*radius;

    if(pow(x,2)+pow(y,2)+pow(z,2)<=radius){OUT=false;}
  }
  array<double,3> res{x,y,z};
  return res;
}
