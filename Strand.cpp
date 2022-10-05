#include "Header.h"
using namespace std;

Strand::Strand(){}
Strand::Strand(array<double, 3> R0,
            map3d<double,double,double,array<double,3>>& linkers,
            std::map<array<double,3>*,vector<Strand*>> linker_to_strand,
            double ell_0,
            double rho,
            bool rho_adjust)
{
  IF(true) { cout << "Strand : creator" << endl; }
  // generator.seed(seed);
  //  the right anchoring point of the strand
  Rleft = R0;
  ell_coordinate_0 = ell_0;
  //if (rho_adjust){rho0 = rho + (0.1 - rho) * (diff(Rright, Rleft) / ell - dsl);}
  rho0 = rho;
}

Strand::~Strand()
{
    //for(auto& p : p_linkers){delete p;} 
}

Strand::Strand(const Strand& strand,map3d<double,double,double,array<double,3>>& linkers,
              std::map<array<double,3>*,vector<Strand*>> linker_to_strand)
{
    Rleft = strand.Rleft;
    rho0 = strand.rho0;
    ell_coordinate_0 = strand.ell_coordinate_0;
    ell = strand.ell;
    V = strand.V;
}

bool Strand::operator<(const Strand &strand_right) const
{
    return ell_coordinate_0 < strand_right.ell_coordinate_0;
}

std::array<double, 3> Strand::get_Rleft() const { return Rleft; }

std::vector<array<double, 3>*> Strand::get_r() const { return p_linkers; }

std::vector<std::vector<double>> Strand::get_rates() const { return rates; }

double Strand::get_ell() const { return ell; }

double Strand::get_ell_coordinate_0() const { return ell_coordinate_0; }

double Strand::get_V() const { return V; }

double Strand::get_total_binding_rates() const { return total_rates; }

void Strand::compute_all_rates()
{
  IF(true) { cout << "Strand : compute all the unbinding rates" << endl; }
  for (auto &rlink : p_linkers)
  {
    vector<double> rates_ell, cum_rates_ell; // rates is just a temporary vector to append the double vector of rates, cum_rates is just the corresponding cumulative array.
    rates_ell.reserve((int)ell);
    // fill the rates vector in which all rates associated to the binding of any crosslinker
    // at any length ell
    double rate_sum_l(0);                                            // initialize the double that just gives the binding rate transition to whatever length
    for (int ELL = 1; ELL < (int)ell; ELL++)
    {
      double li(ELL);
      // rates_ell.push_back(exp(1.5*log(3*ell/(2*Pi*li*(ell-li)))-log(4*Pi)-1.5*(get_square_diff(Rleft,rlink)/li+get_square_diff(rlink,Rright)/(ell-li))+unbound_term)); // push back the rate
      rates_ell.push_back(compute_rate(li,rlink)); // push back the rate
      // ad it to the cumulative vector
      if (cum_rates_ell.size() == 0)
      {
        cum_rates_ell.push_back(rates_ell.back());
      }
      else
      {
        cum_rates_ell.push_back(cum_rates_ell.back() + rates_ell.back());
      }
      rate_sum_l += rates_ell.back();
    }
    rates.push_back(rates_ell);
    cum_rates.push_back(cum_rates_ell);
    if (sum_l_cum_rates.size() == 0)
    {
      sum_l_cum_rates.push_back(rate_sum_l);
    }
    else
    {
      sum_l_cum_rates.push_back(sum_l_cum_rates.back() + rate_sum_l);
    }
  }
  // compute the total binding rate, which is the sum of the rates vector:
  IF(true) { cout << "Strand : compute the total unbinding rate" << endl; }
  total_rates = 0; // make sure it's 0 if there is no bond
  for (auto &ell_r : rates)
  {
    for (auto &rate : ell_r)
    {
      total_rates += rate;
    }
  }
}

void Strand::select_link_length(double &length, array<double, 3>& r_selected) const
{
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////select a crosslinker////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  // check if there is a neighboring crosslinker, it shouldn't be possible
  if (p_linkers.size() == 0)
  {
    throw out_of_range("No crosslinker in the vicinity");
  }
  IF(true) { cout << "Strand : start selecting a link and a length" << endl; }
  IF(true) { cout << "Strand : select a r" << endl; }
  uniform_real_distribution<double> distribution(0, sum_l_cum_rates.back());
  double pick_rate = distribution(generator);
  vector<double> copy_sum_l_cum_rates(sum_l_cum_rates);
  vector<double>::iterator rate_selec = lower_bound(copy_sum_l_cum_rates.begin(), copy_sum_l_cum_rates.end(), pick_rate);
  int rindex(distance(copy_sum_l_cum_rates.begin(), rate_selec));
  /*cout<<"np.array([";
  for(auto& it : copy_sum_l_cum_rates){cout<<it<<",";}cout<<"])";
  cout<<endl<<rindex<<endl;*/
  r_selected = *p_linkers[rindex];
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////select a length/////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  IF(true) { cout << "Strand : select a length" << endl; }
  uniform_real_distribution<double> distribution_ell(0, cum_rates[rindex].back());
  double pick_rate_ell = distribution_ell(generator);
  vector<double> copy_cum_rates_rindex(cum_rates[rindex]);
  vector<double>::iterator rate_ell_selec = lower_bound(copy_cum_rates_rindex.begin(), copy_cum_rates_rindex.end(), pick_rate_ell);
  int ell_index(distance(copy_cum_rates_rindex.begin(), rate_ell_selec) + 1);
  length = (double)ell_index;
  /*cout<<"np.array([";
  for(auto& it :copy_cum_rates_rindex){cout<<it<<",";}cout<<"])";
  cout<<endl<<ell_index<<endl;*/
}

double Strand::compute_rate(double li, array<double,3>* rlinker){return 0;} // dummy function that is overwritten in the child class

array<double,3> Strand::random_in_volume(){return {0,0,0};} // dummy function that is overwritten in the child class

void Strand::get_volume_limit(double& key_0_min,double& key_0_max,
                              double& key_1_min,double& key_1_max,
                              double& key_2_min,double& key_2_max) const{}

void Strand::reset_p_linkers(map3d<double,double,double,array<double,3>>& linkers)
{
  throw bad_function_call();
  for(auto& p : p_linkers){delete p;}
  p_linkers.clear();
  double key_0_min,key_0_max,key_1_min,key_1_max,key_2_min,key_2_max;
  get_volume_limit(key_0_min,key_0_max,key_1_min,key_1_max,key_2_min,key_2_max);
  linkers.cut_slice(key_0_min,key_0_max,key_1_min,key_1_max,key_2_min,key_2_max);
  generate_binding_sites(linkers);
}

void Strand::generate_binding_sites(map3d<double,double,double,array<double,3>>& linkers)
{
  p_linkers.clear();
  //---------------------------------------------------------------
  //-----------------draw a number of crosslinkers ----------------
  //---------------------------------------------------------------
  //poisson_distribution<int> distribution(rho0 * V);
  //int N_crosslinker = distribution(generator);
  //---------------------------------------------------------------
  //_________________store the crosslinkers that are already ______
  //_____________________in the vicinity___________________________
  //---------------------------------------------------------------
  double key_0_min,key_0_max,key_1_min,key_1_max,key_2_min,key_2_max;
  get_volume_limit(key_0_min,key_0_max,key_1_min,key_1_max,key_2_min,key_2_max);
  std::vector<std::array<double, 3>*> in_vicinity(linkers.slice(key_0_min,key_0_max,key_1_min,key_1_max,key_2_min,key_2_max));
  p_linkers.insert(p_linkers.end(),in_vicinity.begin(),in_vicinity.end());
  //---------------------------------------------------------------
  //_____________if there are less crosslinkers than needed :______
  //---------------------------------------------------------------
  /*
  if(N_crosslinker>p_linkers.size())
  {
    // add new ones
    for(int n=0;n<N_crosslinker-p_linkers.size();n++)
    {
      array<double,3> new_linker(random_in_volume());
      //p_linkers.push_back(&new_linker);
      //linkers.add(new_linker[0],new_linker[1],new_linker[2],new_linker);
      p_linkers.push_back(linkers.add_return_address(
                            new_linker[0],
                            new_linker[1],
                            new_linker[2],
                            new_linker));
      //cout<<linkers(new_linker[0],new_linker[1],new_linker[2])[0]<<endl;
    }
  }
  else
  {
    // remove some
    for(int n=0;n<p_linkers.size()-N_crosslinker;n++)
    {
      uniform_int_distribution<int> rand_int(0, p_linkers.size()-1);
      vector<array<double,3>*>::iterator linker_select(p_linkers.begin()+rand_int(generator));
      linkers.remove((*linker_select)->at(0),(*linker_select)->at(1),(*linker_select)->at(2));
      p_linkers.erase(linker_select);
    }
  }
  */
}