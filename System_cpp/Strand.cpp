#include "Header.h"
using namespace std;

Strand::Strand(){}
Strand::Strand(Linker* R0,
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
  IF(true){cout<<"Strand destructor "<<this<<endl;}

}

Strand::Strand(const Strand& strand, Linker* new_left_linker)
{
  IF(true){cout<<"Strand : Copy constructor with left_linker"<<endl;}
  xg = strand.xg;
  yg = strand.yg;
  zg = strand.zg;
  Rleft = new_left_linker;
  rho0 = strand.rho0;
  ell_coordinate_0 = strand.ell_coordinate_0;
  ell = strand.ell;
  V = strand.V;
}

Strand::Strand(const Strand& strand)
{
  IF(true){cout<<"Strand : Copy constructor"<<endl;}
  xg = strand.xg;
  yg = strand.yg;
  zg = strand.zg;
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

void Strand::remove_from_linkers()
{
    for(auto& linker : free_linkers)
    {
      linker->remove_strand(this);
    }
    for(auto& linker : occ_linkers)
    {
      linker->remove_strand(this);
    }
}

void Strand::compute_all_rates()
{
  /*Not needed, because strands never changes !
  IF(true){cout<<"Strand : clean the rates"<<endl;}
  rates.clear();
  cum_rates.clear();
  sum_l_cum_rates.clear();*/
  IF(true) { cout << "Strand : compute all the binding rates" << endl; }
  for (auto &rlink : free_linkers)
  {
    vector<double> rates_ell, cum_rates_ell; // rates is just a temporary vector to append the double vector of rates, cum_rates is just the corresponding cumulative array.
    try{rates_ell.reserve((int)ell);}
    catch(length_error){cout<<(int)ell<<" "<<this<<" "<<Rleft->r()[0]<<" "<<Rleft->r()[1]<<" "<<Rleft->r()[2]<<endl;throw length_error("ell does not have a correct value");}
    // fill the rates vector in which all rates associated to the binding of any crosslinker
    // at any length ell
    double rate_sum_l(0);                                            // initialize the double that just gives the binding rate transition to whatever length
    //IF(true){cout<<"Compute rate for each Ell"<<endl;}
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
    //IF(true){cout<<"compute the comulative array"<<endl;}
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

void Strand::select_link_length(double &length, Linker*& r_selected) const
{
  
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////select a crosslinker////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  // check if there is a neighboring crosslinker, it shouldn't be possible
  if (free_linkers.size() == 0)
  {
    throw out_of_range("No crosslinker in the vicinity");
  }
  IF(true) { cout << "Strand : start selecting a link and a length" << endl; }
  //Check_integrity();
  IF(true) { cout << "Strand : select a r" << endl; }
  uniform_real_distribution<double> distribution(0, sum_l_cum_rates.back());
  double pick_rate = distribution(generator);
  vector<double> copy_sum_l_cum_rates(sum_l_cum_rates);
  vector<double>::iterator rate_selec = lower_bound(copy_sum_l_cum_rates.begin(), copy_sum_l_cum_rates.end(), pick_rate);
  int rindex(distance(copy_sum_l_cum_rates.begin(), rate_selec));
  /*cout<<"np.array([";
  for(auto& it : copy_sum_l_cum_rates){cout<<it<<",";}cout<<"])";
  cout<<endl<<rindex<<endl;*/
  r_selected = free_linkers[rindex];
  //cout<<"rindex = "<<rindex<<" p_linkers size = "<<p_linkers.size()<<" "<<p_linkers[rindex]<<" "<<r_selected<<endl;
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
  /*cout<<(*r_selected)[0]<<" "<<(*r_selected)[1]<<" "<<(*r_selected)[2]<<endl;
  cout<<r_selected<<endl;*/
  /*cout<<"np.array([";
  for(auto& it :copy_cum_rates_rindex){cout<<it<<",";}cout<<"])";
  cout<<endl<<ell_index<<endl;*/
}

void Strand::set_p_linkers(LoopLinkWrap& loop_link)
{
  IF(true){cout<<"set_p_linkers"<<endl;}
  // get the volume limit
  double key_0_min,key_0_max,key_1_min,key_1_max,key_2_min,key_2_max;
  get_volume_limit(key_0_min,key_0_max,key_1_min,key_1_max,key_2_min,key_2_max);
  // select the linkers in the vicinity
  IF(true){cout<<"Strand : slice"<<endl;}
  loop_link.slice_free(key_0_min,key_0_max,
                      key_1_min,key_1_max,
                      key_2_min,key_2_max,
                      free_linkers,occ_linkers);
  // tell the linkers that this strand is around
  for(auto& linker : free_linkers){linker->add_strand(this);}
  for(auto& linker : occ_linkers){linker->add_strand(this);}
}

void Strand::Check_integrity() const
{
  cout<<ell_coordinate_0<<endl;
}

// -----------------------------------------------------------------------------
// -----------------------------accessor----------------------------------------
// -----------------------------------------------------------------------------
Linker* Strand::get_Rleft() const { return Rleft; }

std::vector<Linker*> Strand::get_r() const { return free_linkers; }

std::vector<Linker*> Strand::get_occ_r() const { return occ_linkers; }

std::vector<std::vector<double>> Strand::get_rates() const { return rates; }

double Strand::get_ell() const { return ell; }

double Strand::get_ell_coordinate_0() const { return ell_coordinate_0; }

double Strand::get_V() const { return V; }

double Strand::get_total_binding_rates() const { return total_rates; }