#include "Header.h"
using namespace std;

Strand::Strand(){}
Strand::Strand(Linker* R0,
            double ell_0,
            double rho,
            bool sliding)
{
  IF(true) { cout << "Strand : creator" << endl; }
  // generator.seed(seed);
  //  the right anchoring point of the strand
  Rleft = R0;
  ell_coordinate_0 = ell_0;
  //if (rho_adjust){rho0 = rho + (0.1 - rho) * (diff(Rright, Rleft) / ell - dsl);}
  rho0 = rho;
  slide= sliding;
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
  slide = strand.slide;
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
  slide = strand.slide;
}

void Strand::set_linkers(vector<Linker*> new_free_linkers,vector<Linker*> occupied_linkers)
{
  free_linkers = new_free_linkers;
  occ_linkers=occupied_linkers;
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

void Strand::compute_total_rates()
{
  total_rates = 0.;
  for (auto &rlink : free_linkers) // iterate through every linkers
  {
    for (int ELL = 1; ELL < (int)ell; ELL++)
    {
      total_rates+=compute_binding_rate((double)ELL,rlink);
    }
  }
}

void Strand::compute_cum_rates(vector<double>& sum_l_cum_rates,
                               vector<vector<double>>& cum_rates) const
{
//reserve the correct amount of memory
sum_l_cum_rates.resize(free_linkers.size());
cum_rates.resize(free_linkers.size());
for(auto& cum_rate : cum_rates){cum_rate.resize((int)ell-1);}
// the entry vectors must be empty
int rindex(0);
for (auto &rlink : free_linkers)
  {
    // iterate through each linker
    // and compute a cumulative binding rate vector for each
    // length and.
    int ellindex(0);
    for (int ELL = 1; ELL < (int)ell; ELL++)
    {
      // ad it to the cumulative vector
      if (ellindex == 0)
      {
        cum_rates[rindex][ellindex] = compute_binding_rate((double)ELL,rlink);
      }
      else
      {
        cum_rates[rindex][ellindex] = cum_rates[rindex][ellindex-1]+compute_binding_rate((double)ELL,rlink);
      }
      ellindex++;
    }
    // Add the end of this vector, which is the total probability
    // to bind to this specific linker to sum_l_cum_rates.
    if(rindex==0){sum_l_cum_rates[0] = cum_rates[rindex].back();}
    else{sum_l_cum_rates[rindex] = sum_l_cum_rates[rindex-1]+cum_rates[rindex].back();}
    // add the whole vector to to cum_rates:
    rindex++;
  }
  //cout<<"cum_rate"<<endl;
  //for(auto& rate : cum_rates){for(auto& r : rate){cout<<r<<endl;}}

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
  vector<double> sum_l_cum_rates;
  vector<vector<double>> cum_rates;
  compute_cum_rates(sum_l_cum_rates,cum_rates);
  //Check_integrity();
  IF(true) { cout << "Strand : select a r" << endl; }
  uniform_real_distribution<double> distribution(0, sum_l_cum_rates.back());
  double pick_rate = distribution(generator);
  vector<double>::iterator rate_selec = lower_bound(sum_l_cum_rates.begin(), sum_l_cum_rates.end(), pick_rate);
  int rindex(distance(sum_l_cum_rates.begin(), rate_selec));
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
  IF(true){cout<<"ell selected : "<<length<<endl;}
  IF(true){cout<<"r selected : "<<r_selected->r()[0]<<" "<<r_selected->r()[1]<<" "<<r_selected->r()[2]<<endl;}
}

void Strand::Check_integrity() const
{
  cout<<ell_coordinate_0<<endl;
}

// -----------------------------------------------------------------------------
/*
 /\  _ _ _  _ _ _  _
/~~\(_(_(/__\_\(_)|                     
*/
// -----------------------------------------------------------------------------
Linker* Strand::get_Rleft() const { return Rleft; }

std::vector<Linker*> Strand::get_r() const { return free_linkers; }

std::vector<Linker*> Strand::get_occ_r() const { return occ_linkers; }

double Strand::get_ell() const { return ell; }

double Strand::get_ell_coordinate_0() const { return ell_coordinate_0; }

double Strand::get_V() const { return V; }

double Strand::get_total_binding_rates() const { return total_rates; }