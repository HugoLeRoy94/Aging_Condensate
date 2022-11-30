#include "Header.h"
using namespace std;
Dangling::Dangling() : Strand()
{}
// dsl : distance sur ell : ratio of the two usefull to know how to adjust the concentration of linkers
Dangling::Dangling(Linker* R0,
                  LoopLinkWrap& loop_link,
                  double ell_0,
                  double ell_in,
                  double rho,
                  bool rho_adjust) : Strand(R0,ell_0,rho,rho_adjust)
{
  IF(true) { cout << "Dangling : creator" << endl; }
  ell = ell_in;
  radius = 2*sqrt(ell);
  rho0 = rho;
  //V = 4 / 3 * Pi * pow(2 * ell, 1.5);
  V = pow(radius,3); // linkers are generated into a squared box
  xg = R0->r().at(0);
  yg = R0->r().at(1);
  zg = R0->r().at(2);
  IF(true) { cout << "Dangling : generate the binding sites" << endl; }
  set_p_linkers(loop_link);
  IF(true) { cout << "Dangling : compute the binding rates" << endl; }
  /*
  cout<<"free linkers"<<endl;
  for(auto& link : free_linkers){cout<<link<<endl;}
  cout<<"bounded linkers :"<<endl;
  for(auto& link : occ_linkers){cout<<link<<endl;}
  */
  compute_all_rates();
  IF(true) { cout << "Dangling : add the dangling to loop_link" << endl; }
  loop_link.Insert_Strand(this);
  IF(true) { cout << "Dangling : constructor over." << endl; }
}
Dangling::Dangling(const Dangling &dangling,
                  Linker* new_left_linker,
                  LoopLinkWrap& loop_link) : Strand(dangling,new_left_linker)
{
  IF(true) { cout << "Dangling : copy constructor" << endl; }
  radius = dangling.radius;
  set_p_linkers(loop_link);
  compute_all_rates();
}
Dangling::Dangling(const Dangling &dangling,LoopLinkWrap& loop_link) : Strand(dangling)
{
  IF(true) { cout << "Dangling : copy constructor" << endl; }
  radius = dangling.radius;
  set_p_linkers(loop_link);
  compute_all_rates();
}

// ---------------------------------------------------------------------------
//-----------------------------------accessor---------------------------------
// ---------------------------------------------------------------------------
double Dangling::get_S() const { return ell * log(4 * Pi); }
Linker* Dangling::get_Rright() const{throw out_of_range("dangling");}
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
double Dangling::Omega(double ell) const
{
  return pow(4 * Pi, ell);
}
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//________________________geometric child_dependend computation_____________
array<double, 3> Dangling::random_in_volume()
{
  /*
  bool OUT(true);
  double x(0), y(0), z(0);
  int counter(0);
  while (OUT)
  {
    counter++;
    if (counter > 100000)
    {
      throw runtime_error("cannot find a place to put a crosslinker");
    }
    uniform_real_distribution<double> distribution(-1, 1); // doubles from -1 to 1

    x = distribution(generator) * radius;
    y = distribution(generator) * radius;
    z = distribution(generator) * radius;

    if (pow(x, 2) + pow(y, 2) + pow(z, 2) <= ell)
    {
      OUT = false;
    }
  }*/
  uniform_real_distribution<double> distribution(-1, 1); // doubles from -1 to 1
  double x(distribution(generator) * radius);
  double y(distribution(generator) * radius);
  double z(distribution(generator) * radius);
  array<double, 3> res{x + xg, y + yg, z + zg};
  return res;
}

void Dangling::get_volume_limit(double& key_0_min,double& key_0_max,
                                double& key_1_min,double& key_1_max,
                                double& key_2_min,double& key_2_max) const
{
  key_0_min = xg - radius;
  key_0_max = xg + radius;

  key_1_min = yg - radius;
  key_1_max = yg + radius;

  key_2_min = zg - radius;
  key_2_max = zg + radius;
}
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
double Dangling::compute_rate(double li, Linker* rlinker)
{
  return exp(1.5*log(1.5/(Pi*li))-1.5*get_square_diff(Rleft->r(),rlinker->r())/li);
}