#include "Header.h"
using namespace std;
Dangling::Dangling() : Strand()
{}
// dsl : distance sur ell : ratio of the two usefull to know how to adjust the concentration of linkers
Dangling::Dangling(std::array<double, 3> R0,
                  map3d<double,double,double,array<double,3>>& linkers,
                  double ell_0,
                  double ell_in,
                  double rho,
                  bool rho_adjust) : Strand(R0,linkers,ell_0,rho,rho_adjust)
{
  IF(true) { cout << "Dangling : creator" << endl; }
  ell = ell_in;
  radius = 2*sqrt(ell);
  rho0 = rho;
  //V = 4 / 3 * Pi * pow(2 * ell, 1.5);
  V = pow(radius,3); // linkers are generated into a squared box
  xg = R0[0];
  yg = R0[1];
  zg = R0[2];
  IF(true) { cout << "Dangling : generate the binding sites" << endl; }
  generate_binding_sites(linkers);
  IF(true) { cout << "Dangling : compute the binding rates" << endl; }
  compute_all_rates();
  IF(true) { cout << "Dangling : constructor over." << endl; }
}

Dangling::Dangling(const Dangling &dangling,map3d<double,double,double,array<double,3>>& linkers) : Strand(dangling, linkers)
{
  IF(true) { cout << "Dangling : copy constructor" << endl; }
  V = dangling.V;
  radius = dangling.radius;
  generate_binding_sites(linkers);
  compute_all_rates();
}

// ---------------------------------------------------------------------------
//-----------------------------------accessor---------------------------------
// ---------------------------------------------------------------------------
double Dangling::get_S() const { return ell * log(4 * Pi); }
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
double Dangling::compute_rate(double li, array<double,3>* rlinker)
{
  return exp(1.5*log(1.5/(Pi*li))-1.5*get_square_diff(Rleft,*rlinker)/li);
}