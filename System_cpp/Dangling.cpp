#include "Header.h"
using namespace std;
Dangling::Dangling() : Strand()
{}
// dsl : distance sur ell : ratio of the two usefull to know how to adjust the concentration of linkers
Dangling::Dangling(Linker* R0,
                  double ell_0,
                  double ell_in,
                  double rho) : Strand(R0,ell_0,rho)
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
  IF(true) { cout << "Dangling : constructor over." << endl; }
}
Dangling::Dangling(const Dangling &dangling,
                  Linker* new_left_linker) : Strand(dangling,new_left_linker)
{
  IF(true) { cout << "Dangling : copy constructor" << endl; }
  radius = dangling.radius;
}

Dangling::Dangling(const Dangling &dangling) : Strand(dangling)
{
  IF(true) { cout << "Dangling : copy constructor" << endl; }
  radius = dangling.radius;
}

Strand* Dangling::clone() const{return new Dangling(*this);}
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

pair<unique_ptr<Strand>,unique_ptr<Strand>> Dangling::bind() const
{
  // return a reference to a the two loop that must be created
  // when binding the current loop to a linker randomly choosen
  // in the vicinity of the current loop. The choice is made
  // in accordance with the binding rate of each specific linker
  Linker* linker_selected;
  double length;
  select_link_length(length,linker_selected);
  linker_selected->set_bounded();
  unique_ptr<Loop> left_loop = make_unique<Loop>(Rleft,linker_selected,ell_coordinate_0,ell_coordinate_0+length,rho0);
  unique_ptr<Dangling> dangling = make_unique<Dangling>(linker_selected,ell_coordinate_0+length,ell-length,rho0);
  return  {move(left_loop),move(dangling)};
}
unique_ptr<Strand> Dangling::unbind_from(Strand* left_strand) const
{
  // return a reference to the loop that must be created when
  // unbinding a linker between the current loop (at its right)
  // and "left_loop" that is at its left.
  return make_unique<Dangling>(left_strand->get_Rleft(),
                               left_strand->get_ell_coordinate_0(),
                               ell + left_strand->get_ell(),
                               rho0);
}