#include "Header.h"
using namespace std;
Dangling::Dangling() : Strand()
{}
// dsl : distance sur ell : ratio of the two usefull to know how to adjust the concentration of linkers
Dangling::Dangling(Linker* R0,
                  double ell_0,
                  double ell_in,
                  double rho,
                  bool sliding) : Strand(R0,ell_0,rho,sliding)
{
  IF(true) { cout << "Dangling : creator" << endl; }
  ell = ell_in;
  radius = sqrt(ell/2.);
  rho0 = rho;
  //V = 4 / 3 * Pi * pow(2 * ell, 1.5);
  V = 4./3.*Pi*pow(radius,3); // linkers are generated into a squared box
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
double Dangling::get_S(double dl) const { return (ell+dl) * log(4 * Pi); }
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

void Dangling::get_volume_limit(array<double,3>& main_ax, 
                                array<double,3>& ctr_mass,
                                double& a, double& b) const
{
  ctr_mass={xg,yg,zg};
  a=radius;
  b=radius;
  main_ax = {1.,1.,1.}; // the main ax does not matter when a = b
}
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
double Dangling::compute_binding_rate(double li, Linker* rlinker)const
{
  if(diff(Rleft->r(), rlinker->r())>li)
  {return 0.;}
  return exp(1.5*log(1.5/(Pi*li))-1.5*get_square_diff(Rleft->r(),rlinker->r())/li)/ell;
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
  //linker_selected->set_bounded();
  unique_ptr<Loop> left_loop = make_unique<Loop>(Rleft,linker_selected,ell_coordinate_0,ell_coordinate_0+length,rho0,slide);
  unique_ptr<Dangling> dangling = make_unique<Dangling>(linker_selected,ell_coordinate_0+length,ell-length,rho0,slide);
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
                               rho0,
                               slide);
}
unique_ptr<Strand> Dangling::do_slide(double dl,bool right) const
{
  return make_unique<Dangling>(Rleft,ell_coordinate_0+dl,ell-dl,rho0,slide);
}
double Dangling::get_ell_coordinate_1() const{return ell_coordinate_0+ell;}