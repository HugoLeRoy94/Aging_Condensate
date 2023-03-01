#include "Header.h"
using namespace std;
int Linker::counter = 0.;
Loop::Loop(Linker* R0,
           Linker* R1,
           double ell_0,
           double ell_1,
           double rho,
           bool sliding) : Strand(R0,ell_0,rho,sliding)
{
  IF(true) { cout << "Loop : creator" << endl; }
  Rright = R1;
  ell_coordinate_1 = ell_1;
  ell = ell_coordinate_1 - ell_coordinate_0;
    // Compute the volume covered by the polymer
  IF(true) { cout << "Loop : compute the volume covered by the polymer" << endl; }
  a=0;
  b=0;
  V=0;
  if (ell < 2.)
  {
    V = 0.;
  }
  
  //else if (diff(Rleft->r(), Rright->r()) < 0.1 * ell)
  else if(diff(Rleft->r(), Rright->r())/2.<sqrt(ell/2.))
  {
    //V = 4 / 3 * Pi * pow(2*ell, 1.5);
    a = sqrt(ell/2.);
    b = sqrt(ell/2.);
    V = 4/3*Pi*a*b*b;
  }
  else
  {
    //V = Pi / 6 * diff(Rleft, Rright)*1.5 * ell*2;
    a = norm(Minus(Rleft->r(), Rright->r()))/2.;
    b = sqrt(ell/2.);
    V = 4/3*Pi*a*b*b;
  }
  xg = 0.5 * (Rright->r().at(0) + Rleft->r().at(0));
  yg = 0.5 * (Rright->r().at(1) + Rleft->r().at(1));
  zg = 0.5 * (Rright->r().at(2) + Rleft->r().at(2));
  unbound_term = 1.5 * get_square_diff(Rleft->r(), Rright->r()) / ell; // intermediate fastener computation
  IF(true){cout<<"Loop::Constructor : OVER"<<endl;}
}

Loop::Loop(const Loop &loop,
            Linker* new_left_linker,
            Linker* new_right_linker) : Strand(loop,new_left_linker)
{
  IF(true){cout<<"Loop constructor with new link"<<endl;}
  Rright = new_right_linker;
  ell_coordinate_1 = loop.ell_coordinate_1;
  a = loop.a;
  b = loop.b;
  unbound_term = loop.unbound_term;
  //compute_all_rates();
}

Loop::Loop(const Loop &loop) : Strand(loop)
{
  Rright = loop.Rright;
  ell_coordinate_1 = loop.ell_coordinate_1;
  a = loop.a;
  b = loop.b;
  unbound_term = loop.unbound_term;
  //compute_all_rates();
}

Strand* Loop::clone() const{return new Loop(*this);}

double Loop::Omega(Linker* r1, Linker* r2, double ell) const
{
  // return pow(4*Pi,ell)*pow(3/(2*Pi*ell),1.5)*exp(-3/2*(get_square_diff(Rleft,Rright))/ell);
  return pow(3 / (2 * Pi * ell), 1.5) * exp(-3 / 2 * get_square_diff(r1->r(), r2->r()) / ell);
}

array<double, 3> Loop::random_in_volume()
{
 double theta(atan2(yg, xg)), phi(atan2(xg, zg) - Pi / 2.);
  uniform_real_distribution<double> distribution(-1, 1); // doubles from -1 to 1
  double x(distribution(generator) * a);
  double y(distribution(generator) * b);
  double z(distribution(generator) * b);
  array<double, 3> res{cos(phi) * (cos(theta) * x - sin(theta) * y + sin(phi) * z),
                       sin(theta) * x + cos(theta) * y,
                       -sin(phi) * (cos(theta) * x - sin(theta) * y) + cos(phi) * z};
  res[0] += xg;
  res[1] += yg;
  res[2] += zg;
  // cout<<res[0]<<" "<<res[1]<<" "<<res[2]<<endl;
  return res;
}

void Loop::get_volume_limit(array<double,3>& main_ax, 
                            array<double,3>& ctr_mass,
                            double& a_in, double& b_in) const
{
  IF(true){cout<<"loop : get_volume_limit"<<endl;}
  ctr_mass = {xg,yg,zg};
  //a = norm(Minus(Rleft->r(),Rright->r()))*0.5;
  //b = sqrt(ell/2);
  a_in = a;
  b_in = b;
  main_ax = Minus(Rleft->r(),Rright->r());
}
                       

pair<unique_ptr<Strand>,unique_ptr<Strand>> Loop::bind() const
{
  // return a reference to a the two loop that must be created
  // when binding the current loop to a linker randomly choosen
  // in the vicinity of the current loop. The choice is made
  // in accordance with the binding rate of each specific linker
  Linker* linker_selected;
  double length;
  select_link_length(length,linker_selected);
  unique_ptr<Loop> left_loop =make_unique<Loop>(Rleft,linker_selected,ell_coordinate_0,ell_coordinate_0+length,rho0,slide);
  unique_ptr<Loop> right_loop = make_unique<Loop>(linker_selected,Rright,ell_coordinate_0+length,ell_coordinate_1,rho0,slide);
  return  {move(left_loop),move(right_loop)};
}
unique_ptr<Strand> Loop::unbind_from(Strand* left_strand) const
{
  // return a reference to the loop that must be created when
  // unbinding a linker between the current loop (at its right)
  // and "left_loop" that is at its left.
  return make_unique<Loop>(left_strand->get_Rleft(),
                          Rright,
                          left_strand->get_ell_coordinate_0(),
                          ell_coordinate_1,
                          rho0,
                          slide);
}
unique_ptr<Strand> Loop::do_slide(double dl,bool right) const
{
  // simply add dl to the right linker if right otherwise add dl to the left one
  if(right){
    return make_unique<Loop>(Rleft,Rright,ell_coordinate_0,ell_coordinate_1+dl,rho0,slide);
    }
  else{
    return make_unique<Loop>(Rleft,Rright,ell_coordinate_0+dl,ell_coordinate_1,rho0,slide);
  }
}
double Loop::compute_binding_rate(double li, Linker* linker) const
{
  if(diff(Rleft->r(), linker->r())>li or diff(linker->r(), Rright->r())>ell - li)
  {return 0.;}
  return 1/ell*exp(1.5 * log(3 * ell / (2 * Pi * li * (ell - li))) - 1.5 * (get_square_diff(Rleft->r(), linker->r()) / li + get_square_diff(linker->r(), Rright->r()) / (ell - li)) + unbound_term);
}
void Loop::Check_integrity() const
{
  cout<<ell_coordinate_0<<" "<<ell_coordinate_1<<endl;
}

// ---------------------------------------------------------------------------
//-----------------------------------accessor---------------------------------
// ---------------------------------------------------------------------------
Linker* Loop::get_Rright() const { return Rright; }

double Loop::get_theta() const { return atan2(0.5 * (Rright->r().at(1) - Rleft->r().at(1)), 0.5 * (Rright->r().at(0) - Rleft->r().at(0))); }

double Loop::get_phi() const { return atan2(0.5 * (Rright->r().at(0) - Rleft->r().at(0)), 0.5 * (Rright->r().at(2) - Rleft->r().at(2))) - Pi / 2.; }

array<double, 3> Loop::get_Rg() const { return {0.5 * (Rright->r().at(0) + Rleft->r().at(0)), 0.5 * (Rright->r().at(1) + Rleft->r().at(1)), 0.5 * (Rright->r().at(2) + Rleft->r().at(2))}; }

double Loop::get_ell_coordinate_1() const { return ell_coordinate_1; }

double Loop::get_S(double dl) const
{
  return log(pow(3 / (2 * Pi * (ell+dl)), 1.5)) - 3 / 2 * get_square_diff(Rleft->r(), Rright->r()) / (ell+dl);
}
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------