#include "Header.h"
using namespace std;

Loop::Loop(array<double,3> R0,array<double,3> R1, double ell_loop){
  // the right anchoring point of the loop
  Rleft=R0;
  Rright=R1;
  ell = ell_loop;
  // create r by randmoly generating positions
  generate_binding_sites();
}
Loop::~Loop(){}
std::array<double,3> Loop::get_Rright() const{
  return Rright;
}
std::vector<array<double,3>> Loop::get_r() const{
  return r;
}
double Loop::get_ell() const{
  return ell;
}
double Loop::compute_binding_rate(int r_index, double ell){
  return 0;
}
void Loop::select_link_length(double& length, array<double,3>& r_selected){

}
void Loop::generate_binding_sites(){

}
