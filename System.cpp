#include "Header.h"

using namespace std;


System::System(){}
System::System(double ell_tot,double distance_anchor,double rho0,double temperature){
  int seed(1204985);
  srand(seed);
  generator.seed(seed);
  // ---------------------------------------------------------------------------
  // -----------------------constant of the simulation--------------------------
  ell = ell_tot;
  D = distance_anchor;
  rho = rho0;
  kBT = temperature;
  // ---------------------------------------------------------------------------
  //-----------------------------initialize loops-------------------------------
  // *reserve* makes the vector to never reallocate memory, and fasten the push_backs
  loops.reserve(ell_tot); // maximum number of loop that can be inserted
  array<double,3> Rf={D,0.,0.},R0={0.,0.,0.};
  loops.push_back(new Loop(R0,Rf,ell,rho0)); // create the first loop.
}
System::~System(){
  for(auto& it : loops){delete it;}
  loops.clear();
}
// -----------------------------------------------------------------------------
// ----------------------------Main function------------------------------------
// -----------------------------------------------------------------------------
double System::evolve(){
  // cumpute the cumulative transition rates for each loop
  // !!!!!!!!!!!! only binding move implemented!!!!!!!!!!!!
  // -----------------------------------------------------------------------------
  // -----------------------------------------------------------------------------
  vector<double> cum_rates(loops.size());
  cum_rates[0]=loops[0]->get_total_binding_rates();
  for(int n=1; n<loops.size();n++){
    cum_rates[n] = cum_rates[n-1]+loops[n]->get_total_binding_rates();//+exp(1/kBT);
  }
  // -----------------------------------------------------------------------------
  // -----------------------------------------------------------------------------
  // pick a random loop
  uniform_real_distribution<double> distribution(0,cum_rates.back());
  auto rate_selec = lower_bound(cum_rates.begin(),cum_rates.end(),distribution(generator));
  Loop* loop_selec(loops[distance(cum_rates.begin(),rate_selec)]);
  // ask the loop for a random linker, and a random length
  double length;
  double time;
  array<double,3> r_selected;
  loop_selec->select_link_length(length,r_selected,time);
  cout<<length<<" "<<diff(loop_selec->get_Rleft(),r_selected)<<endl;
  // create the two new loops
  Loop* l1 = new Loop(loop_selec->get_Rleft(),r_selected,length,rho);
  Loop* l2 = new Loop(r_selected,loop_selec->get_Rright(),loop_selec->get_ell()-length,rho);
  // delete the loop
  loops.erase(loops.begin()+distance(cum_rates.begin(),rate_selec));
  delete loop_selec;
  // add the newly created loop to the vector of loops
  loops.push_back(l1);
  loops.push_back(l2);
  return time;
}
// -----------------------------------------------------------------------------
// -----------------------------accessor----------------------------------------
// -----------------------------------------------------------------------------
int System::get_N() const{return loops.size();}

void System::get_R(double* R, int size) const{
  if(size!=3*(loops.size()+1)){
    throw invalid_argument("invalid size in System::get_R");}
  // fill the vector R with the value of R of each loops
  R[0] = 0.;
  R[1] = 0.;
  R[2] = 0.;
  for(int n=1;n<loops.size()+1;n++){
    R[3*n] = loops[n-1]->get_Rright()[0];
    R[3*n+1] = loops[n-1]->get_Rright()[1];
    R[3*n+2] = loops[n-1]->get_Rright()[2];
  }
}

void System::get_ell(double* ells,int size) const{
  if(size!=loops.size()){throw invalid_argument("invalid size in System::get_ell");}
  for(int n=0;n<size;n++){
    ells[n] = loops[n]->get_ell();
  }
}

void System::get_r(double* r, int size) const{
  if(size!=get_r_size()){throw invalid_argument("invalid size in System::get_r");}
  int n(0);
  for(auto& loop :loops ){
    for(auto& link : loop->get_r()){
      r[3*n] = link[0];
      r[3*n+1] = link[1];
      r[3*n+2] = link[2];
      n++;
    }
  }
}

int System::get_r_size()const{
  int size(0);
  for(auto& it :loops){
    size+=it->get_r().size();
  }
  return 3*size;
}

void System::Print_Loop_positions(){
  for(auto& loop : loops){
    cout<<"theta Phi "<<loop->get_theta()<<" "<<loop->get_phi()<<endl;
    cout<<"xg,yg,zg "<<loop->get_Rg()[0]<<" "<<loop->get_Rg()[1]<<" "<<loop->get_Rg()[2]<<endl;
    cout<< "volume "<<loop->get_V()<<endl;
    cout<<"ell "<<loop->get_ell()<<endl;
    cout<< "anchoring points. left :"<<loop->get_Rleft()[0]<<" "<<
                                       loop->get_Rleft()[1]<<" "<<
                                       loop->get_Rleft()[2]<<
                                       " right : "<<
                                       loop->get_Rright()[0]<<" "<<
                                       loop->get_Rright()[1]<<" "<<
                                       loop->get_Rright()[2]<<endl;
    cout<<"number of crosslinkers of this loop :"<<loop->get_r().size()<<endl<<endl;
  }
}
