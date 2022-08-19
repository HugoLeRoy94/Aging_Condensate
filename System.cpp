#include "Header.h"

using namespace std;

// -----------------------------------------------------------------------------
// -----------------------------creator/destructor------------------------------
// -----------------------------------------------------------------------------
System::System(){}
System::System(double ell_tot,double distance_anchor,double rho0,double temperature,int seed,bool adjust) : distrib(1,10000000){
  //srand(seed);
  generator.seed(seed);
  // ---------------------------------------------------------------------------
  // -----------------------constant of the simulation--------------------------
  ell = ell_tot;
  D = distance_anchor;
  anchor = {D,0.,0.};
  rho = rho0;
  rho_adjust = adjust;
  kBT = temperature;
  // ---------------------------------------------------------------------------
  //-----------------------------initialize loops-------------------------------
  // *reserve* makes the vector to never reallocate memory, and fasten the push_backs
  array<double,3> Rf={D,0.,0.},R0={0.,0.,0.};
  loops.insert(loops.begin(),Loop(R0,Rf,ell,rho,D/ell,rho_adjust)); // create the first loop.
  //loops.push_back(new Loop(ell,rho0,distrib(generator))); // create the first loop.
}
System::~System(){
  //for(auto& it : loops){delete it;}
  loops.clear();
}
// -----------------------------------------------------------------------------
// ----------------------------Main function------------------------------------
// -----------------------------------------------------------------------------
double System::evolve(bool* bind){
  // -----------------------------------------------------------------------------
  // -----------------------------------------------------------------------------
  // compute the cumulative transition rates for each loop
  vector<double> cum_rates(loops.size()+1); // the +1 is for removing a bond
  IF(true){cout<<"System : Start computing the cumulative probability array"<<endl;}
  int n(0);
  for(auto& it : loops){
    if( n ==0 ){ cum_rates[n] = it.get_total_binding_rates();}
    else{cum_rates[n] = cum_rates[n-1]+it.get_total_binding_rates();}
    n++;
  }
  cum_rates.back() = (loops.size()-1)*exp(-1/kBT)/(1-exp(-1/kBT))+cum_rates.rbegin()[1];
  if(cum_rates.back() == 0){
    // if the system is blocked, we remake a loop with the possibility of having a crosslinker this time
    return reset_crosslinkers();
  }
  //for(auto& loop : loops){for(auto& it : loop.get_rates()){for(auto& it2 : it){cout<<it2<<" ";}}cout<<endl;}
  // pick a random process
  IF(true){cout<<"System : draw a random process"<<endl;}
  uniform_real_distribution<double> distribution(0,cum_rates.back());
  double pick_rate = distribution(generator);
    // becareful : if the rate_selec number is higher than cum_rates.back()  lower_bound returns cum_rates.back()
  vector<double>::iterator rate_selec = lower_bound(cum_rates.begin(),cum_rates.end(),pick_rate);
  if(rate_selec != cum_rates.end()-1){
    IF(true){cout<<"System : add a bond"<<endl;}
    add_bond(cum_rates,rate_selec);
    *bind = true;
    //catch(out_of_range& e){for(auto& it : cum_rates){cout<<it<<endl;}throw;}
  }
  else{
    // select a random loop:
    IF(true){cout<<"System : remove a bond"<<endl;}
    uniform_int_distribution<int> distribution(0,loops.size()-2);
    set<Loop>::iterator loop_selec_left(loops.begin());
    set<Loop>::iterator loop_selec_right(loops.begin());
    int index(distribution(generator));
    advance(loop_selec_left,index);
    advance(loop_selec_right,index+1);
    // unbind
    unbind_loop(loop_selec_left, loop_selec_right);
    *bind = false;
  }
  //cout<<cum_rates.back()<<endl;
  return draw_time(cum_rates.back());
}
// -----------------------------------------------------------------------------
//------------------------------Remove a bond-----------------------------------
// -----------------------------------------------------------------------------
void System::unbind_loop(set<Loop>::iterator& loop_selec_left,set<Loop>::iterator& loop_selec_right){
  IF(true){cout<<"unbind loop"<<endl;}
  if(loops.size()==1){return;} // skip the unbinding if it's the last loop
  // create a new loop that is the combination of both inputs
  Loop loop = Loop(loop_selec_left->get_Rleft(),
                   loop_selec_right->get_Rright(),
                   loop_selec_left->get_ell()+loop_selec_right->get_ell(),
                   rho,
                   D/ell,rho_adjust
                   );
  loops.erase(loop_selec_right);
  loops.insert(loop_selec_left,loop);
  loops.erase(loop_selec_left);
  //time = draw_time(exp(-1/kBT)/(1+exp(-1/kBT)));
  //return time;
}
// -----------------------------------------------------------------------------
//------------------------------Add a bond--------------------------------------
// -----------------------------------------------------------------------------
void System::add_bond(vector<double>& cum_rates, vector<double>::iterator& rate_selec){
  IF(true){cout<<"System : select the associated loop"<<endl;}
  set<Loop>::iterator loop_selec(next(loops.begin(),distance(cum_rates.begin(),rate_selec)));
  // ask the loop for a random linker, and a random length
  double length;
  array<double,3> r_selected;
  IF(true){cout<<"System : select a length and a r"<<endl;}
  loop_selec->select_link_length(length,r_selected);
  /*catch(const exception& e){
    cout<<e.what()<<"\n print all the crosslinkers r \n";
    cout<<"loops size :"<<loops.size()<<endl;
    for(auto& loop : loops){cout<<"loop.r.size = " << loop.get_r().size()<<endl;}//;for(auto& it : loop->get_r()){
        //cout<<it[0]<<" "<<it[1]<<" "<<it[2]<<endl;
      //}
    throw;
  }*/
  //cout<<length<<" "<<diff(loop_selec->get_Rleft(),r_selected)<<endl;
  // create the two new loops
  IF(true){cout<<"System : create two new loops"<<endl;}
  Loop l1 = Loop(loop_selec->get_Rleft(),r_selected,length,rho,D/ell,rho_adjust);
  Loop l2 = Loop(r_selected,loop_selec->get_Rright(),loop_selec->get_ell()-length,rho,D/ell,rho_adjust);
  // delete the loop
  IF(true){cout<<"System : delete the old loop"<<endl;}
  //loops.erase(loops.begin()+distance(cum_rates.begin(),rate_selec));
  //delete loop_selec;
  set<Loop>::iterator hint = loops.erase(loop_selec);
  /*cout<<"////////////////////////// two new loops /////////////////////"<<endl;
  cout<<l1.get_ell()<<" "<<diff(l1.get_Rleft(),l1.get_Rright())<<endl;
  cout<<l2.get_ell()<<" "<<diff(l2.get_Rleft(),l2.get_Rright())<<endl;
  cout<<"//////////////////////////////////////////////////////////////"<<endl;*/
  // add the newly created loop to the vector of loops
  IF(true){cout<<"System : add the two newly created loop"<<endl;}
  loops.insert(hint,l1);
  loops.insert(hint,l2);
}
// -----------------------------------------------------------------------------
//-------reset the crosslinkers if no adding or removing is possible------------
// -----------------------------------------------------------------------------
double System::reset_crosslinkers(){
  // it basically consist in remaking every loops
  set<Loop> new_loops;
  for(auto& it : loops){
    Loop newloop(it);
    new_loops.insert(new_loops.end(),newloop);
  }
  loops = new_loops;
  return 0.;
}
// -----------------------------------------------------------------------------
// -----------------------------accessor----------------------------------------
// -----------------------------------------------------------------------------
int System::get_N() const{return loops.size();}

void System::get_R(double* R, int size) const{
  if(size!=3*(loops.size()+1)){
    throw invalid_argument("invalid size in System::get_R");}
  // We construct an vector of the x,y,z coordinates to sort them according to x
  //for( auto& it : R_to_sort){for(auto& it2 : it){cout<<it2<<" ";}cout<<endl;}
  // fill the R array
  R[0] = 0.;
  R[1] = 0.;
  R[2] = 0.;
  int n(3);
  for(auto& it : loops){
    R[n] = it.get_Rright()[0];
    R[n+1] = it.get_Rright()[1];
    R[n+2] = it.get_Rright()[2];
    n+=3;
  }
}

void System::get_ell(double* ells,int size) const{
  if(size!=loops.size()){throw invalid_argument("invalid size in System::get_ell");}
  int n(0);
  for(auto& it : loops){
    ells[n] = it.get_ell();
    n++;
  }
}

void System::get_r(double* r, int size) const{
  if(size!=get_r_size()){throw invalid_argument("invalid size in System::get_r");}
  int n(0);
  for(auto& loop :loops ){
    for(auto& link : loop.get_r()){
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
    size+=it.get_r().size();
  }
  return 3*size;
}

double System::get_F()const{
  double F(-get_N());
  for(auto& it : loops){
    F+= - kBT*it.get_S();
  }
  return F;
}

void System::Print_Loop_positions() const{
  for(auto& loop : loops){
    cout<<"theta Phi "<<loop.get_theta()<<" "<<loop.get_phi()<<endl;
    cout<<"xg,yg,zg "<<loop.get_Rg()[0]<<" "<<loop.get_Rg()[1]<<" "<<loop.get_Rg()[2]<<endl;
    cout<< "volume "<<loop.get_V()<<endl;
    cout<<"ell "<<loop.get_ell()<<endl;
    cout<< "anchoring points. left :"<<loop.get_Rleft()[0]<<" "<<
                                       loop.get_Rleft()[1]<<" "<<
                                       loop.get_Rleft()[2]<<
                                       " right : "<<
                                       loop.get_Rright()[0]<<" "<<
                                       loop.get_Rright()[1]<<" "<<
                                       loop.get_Rright()[2]<<endl;
    cout<<"number of crosslinkers of this loop :"<<loop.get_r().size()<<endl<<endl;
  }
}

double System::draw_time(double rate) const{
  uniform_real_distribution<double> distrib;
  double xi(distrib(generator));
  return -log(1-xi)/rate;
}
