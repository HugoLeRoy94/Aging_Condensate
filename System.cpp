#include "Header.h"

using namespace std;

// -----------------------------------------------------------------------------
// -----------------------------creator/destructor------------------------------
// -----------------------------------------------------------------------------
System::System(){}
System::System(double ell_tot,double distance_anchor,double rho0,double temperature,int seed) : distrib(1,10000000){
  //srand(seed);
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
  loops.push_back(new Loop(R0,Rf,ell,rho0,distrib(generator))); // create the first loop.
  //loops.push_back(new Loop(ell,rho0,distrib(generator))); // create the first loop.
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
  vector<double> cum_rates(loops.size()+1); // the +1 is for removing a bond
  IF(true){cout<<"System : Start computing the cumulative probability array"<<endl;}
  cum_rates[0]=loops[0]->get_total_binding_rates();
  for(int n=1; n<loops.size();n++){
    cum_rates[n] = cum_rates[n-1]+loops[n]->get_total_binding_rates();//+exp(1/kBT);
  }
  cum_rates.back() = loops.size()*exp(-1/kBT)+cum_rates.rbegin()[1];
  // pick a random process
  IF(true){cout<<"System : draw a random process"<<endl;}
  uniform_real_distribution<double> distribution(0,cum_rates.back());
    // becareful : if the rate_selec number is higher than cum_rates.back()  lower_bound returns cum_rates.back()
  vector<double>::iterator rate_selec = lower_bound(cum_rates.begin(),cum_rates.end(),distribution(generator));
  int time(0);
  if(*rate_selec != cum_rates.back()){
    time = add_bond(cum_rates,rate_selec);}
  else{
    // select a random loop:
    uniform_int_distribution<int> distribution(0,loops.size()-2);
    vector<Loop*>::iterator loop_selec_left(loops.begin());
    vector<Loop*>::iterator loop_selec_right(loops.begin());
    int index(distribution(generator));
    advance(loop_selec_left,index);
    advance(loop_selec_right,index+1);
    // unbind
    time = unbind_loop(loop_selec_left, loop_select_right);
  }
  return time;
}
// -----------------------------------------------------------------------------
//------------------------------Remove a bond-----------------------------------
// -----------------------------------------------------------------------------
double System::unbind_loop(vector<Loop*>::iterator& loop_selec_left,vector<Loop*>::iterator& loop_selec_left){
  double time(0);
  IF(true){cout<<"unbind loop"<<endl;}
  if(loops.size()==1){return time;} // skip the unbinding if it's the last loop
  // the unbinding of a loop

  return time;
}
// -----------------------------------------------------------------------------
//------------------------------Add a bond--------------------------------------
// -----------------------------------------------------------------------------
double System::add_bond(vector<double>& cum_rates, vector<double>::iterator& rate_selec){
  IF(true){cout<<"System : select the associated loop"<<endl;}
  Loop* loop_selec(loops[distance(cum_rates.begin(),rate_selec)]);
  // ask the loop for a random linker, and a random length
  double length;
  double time;
  array<double,3> r_selected;
  IF(true){cout<<"System : select a length and a r"<<endl;}
  try{loop_selec->select_link_length(length,r_selected,time);}
  catch(const exception& e){
    cout<<e.what()<<"\n print all the crosslinkers r \n";
    cout<<"loops size :"<<loops.size()<<endl;
    for(auto& loop : loops){cout<<"loop.r.size = " << loop->get_r().size()<<endl;}//;for(auto& it : loop->get_r()){
        //cout<<it[0]<<" "<<it[1]<<" "<<it[2]<<endl;
      //}
    throw;
  }
  //cout<<length<<" "<<diff(loop_selec->get_Rleft(),r_selected)<<endl;
  // create the two new loops
  IF(true){cout<<"System : create two new loops"<<endl;}
  Loop* l1 = new Loop(loop_selec->get_Rleft(),r_selected,length,rho,distrib(generator));
  Loop* l2 = new Loop(r_selected,loop_selec->get_Rright(),loop_selec->get_ell()-length,rho,distrib(generator));
  // delete the loop
  IF(true){cout<<"System : delete the old loop"<<endl;}
  loops.erase(loops.begin()+distance(cum_rates.begin(),rate_selec));
  delete loop_selec;
  // add the newly created loop to the vector of loops
  IF(true){cout<<"System : add the two newly created loop"<<endl;}
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
  // We construct an vector of the x,y,z coordinates to sort them according to x
  vector<array<double,3>> R_to_sort(loops.size()+1);
  // fill the vector R with the value of R of each loops
  R_to_sort[0][0] = 0.;
  R_to_sort[0][1] = 0.;
  R_to_sort[0][2] = 0.;
  for(int n=1;n<loops.size()+1;n++){
    R_to_sort[n][0] = loops[n-1]->get_Rright()[0];
    R_to_sort[n][1] = loops[n-1]->get_Rright()[1];
    R_to_sort[n][2] = loops[n-1]->get_Rright()[2];
  }
  // sort according to the value in 0
  sort(R_to_sort.begin(),R_to_sort.end(),customLess);
  //for( auto& it : R_to_sort){for(auto& it2 : it){cout<<it2<<" ";}cout<<endl;}
  // fill the R array
  for(int n =0;n<loops.size()+1;n++){
    for(int k=0;k<3;k++){R[3*n+k] = R_to_sort[n][k];}
  }
}

void System::get_ell(double* ells,int size) const{
  if(size!=loops.size()){throw invalid_argument("invalid size in System::get_ell");}
  vector<pair<double,double>> ells_R_to_sort(size);
  for(int n=0;n<size;n++){
    ells_R_to_sort[n].first = loops[n]->get_ell();
    ells_R_to_sort[n].second = loops[n]->get_Rright()[0];
  }
  sort(ells_R_to_sort.begin(),ells_R_to_sort.end(),customLess2);
  //for(auto& it : ells_R_to_sort){cout<<it.first<<" "<<it.second<<endl;}
  for(int n = 0;n<size;n++){
    ells[n] = ells_R_to_sort[n].first;
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
