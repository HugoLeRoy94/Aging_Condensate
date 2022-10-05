#include "Header.h"

using namespace std;

// -----------------------------------------------------------------------------
// -----------------------------creator/destructor------------------------------
// -----------------------------------------------------------------------------
System::System(double ell_tot, double rho0, double temperature, int seed, bool adjust) : distrib(1, 10000000)
{
  IF(true) { cout << "System : creator" << endl; }
  // srand(seed);
  generator.seed(seed);
  // ---------------------------------------------------------------------------
  // -----------------------constant of the simulation--------------------------
  ell = ell_tot;
  rho = rho0;
  rho_adjust = adjust;
  kBT = temperature;
  array<double, 3> R0 = {0, 0, 0};
  dangling = Dangling(R0,linkers, 0., ell, rho, rho_adjust); // dummy dangling that helps generate crosslinkers but has none initially
  // ---------------------------------------------------------------------------
  //-----------------------------initialize crosslinkers------------------------
  generate_crosslinkers();
  // ---------------------------------------------------------------------------
  //-----------------------------initialize dangling----------------------------
  IF(true){cout<< "System : create dangling" << endl;}
  dangling = Dangling(R0,linkers, 0., ell, rho, rho_adjust);
  IF(true) { cout << "System : created" << endl; }
}

System::~System()
{
  loops.clear();
}
// -----------------------------------------------------------------------------
// ----------------------------Main function------------------------------------
// -----------------------------------------------------------------------------

double System::evolve(bool *bind)
{
  IF(true) { cout << "System : start evolve" << endl; }
  // -----------------------------------------------------------------------------
  // -----------------------------------------------------------------------------
  // compute the cumulative transition rates for each loop
  vector<double> cum_rates(loops.size() + 2); // the +1 is for removing a bond and +2 for binding the dangling end
  IF(true) { cout << "System : Start computing the cumulative probability array" << endl; }
  cum_rates[0] = (loops.size()) * exp(-1 / kBT) / (1 - exp(-1 / kBT));
  cum_rates[1] = dangling.get_total_binding_rates() + cum_rates[0];
  int n(2);
  for (auto &it : loops)
  {
    cum_rates[n] = cum_rates[n - 1] + it.get_total_binding_rates();
    n++;
  }
  if (cum_rates.back() == 0)
  {
    // if the system is blocked, we remake a loop with the possibility of having a crosslinker this time
    reset_crosslinkers();
    return 0;
  }
  // pick a random process
  IF(true) { cout << "System : draw a random process" << endl; }
  uniform_real_distribution<double> distribution(0, cum_rates.back());
  double pick_rate = distribution(generator);
  // becareful : if the rate_selec number is higher than cum_rates.back()  lower_bound returns cum_rates.back()
  vector<double>::iterator rate_selec = lower_bound(cum_rates.begin(), cum_rates.end(), pick_rate);
  if (rate_selec - cum_rates.begin() != 0)
  {
    IF(true) { cout << "System : add a bond" << endl; }
    if (rate_selec - cum_rates.begin() != 1)
    {
      IF(true) { cout << "System : add a bond to a loop" << endl; }
      add_bond(cum_rates, rate_selec);
    }
    else
    {
      IF(true) { cout << "System : add a bond to the polymer end" << endl; }
      add_bond_to_dangling();
    }
    *bind = true;
  }
  else
  {
    // select a random loop:
    IF(true) { cout << "System : remove a bond" << endl; }
    uniform_int_distribution<int> distribution(0, loops.size() - 1);
    set<Loop>::iterator loop_selec_left(loops.begin());
    int index(distribution(generator));
    advance(loop_selec_left, index);
    if (index < loops.size() - 1)
    {
      IF(true) { cout << "System : a bond between two loop" << endl; }
      set<Loop>::iterator loop_selec_right(loops.begin());
      advance(loop_selec_right, index + 1);
      // unbind
      unbind_loop(loop_selec_left, loop_selec_right);
    }
    else
    {
      IF(true) { cout << "System : remove the extremety one" << endl; }
      unbind_extremity(loop_selec_left);
    }
    *bind = false;
  }
  // cout<<cum_rates.back()<<endl;
  //check_loops_integrity();
  return draw_time(cum_rates.back());
}
// -----------------------------------------------------------------------------
//------------------------------Remove a bond-----------------------------------
// -----------------------------------------------------------------------------

void System::unbind_loop(set<Loop>::iterator &loop_selec_left, set<Loop>::iterator &loop_selec_right)
{
  IF(true) { cout << "System : unbind loop" << endl; }
  // create a new loop that is the combination of both inputs
  Loop loop = Loop(loop_selec_left->get_Rleft(),
                   loop_selec_right->get_Rright(),
                   linkers,
                   loop_selec_left->get_ell_coordinate_0(),
                   loop_selec_right->get_ell_coordinate_1(),
                   rho,
                   rho_adjust);

  //cout << loop_selec_left->get_ell()<<" "<<loop_selec_right->get_ell()<<" "<<loop.get_ell()<<endl;
  linkers.add(loop_selec_right->get_Rleft()[0],
              loop_selec_right->get_Rleft()[1],
              loop_selec_right->get_Rleft()[2],
              loop_selec_right->get_Rleft());
  loops.erase(loop_selec_right);
  loops.erase(loop_selec_left);
  loops.insert(loop);
  // time = draw_time(exp(-1/kBT)/(1+exp(-1/kBT)));
  // return time;
}

void System::unbind_extremity(set<Loop>::iterator &loop_selec_left)
{
  IF(true) { cout << "System : unbind the extremity" << endl; }
  dangling = Dangling(loop_selec_left->get_Rleft(),
                      linkers,
                      loop_selec_left->get_ell_coordinate_0(),
                      loop_selec_left->get_ell() + dangling.get_ell(),
                      rho,
                      rho_adjust);
  // ad the linkers from which we just unbind to the list
  linkers.add(loop_selec_left->get_Rright()[0],
            loop_selec_left->get_Rright()[1],
            loop_selec_left->get_Rright()[2],
            loop_selec_left->get_Rright());
  loops.erase(loop_selec_left);
}

// -----------------------------------------------------------------------------
//------------------------------Add a bond--------------------------------------
// -----------------------------------------------------------------------------

void System::add_bond(vector<double> &cum_rates, vector<double>::iterator &rate_selec)
{
  IF(true) { cout << "System : select the associated loop" << endl; }
  set<Loop>::iterator loop_selec(next(loops.begin(), distance(cum_rates.begin(), rate_selec)-2));
  // ask the loop for a random linker, and a random length
  double length;
  array<double, 3> r_selected;
  IF(true) { cout << "System : select a length and a r" << endl; }
  loop_selec->select_link_length(length, r_selected);
  linkers.remove(r_selected[0],r_selected[1],r_selected[2]);
  /*catch(const exception& e){
    cout<<e.what()<<"\n print all the crosslinkers r \n";
    cout<<"loops size :"<<loops.size()<<endl;
    for(auto& loop : loops){cout<<"loop.r.size = " << loop.get_r().size()<<endl;}//;for(auto& it : loop->get_r()){
        //cout<<it[0]<<" "<<it[1]<<" "<<it[2]<<endl;
      //}
    throw;
  }*/
  // cout<<length<<" "<<diff(loop_selec->get_Rleft(),r_selected)<<endl;
  //  create the two new loops
  IF(true) { cout << "System : create two new loops" << endl; }
  Loop l1 = Loop(loop_selec->get_Rleft(),
                 r_selected,
                 linkers,
                 loop_selec->get_ell_coordinate_0(),
                 loop_selec->get_ell_coordinate_0() + length,
                 rho,
                 rho_adjust);
  Loop l2 = Loop(r_selected,
                 loop_selec->get_Rright(),
                 linkers,
                 loop_selec->get_ell_coordinate_0() + length,
                 loop_selec->get_ell_coordinate_1(),
                 rho,
                 rho_adjust);
  //---------------------------------------------------------------------------------
  //------------------------------Actualize the linkers------------------------------
  // delete the loop
  IF(true) { cout << "System : delete the old loop" << endl; }
  // loops.erase(loops.begin()+distance(cum_rates.begin(),rate_selec));
  // delete loop_selec;
  set<Loop>::iterator hint = loops.erase(loop_selec);
  // add the newly created loop to the vector of loops
  IF(true) { cout << "System : add the two newly created loops" << endl; }
  loops.insert(hint,l1);
  loops.insert(hint,l2);
}

void System::add_bond_to_dangling()
{
  double length;
  array<double, 3> r_selected;
  IF(true) { cout << "System : select a length and r in dangling" << endl; }
  dangling.select_link_length(length, r_selected);
  linkers.remove(r_selected[0],r_selected[1],r_selected[2]);
  IF(true) { cout << "System : create a new loop and a new dangling" << endl; }
  Loop l = Loop(dangling.get_Rleft(),
                r_selected,
                linkers,
                dangling.get_ell_coordinate_0(),
                dangling.get_ell_coordinate_0() + length,
                rho,
                rho_adjust);
  dangling = Dangling(r_selected,
                      linkers,
                      dangling.get_ell_coordinate_0() + length,
                      dangling.get_ell() - length,
                      rho,
                      rho_adjust);
  loops.insert(loops.end(),l);
}
// -----------------------------------------------------------------------------
//-------reset the crosslinkers if no adding or removing is possible------------
// -----------------------------------------------------------------------------
void System::reset_crosslinkers()
{
  
  linkers.clear();
  generate_crosslinkers();
  // it basically consist in remaking every loops
  set<Loop> new_loops;
  for (auto &it : loops)
  {
    Loop newloop(it,linkers);
    new_loops.insert(new_loops.end(), newloop);
  }
  loops = new_loops;
  dangling = Dangling(dangling);

}
// -----------------------------------------------------------------------------
// -----------------------------accessor----------------------------------------
// -----------------------------------------------------------------------------
int System::get_N() const { return loops.size(); }

void System::get_R(double *R, int size) const
{
  if (size != 3 * (loops.size() + 1))
  {
    throw invalid_argument("invalid size in System::get_R");
  }
  // We construct an vector of the x,y,z coordinates to sort them according to x
  // for( auto& it : R_to_sort){for(auto& it2 : it){cout<<it2<<" ";}cout<<endl;}
  // fill the R array
  R[0] = 0.;
  R[1] = 0.;
  R[2] = 0.;
  int n(3);
  for (auto &it : loops)
  {
    R[n] = it.get_Rright()[0];
    R[n + 1] = it.get_Rright()[1];
    R[n + 2] = it.get_Rright()[2];
    n += 3;
  }
}

void System::get_ell_coordinates(double* ell_coordinate,int size)const 
{
  if(size != loops.size()+1){throw invalid_argument("invalid size in System::get_R");}
int n(0);
  for(auto& loop : loops){
    ell_coordinate[n] = loop.get_ell_coordinate_0();
    n++;
  }
  ell_coordinate[n] = dangling.get_ell_coordinate_0();
}

void System::get_ell(double *ells, int size) const
{
  if (size != loops.size())
  {
    throw invalid_argument("invalid size in System::get_ell");
  }
  int n(0);
  for (auto &it : loops)
  {
    ells[n] = it.get_ell();
    n++;
  }
}

void System::get_r(double *r, int size) const
{
  if (size != get_r_size())
  {
    throw invalid_argument("invalid size in System::get_r");
  }
  int n(0);
  for (auto &loop : loops)
  {
    for (auto &link : loop.get_r())
    {
      r[3 * n] = (*link)[0];
      r[3 * n + 1] = (*link)[1];
      r[3 * n + 2] = (*link)[2];
      n++;
    }
  }
  for (auto &link : dangling.get_r())
  {
    r[3 * n] = (*link)[0];
    r[3 * n + 1] = (*link)[1];
    r[3 * n + 2] = (*link)[2];
    n++;
  }
}

int System::get_r_size() const
{
  int size(0);
  for (auto &it : loops)
  {
    size += it.get_r().size();
  }
  size += dangling.get_r().size();
  return 3 * size;
}

double System::get_F() const
{
  double F(-get_N());
  for (auto &it : loops)
  {
    F += -kBT * it.get_S();
  }
  F += dangling.get_S();
  return F;
}

void System::Print_Loop_positions() const
{
  for (auto &loop : loops)
  {
    cout << "theta Phi " << loop.get_theta() << " " << loop.get_phi() << endl;
    cout << "xg,yg,zg " << loop.get_Rg()[0] << " " << loop.get_Rg()[1] << " " << loop.get_Rg()[2] << endl;
    cout << "volume " << loop.get_V() << endl;
    cout << "ell " << loop.get_ell() << endl;
    cout << "anchoring points. left :" << loop.get_Rleft()[0] << " " << loop.get_Rleft()[1] << " " << loop.get_Rleft()[2] << " right : " << loop.get_Rright()[0] << " " << loop.get_Rright()[1] << " " << loop.get_Rright()[2] << endl;
    cout << "number of crosslinkers of this loop :" << loop.get_r().size() << endl
         << endl;
  }
}

double System::draw_time(double rate) const
{
  uniform_real_distribution<double> distrib;
  double xi(distrib(generator));
  return -log(1 - xi) / rate;
}

void System::check_loops_integrity() const
{
  double ell_tot(0);
  for(auto& loop:loops){ell_tot+=loop.get_ell();}
  ell_tot+=dangling.get_ell();
  if(ell_tot != ell){cout<< "ell tot = " << ell_tot<<endl;throw range_error("ell is not equal to the sum of the ells");}
  for(auto& loop : loops)
  {
    try
    {
      cout<<loop.get_r().size()<<endl;
    }
    catch(const exception& e){throw e;}
    for(int r=0;r<loop.get_r().size();r++)
    {
      try
      {
        loop.get_r().at(r);
      }
      catch(const exception& e){throw e;}
    }
  }
  cout<<" r vector i working"<<endl;
}

void System::get_r_system(double *r, int size)const
{
  if (size != get_r_system_size())
  {
    throw invalid_argument("invalid size in System::get_r");
  }
  int n(0);
  for(auto& it : linkers.underlying_array())
  {
    for(auto& it2 : it.second)
    {
      for(auto& it3 : it2.second)
      {
        r[3*n] = it3.second[0];
        r[3*n+1] = it3.second[1];
        r[3*n+2] = it3.second[2];
        n++;
      }
    }
  }
}

int System::get_r_system_size() const
{
  return 3*linkers.get_number_of_elements();
}

void System::print_random_stuff() const{
  cout<<"dangling linkers : "<<endl;
  for(auto& it : dangling.get_r()){
    cout<<(*it)[0]<< " " << (*it)[1] << " "<<(*it)[2]<<endl;
  } 
  cout<<"loops linkers : "<<endl;
  for(auto& it : loops){for(auto& link : it.get_r()){
    cout<<(*link)[0]<< " " << (*link)[1] << " "<<(*link)[2]<<endl;
  }}
}

void System::generate_crosslinkers(){
  set<double> xminl,xmaxl,yminl,ymaxl,zminl,zmaxl;
  for(auto& loop : loops){
    double ximin,ximax,yimin,yimax,zimin,zimax;
    loop.get_volume_limit(ximin,ximax,yimin,yimax,zimin,zimax);
    xminl.insert(ximin);xmaxl.insert(ximax);yminl.insert(yimin);ymaxl.insert(yimax);zminl.insert(zimin);zmaxl.insert(zimax);
  }
  double dxmin,dxmax,dymin,dymax,dzmin,dzmax;
  dangling.get_volume_limit(dxmin,dxmax,dymin,dymax,dzmin,dzmax);
  xminl.insert(dxmin);xmaxl.insert(dxmax);yminl.insert(dymin);ymaxl.insert(dymax);zmaxl.insert(dzmax);zminl.insert(dzmin);
  double  xmin(*(xminl.begin())),xmax(*(xmaxl.rbegin()));
  double ymin(*(yminl.begin())),ymax(*(ymaxl.rbegin()));
  double zmin(*(zminl.begin())),zmax(*(zmaxl.rbegin()));
  double V((xmax-xmin)*(ymax-ymin)*(zmax-zmin));
  /*cout<<V<<endl;
  cout<<xmin<<" "<<xmax<<endl;
  cout<<ymin<<" "<<ymax<<endl;
  cout<<zmin<<" "<<zmax<<endl;*/
  poisson_distribution<int> distribution(rho * V);
  int N_crosslinker = distribution(generator);
  for(int n =0; n<N_crosslinker;n++)
  {
    uniform_real_distribution<double> distribution(-1, 1); // doubles from -1 to 1
    double x(distribution(generator) * (xmax-xmin)+xmin);
    double y(distribution(generator) * (ymax-ymin)+ymin);
    double z(distribution(generator) * (zmax-zmin)+zmin);
    array<double,3> link({x,y,z});
    linkers.add(x,y,z,link);
  }
}