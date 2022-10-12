#include "Header.h"

using namespace std;

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
    Linker* R0 = new Linker({0, 0, 0});
    Strand* dummy_dangling = new Dangling(R0,loop_link, 0., ell, rho, rho_adjust); // dummy dangling that helps generate crosslinkers but has none initially
    loop_link.Insert_Strand(dummy_dangling);
    // ---------------------------------------------------------------------------
    //-----------------------------initialize crosslinkers------------------------
    generate_crosslinkers();
    loop_link.Remove_Strand(dummy_dangling);
    // ---------------------------------------------------------------------------
    //-----------------------------initialize dangling----------------------------
    IF(true){cout<< "System : create dangling" << endl;}
    loop_link.Insert_Strand(new Dangling(R0,loop_link, 0., ell, rho, rho_adjust));
    //print_random_stuff();
    //for(auto& it : linker_to_strand){for(auto& it2 : it.second){cout<<it2->get_Rleft()[0]<<" "<<it2->get_Rleft()[1]<<" "<<it2->get_Rleft()[2]<<endl;}}
    IF(true) { cout << "System : created" << endl; }
}

System::~System()
{
  loop_link.delete_pointers();
}

double System::evolve(bool *bind)
{
  
  IF(true) { cout << "System : start evolve" << endl; }
  IF(true){check_loops_integrity();}
  IF(true){cout<< "System : check integrity of every loops"<<endl;}
  // -----------------------------------------------------------------------------
  // -----------------------------------------------------------------------------
  // compute the cumulative transition rates for each loop
  vector<double> cum_rates(loop_link.get_strand_size() + 1); // the +1 is for removing a bond and +2 for binding the dangling end
  IF(true) { cout << "System : Start computing the cumulative probability array" << endl; }
  cum_rates[0] = loop_link.get_strand_size()  * exp(-1 / kBT) / (1 - exp(-1 / kBT));
  int n(1);
  for (auto &it : loop_link.get_loops())
  {
    cum_rates[n] = cum_rates[n - 1] + it->get_total_binding_rates();
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
    IF(true) { cout << "System : add a bond to a loop" << endl; }
    add_bond(cum_rates, rate_selec);
    *bind = true;
  }
  else
  {
    // select a random loop:
    IF(true) { cout << "System : remove a bond" << endl; }
    // select the index of the left bond to remove
    // the left bond index maximum is loop.size -2
    uniform_int_distribution<int> distribution(0, loop_link.get_strand_size() - 2);
    int index(distribution(generator));
    set<Strand*>::iterator loop_selec_left(loop_link.get_loop(index));
    set<Strand*>::iterator loop_selec_right(loop_link.get_loop(index+1));
    // unbind
    unbind_loop(loop_selec_left, loop_selec_right);
    *bind = false;
  }
  return draw_time(cum_rates.back());
}

void System::add_bond(vector<double> &cum_rates, vector<double>::iterator &rate_selec)
{
    IF(true) { cout << "System : select the associated loop" << endl; }
    // cum_rates.begin() is unbinding
    set<Strand*>::iterator loop_selec(loop_link.get_loop(distance(cum_rates.begin()+1, rate_selec))); 
    // ask the loop for a random linker, and a random length
    double length;
    Linker* r_selected;
    IF(true) { cout << "System : select a length and a r" << endl; }
    // select a link to bind to
    (*loop_selec)->select_link_length(length, r_selected);
    //  create the two new loops
    IF(true) { cout << "System : create two new loops" << endl; }
    Loop* loop_left = new Loop((*loop_selec)->get_Rleft(),
                    r_selected,
                    loop_link,
                    (*loop_selec)->get_ell_coordinate_0(),
                    (*loop_selec)->get_ell_coordinate_0() + length,
                    rho,
                    rho_adjust);
    Strand* loop_right;
    try
    {// if loop_selec isn't dangling :
        (*loop_selec)->get_Rright();
        // then the right side is a loop
        loop_right = new Loop(r_selected,
                            (*loop_selec)->get_Rright(),
                            loop_link,
                            (*loop_selec)->get_ell_coordinate_0() + length,
                            (reinterpret_cast<Loop*>(*loop_selec))->get_ell_coordinate_1(),
                            rho,
                            rho_adjust);
    }
    catch(out_of_range oor)
    {
        // then the right side is a Dangling
        loop_right = new Dangling(r_selected,
                                loop_link,
                                (*loop_selec)->get_ell_coordinate_0() + length,
                                reinterpret_cast<Dangling*>(*loop_selec)->get_ell() - length,
                                rho,
                                rho_adjust);
    }
    //---------------------------------------------------------------------------------
    //------------------------------Actualize the linkers------------------------------
    // delete the loop
    IF(true) { cout << "System : delete the old loop" << endl; }
    loop_link.Remove_Strand((*loop_selec));
    IF(true) { cout << "System : add the two newly created loops" << endl; }
    loop_link.Insert_Strand(loop_left);
    loop_link.Insert_Strand(loop_right);
}

void System::unbind_loop(set<Strand*>::iterator &loop_selec_left, set<Strand*>::iterator &loop_selec_right)
{
    IF(true) { cout << "System : unbind loop" << endl; }
    // create a new loop that is the combination of both inputs
    Strand* loop;
    try
    {   
        (*loop_selec_right)->get_Rright();
        loop = new Loop((*loop_selec_left)->get_Rleft(),
                        (*loop_selec_right)->get_Rright(),
                        loop_link,
                        (*loop_selec_left)->get_ell_coordinate_0(),
                        reinterpret_cast<Loop*>(*loop_selec_right)->get_ell_coordinate_1(),
                        rho,
                        rho_adjust);
    }
    catch(out_of_range oor)
    {
        loop = new Dangling((*loop_selec_left)->get_Rleft(),
                        loop_link,
                        (*loop_selec_left)->get_ell_coordinate_0(),
                        reinterpret_cast<Loop*>(*loop_selec_right)->get_ell(),
                        rho,
                        rho_adjust);
    }
    IF(true){cout<<"System : remove the old strand from loop_link"<<endl;}
    loop_link.Remove_Strand(*loop_selec_left);
    loop_link.Remove_Strand(*loop_selec_right);
    IF(true){cout<<"System : Insert the newly created loop"<<endl;}
    loop_link.Insert_Strand(loop);
}

double System::draw_time(double rate) const
{
  uniform_real_distribution<double> distrib;
  double xi(distrib(generator));
  return -log(1 - xi) / rate;
}

void System::reset_crosslinkers()
{
    // clear the container which clear also the loops from their linkers
    loop_link.delete_linkers();
    // create a new set of free crosslinkers
    generate_crosslinkers();
    // set all the crosslinker into the linker_to_strand
    // which also add the bound extremities
    reset_loops();
}

void System::generate_crosslinkers(){
  set<double> xminl,xmaxl,yminl,ymaxl,zminl,zmaxl;
  for(auto& loop : loop_link.get_loops()){
    double ximin,ximax,yimin,yimax,zimin,zimax;
    loop->get_volume_limit(ximin,ximax,yimin,yimax,zimin,zimax);
    xminl.insert(ximin);xmaxl.insert(ximax);yminl.insert(yimin);ymaxl.insert(yimax);zminl.insert(zimin);zmaxl.insert(zimax);
  }
  double  xmin(*(xminl.begin())),xmax(*(xmaxl.rbegin()));
  double ymin(*(yminl.begin())),ymax(*(ymaxl.rbegin()));
  double zmin(*(zminl.begin())),zmax(*(zmaxl.rbegin()));
  double V((xmax-xmin)*(ymax-ymin)*(zmax-zmin));
  poisson_distribution<int> distribution(rho * V);
  int N_crosslinker = distribution(generator);
  for(int n =0; n<N_crosslinker;n++)
  {
    uniform_real_distribution<double> distribution(0, 1); // doubles from -1 to 1
    double x(distribution(generator) * (xmax-xmin)+xmin);
    double y(distribution(generator) * (ymax-ymin)+ymin);
    double z(distribution(generator) * (zmax-zmin)+zmin);
    loop_link.create_new_free_linker(x,y,z);
  }
}

void System::reset_loops()
{
    set<Strand*,LessLoop> new_loops;
    for(auto& loop : loop_link.get_loops())
    {
        try
        {
            loop->get_Rright(); // check to know if it's dangling
            new_loops.insert(new Loop(*reinterpret_cast<Loop*>(loop),loop_link));
        }
        catch(out_of_range oor)
        {
            new_loops.insert(new Dangling(*reinterpret_cast<Dangling*>(loop),loop_link));
        }
    }
    loop_link.delete_loops(),
    loop_link.reset_loops(new_loops);
}
