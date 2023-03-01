#include "../Shared_Objects/Header.h"
#include "Gillespie.h"

using namespace std;

Gillespie::Gillespie(double ell_tot, double rho0, double BindingEnergy,double k_diff, int seed,bool sliding, int Nlinker,int dimension) : distrib(1, 10000000)
{
    IF(true) { cout << "Gillespie : creator" << endl; }
    // srand(seed);
    generator.seed(seed);
    // ---------------------------------------------------------------------------
    // -----------------------constant of the simulation--------------------------
    ell = ell_tot;
    rho = rho0;
    slide=sliding;
    LoopLinkWrap::dimension=dimension;
    N_linker_max = Nlinker;
    Linker::counter = 0;
    binding_energy = BindingEnergy;
    kdiff = k_diff;
    loop_link.create_new_occupied_linker(0.,0.,0.);
    Linker* R0 = loop_link.get_linkers().at({0.,0.,0.});
    Dangling dummy_dangling(R0, 0., ell, rho,slide); // dummy dangling that helps generate crosslinkers but has none initially
    Strand* dummy_strand(loop_link.Create_Strand(dummy_dangling));
    // ---------------------------------------------------------------------------
    //-----------------------------initialize crosslinkers------------------------
    set<array<double,3>> linkers(generate_crosslinkers(0));
  
    for(auto& linker : linkers)
    {
      loop_link.create_new_free_linker(linker.at(0),linker.at(1),linker.at(2));// the linker is free by default
    }
    loop_link.Remove_Strand(dummy_strand);
    // ---------------------------------------------------------------------------
    //-----------------------------initialize dangling----------------------------
    IF(true){cout<< "Gillespie : create dangling" << endl;}
    loop_link.Create_Strand(Dangling(R0, 0., ell, rho,slide));
    //print_random_stuff();
    //for(auto& it : linker_to_strand){for(auto& it2 : it.second){cout<<it2->get_Rleft()[0]<<" "<<it2->get_Rleft()[1]<<" "<<it2->get_Rleft()[2]<<endl;}}
    IF(true) { cout << "Gillespie : created" << endl; }
}

Gillespie::~Gillespie()
{
  loop_link.delete_pointers();
  Linker::counter = 0;
}

void Gillespie::compute_cum_rates(vector<double>& cum_rates) const
{
  IF(true) { cout << "Gillespie : Start computing the cumulative probability array" << endl; }
  cum_rates[0] = (loop_link.get_strand_size()-1)  * exp(binding_energy);
  cum_rates[1] = kdiff*loop_link.get_N_free_linker()+cum_rates[0];
  int n(2);
  for (auto &it : loop_link.get_strands())
  {
    cum_rates[n] = cum_rates[n - 1] + it->get_total_binding_rates();
    n++;
  }
  if(slide)
  {
  // iterate over each loop and the next one (skip the last)
  //-------------------------------------------------------------------
  for(set<Strand*>::iterator it = loop_link.get_strands().begin();
      it!=prev(loop_link.get_strands().end());
      it++)
  {    

    auto next_strand = next(it);
    cum_rates[n] = cum_rates[n - 1] + get_slide_rate(*it,*next_strand,1)+//slide right
                                      get_slide_rate(*it,*next_strand,-1);//slide left
    n++;
  }
  }
  //-------------------------------------------------------------------
  if (cum_rates.back() == 0)
  {
    IF(true){
    cout<<"cumulative rate are 0, let's output the number of linkers in the loops :" <<endl; 
    for(auto& loop : loop_link.get_strands())
    {cout<<loop->get_r().size()<<endl;}}
    throw invalid_argument("cum_rates.back()=0 no available process"); 
  }
  //for(auto& rate : cum_rates){cout<<rate<<endl;}
}

int Gillespie::pick_random_process(vector<double>& cum_rates) const
{
   IF(true) { cout << "Gillespie : draw a random process" << endl; }
  uniform_real_distribution<double> distribution(0, cum_rates.back());
  double pick_rate = distribution(generator);
  // becareful : if the rate_selec number is higher than cum_rates.back()  lower_bound returns cum_rates.back()
  vector<double>::iterator rate_selec = lower_bound(cum_rates.begin(), cum_rates.end(), pick_rate);
  //cout<<"size of the cumulative rate array, and the distance from the beginning  of the rate selected "<<cum_rates.size()<<" "<<distance(cum_rates.begin(),rate_selec)<<endl;
  if(rate_selec==cum_rates.end()){
    cout<<"pick_rate is too high"<<pick_rate<<endl;
    for(auto& rate : cum_rates){cout<<rate<<endl;}
    exit(0);
    }
  return distance(cum_rates.begin(),rate_selec);
}

double Gillespie::evolve(int *bind)
{
  IF(true){cout<<"-------------------------------------------------"<<endl;}
  IF(true) { cout << "Gillespie : start evolve" << endl; }
  IF(true){cout<<"-------------------------------------------------"<<endl;}
  //IF(true){check_loops_integrity();}
  // -----------------------------------------------------------------------------
  // -----------------------------------------------------------------------------
  // compute the cumulative transition rates for each loop
  vector<double> cum_rates;
  if(slide){
  // the +1 is for removing a bond
  // the -1 is for the (0,0,0) linker that cannot be slide
  // +1 is for diffusion of free linkers
    cum_rates.resize(loop_link.get_strand_size()*2 + 1-1+1,0);}
  else{
    // the +1 is for removing a bond
    // the +1 is for diffusion
    cum_rates.resize(loop_link.get_strand_size() + 1 +1 ,0);}
  // -----------------------------------------------------------------------------
  // -----------------------------------------------------------------------------
  try{compute_cum_rates(cum_rates);}
  catch(invalid_argument& e){return 0.;}
  // -----------------------------------------------------------------------------
  // -----------------------------------------------------------------------------
  // pick a random process
  int rate_select(pick_random_process(cum_rates));
  //cout<<"index of the rate selected = "<<rate_select<<endl;
  // Exectute the process
  if (rate_select == 0)
  {
        // Unbind a loop
    IF(true) { cout << "Gillespie : remove a bond" << endl; }
    // unbind
    unbind_random_loop();
    *bind = 0;
  }
  else if(rate_select == 1)
  {
    IF(true){cout<<"Gillespie : move the linkers"<<endl;}
    move_random_free_linkers();
    *bind = 1;
  }
  else if(rate_select>=loop_link.get_strand_size()+2)
  {
    // slide
    IF(true){cout<<"Gillespie : slide a bond"<<endl;}
    int loop_index_left(rate_select-loop_link.get_strand_size()-2);
    slide_bond(loop_index_left);
    *bind=2;
  }
  else
  {
    // add a linker to a strand or slide it
    IF(true) { cout << "Gillespie : add a bond" << endl; }
    int loop_index(rate_select-2); // rate_select = 1 => first bond

    add_bond(loop_index);
    *bind = 3;
  }
  // remake the strands whose rates have been modified by the event.
  //IF(true){check_loops_integrity();}
  IF(true){cout<<"output the value of the rates"<<endl;}
  IF(true){for(auto& rates : cum_rates){cout<<rates<<endl;}}
  return draw_time(cum_rates.back());
}

void Gillespie::add_bond(int loop_index)
{
    IF(true) { cout << "Gillespie : select the associated loop" << endl;}
    // cum_rates.begin() is unbinding
    set<Strand*,LessLoop>::iterator loop_selec(loop_link.get_strand(loop_index));
    IF(true){cout<<"the affected loop id is : "<<*loop_selec<<endl;}
    // ask the loop for a random linker, and a random length
    IF(true) { cout << "Gillespie : select a length and a r" << endl; }
    // select a link to bind to and ask the selected loop to return two new loops.
    pair<unique_ptr<Strand>,unique_ptr<Strand>> new_strands((*loop_selec)->bind());
    //new_strands.first->get_Rright()->set_bounded(); // set the linker to bounded
    loop_link.set_occupied(new_strands.first->get_Rright());
    // Create the pointers in the wrapper
    IF(true){cout<< "Gillespie : add_bond : build the list of strands that are affected"<<endl;}
    Strand* strand_left = loop_link.Create_Strand(*new_strands.first);
    Strand* strand_right = loop_link.Create_Strand(*new_strands.second);
    IF(true) { cout << "Gillespie : add_bond : delete the old loop" << endl; }
    // delete the loop
    loop_link.Remove_Strand((*loop_selec));
    // actualize the rates of the vicinity of the move
    // 1) store the affected strands
    set<Strand*,LessLoop> strands_affected = strand_left->get_Rright()->get_strands();
    // 2) remove those that have just been created
    strands_affected.erase(strand_right);//do not remake the strand that has just been made
    strands_affected.erase(strand_left);
    // 3) actualize the remaining strands
    loop_link.remake_strands(strands_affected);
}

void Gillespie::unbind_random_loop()
{
  // select the index of the left bond to remove
  // last index is loop.size-1
  // the left bond index maximum is loop.size -2
  uniform_int_distribution<int> distribution(0, loop_link.get_strand_size() - 2);
  int index(distribution(generator));
  Strand* loop_selec_left(*loop_link.get_strand(index));
  Strand* loop_selec_right(*loop_link.get_strand(index+1));
  IF(true) { cout << "Gillespie : unbind loop from the loop_left : "<<loop_selec_left<<" and the right : "<<loop_selec_right << endl; }

  // set the linker that was bound to unbound
  //loop_selec_left->get_Rright()->set_free();
  loop_link.set_free(loop_selec_left->get_Rright());
  // create a new loop that is the combination of both inputs

  Strand* loop(loop_link.Create_Strand(*loop_selec_right->unbind_from(loop_selec_left)));
  IF(true){cout<<"Gillespie : remove the old strand from loop_link"<<endl;}
  // save the reference of the linker to actualize the linker's state
  Linker* freed((loop_selec_right)->get_Rleft());
  // delete the loop before actualize vicinity.
  loop_link.Remove_Strand(loop_selec_left);
  loop_link.Remove_Strand(loop_selec_right);
  IF(true){cout<<"Gillespie : unbind_loop actualize vicinity"<<endl;}
  // 1) access all the affected strands in the neighboring
  set<Strand*,LessLoop> strands_affected = freed->get_strands();
  // 2) remove those that have just been created
  strands_affected.erase(loop);
  // 3) recompute the rates
  loop_link.remake_strands(strands_affected);
}

void Gillespie::slide_bond(int left_loop_index)
{ 
  double dl(choose_dl(left_loop_index));
  IF(true){cout<<"slide a bond by a dl = "<<dl<<endl;}
  Strand* left_strand((*loop_link.get_strand(left_loop_index)));
  Strand* right_strand((*loop_link.get_strand(left_loop_index+1)));
  IF(true){
    cout<<"number of strands = "<<loop_link.get_strands().size()<<endl;
    cout<<"slide the linker between the index : "<<left_loop_index<<" and the index : "<<left_loop_index+1<<endl;
    cout<<"address of the strands selected "<<left_strand<<" "<<right_strand<<endl;
    cout<<"with the coordinates : "<<left_strand->get_ell_coordinate_0()<<" "<<right_strand->get_ell_coordinate_0()<<endl;
  }
  
  Strand* new_left_strand = loop_link.Create_Strand(*left_strand->do_slide(dl,true));
  Strand* new_right_strand = loop_link.Create_Strand(*right_strand->do_slide(dl,false));
  loop_link.Remove_Strand(left_strand);
  loop_link.Remove_Strand(right_strand);
}

void Gillespie::move_random_free_linkers()
{
  //LoopLinkWrap new_loop_link;
  // move the linkers
  IF(true){cout<<"Gillespie : Move_linkers : diffuse linkers"<<endl;}
  IF(true){cout<<"Gillespie : move_random_free_linker : select a Linker"<<endl;}
  Linker* moved_linker(loop_link.diffuse_random_free_linker());
  // 1) access all the affected strands in the neighboring
  IF(true){cout<<"Gillespie : move_random_free_linker : move the linker"<<endl;}
  set<Strand*,LessLoop> strands_affected = moved_linker->get_strands();    
  // 3) recompute the rates
  IF(true){cout<<"Gillespie : move_random_free_linker : remake the affected strands"<<endl;}
  loop_link.remake_strands(strands_affected);
  // check if the linkers still belong to a loop:
  // check if the linker belong somewhere
  if(moved_linker->get_strands().size()==0){
    //remake a linker in the vicinity of one of the strand
    reset_crosslinkers();
  }
}

double Gillespie::choose_dl(int left_loop_index)
{
  set<Strand*,LessLoop>::iterator left_loop(loop_link.get_strand(left_loop_index));
  set<Strand*,LessLoop>::iterator right_loop(loop_link.get_strand(left_loop_index+1));
  uniform_real_distribution<double> distribution(0,get_slide_rate(*left_loop,*right_loop,1)+
                                                   get_slide_rate(*left_loop,*right_loop,-1));
  double pick_rate = distribution(generator);
  if(pick_rate<=get_slide_rate(*left_loop,*right_loop,+1)){
    return +1;
  }
  else{ return -1;}
}

double Gillespie::draw_time(double rate) const
{
  uniform_real_distribution<double> distrib;
  double xi(distrib(generator));
  return -log(1 - xi) / rate;
}

void Gillespie::reset_crosslinkers()
{
  IF(true){cout<<"------------------------------------------------------"<<endl;}
  IF(true){cout<<"reset all the crosslinkers"<<endl;}
  IF(true){cout<<"------------------------------------------------------"<<endl;}
  IF(true){check_loops_integrity();}
  LoopLinkWrap new_loop_link;
  // create a new set of occupied crosslinkers
  set<array<double,3>> occ_linkers_to_remake;
  // save the occupied linkers that stay
  for(auto& strand : loop_link.get_strands()){
    occ_linkers_to_remake.insert(strand->get_Rleft()->r());
    }
  // remake them
  for(auto& linker : occ_linkers_to_remake){
    new_loop_link.create_new_occupied_linker(linker.at(0),linker.at(1),linker.at(2));
  }
  // delete the linkers before
  // generate a bunch of free linkers
  set<array<double,3>> free_linkers_to_remake(generate_crosslinkers(loop_link.get_N_free_linker()));
  IF(true){cout<<"Number of free linkers to remake : "<<free_linkers_to_remake.size()<<endl;}
  // new create the associated linkers. carefull, they are free !
  for(auto& linker : free_linkers_to_remake){
    new_loop_link.create_new_free_linker(linker.at(0),linker.at(1),linker.at(2));
  }
  // set all the crosslinker into the linker_to_strand
  // which also add the bound extremities
  reset_loops(new_loop_link);
  loop_link.delete_linkers();
  loop_link = new_loop_link;
  IF(true){check_loops_integrity();}
}
/*
set<array<double,3>> Gillespie::generate_crosslinkers(int N_linker_already){
  //This function generates crosslinker at random position within a sphere of given radius
  IF(true){cout<<"Gillespie : generate crosslinkers"<<endl;}
  // create a set with all the limits
  set<double> xminl,xmaxl,yminl,ymaxl,zminl,zmaxl;
  
  for(auto& strand : loop_link.get_strands()){
    double ximin,ximax,yimin,yimax,zimin,zimax;
    //try{cout<<loop->get_Rright()->r().at(0)<<" "<<loop->get_Rright()->r().at(1)<<" "<<loop->get_Rright()->r().at(2)<<endl;}
    //catch(out_of_range oor){cout<<endl;}
    strand->get_volume_limit(ximin,ximax,yimin,yimax,zimin,zimax);
    xminl.insert(ximin);xmaxl.insert(ximax);yminl.insert(yimin);ymaxl.insert(yimax);zminl.insert(zimin);zmaxl.insert(zimax);
  }
  // select the min and max value
  double  xmin(*(xminl.begin())),xmax(*(xmaxl.rbegin()));
  double ymin(*(yminl.begin())),ymax(*(ymaxl.rbegin()));
  double zmin(*(zminl.begin())),zmax(*(zmaxl.rbegin()));
  // compute the volume and draw a number of linkers
  double V((xmax-xmin)*(ymax-ymin)*(zmax-zmin));
  poisson_distribution<int> distribution(rho * V);
  int N_crosslinker;
  if(N_linker_max>0)
  {
    N_crosslinker = N_linker_max-Linker::counter;
    //if(Linker::counter<=N_linker_max){N_crosslinker= 1;}
    //else{N_crosslinker = 0;}
  }
  else
  {
    N_crosslinker = max(0,distribution(generator)-N_linker_already);
  }
  // add all the occupied linkers of the strands :
  set<array<double,3>> res; // the result is a set of coordinates
  // draw a set of random position and create a linker at this position
  for(int n =0; n<N_crosslinker;n++)
  {
    uniform_real_distribution<double> distribution(0, 1); // doubles from 0 to 1 included
    double x(distribution(generator) * (xmax-xmin)+xmin);
    double y(distribution(generator) * (ymax-ymin)+ymin);
    double z(distribution(generator) * (zmax-zmin)+zmin);
    res.insert({x,y,z});
  }
  return res;
}
*/

set<array<double,3>> Gillespie::generate_crosslinkers(bool remake){
// loop over all the strand, and generate crosslinkers within their ellipse.
  set<array<double,3>> res;
  int N_linker_to_make(0.);
  // if we generate a fixed number of linkers
  if(N_linker_max>0)
    {
      Linker::counter = loop_link.get_linker_size()-loop_link.get_N_free_linker();// number of occupied linkers      
      N_linker_to_make = N_linker_max-Linker::counter; // total number of linkers that has to be added
      //cout<<"counter : "<<Linker::counter<<endl;
      //cout<<"N_linker_max :"<<N_linker_max<<endl;      
      // compute the probability to place N linkers propto the volume of each strands
      double total_volume(0.);
      for(auto& strand : loop_link.get_strands()){total_volume+=strand->get_V();}
      for(auto& strand : loop_link.get_strands())
      {
        double a,b;
        array<double,3> ctr_mass,main_ax;
        // number of linker to add to this specific strand
        int Nlinkers_strand(round(strand->get_V()/total_volume * N_linker_to_make));
        //cout<<strand->get_V()/total_volume * N_linker_to_make<<endl;
        //cout<<Nlinkers_strand<<endl;
        strand->get_volume_limit(main_ax,ctr_mass,a,b);
        generate_point_in_ellipse(main_ax,ctr_mass,a,b,res,Nlinkers_strand);
      }
    }
  else
  {
  for(auto & strand : loop_link.get_strands())
  {
    double a,b;
    array<double,3> ctr_mass,main_ax;
    poisson_distribution<int> distribution(rho * strand->get_V());
    cout<<"volume of the strand "<<strand->get_V()<<endl;
    
    if(remake){N_linker_to_make =N_linker_to_make = max(0.,distribution(generator)-(double)strand->get_occ_r().size());}
    else{N_linker_to_make = max(0.,distribution(generator)-(double)strand->get_r().size()-(double)strand->get_occ_r().size());}
    
    strand->get_volume_limit(main_ax,ctr_mass,a,b);
    generate_point_in_ellipse(main_ax,ctr_mass,a,b,res,N_linker_to_make);
  }
  }
  // transform res depending on the dimension
  set<array<double,3>> dimensional_res;
  if (LoopLinkWrap::dimension == 2){for(auto& xyz: res){dimensional_res.insert({xyz[0],xyz[1],0.});}}
  else if(LoopLinkWrap::dimension==1){for(auto& xyz: res){dimensional_res.insert({xyz[0],0.,0.});}}
  else{return res;}
  return dimensional_res;
}
/*
set<array<double,3>> Gillespie::generate_crosslinkers(int N_linker_already){
  set<array<double,3>> res; // the result is a set of coordinates
  for(int n =1; n<rho*2*sqrt(ell);n++){
    double x(n/rho),y(0),z(0);
    res.insert({n/rho,0,0});
  }
  return res;
}
*/

double Gillespie::compute_slide_S(Strand* left_strand, Strand* right_strand,double dl) const
{
  // left strand is always a loop
  if(
    // if the slided linker touches its left neighbor
    left_strand->get_ell_coordinate_0()-(left_strand->get_ell_coordinate_1()+dl)<=1 or 
    // uf the slided linker touches its right neighbor
    (right_strand->get_ell_coordinate_0()+dl)-right_strand->get_ell_coordinate_1()<=1 or
    // if the left strand is not long enough anymore:
    left_strand->get_ell_coordinate_0()-(left_strand->get_ell_coordinate_1()+dl)<diff(left_strand->get_Rleft()->r(),left_strand->get_Rright()->r())
    )
  {
    throw invalid_argument("the two left linkers overlap");
  }
  try{right_strand->get_Rright();}
  catch(out_of_range oor){return left_strand->get_S(-dl)+right_strand->get_S(dl)-//slide right
                                  (left_strand->get_S(0)+right_strand->get_S(0));}
  // if the right strand is not long enough anymore:
  if((right_strand->get_ell_coordinate_0()+dl)-right_strand->get_ell_coordinate_1()<diff(right_strand->get_Rleft()->r(),right_strand->get_Rright()->r()))
  {return 0.;}
  return left_strand->get_S(-dl)+right_strand->get_S(dl)-//slide right
         (left_strand->get_S(0)+right_strand->get_S(0));

}

double Gillespie::get_slide_rate(Strand* left_strand, Strand* right_strand,double dl) const
{
  // This function simply returns the rate of sliding
  // if it is a forbiden move it returns 0;
  try{
    return exp(0.5*compute_slide_S(left_strand,right_strand,+1));
  }
  catch(invalid_argument ia){return 0.;}
}

void Gillespie::reset_loops(LoopLinkWrap& new_loop_link)
{
  // recreate all the loops while changing the extremities to bounded
  IF(true){cout<<"Gillespie : reset loops"<<endl;}
  for(auto& strand : loop_link.get_strands())
  {
      try
      {
          strand->get_Rright(); // check to know if it's dangling
          // access the new linker via the coordinate of the old one : should always be valid
          Linker* const new_linker_right = new_loop_link.get_linkers().at({strand->get_Rright()->r().at(0),
                                                                         strand->get_Rright()->r().at(1),
                                                                         strand->get_Rright()->r().at(2)});
          Linker* const new_linker_left = new_loop_link.get_linkers().at({strand->get_Rleft()->r().at(0),
                                                                        strand->get_Rleft()->r().at(1),
                                                                        strand->get_Rleft()->r().at(2)});
          
          // seems that the new_linkers are invalid ?
          Strand* newloop(new_loop_link.Create_Strand(Loop(*reinterpret_cast<Loop*>(strand),
                                                            new_linker_left,new_linker_right)));
      }
      catch(out_of_range oor)
      {
        // access the new linker via the coordinate of the old one : should always be valid          
          Linker* const new_linker_left = new_loop_link.get_linkers().at({strand->get_Rleft()->r().at(0),
                                                          strand->get_Rleft()->r().at(1),
                                                          strand->get_Rleft()->r().at(2)});
          
          Strand* newloop(new_loop_link.Create_Strand(Dangling(*reinterpret_cast<Dangling*>(strand),
                                                            new_linker_left)));
      }
  }
  IF(true){cout<<"Gillespie : copy finished : delete old strands"<<endl;}
  loop_link.delete_strands();
}
/*
set<Strand*,LessLoop> Gillespie::get_vicinity(Linker* modified_linker,set<Strand*,LessLoop> strand_created)
{
  set<Strand*,LessLoop>  modified(modified_linker->get_strands());//.begin(),modified_linker->get_strands().end());
  for(auto& it : modified_linker->get_strands()){modified.insert(it);}
  set<Strand*,LessLoop> not_to_modify(strand_created.begin(),strand_created.end());
  set<Strand*,LessLoop>result0;

  set_difference(modified.begin(),modified.end(),
                not_to_modify.begin(),not_to_modify.end(),
                inserter(result0,result0.end()));
  IF(true){cout<<"strands that are reset : "<<endl;for(auto& strand : result0){cout<<strand<<endl;}}
  return result0;
}
*/