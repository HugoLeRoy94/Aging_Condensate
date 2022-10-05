#ifndef System_h
#define System_h
class System{
public:
  System(double ell_tot,double rho0,double temperature,int seed,bool rho_adjust); // constructor with the Parameters
  ~System(); // destructor that have to delete all the loops
  double evolve(bool* bind); // make the system evolve and return the time increment.
  int get_N() const; // return N the number of loop.
  void get_R(double* R, int size) const; // return the position of the anchored points
  void get_ell_coordinates(double* ell_coordinate,int size)const; // get the curvilinear coordinates of the links.
  void get_ell(double* ells, int size) const;// return the list of length of the loops.
  void get_r(double* r,int size) const;
  void get_r_system(double* r, int size)const;
  double get_S() const;
  double get_F() const;
  int get_r_size()const;
  int get_r_system_size() const;
  void Print_Loop_positions() const;
  void print_random_stuff() const;
  void reset_crosslinkers();
private:
  void add_bond(std::vector<double>& cum_rates, std::vector<double>::iterator& rates_selec);
  void add_bond_to_dangling();
  void unbind_loop(std::set<Loop*>::iterator& loop_select_left, std::set<Loop*>::iterator& loop_select_right);
  void unbind_extremity(std::set<Loop*>::iterator& loop_selec_left);
  void generate_crosslinkers();
  double draw_time(double rate) const;
  void check_loops_integrity() const;
  double get_key(std::array<double,3>* p);
  //std::mt19937_64 generator;
  std::uniform_int_distribution<int> distrib;
  std::set<Loop*> loops;
  double ell,D,rho,kBT;
  bool rho_adjust;
  Dangling dangling;
  map3d<double,double,double,std::array<double,3>> linkers;
  std::map<array<double,3>*,vector<Strand*>> linker_to_strand;
};
#endif
