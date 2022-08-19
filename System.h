#ifndef System_h
#define System_h
class System{
public:
  System(); // constructor
  System(double ell_tot,double distance_anchor,double rho0,double temperature,int seed,bool rho_adjust); // constructor with the Parameters
  ~System(); // destructor that have to delete all the loops
  double evolve(bool* bind); // make the system evolve and return the time increment.
  int get_N() const; // return N the number of loop.
  void get_R(double* R, int size) const; // return the position of the anchored points
  void get_ell(double* ells, int size) const;// return the list of length of the loops.
  void get_r(double* r,int size) const;
  double get_S() const;
  double get_F() const;
  int get_r_size()const;
  void Print_Loop_positions() const;
private:
  void add_bond(std::vector<double>& cum_rates, std::vector<double>::iterator& rates_selec);
  void unbind_loop(std::set<Loop>::iterator& loop_select_left, std::set<Loop>::iterator& loop_select_right);
  double reset_crosslinkers();
  double draw_time(double rate) const;
  //std::mt19937_64 generator;
  std::uniform_int_distribution<int> distrib;
  std::set<Loop> loops;
  double ell,D,rho,kBT;
  bool rho_adjust;
};
#endif
