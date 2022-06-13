#ifndef System_h
#define System_h
class System{
public:
  System(); // constructor
  System(double ell_tot,double distance_anchor,double rho0,double temperature,int seed); // constructor with the Parameters
  ~System(); // destructor that have to delete all the loops
  double evolve(); // make the system evolve and return the time increment.
  int get_N() const; // return N the number of loop.
  void get_R(double* R, int size) const; // return the position of the anchored points
  void get_ell(double* ells, int size) const;// return the list of length of the loops.
  void get_r(double* r,int size) const;
  int get_r_size()const;
  void Print_Loop_positions();
private:
  double add_bond(std::vector<double>& cum_rates, std::vector<double>::iterator& rates_selec);
  double unbind_loop(std::set<Loop>::iterator& loop_select_left, std::set<Loop>::iterator& loop_select_right);
  double reset_crosslinkers();
  //std::mt19937_64 generator;
  std::uniform_int_distribution<int> distrib;
  std::set<Loop> loops;
  double ell,D,rho,kBT;
};
#endif
