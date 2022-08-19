#ifndef Loop_h
#define Loop_h
class Loop{
public:
  //Loop(double ell_loop,double rho);
  Loop(std::array<double,3> R0,std::array<double,3> R1, double ell_loop, double rho,double dsl,bool rho_adjust);
  ~Loop();
  Loop(const Loop& loop);
  bool operator<(const Loop& otherloop) const;
  // main function for the evolution
  void select_link_length(double& length, std::array<double,3>& r_selected) const;
  // Accessors :
  std::array<double,3> get_Rright() const;
  std::array<double,3> get_Rleft() const;
  std::vector<std::array<double,3>> get_r() const;
  double get_ell() const;
  std::vector<std::vector<double>> get_rates() const;
  double get_total_binding_rates() const;
  double get_theta()const; // angle with x axis
  double get_phi() const; // angle with z axis
  std::array<double,3> get_Rg() const;
  double get_V()const;
  double get_S()const;
private:
  std::array<double,3> Rright,Rleft; // Position of the right and left anchor
  std::vector<std::array<double,3>> r; // position of all crosslinkers
  std::vector<std::vector<double>> rates,cum_rates; // rate of binding at any linkers for every length
  std::vector<double> sum_l_cum_rates; // rate of binding at any linkers
  double ell,V; // size of the polymer and volume it can occupy
  double a,b; // small and large axis
  double rho0; // volume fraction (initial of crosslinkers)
  double total_rates; // total binding rates to crosslinkers

  // returns a random position in an ellipse
  std::array<double,3> random_in_ellipse(double a,double b,double c,double xg,double yg,double zg);
  // number of configuration of a polymer bound in r1 and r2 and length ell
  double Omega(std::array<double,3> r1,std::array<double,3> r2,double ell) const;
  // inner function to compute all rates of the loop
  void compute_all_rates();
  // use random_in_ellipse to generate  a number of linkers
  void generate_binding_sites();

};
#endif
