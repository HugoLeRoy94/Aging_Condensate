#ifndef Dangling_h
#define Dangling_h
class Dangling{
  /*
  Dangling is the extremity of the polymer. The polymer is always bound at its
  left and free at its right. This class is very similar to Loop, except that
  linkers are in a sphere.
  */
public:
  Dangling(std::array<double,3> R0,double ell, double rho,double dsl,bool rho_adjust);
  Dangling(const Dangling& dangling);
  ~Dangling();
  void select_link_length(double& length, std::array<double,3>& r_selected) const;
  // Accessors :
  std::array<double,3> get_Rleft() const;
  std::vector<std::array<double,3>> get_r() const;
  double get_ell() const;
  std::vector<std::vector<double>> get_rates() const;
  double get_total_binding_rates() const;
  double get_V()const;
  double get_S()const; // entropy of the polymer.

private:

  std::array<double,3> Rleft; // Position of the left anchor
  std::vector<std::array<double,3>> r; // position of all crosslinkers
  std::vector<std::vector<double>> rates,cum_rates; // rate of binding at any linkers for every length
  std::vector<double> sum_l_cum_rates; // rate of binding at any linkers
  double ell,V; // size of the polymer and volume it can occupy
  double rho0; // volume fraction (initial of crosslinkers)
  double total_rates; // total binding rates to crosslinkers
  double radius;

  // returns a random position in a sphere
  std::array<double,3> random_in_sphere(double xg,double yg,double zg);
  // number of configuration of a polymer bound in r1 and length ell
  double Omega(double ell) const;
  // inner function to compute all rates of the loop
  void compute_all_rates();
  // use random_in_sphere to generate  a number of linkers
  void generate_binding_sites();

};
#endif
