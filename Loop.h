#ifndef Loop_h
#define Loop_h
class Loop : public Strand
{
public:
  // Loop(double ell_loop,double rho);
  Loop(std::array<double, 3> R0,
       std::array<double, 3> R1,
       map3d<double,double,double,std::array<double,3>>& linkers,
       std::map<array<double,3>*,vector<Strand*>> linker_to_strand,
       double ell_coordinate_0,
       double ell_coordinate_1,
       double rho,
       bool rho_adjust);
  // main function for the evolution
  Loop(const Loop &loop,map3d<double,double,double,std::array<double,3>>& linkers);
  // Accessors :
  std::array<double, 3> get_Rright() const;
  
  double get_ell_coordinate_1() const;
  double get_theta() const; // angle with x axis
  double get_phi() const;   // angle with z axis
  std::array<double, 3> get_Rg() const;

  double get_S() const;
  void get_volume_limit(double& key_0_min,double& key_0_max,
                        double& key_1_min,double& key_1_max,
                        double& key_2_min,double& key_2_max) const override;

private:
  std::array<double, 3> Rright;
  double unbound_term;
  double ell_coordinate_1;         // curvilinear coordinate of the linkers along the polymer
  double a,b,c,xg,yg,zg;                                       // small and large axis
  // returns a random position in an ellipse
  std::array<double, 3> random_in_volume() override;
  // number of configuration of a polymer bound in r1 and r2 and length ell
  double Omega(std::array<double, 3> r1, std::array<double, 3> r2, double ell) const;
  // inner function to compute all rates of the loop
  double compute_rate(double li, std::array<double,3>* rlinker) override;
};
#endif
