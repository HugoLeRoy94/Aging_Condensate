#ifndef Loop_h
#define Loop_h
class Loop : public Strand
{
public:
  // Loop(double ell_loop,double rho);
  Loop(Linker* R0,
       Linker* R1,
       double ell_coordinate_0,
       double ell_coordinate_1,
       double rho,
       bool sliding);
  // main function for the evolution
  Loop(const Loop &loop);
  Loop(const Loop &loop,
            Linker* new_left_linker,
            Linker* new_rightlinker); // overloaded constructor with new linkers
  // Accessors :
  Linker* get_Rright() const;
  
  double get_ell_coordinate_1() const override;
  double get_theta() const; // angle with x axis
  double get_phi() const;   // angle with z axis
  std::array<double, 3> get_Rg() const;
  double get_S(double dl=0) const override;
  void get_volume_limit(std::array<double,3>& main_ax, std::array<double,3>& ctr_mass,double& a, double& b) const override;
  void Check_integrity() const;
  std::unique_ptr<Strand> unbind_from(Strand* left_strand) const override;
  std::pair<std::unique_ptr<Strand>,std::unique_ptr<Strand>> bind() const override;
  std::unique_ptr<Strand> do_slide(double dl,bool right) const override;
private:
  Strand* clone() const override;
  Linker* Rright;
  double unbound_term;
  double ell_coordinate_1;         // curvilinear coordinate of the linkers along the polymer
  double a,b;                                       // small and large axis
  // returns a random position in an ellipse
  std::array<double, 3> random_in_volume() override;
  // number of configuration of a polymer bound in r1 and r2 and length ell
  double Omega(Linker* r1, Linker* r2, double ell) const;
  // inner function to compute all rates of the loop
  double compute_binding_rate(double li, Linker* linker) const override;
  
};
#endif
