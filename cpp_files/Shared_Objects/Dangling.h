#ifndef Dangling_h
#define Dangling_h
class Dangling : public Strand
{
  /*
  Dangling is the extremity of the polymer. The polymer is always bound at its
  left and free at its right. This class is very similar to Loop, except that
  linkers are in a sphere.
  */
public:
  Dangling();
  Dangling(Linker* R0,
          double ell_0,  // ell_0 is the coordinate
          double ell_in, // this is the remaining length
          double rho,
          bool sliding);
  Dangling(const Dangling& dangling);
  Dangling(const Dangling& dangling,
            Linker* new_left_linker);
  double get_S(double dl=0)const override; // entropy of the polymer.

  void get_volume_limit(std::array<double,3>& main_ax, std::array<double,3>& ctr_mass,double& a, double& b) const override;
  Linker* get_Rright() const;
  std::unique_ptr<Strand> unbind_from(Strand* left_strand) const override;
  std::pair<std::unique_ptr<Strand>,std::unique_ptr<Strand>> bind() const override;
  std::unique_ptr<Strand> do_slide(double dl,bool right) const override;
  double get_ell_coordinate_1() const override;
private:
  Strand* clone() const override;
  double radius;
  // returns a random position in a sphere
  std::array<double,3> random_in_volume() override;
  // number of configuration of a polymer bound in r1 and length ell
  double Omega(double ell) const;
  // use random_in_sphere to generate  a number of linkers
   // build le vector p_linkers from the overall map of the gillespie:
  double compute_binding_rate(double li, Linker* rlinker) const override;
  
};
#endif
