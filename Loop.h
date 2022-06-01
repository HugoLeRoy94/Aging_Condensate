#ifndef Loop_h
#define Loop_h
class Loop{
public:
  Loop(std::array<double,3> R0,std::array<double,3> R1, double ell_loop);
  ~Loop();
  std::array<double,3> get_Rright() const;
  std::vector<std::array<double,3>> get_r() const;
  double get_ell() const;
private:
  std::array<double,3> Rright,Rleft;
  std::vector<std::array<double,3>> r;
  double ell;

  double compute_binding_rate(int r_index, double ell);
  void select_link_length(double& length, std::array<double,3>& r_selected);
  void generate_binding_sites();
};
#endif
