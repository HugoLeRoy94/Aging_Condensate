#ifndef Loop_h
#define Loop_h
class Loop{
public:
  Loop(double R, double ell);
  ~Loop();
  std::array<double,3> get_R() const;
  std::vector<array<double,3>> get_r() const;
  double get_ell() const;
private:
  std::array<double,3> R;
  std::vector<array<double,3>> r;
  double ell_loop;

  double compute_binding_rate(int r_index, double ell);
  void select_link_length(double& length, array<double,3>& r_selected){

  }

};
#endif
