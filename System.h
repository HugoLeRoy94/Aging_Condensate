#ifndef System_h
#define System_h
class System{
public:
  System(); // constructor
  System(double ell_tot,double D,double rho0,double kbT) // constructor with the Parameters
  ~System(); // destructor that have to delete all the loops
  double evolve(); // make the system evolve and return the time increment.
  int get_N() const; // return N the number of loop.
  void get_R(double* R, int size) const; // return the position of the anchored points
  void get_ell(double* ells, int size) const;// return the list of length of the loops.

private:
  std::vector<loop*> loops;
};
#endif
