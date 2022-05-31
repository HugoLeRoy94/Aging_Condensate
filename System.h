#ifndef System_h
#define System_h
class System{
public:
  System(); // constructor
  ~System(); // destructor
  double evolve(); // make the system evolve and return the time increment.
  int get_N() const; // return N the number of loop.
  void get_R(double* R, int size) const; // return the position of the anchored points
  void get_ell(double* ells, int size) const;// return the list of length of the loops.

private:
  std::vector<loop*> loops;
};
#endif
