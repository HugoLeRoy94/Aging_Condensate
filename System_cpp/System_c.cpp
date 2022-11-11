#include "Header.h"
extern "C"
{
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void* create_system(double ell_tot,double rho0,double temperature,int seed,bool adjust)
  {
    return new System(ell_tot, rho0, temperature,seed,adjust);
  }
  
  void* CopySystem(void* ptr)
  {
    System* system = reinterpret_cast<System* >(ptr);
    return new(std::nothrow) System(*system);
  }
  
  void get_R(void* ptr, double* R, int size){
    System* system = reinterpret_cast<System*>(ptr);
    system->get_R(R,size);
  }
  
  void get_ell(void* ptr, double* ell, int size){
    System* system = reinterpret_cast<System*>(ptr);
    system->get_ell(ell,size);
  }
  
  double get_F(void* ptr){
    System* system = reinterpret_cast<System*>(ptr);
    return system->get_F();
  }
  
  int get_N_strand(void* ptr){
    System* system = reinterpret_cast<System*>(ptr);
    return system->get_N_strand();
  }
  
  double evolve(void* ptr, bool* bind)
  {
    System* system = reinterpret_cast<System*>(ptr);
    return system->evolve(bind);
  }
  
  int get_r_size(void* ptr){
    System* system = reinterpret_cast<System*>(ptr);
    return system->get_r_size();
  }
  
  int get_r_system_size(void* ptr){
    System* system = reinterpret_cast<System*>(ptr);
    return system->get_r_system_size();
  }
  
  void get_r(void* ptr, double* r, int size){
    System* system = reinterpret_cast<System*>(ptr);
    system->get_r(r,size);
  }
  
  void get_r_system(void* ptr, double* r, int size){
    System* system = reinterpret_cast<System*>(ptr);
    system->get_r_system(r,size);
  }
  
  void Print_Loop_positions(void* ptr){
    System* system = reinterpret_cast<System*>(ptr);
    system->Print_Loop_positions();
  }
  
  void print_random_stuff(void* ptr){
    System* system = reinterpret_cast<System*>(ptr);
    system->print_random_stuff();
  }
  
  void reset_crosslinkers(void* ptr){
    System* system = reinterpret_cast<System*>(ptr);
    system->reset_crosslinkers();
  }
  
  void get_ell_coordinates(void* ptr, double* ell_coordinates, int size){
    System* system = reinterpret_cast<System*>(ptr);
    system->get_ell_coordinates(ell_coordinates, size);
  }
  
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------

}
