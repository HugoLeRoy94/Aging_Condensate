#include "Header.h"
extern "C"
{
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  void* create_gillespie(double ell_tot,double rho0,double BindingEnergy,double k_diff,int seed,bool slide,int Nlinker,int dimension)
  {
    return new Gillespie(ell_tot, rho0, BindingEnergy,k_diff,seed,slide,Nlinker,dimension);
  }
  void* CopyGillespie(void* ptr)
  {
    Gillespie* gillespie = reinterpret_cast<Gillespie* >(ptr);
    return new(std::nothrow) Gillespie(*gillespie);
  }
  
  void get_R(void* ptr, double* R, int size){
    Gillespie* gillespie = reinterpret_cast<Gillespie*>(ptr);
    gillespie->get_R(R,size);
  }
  
  void get_ell(void* ptr, double* ell, int size){
    Gillespie* gillespie = reinterpret_cast<Gillespie*>(ptr);
    gillespie->get_ell(ell,size);
  }
  
  double get_F(void* ptr){
    Gillespie* gillespie = reinterpret_cast<Gillespie*>(ptr);
    return gillespie->get_F();
  }
  
  int get_N_strand(void* ptr){
    Gillespie* gillespie = reinterpret_cast<Gillespie*>(ptr);
    return gillespie->get_N_strand();
  }
  
  double evolve(void* ptr, int* bind)
  {
    Gillespie* gillespie = reinterpret_cast<Gillespie*>(ptr);
    return gillespie->evolve(bind);
  }
  void delete_Gillespie(void* ptr)
  {
    Gillespie* gillespie = reinterpret_cast<Gillespie*>(ptr);
    delete gillespie;
  }
  int get_r_size(void* ptr){
    Gillespie* gillespie = reinterpret_cast<Gillespie*>(ptr);
    return gillespie->get_r_size();
  }
  
  int get_r_gillespie_size(void* ptr){
    Gillespie* gillespie = reinterpret_cast<Gillespie*>(ptr);
    return gillespie->get_r_gillespie_size();
  }
  
  void get_r(void* ptr, double* r, int size){
    Gillespie* gillespie = reinterpret_cast<Gillespie*>(ptr);
    gillespie->get_r(r,size);
  }
  
  void get_r_gillespie(void* ptr, double* r, int size){
    Gillespie* gillespie = reinterpret_cast<Gillespie*>(ptr);
    gillespie->get_r_gillespie(r,size);
  }
  
  void Print_Loop_positions(void* ptr){
    Gillespie* gillespie = reinterpret_cast<Gillespie*>(ptr);
    gillespie->Print_Loop_positions();
  }
  
  void print_random_stuff(void* ptr){
    Gillespie* gillespie = reinterpret_cast<Gillespie*>(ptr);
    gillespie->print_random_stuff();
  }
  
  void reset_crosslinkers(void* ptr){
    Gillespie* gillespie = reinterpret_cast<Gillespie*>(ptr);
    gillespie->reset_crosslinkers();
  }
  
  void get_ell_coordinates(void* ptr, double* ell_coordinates, int size){
    Gillespie* gillespie = reinterpret_cast<Gillespie*>(ptr);
    gillespie->get_ell_coordinates(ell_coordinates, size);
  }
  
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------

}
