#ifndef Gillespie_h
#define Gillespie_h
class Gillespie{
public:
    Gillespie(double ell_tot,double rho0,double BindingEnergy,double k_diff,int seed,bool sliding=0,int Nlinker=0,int dimension=3); // constructor with the Parameters
    ~Gillespie(); // destructor that have to delete all the loops
    double evolve(int* bind); // make the gillespie evolve and return the time increment.
    void reset_crosslinkers();
    // -----------------------------------------------------------------------------
    // -----------------------------accessor----------------------------------------
    // -----------------------------------------------------------------------------
    int get_N_strand() const; // return N the number of loop.
    void get_R(double* R, int size) const; // return the position of the anchored points
    void get_ell_coordinates(double* ell_coordinate,int size)const; // get the curvilinear coordinates of the links.
    void get_ell(double* ells, int size) const;// return the list of length of the loops.
    void get_r(double* r,int size) const;
    void get_r_gillespie(double* r, int size)const;
    double get_S() const;
    double get_F() const;
    int get_r_size()const;
    int get_r_gillespie_size() const;
    void Print_Loop_positions() const;
    void print_random_stuff() const;
    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------

private:
    void compute_cum_rates(std::vector<double>& cum_rates) const;
    int pick_random_process(std::vector<double>& cum_rates) const;

    void add_bond(int loop_selected);
    void unbind_random_loop();
    void slide_bond(int left_loop);
    void move_random_free_linkers();

    double choose_dl(int left_loop_index);
    double compute_slide_S(Strand* left_strand, Strand* right_strand,double dl) const;
    double get_slide_rate(Strand* left_strand, Strand* right_strand,double dl) const;

    std::set<std::array<double,3>>  generate_crosslinkers(bool remake);
    void reset_loops(LoopLinkWrap& new_loop_link);
    double draw_time(double rate) const;

    void check_loops_integrity();
    //std::mt19937_64 generator;
    std::uniform_int_distribution<int> distrib;
    double ell,D,rho,binding_energy,kdiff;
    bool slide;
    int N_linker_max;    

    LoopLinkWrap loop_link;

};
#endif