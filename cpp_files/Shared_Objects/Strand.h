#ifndef Strand_h
#define Strand_h
class Accessor;
class Strand
{
    //friend Strand* LoopLinkWrap::Create_Strand(const Strand& new_strand);
    friend class Accessor;
    public:
        // operators
        Strand();
        Strand(Linker* R0,
                double ell_coordinate_0,
                double rho,
                bool sliding);
        virtual ~Strand();
        Strand(const Strand& strand, Linker* new_left_linker);
        Strand(const Strand& strand);
        bool operator<(const Strand &otherstrand) const;
        // do the moves
        virtual std::unique_ptr<Strand> unbind_from(Strand* left_strand) const=0;
        virtual std::pair<std::unique_ptr<Strand>,std::unique_ptr<Strand>> bind() const = 0;
        virtual std::unique_ptr<Strand> do_slide(double dl,bool right) const = 0;
        // Accessors :
        Linker* get_Rleft() const;
        virtual Linker* get_Rright() const =0;
        std::vector<Linker*> get_r() const;
        std::vector<Linker*> get_occ_r() const;
        double get_ell() const;
        double get_ell_coordinate_0() const;
        virtual double get_ell_coordinate_1() const = 0; 
        double get_total_binding_rates() const;
        double get_V() const;
        virtual double get_S(double dl=0) const =0;
        // surrounding linkers interaction
        void set_linkers(std::vector<Linker*> new_free_linkers,std::vector<Linker*> occupied_linkers);
        virtual void get_volume_limit(std::array<double,3>& main_ax, std::array<double,3>& ctr_mass,double& a, double& b) const =0;
        void remove_from_linkers();
        virtual void Check_integrity() const;
    protected:
        // compute the rates
        //void compute_all_rates();
        void compute_total_rates();
        void compute_cum_rates(std::vector<double>& sum_l_cum_rates, 
                               std::vector<std::vector<double>>& cum_rates) const;

        void select_link_length(double &length, Linker*& r_selected) const;
        virtual Strand* clone() const = 0;
        virtual std::array<double,3> random_in_volume() = 0;
        Linker* Rleft;               // Position of the right and left anchor
        std::vector<Linker*> free_linkers,occ_linkers;     // position of all crosslinkers
        //std::vector<std::vector<double>> rates, cum_rates; // rate of binding at any linkers for every length
        //std::vector<double> sum_l_cum_rates;               // rate of binding at any linkers
        double ell,V,xg,yg,zg;                             // size of the polymer and volume it can occupy
        double ell_coordinate_0;         // curvilinear coordinate of the linkers along the polymer
        double rho0;                                       // volume fraction (initial of crosslinkers)
        double total_rates;                                // total binding rates to crosslinkers
        bool slide;
        double slide_left,slide_right;

        void generate_binding_sites(LoopLinkWrap& loop_link);
        // inner function to compute all rates of the loop
        // use random_in_ellipse to generate  a number of linkers
        virtual double compute_binding_rate(double li, Linker* rlinker) const =0;

    //friend class LoopLinkWrap;
};
#endif