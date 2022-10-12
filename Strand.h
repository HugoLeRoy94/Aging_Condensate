#ifndef Strand_h
#define Strand_h
class Strand
{
    public:
        Strand();
        Strand(Linker* R0,
                double ell_coordinate_0,
                double rho,
                bool rho_adjust);
        virtual ~Strand();
        Strand(const Strand& strand);
        bool operator<(const Strand &otherstrand) const;
        // main function for the evolution
        void select_link_length(double &length, Linker*& r_selected) const;
        // Accessors :
        Linker* get_Rleft() const;
        virtual Linker* get_Rright() const =0;
        std::vector<Linker*> get_r() const;
        double get_ell() const;
        double get_ell_coordinate_0() const;
        std::vector<std::vector<double>> get_rates() const;
        double get_total_binding_rates() const;
        double get_V() const;
        virtual double get_S() const =0;
        void set_p_linkers(LoopLinkWrap& loop_link);
        virtual void get_volume_limit(double& key_0_min,double& key_0_max,
                                      double& key_1_min,double& key_1_max,
                                      double& key_2_min,double& key_2_max) const =0;
        virtual void Check_integrity() const;
    protected:
        virtual std::array<double,3> random_in_volume() = 0;
        Linker* Rleft;               // Position of the right and left anchor
        std::vector<Linker*> free_linkers,occ_linkers;     // position of all crosslinkers
        std::vector<std::vector<double>> rates, cum_rates; // rate of binding at any linkers for every length
        std::vector<double> sum_l_cum_rates;               // rate of binding at any linkers
        double ell, V;                                     // size of the polymer and volume it can occupy
        double ell_coordinate_0;         // curvilinear coordinate of the linkers along the polymer
        double rho0;                                       // volume fraction (initial of crosslinkers)
        double total_rates;                                // total binding rates to crosslinkers

        void generate_binding_sites(LoopLinkWrap& loop_link);
        // inner function to compute all rates of the loop
        void compute_all_rates();
        // use random_in_ellipse to generate  a number of linkers
        virtual double compute_rate(double li, Linker* rlinker)=0;
};
#endif