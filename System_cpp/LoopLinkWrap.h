#ifndef LoopLinkWrap_h
#define LoopLinkWrap_h
//class LoopLinkWrap;
class LoopLinkWrap
{
    private:
        map3dLink linkers;
        std::set<Strand*,LessLoop> strands;

        void set_p_linkers(Strand* newly_created_strand);
    public:
        ~LoopLinkWrap();
        /*
            (~_|_ _ _  _  _| _
            _) | | (_|| |(_|_\                                                              
        */
        void reset_strands(std::set<Strand*,LessLoop> new_strands);

        Strand* Create_Strand(const Strand& new_strand);

        void Remove_Strand(Strand* strand_to_remove);
        
        void delete_strands();

        std::set<Strand*,LessLoop>::iterator get_strand(int distance);

        std::set<Strand*,LessLoop> const& get_strands() const;

        int get_strand_size() const;
        
       /*
            | . _ |  _  _ _
            |_|| ||<(/_| _\                    
       */
                
        //void add_linker(Linker* to_add);

        void create_new_free_linker(double x,double y, double z);
        
        void create_new_occupied_linker(double x,double y,double z);
                
        void delete_linkers();
        
        int get_linker_size() const;
        
        std::map<double,std::map<double,std::map<double,Linker*>>>const &  get_linkers() const;
        
        map3dLink get_linkers3d() const;

        void diffuse_linkers();
        
        /*
        |\/|. _ _ _ || _  _  _  _     _
        |  ||_\(_(/_||(_|| |(/_(_)|_|_\
        */

        void slice_free(double key_0_min, double key_0_max, 
                                        double key_1_min,double key_1_max,
                                        double key_2_min, double key_2_max,
                                        std::vector<Linker*>& free_linkers,
                                        std::vector<Linker*>& occ_linkers) const;

        void remake_strands(std::set<Strand*,LessLoop> to_remake);

        void delete_pointers();
};
class Accessor
{
    private:
        static void compute_rates(Strand* strand_to_compute_rates);
        static Strand* clone(const Strand& strand_to_clone);
    friend Strand* LoopLinkWrap::Create_Strand(const Strand& new_strand);
};
#endif