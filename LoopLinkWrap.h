#ifndef LoopLinkWrap_h
#define LoopLinkWrap_h
class LoopLinkWrap
{
    private:
        map3dLink linkers;
        std::set<Strand*,LessLoop> strands;
    public:

        void set_loops(std::set<Strand*,LessLoop> new_loops);
       
        void Insert_Strand(Strand* new_loop);
        
        void Remove_Strand(Strand* loop_to_remove);
        
        void set_linker_bounded(Linker* linker);
        
        void set_linker_free(Linker* linker_to_add);
        
        void create_new_free_linker(double x,double y, double z);
        
        std::set<Strand*,LessLoop>::iterator get_loop(int distance);
        
        void delete_pointers();
        
        void delete_loops();
        
        void delete_linkers();
        
        int get_strand_size() const;

        int get_linker_size() const;

        std::set<Strand*,LessLoop> get_loops() const;
        
        std::map<double,std::map<double,std::map<double,Linker*>>>  get_linkers() const;

        void slice_free(double key_0_min, double key_0_max, 
                                        double key_1_min,double key_1_max,
                                        double key_2_min, double key_2_max,
                                        std::vector<Linker*>& free_linkers,
                                        std::vector<Linker*>& occ_linkers) const;
        
        void actualize_vicinity(std::set<Strand*,LessLoop> to_remake);


};
#endif