#ifndef Linker_h
#define Linker_h
class Linker
{
    private:
        std::array<double,3> R;
        bool free;
        std::set<Strand*,LessLoop> strands;
    public:
        Linker(std::array<double,3> r_c);
        std::array<double,3> r() const;
        void set_free(); // also 
        void set_bounded();
        bool is_free() const;
        void add_strand(Strand* strand);
        void remove_strand(Strand* strand);
        std::set<Strand*,LessLoop> get_strands() const;
};
class map3dLink : public map3d<double,double,double,Linker*>
{
    public:
        void delete_pointers(){
            for(auto& slice1 : m){
            for(auto& slice2 : slice1.second){
                for(auto& value: slice2.second){
                    delete value.second;
                }
            }
            }}
        void slice_free(double key_0_min, double key_0_max, 
                            double key_1_min,double key_1_max,
                            double key_2_min, double key_2_max,
                            std::vector<Linker*>& free_linkers,
                            std::vector<Linker*>& occ_linkers) const
        {
            /*This function returns the list of key that are between 
            [a,b] for the first key, [c,d] for the second, [e,f] for
            the third*/
            /*std::map<T0,std::map<T1,std::map<T2,T>>>::iterator*/ 
            auto key_0_min_it(m.lower_bound(key_0_min));
            //std::map<T0,std::map<T1,std::map<T2,T>>>::iterator 
            auto key_0_max_it(m.upper_bound(key_0_max));
            for(auto& it = key_0_min_it;it!=key_0_max_it;++it)
            {
            //std::map<<T1,std::map<T2,T>>::iterator 
            auto key_1_min_it(m.at(it->first).lower_bound(key_1_min)); 
            //std::map<<T1,std::map<T2,T>>::iterator 
            auto key_1_max_it(m.at(it->first).upper_bound(key_1_max));
            for(auto& it2 = key_1_min_it; it2!=key_1_max_it;++it2)
            {
                //std::map<T2,T>::iterator 
                auto key_2_min_it(m.at(it->first).at(it2->first).lower_bound(key_2_min));
                //std::map<T2,T>>::iterator 
                auto key_2_max_it(m.at(it->first).at(it2->first).upper_bound(key_2_max));
                for(auto& it3 = key_2_min_it; it3!=key_2_max_it;++it3)
                {
                if(m.at(it->first).at(it2->first).at(it3->first)->is_free())
                {
                    free_linkers.push_back(m.at(it->first).at(it2->first).at(it3->first));
                }
                else{occ_linkers.push_back(m.at(it->first).at(it2->first).at(it3->first));}
                }
            }
            }
        }
};
#endif