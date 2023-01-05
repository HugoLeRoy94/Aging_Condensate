#include "Header.h"

using namespace std;

void Accessor::compute_rates(Strand* strand_to_compute_rates)
{
    strand_to_compute_rates->compute_total_rates();
}
Strand* Accessor::clone(const Strand& strand_to_clone)
{
    return strand_to_clone.clone();
}

LoopLinkWrap::~LoopLinkWrap()
{
    //delete_pointers();
}

void LoopLinkWrap::set_p_linkers(Strand* newly_created_strand)
{
  IF(true){cout<<"set_p_linkers"<<endl;}
  // get the volume limit
  double key_0_min,key_0_max,key_1_min,key_1_max,key_2_min,key_2_max;
  newly_created_strand->get_volume_limit(key_0_min,key_0_max,key_1_min,key_1_max,key_2_min,key_2_max);
  // select the linkers in the vicinity
  vector<Linker*> free_linkers,occ_linkers;
  IF(true){cout<<"Strand : slice"<<endl;}
  slice_free(key_0_min,key_0_max,
             key_1_min,key_1_max,
             key_2_min,key_2_max,
             free_linkers,occ_linkers);
  // tell the linkers that this strand is around
  for(auto& linker : free_linkers){linker->add_strand(newly_created_strand);}
  for(auto& linker : occ_linkers){linker->add_strand(newly_created_strand);}
  newly_created_strand->set_linkers(free_linkers,occ_linkers);
}

/*
    (~_|_ _ _  _  _| _
    _) | | (_|| |(_|_\                                                              
*/

void LoopLinkWrap::reset_strands(set<Strand*,LessLoop> new_strands)
{
    strands = new_strands;
}

Strand* LoopLinkWrap::Create_Strand(const Strand& new_strand)
{
    //Strand* new_created_strand = new_strand.clone();
    Strand* new_created_strand = Accessor::clone(new_strand);
    strands.insert(new_created_strand);
    set_p_linkers(new_created_strand);
    //new_created_strand->compute_all_rates();
    Accessor::compute_rates(new_created_strand);
    return new_created_strand;
}

void LoopLinkWrap::Remove_Strand(Strand* strand_to_remove)
{
    IF(true){cout<<"LoopLinkWrapper : start removing a strand"<<endl;}
    strands.erase(strand_to_remove);
    strand_to_remove->remove_from_linkers();
    delete strand_to_remove;
}

set<Strand*,LessLoop>::iterator LoopLinkWrap::get_strand(int distance)
{return next(strands.begin(),distance);}

void LoopLinkWrap::delete_strands()
{for(auto& strand : strands){delete strand;}strands.clear();}

int LoopLinkWrap::get_strand_size() const {return strands.size();}

set<Strand*,LessLoop> const & LoopLinkWrap::get_strands() const {return strands;}


/*
    | . _ |  _  _ _
    |_|| ||<(/_| _\                    
*/



void LoopLinkWrap::create_new_free_linker(double x,double y, double z)
{
    linkers.add(x,y,z,new Linker({x,y,z}));
}

void LoopLinkWrap::create_new_occupied_linker(double x, double y, double z)
{
    Linker* link = new Linker({x,y,z});
    link->set_bounded();
    linkers.add(x,y,z,link);
}

void LoopLinkWrap::delete_linkers()
{
    // Simply delete all the pointer in linkers and initialize a new empty map
    IF(true){cout<<"LoopLinkWrap : delete the linkers pointer"<<endl;}
    // doesn t delete them from the strands objects
    linkers.delete_pointers();
    map3dLink newlinkers; // make a new empty map3d
    linkers = newlinkers;
}

map<double,map<double,map<double,Linker*>>> const &LoopLinkWrap::get_linkers() const
{
    return linkers.underlying_array();
}

map3dLink LoopLinkWrap::get_linkers3d() const { return linkers;}

int LoopLinkWrap::get_linker_size() const{
    return linkers.get_number_of_elements();
}

//void LoopLinkWrap::add_linker(Linker* to_add){
//    linkers.add(to_add->r()[0],
//                to_add->r()[1],
//                to_add->r()[2],
//                to_add);
//}

void LoopLinkWrap::diffuse_linkers(){
    // make a new map3d
    map3dLink new_map3d;
for(auto& slice1 : linkers.underlying_array()){
    for(auto& slice2 : slice1.second){
        for(auto& linker: slice2.second){
            if(linker.second->is_free()){
            linker.second->diffuse();
            new_map3d.add(linker.second->r()[0],
                          linker.second->r()[1],
                          linker.second->r()[2],
                          linker.second);}
            else{
            new_map3d.add(linker.second->r()[0],
                          linker.second->r()[1],
                          linker.second->r()[2],
                          linker.second);}
            }
        }
    }
    linkers = new_map3d;
}
/*
    |\/|. _ _ _ || _  _  _  _     _
    |  ||_\(_(/_||(_|| |(/_(_)|_|_\
*/
void LoopLinkWrap::delete_pointers()
{
    IF(true){cout<<"delete_pointers"<<endl;}
    delete_strands();
    delete_linkers();
}

void LoopLinkWrap::slice_free(double key_0_min, double key_0_max, 
                                        double key_1_min,double key_1_max,
                                        double key_2_min, double key_2_max,
                                        std::vector<Linker*>& free_linkers,
                                        std::vector<Linker*>& occ_linkers) const
{
    return linkers.slice_free(key_0_min,key_0_max,
                        key_1_min,key_1_max,
                        key_2_min,key_2_max,
                        free_linkers,
                        occ_linkers);
}

void LoopLinkWrap::remake_strands(set<Strand*,LessLoop> to_remake)
{
    // this function just remake the list of strands given.
    // Recreating has only one function : recomputing the rates.
    IF(true){cout<<"LoopLinkWrap : actualize_vicinity"<<endl;}
    for(auto& strand : to_remake)
    {
        Create_Strand(*strand);
    }
    for(auto& strand : to_remake){Remove_Strand(strand);}
    IF(true){cout<<"LoopLinkWrap : finished actualizing vicinity"<<endl;}
}