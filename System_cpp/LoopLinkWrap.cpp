#include "Header.h"

using namespace std;
void LoopLinkWrap::set_loops(set<Strand*,LessLoop> new_loops)
{
    strands = new_loops;
}

void LoopLinkWrap::Insert_Strand(Strand* new_loop)
{
    strands.insert(new_loop);
}

void LoopLinkWrap::Remove_Strand(Strand* loop_to_remove)
{
    IF(true){cout<<"LoopLinkWrapper : start removing a strand"<<endl;}
    strands.erase(loop_to_remove);
    loop_to_remove->remove_from_linkers();
    delete loop_to_remove;
}

void LoopLinkWrap::set_linker_bounded(Linker* linker)
{
    linker->set_bounded();
}

void LoopLinkWrap::set_linker_free(Linker* linker)
{
    linker->set_free();
}

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

set<Strand*,LessLoop>::iterator LoopLinkWrap::get_loop(int distance)
{return next(strands.begin(),distance);}

void LoopLinkWrap::delete_pointers()
{
    IF(true){cout<<"delete_pointers"<<endl;}
    delete_loops();
    delete_linkers();
}

void LoopLinkWrap::delete_loops()
{for(auto& strand : strands){delete strand;}strands.clear();}

void LoopLinkWrap::delete_linkers()
{
    // Simply delete all the pointer in linkers and initialize a new empty map
    IF(true){cout<<"LoopLinkWrap : delete the linkers pointer"<<endl;}
    // doesn t delete them from the strands objects
    linkers.delete_pointers();
    map3dLink newlinkers; // make a new empty map3d
    linkers = newlinkers;
}

int LoopLinkWrap::get_strand_size() const {return strands.size();}

set<Strand*,LessLoop> LoopLinkWrap::get_loops() const {return strands;}

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

void LoopLinkWrap::actualize_vicinity(set<Strand*,LessLoop> to_remake)
{
    IF(true){cout<<"LoopLinkWrap : actualize_vicinity"<<endl;}
    for(auto& strand : to_remake)
    {
        try
        {
            strand->get_Rright();
            Insert_Strand(new Loop(*reinterpret_cast<Loop*>(strand),*this));
        }
        catch(out_of_range oor)
        {
            Insert_Strand(new Dangling(*reinterpret_cast<Dangling*>(strand),*this));
        }
    }
    for(auto& strand : to_remake){Remove_Strand(strand);}
    IF(true){cout<<"LoopLinkWrap : finished actualizing vicinity"<<endl;}
}

map<double,map<double,map<double,Linker*>>> LoopLinkWrap::get_linkers() const
{
    return linkers.underlying_array();
}
map3dLink LoopLinkWrap::get_linkers3d() const{ return linkers;}
int LoopLinkWrap::get_linker_size() const{
    return linkers.get_number_of_elements();
}