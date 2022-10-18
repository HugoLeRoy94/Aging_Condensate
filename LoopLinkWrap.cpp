#include "Header.h"

using namespace std;

/*void LoopLinkWrap::set_linker_to_strand()
{
    linker_to_strand.clear();
    for(auto& loop : loops)
    {
        for(auto& r : loop->get_r())
        {
            linker_to_strand[r].insert(loop);
        }
    }
    // need to also set the bound extremities to linker_to_strand
    // start by creating a map3d of the extremities:
    map3d<double,double,double,std::array<double,3>*> connected_linkers;
    for(auto& loop:loops)
    {
        array<double,3>* left = new array<double,3>(loop->get_Rleft());
        connected_linkers.add(left->at(0),left->at(1),left->at(2),left);
    }
    // now, for each loop, select the ones that are in the vincinity 
    // of a given loop
    for(auto& loop:loops)
    {
        double k0,k1,k2,k3,k4,k5;
        loop->get_volume_limit(k0,k1,k2,k3,k4,k5);
        vector<array<double,3>*> r_in_vicinity(connected_linkers.slice(k0,k1,k2,k3,k4,k5));
        for(auto& r : r_in_vicinity)
        {
            linker_to_strand[r].insert(loop);
        }
    }

}*/

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
    // tell all the concerned linkers that this strand no longer exists
    //for(auto& linker : loop_to_remove->get_r())
    //{linker->remove_strand(loop_to_remove);}
    //for(auto& linker : loop_to_remove->get_occ_r())
    //{linker->remove_strand(loop_to_remove);}
    // Remove and delete
    strands.erase(loop_to_remove);
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

set<Strand*,LessLoop>::iterator LoopLinkWrap::get_loop(int distance)
{return next(strands.begin(),distance);}

void LoopLinkWrap::delete_pointers()
{
    delete_loops();
    delete_linkers();
}

void LoopLinkWrap::delete_loops()
{for(auto& strand : strands){Remove_Strand(strand);}strands.clear();}

void LoopLinkWrap::delete_linkers()
{
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

int LoopLinkWrap::get_linker_size() const{
    return linkers.get_number_of_elements();
}