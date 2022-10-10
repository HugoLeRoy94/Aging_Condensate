#include "Header.h"
using namespace std;
class CustomContainer
{
    private:
        map3d<double,double,double,std::array<double,3>> linkers;
        set<Strand*> loops;
        map<Strand*,set<array<double,3>*>> strand_to_linker;
        map<array<double,3>*,set<Strand*>> linker_to_strand; // keep them all all the time
    public:
        Insert_Strand(Strand* new_loop)
        {
            loops.insert(new_loop);
            for(auto& r new_loop->get_r())
            {
                strand_to_linker[new_loop].insert(r);
                linker_to_strand[r].insert(new_loop);
            }
        };
        Remove_Strand(Strand* loop_to_remove)
        {
            // Remove and delete
            loops.erase(loop_to_remove);
            for(auto& r : loop_to_remove->get_r())
            {
                strand_to_linker[r].erase(loop_to_remove);
            }
            linker_to_strand.erase(loop_to_remove);
            delete loop_to_remove;
        };
        Remove_linker(array<double,3>* linker_to_remove)
        {
            for(auto& strand : linker_to_strand[linker_to_remove])
            {
                strand_to_linker[strand].erase(linker_to_remove);
            }
            linker_to_strand.erase(linker_to_remove);
            linkers[linker_to_remove->at(0)][linker_to_remove->at(1)].erase(linker_to_remove->at(2));
        };
        Add_linker(array<double,3> linker_to_add)
        {
            linkers.add(linker_to_add[0],linker_to_add[1],linker_to_add[2],linker_to_add);

        };
};