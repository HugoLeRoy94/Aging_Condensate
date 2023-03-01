#include "Header.h"

using namespace std;

Linker::Linker(std::array<double,3> r_c){R = r_c;free=true; Linker::counter++;}
Linker::~Linker(){Linker::counter--;}

array<double,3> Linker::r() const{
    if(LoopLinkWrap::dimension == 3){
        return R;}
    else if(LoopLinkWrap::dimension == 2){
        return {R[0],R[1],0.};
    }
    else if(LoopLinkWrap::dimension == 1){
        return {R[0],0.,0.};
    }
    else{throw invalid_argument("invalid dimension value");}
}

void Linker::set_free(){free = true;}

void Linker::set_bounded(){free=false;}

bool Linker::is_free() const{return free;}

void Linker::add_strand(Strand* strand){strands.insert(strand);}

void Linker::remove_strand(Strand* strand){strands.erase(strand);}

set<Strand*,LessLoop> Linker::get_strands() const
{
    return strands;
}
void Linker::print_position(string end)const
{
    cout<<R[0]<<" "<<R[1]<<" "<<R[2]<<end;
}
void Linker::diffuse()
{
    // chose a direction to make the move.
    normal_distribution<double> distribution(0.,1.);
    double dx(distribution(generator));
    double dy(distribution(generator));
    double dz(distribution(generator));
    double norm = sqrt(dx*dx+dy*dy+dz*dz);
    R[0]+=dx/norm;
    R[1]+=dy/norm;
    R[2]+=dz/norm;

}