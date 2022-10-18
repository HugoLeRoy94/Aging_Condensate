#include "Header.h"

using namespace std;

Linker::Linker(std::array<double,3> r_c){R = r_c;free=true;}

array<double,3> Linker::r() const{return R;}

void Linker::set_free(){free = true;}

void Linker::set_bounded(){free=false;}

bool Linker::is_free() const{return free;}

void Linker::add_strand(Strand* strand){strands.insert(strand);}

void Linker::remove_strand(Strand* strand){strands.erase(strand);}

set<Strand*,LessLoop> Linker::get_strands() const
{
    return strands;
}