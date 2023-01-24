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
LoopLinkWrap::LoopLinkWrap(){Nfree_linker=0.;}
LoopLinkWrap::~LoopLinkWrap()
{
    //delete_pointers();
}

void LoopLinkWrap::set_p_linkers(Strand* newly_created_strand)
{
  IF(true){cout<<"set_p_linkers"<<endl;}
  // get the volume limit
  array<double,3> main_ax,ctr_mass;
  double a,b;
  newly_created_strand->get_volume_limit(main_ax,ctr_mass,a,b);
  // select the linkers in the vicinity
  std::vector<Linker*> free_linkers,occ_linkers;
  IF(true){cout<<"Strand : slice"<<endl;}
  get_in_ellipse(ctr_mass,main_ax,a,b,free_linkers,occ_linkers);
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

void LoopLinkWrap::set_free(Linker* link)
{
    Nfree_linker++;
    link->set_free();
}

void LoopLinkWrap::set_occupied(Linker* link)
{
    Nfree_linker--;
    link->set_bounded();
}

void LoopLinkWrap::create_new_free_linker(double x,double y, double z)
{
    linkers[{x,y,z}] = new Linker({x,y,z});
    Nfree_linker++;
}

int LoopLinkWrap::get_N_free_linker() const
{
    return Nfree_linker;
}

void LoopLinkWrap::create_new_occupied_linker(double x, double y, double z)
{
    Linker* link = new Linker({x,y,z});
    link->set_bounded();
    linkers[{x,y,z}] = link;
}

void LoopLinkWrap::delete_linkers()
{
    // Simply delete all the pointer in linkers and initialize a new empty map
    IF(true){cout<<"LoopLinkWrap : delete the linkers pointer"<<endl;}
    // doesn t delete them from the strands objects
    Nfree_linker = 0.;
    for(auto& linker : linkers){
        delete linker.second;
    }
    map<array<double,3>,Linker*> newlinkers; // make a new empty map3d
    linkers = newlinkers;
}

map<array<double,3>,Linker*> const &LoopLinkWrap::get_linkers() const
{
    return linkers;
}

Linker* LoopLinkWrap::get_random_free_linker() const
{
    int step(0);
    while(step<pow(10.,17.))
    {
        auto item(linkers.begin());
        std::uniform_int_distribution<int> distribution(0,linkers.size()-1);
        std::advance(item, distribution(generator));
        if((*item).second->is_free())
        {
            return (*item).second;
        }
    }
    std::cout<<"no free linker to draw Nfree_linker certainly wrong"<<std::endl;
    throw std::out_of_range("Nfree_wrong ?");
}

int LoopLinkWrap::get_linker_size() const{
    return linkers.size();
}

//void LoopLinkWrap::add_linker(Linker* to_add){
//    linkers.add(to_add->r()[0],
//                to_add->r()[1],
//                to_add->r()[2],
//                to_add);
//}
Linker* LoopLinkWrap::diffuse_random_free_linker(){
    //cout<<"looplink::diffuse_random_free_linker : Nfree = "<<Nfree_linker<<endl;
    //cout<<"looplink::diffuse_random_free_linker : select a linker"<<endl;

    Linker* random_link(get_random_free_linker());
    //cout<<"looplink::diffuse_random_free_linker : remove it from the map"<<endl;
    linkers.erase({random_link->r()[0],random_link->r()[1],random_link->r()[2]});
    //cout<<"looplink::diffuse_random_free_linker : move the linker"<<endl;
    random_link->diffuse();
    //cout<<"looplink::diffuse_random_free_linker : add the linker into the map with new key"<<endl;
    linkers[{random_link->r()[0],random_link->r()[1],random_link->r()[2]}] = random_link;
    //cout<<"looplink::diffuse_random_free_linker : return the linker"<<endl;
    return random_link;
}
/*
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
*/
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

void LoopLinkWrap::get_in_ellipse(  std::array<double,3> ctr_mass,
                                    std::array<double,3> main_ax,
                                    double a,
                                    double b,
                                    std::vector<Linker*>& free_linkers,
                                    std::vector<Linker*>& occ_linkers)const
{
    for(auto& linker: linkers)
    {  
        array<double,3> r(rotate_point(linker.second->r(),main_ax,ctr_mass));
        if(pow(r[0]/a,2)+pow(r[1]/b,2)+pow(r[2]/b,2) <= 1)
        {
            if(linker.second->is_free()){
                free_linkers.push_back(linker.second);}
            else{
                occ_linkers.push_back(linker.second);}
        }
    }
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