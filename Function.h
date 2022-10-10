#ifndef Function_h
#define Function_h
extern double Pi;
//extern std::default_random_engine generator;
extern std::mt19937 generator;
extern std::array<double,3> anchor;
double get_square_diff(std::array<double,3> v1,std::array<double,3> v2);
double diff(std::array<double,3> v1,std::array<double,3> v2);
struct {
        bool operator()(std::array<double,3> a, std::array<double,3> b) const { return a[0] < b[0]; }
    } customLess;
struct {
        bool operator()(std::pair<double,double> a, std::pair<double,double> b) const { return a.second < b.second; }
    } customLess2;

template <typename T0, typename T1, typename T2, typename T>
class map3d {
    std::map<T0,std::map<T1,std::map<T2,T>>> m;
    T0 size_x; // use to optimize the way we arange things.
    T1 size_y;
    T2 size_z;
public:
    map3d() = default;
    map3d(map3d const&) = default;
    //void set_sizes(T0 x,T1 y,T2 z){size_x = x;size_y=y;size_z=z;}
    constexpr T& operator()(T0 x, T1 y, T2 z) { // C++23 required
        return m[x][y][z];
    }
  
    void print()const{for(auto& it : m){for(auto& it2 : it.second){for(auto& it3 : it2.second){std::cout<<it3.second[0]<<" "<<it3.second[1]<<" "<<it3.second[2]<<std::endl;}}}}
    auto& underlying_array() const { return m; }
    std::vector<T*> slice(T0 key_0_min, T0 key_0_max, 
                        T1 key_1_min,T1 key_1_max,
                        T2 key_2_min, T2 key_2_max)
    {
      /*This function returns the list of key that are between 
      [a,b] for the first key, [c,d] for the second, [e,f] for
       the third*/
      std::vector<T*> res;
      /*std::map<T0,std::map<T1,std::map<T2,T>>>::iterator*/ 
      auto key_0_min_it(m.lower_bound(key_0_min));
      //std::map<T0,std::map<T1,std::map<T2,T>>>::iterator 
      auto key_0_max_it(m.upper_bound(key_0_max));
      for(auto& it = key_0_min_it;it!=key_0_max_it;++it)
      {
        //std::map<<T1,std::map<T2,T>>::iterator 
        auto key_1_min_it(m[it->first].lower_bound(key_1_min)); 
        //std::map<<T1,std::map<T2,T>>::iterator 
        auto key_1_max_it(m[it->first].upper_bound(key_1_max));
        for(auto& it2 = key_1_min_it; it2!=key_1_max_it;++it2)
        {
          //std::map<T2,T>::iterator 
          auto key_2_min_it(m[it->first][it2->first].lower_bound(key_2_min));
          //std::map<T2,T>>::iterator 
          auto key_2_max_it(m[it->first][it2->first].upper_bound(key_2_max));
          for(auto& it3 = key_2_min_it; it3!=key_2_max_it;++it3)
          {
            res.push_back(&(m[it->first][it2->first][it3->first]));
          }
        }
      }
      return res;
    }
    void cut_slice(T0 key_0_min, T0 key_0_max, 
                        T1 key_1_min,T1 key_1_max,
                        T2 key_2_min, T2 key_2_max)
    {
      auto key_0_min_it(m.lower_bound(key_0_min));
      //std::map<T0,std::map<T1,std::map<T2,T>>>::iterator 
      auto key_0_max_it(m.upper_bound(key_0_max));
      for(auto& it = key_0_min_it;it!=key_0_max_it;++it)
      {
        //std::map<<T1,std::map<T2,T>>::iterator 
        auto key_1_min_it(m[it->first].lower_bound(key_1_min)); 
        //std::map<<T1,std::map<T2,T>>::iterator 
        auto key_1_max_it(m[it->first].upper_bound(key_1_max));
        for(auto& it2 = key_1_min_it; it2!=key_1_max_it;++it2)
        {
          //std::map<T2,T>::iterator 
          auto key_2_min_it(m[it->first][it2->first].lower_bound(key_2_min));
          //std::map<T2,T>>::iterator 
          auto key_2_max_it(m[it->first][it2->first].upper_bound(key_2_max));
          for(auto& it3 = key_2_min_it; it3!=key_2_max_it;)
          {
            it3 = m[it->first][it2->first].erase(it3);
          }
        }
      }
    }
    void clear()
    {
      /*
      for(auto& it : m)
      {
        for(auto& it2 : it.second)
        {
          it2.second.clear();
        }
        it.second.clear();
      }
      */
      m.clear();
    }
    T is_in(T0 x_i,T1 y_i,T2 z_i)
    {
        try{m.at(x_i).at(y_i).at(z_i);}
        catch(std::out_of_range){throw;}
    }
    void add(T0 key1, T1 key2, T2 key3, T value)
    {
        m[key1][key2][key3] = value;
    }
    T* add_return_address(T0 key1, T1 key2, T2 key3, T value)
    {
      m[key1][key2][key3] = value;
      return &m[key1][key2][key3];
    }
    void remove(T0 key1, T1 key2, T2 key3)
    {
        m[key1][key2].erase(key3);
    }
    int get_number_of_elements() const
    {
      int size(0);
      for(auto& it : m){
        for(auto& it2 : it.second){
          size+=it2.second.size();
        }}
      return size;
    }
};
#endif
