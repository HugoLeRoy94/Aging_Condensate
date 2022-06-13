#ifndef Function_h
#define Function_h
extern double Pi;
//extern std::default_random_engine generator;
double get_square_diff(std::array<double,3> v1,std::array<double,3> v2);
double diff(std::array<double,3> v1,std::array<double,3> v2);
struct {
        bool operator()(std::array<double,3> a, std::array<double,3> b) const { return a[0] < b[0]; }
    } customLess;
struct {
        bool operator()(std::pair<double,double> a, std::pair<double,double> b) const { return a.second < b.second; }
    } customLess2;
#endif
