#ifndef Function_h
#define Function_h
extern double Pi;
extern std::default_random_engine generator;
double get_square_diff(std::array<double,3> v1,std::array<double,3> v2);
double diff(std::array<double,3> v1,std::array<double,3> v2);
std::array<double,3> random_in_ellipse(double a,double b,double c,double xg,double yg,double zg);
#endif
