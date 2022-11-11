#include "Header.h"
using namespace std;
double Pi = acos(-1);
//default_random_engine generator;
mt19937 generator;
array<double,3> anchor;
double get_square_diff(array<double,3> v1,array<double,3> v2){
  double res(0);
  for(int n=0;n<3;n++){
    res+=pow(v1[n]-v2[n],2);
  }
  return res;
}
double diff(array<double,3> v1,array<double,3> v2){
  double res(0);
  for(int n=0;n<3;n++){
    res+=pow(v1[n]-v2[n],2);
  }
  return sqrt(res);
}