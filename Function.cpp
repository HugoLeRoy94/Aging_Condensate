#include "Header.h"
using namespace std;
double get_square(array<double,3> v1){
  double res;
  for(auto it:v1){
    res+=it*it;
  }
  return res;
};
