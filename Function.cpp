#include "Header.h"
using namespace std;
double Pi = acos(-1);
default_random_engine generator;
double get_square_diff(array<double,3> v1,array<double,3> v2){
  double res;
  for(auto it:v1){
    res+=it*it;
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
array<double,3> random_in_ellipse(double a,double b,double c,double xg,double yg,double zg){
  /*
    draw three random number that  must lie within an ellipse of
    revolution of big axe a and two small axes b.
  */
  // draw a x,y,z x in [-a,a] and y,z in [-b,b]
  bool OUT(true);
  double x(0),y(0),z(0);
  while(OUT){
    uniform_real_distribution<double> distribution(-1,1); //doubles from -1 to 1
    /*
    int randmax(100000);
    int xint = rand() % 2*randmax;
    int yint = rand() % 2*randmax;
    int zint = rand() % 2*randmax;
    xint-= randmax;
    yint-= randmax;
    zint-= randmax;
    */

    x = distribution(generator)*a;
    y = distribution(generator)*b;
    z = distribution(generator)*c;

    /*x = (double)xint/(double)randmax*a;
    y = (double)yint/randmax*b;
    z = (double)zint/randmax*c;*/
    //cout<<"x,y,z"<<x<<","<<y<<","<<z<<endl;
    //cout<<pow(x/a,2)+pow(y/b,2)+pow(z/c,2)<<endl;
    if(pow(x/a,2)+pow(y/b,2)+pow(z/c,2)<=1){OUT=false;}
  }
  double theta(atan2(yg,xg)),phi(atan2(xg,zg)-Pi/2.);
  //cout<<"xg yg zg ="<<xg<<" "<<yg<<" "<<zg<<endl;
  //cout<< "theta phi ="<<theta<<" "<<phi<<endl;
  //cout<<x<<" "<<y<<" "<<z<<endl;
  //cout<<x<<" "<<y<<" "<<z<<endl;
  array<double,3> res{cos(phi)*(cos(theta)*x-sin(theta)*y+sin(phi)*z),
                      sin(theta)*x+cos(theta)*y,
                      -sin(phi)*(cos(theta)*x-sin(theta)*y)+cos(phi)*z};
  res[0]+=xg;
  res[1]+=yg;
  res[2]+=zg;
  //cout<<res[0]<<" "<<res[1]<<" "<<res[2]<<endl;
  return res;
}
