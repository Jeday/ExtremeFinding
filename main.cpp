#include <iostream>



class funct{
  public:
  double (*f)(double,double);
  funct(double (*_f)(double,double)){
      f = _f;
  }

  double operator()( double x, double  y){
      return f(x,y);
  }

};

double curve (double x,double y) {
    return 1/(1+x*x)+1/(1+y*y);
}


std::pair<double,double> annealing(funct F){
    std::pair<double,double> point = {0.0,0.0};


}


int main()
{
    funct f(curve);


    return 0;
}
