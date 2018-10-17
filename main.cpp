#include <iostream>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <cstdint>
#include <random>

#define EPS 0.00001
#define MIN_P -100
#define MAX_P 100

std::default_random_engine generator;
std::gamma_distribution<double> distribution(2.0,2.0);

class point {
public:
	double x; 
	double y;

	point(double _x, double _y) {
		x = _x;
		y = _y;
	}

	point() {
		x = 0;
		y = 0;
	}

	double distance_to(point p) {
		return std::sqrt((x - p.x) * (x - p.x) + (y - p.y) * (y - p.y));
	}
};

class funct{
  public:
  double (*f)(double,double);
  funct(double (*_f)(double,double)){
      f = _f;
  }

  double operator()(double x, double  y){
      return f(x,y);
  }

  double operator()(point a) {
	  return f(a.x, a.y);
  }

};

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}


double curve (double x,double y) {
    return x*x + y*y;
}

bool cmp(const point& a, const point&b) {
    return curve(a.x, a.y) < curve(b.x, b.y);
}

bool negcmp(const point& a, const point&b) {
    return curve(a.x, a.y) > curve(b.x, b.y);
}

bool abscmp(const point& a, const point&b) {
    return std::abs(curve(a.x, a.y)) < std::abs(curve(b.x, b.y));
}


class genetics{
   int iter = 100;
   double (*crv)(double,double);

   double offset(double dist) {
       return fRand(-dist, dist);
   }

    bool abscmp(const point& a, const point&b) {
       return std::abs(crv(a.x, a.y)) < std::abs(crv(b.x, b.y));
   }





   point * genetic(bool f ) {
       std::vector<point> population(105);
       for (int i = 0; i < 105; ++i) {
           double x = fRand(MIN_P, MAX_P);
           double y = fRand(MIN_P, MAX_P);
           population[i] = point(x, y);
       }

       f ? std::sort(population.begin(), population.end(), [this](const point& a, const point&b)->bool{ return this->crv(a.x, a.y) < this->crv(b.x, b.y);}) :
           std::sort(population.begin(), population.end(), [this](const point& a, const point&b)->bool{ return this->crv(a.x, a.y) > this->crv(b.x, b.y);});
       point alpha_male = *population.rbegin();
       do {

           auto b = population.rbegin();
           std::vector<point> selection(15);
           for (int i = 0; i < 15; ++i) {
               selection[i] = *b;
               ++b;
           }

           int k = 0;
           for (int i = 0; i < 15; ++i) {
               for (int j = i + 1; j < 15; ++j) {
                   double dist = selection[i].distance_to(selection[j]);
                   population[k] = point((selection[i].x+selection[j].x)/2 + offset(dist)/4,(selection[i].y + selection[j].y)/2 + offset(dist)/4);
                   k++;
               }
           }
           f ? std::sort(population.begin(), population.end(), [this](const point& a, const point&b)->bool{ return this->crv(a.x, a.y) < this->crv(b.x, b.y);}) :
               std::sort(population.begin(), population.end(), [this](const point& a, const point&b)->bool{ return this->crv(a.x, a.y) > this->crv(b.x, b.y);});
          std::cout << population.rbegin()->x << "," << population.rbegin()->y << std::endl;
           --iter;

           if (population.rbegin()->distance_to(alpha_male) < EPS ) {
               alpha_male = *population.rbegin();
               break;
           }
           else if(iter <0){
               return nullptr;

           }
           else
               alpha_male = *population.rbegin();

       } while (true);

       return new point(alpha_male);
   }

public:

   genetics(double (*f)(double,double),int iterations = 100 ){
        crv = f;
        iter = iterations;
   }

   point * operator()(){
       point * res_max = genetic(true);
       point * res_min = genetic(false);
       if(res_max){
           if(res_min)
               res_max = abscmp(*res_min,*res_max) ? res_max :res_min;

       }
       else res_max = res_min;

       if(!res_max)
           return nullptr;
       return res_max;
   }

};


point offset(point p){
    double angle = fRand(0,2*3.14159265359);
    double dist = fRand(0,5);
    return point(p.x+dist*std::cos(angle),p.y+dist*std::sin(angle));

}


point annealing(double (*f)(double,double),double T = 25){\
    point p = point(0,0);
    double fp = f(p.x,p.y);
    while(T>0){
        point np = offset(p);
        double fnp = f(np.x,np.y);
        if(fp < fnp || fRand(0,1)<=std::exp((fnp-fp)/T)){
            p = np;
            fp = fnp;
        }
        t-=0.01;
    }
    return p;


}










int main()
{

	std::srand(time(NULL));
    genetics g(curve);
    point * res = g();
    std::cout << res->x << ',' << res->y << std::endl;
    std::cout << curve(res->x, res->y) << std::endl;
    point p = annealing(curve);

	std::system("pause");
    return 0;
}
