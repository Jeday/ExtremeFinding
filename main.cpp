#include <iostream>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>

#define EPS 0.00001
#define MIN_P -100
#define MAX_P 100

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
  };

};


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

std::pair<double,double> annealing(funct F){
    std::pair<double,double> point = {0.0,0.0};
	return { 0,0 };

}

double fRand(double fMin, double fMax)
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}

double offset(double dist) {
	return fRand(-dist, dist);
}

int iter = 100;

point genetic(bool(*comp)(const point &, const point &)) {
	std::vector<point> population(105);
	for (int i = 0; i < 105; ++i) {
		double x = fRand(MIN_P, MAX_P);
		double y = fRand(MIN_P, MAX_P);
		population[i] = point(x, y);
	}
	std::sort(population.begin(), population.end(), comp);
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
		std::sort(population.begin(), population.end(), comp);
		std::cout << population.rbegin()->x << "," << population.rbegin()->y << std::endl;
		--iter;

		if (population.rbegin()->distance_to(alpha_male) < EPS || iter == 0) {
			alpha_male = *population.rbegin();
			break;
		}
		else 
			alpha_male = *population.rbegin();

	} while (true);

	return alpha_male;
}



int main()
{
	std::srand(time(NULL));
    funct f(curve);
	point res_max = genetic(cmp);
	point res_min = genetic(negcmp);
	point res = abscmp(res_min, res_max) ? res_max : res_min;
	std::cout << res.x << ',' << res.y << std::endl;
	std::cout << curve(res.x, res.y) << std::endl;
	std::system("pause");
    return 0;
}
