
#include <iostream>
#include <cmath>

const double EulerConstant = std::exp(1.0);

//function declarations

double function(double x);
void fixedPoint(double p0, double tol, double max_iter);

double function(double x){
    // g(x)=x-0.25f(x) as the fixed-point function where f(x)=e^x+2^(-x)+2cos(x)-6.
      double efun = pow(EulerConstant, x);
      double xpow = pow(2, -x);
      double cosin = 2*cos(x);
      double constant = -6;
     // double cube = pow(0, 3);
     // double cosin = cos(x);
      //double square = pow(x, 2);
      //double squareRoot = sqrt(x);
    //  double variable = 0;
     
      double fx = efun + xpow + cosin + constant; //e^x + 2^-x + 2cosx - 6
      //std::cout << result << std::endl;
      double gx = x - 0.25*(fx);
      return gx;
}

void fixedPoint(double p0, double tol, double max_iter){
    
    for (int i = 1; i < max_iter; i++){
        double p = function(p0); //computes pi
        std::cout << "P" << i << " = " << p << std::endl;
        double absoluteValue = abs(p-p0);
    
        if (absoluteValue < tol){ //successful procedure
            std::cout << "Finished in " << i << " iterations" << std::endl;
            return;
            
        }
        p0 = p; //update p0
    }
    std::cout << "Method Failed"<< std::endl;
}

int main() {
    double p0 = 1.0;
    double tol = pow(10.0, -5.0);
    
    //iterations
    double max_iter = 100.0;
    
    fixedPoint(p0, tol, max_iter);
}
