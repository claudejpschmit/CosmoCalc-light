#pragma once

/**
 * @file This File contains a template for an integration function as well
 * as a variety of numerical integration methods. To add additional 
 * integration methods, one simply has to add additional methods
 * classes below.
 * source: rosettacode.org/wiki/Numerical_integration#C.2B.2B
 */

/**
 * The integration function requires an integration method and returns
 * the integral of the function f between boundaries a & b.
 *
 * @param f is the function to be integrated. eg. double f(double x){...}
 * @param a is the lower bound of the integral.
 * @param b is the upper bound of the integral.
 * @param steps denotes the number of steps to be taken between a & b.
 * @param m is an integration method according to which the integral is
 *          evaluated.
 * @details An example to use this function would be:
 *          double s = integrate(f, 0.0, 1.0, 100, simpson());
 */

#include <iostream>
#include <cmath>
#include <vector>


using namespace std;

template<typename Method, typename Function, typename Double> 
double integrate(Function f, Double a, Double b, int steps, Method m)
{
    double s = 0;
    double h = (b-a)/steps;
    for (int i = 0; i < steps; ++i)
        s += m(f, a + h*i, h);
    return h*s;
}


/**
 * This is a very simple simpson integration function.
 */
template<typename Function>
double integrate_simps(Function f, double a, double b, int steps)
{
    int extrastep = 0;
    if (steps % 2 == 1){
        //cout << "ERROR in integrator" << endl;
        //return 0;
        extrastep = 1;

    }
    double stepsize = (b-a)/(double)(steps+extrastep);
    double res = f(a) + f(b);

    for (int n = 1; n < steps; ++n){
        res += 4* f(a+n*stepsize);
        ++n;
    }
    for (int n = 2; n < steps - 1; ++n){
        res += 2 * f(a+n*stepsize);
        ++n;
    }
    return res * stepsize / 3.0;
}

/**
 * This class contains the simpson integration method, any other methods
 * are to be added in the same manner.
 */
class simpson
{
    public:
        template<typename Function, typename Double>
            double operator()(Function f, Double x, Double h) const
            {
                return (f(x) + 4*f(x+h/2) + f(x+h))/6;
            }
};
