#include "LevinIntegrator.hpp"
#include "armadillo"
#include <boost/math/special_functions/bessel.hpp>
#include <iostream>

using namespace arma;

Levin::Levin(double a, double b)
    :
        a(a),
        b(b)
{
    //d = (a+b)/2.0 + 0.00000000000001;
    d = (a+b)/2.0;
    Basis = new CosmoBasis();
}

Levin::~Levin()
{
    delete Basis;
}

double Levin::integrate_singleJ(double (*f)(double), double r, double  nu, int n)
{
    mat matrix, rhs, c;
    matrix = randu<mat>(2*n,2*n);
    rhs = randu<mat>(2*n,1);
 
    for (int i = 0; i < n; i++)
    {
        //double x = a + (i-1.0)*(b-a)/(double)(n-1.0);
        
        //Using Chebychev points
        double x = 0.5*(a+b) + 0.5*(b-a)*cos((double)(2*(i+1)-1)*pi/(double)(2.0*n));   
        
        rhs(i,0) = 0;
        rhs(i+n,0) = f(x); 

        for (int k = 0; k < n; k++)
        {
            matrix(i, k) = up(k,x) + (double)(nu-1) * u(k,x)/x;
            matrix(i, k + n) = r * u(k,x);
            matrix(i + n, k) = -r * u(k,x);
            matrix(i + n, k + n) = up(k,x) - (double)nu * u(k,x)/x;
        }
    }
    matrix = matrix.i();
    c = randu<mat>(2*n,1);
    c = matrix * rhs;

    double sum1 = 0;
    double sum2 = 0;
    double sum3 = 0;
    double sum4 = 0;
    
    for (int k = 0; k < n; k++)
    {
        sum1+=c(k,0) * u(k,b);
        sum2+=c(k,0) * u(k,a);
        sum3+=c(k+n,0) * u(k,b);
        sum4+=c(k+n,0) * u(k,a);
    }
    double j1 = bessel_J(nu-1,r*b);
    double j2 = bessel_J(nu-1,r*a);
    double j3 = bessel_J(nu,r*b);
    double j4 = bessel_J(nu,r*a);
    
    return sum1 * j1 - sum2 * j2 + sum3 * j3 - sum4 * j4;
}

double Levin::integrate_doubleJ(double (*f)(double), double r, double nu, int n)
{
    mat matrix, rhs, c;
    matrix = randu<mat>(3*n,3*n);
    rhs = randu<mat>(3*n,1);
   
    for (int i = 0; i < n; i++)
    {
        //double x = a + (i-1.0)*(b-a)/(double)(n-1.0);
        
        //Using Chebychev points
        double x = 0.5*(a+b) + 0.5*(b-a)*cos((double)(2*(i+1)-1)*pi/(double)(2.0*n));        
        
        rhs(i,0) = 0;
        rhs(i+n,0) = 0;
        rhs(i+2*n,0) = f(x); 

        for (int k = 0; k < n; k++)
        {
            matrix(i, k) = up(k,x) + 2*(nu-1.0) * u(k,x)/x;
            matrix(i, k + n) = r * u(k,x);
            matrix(i, k + 2*n) = 0;
            
            matrix(i + n, k) = -2*r * u(k,x);
            matrix(i + n, k + n) = up(k,x) - u(k,x)/x;
            matrix(i + n, k + 2*n) = 2*r * u(k,x);
            
            matrix(i + 2*n, k) = 0;
            matrix(i + 2*n, k + n) = -r * u(k,x);
            matrix(i + 2*n, k + 2*n) = up(k,x)-2.0*nu*u(k,x)/x;
        }
    }
    matrix = matrix.i();
    c = randu<mat>(3*n,1);
    c = matrix * rhs;

    double sum1 = 0;
    double sum2 = 0;
    double sum3 = 0;
    double sum4 = 0;
    double sum5 = 0;
    double sum6 = 0;

    for (int k = 0; k < n; k++)
    {
        sum1+=c(k,0) * u(k,b);
        sum2+=c(k,0) * u(k,a);
        sum3+=c(k+n,0) * u(k,b);
        sum4+=c(k+n,0) * u(k,a);
        sum5+=c(k+2*n,0) * u(k,b);
        sum6+=c(k+2*n,0) * u(k,a);
    }
    double j1 = bessel_J(nu-1,r*b);
    double j2 = bessel_J(nu-1,r*a);
    double j3 = bessel_J(nu,r*b);
    double j4 = bessel_J(nu,r*a);
    
    double result = sum1 * j1 * j1 - sum2 * j2 * j2 +\
                    sum3 * j1 * j3 - sum4 * j2 * j4 +\
                    sum5 * j3 * j3 - sum6 * j4 * j4;
    return result;
}

double Levin::integrate_2sphj_1r(double (*f)(double), double r, int l, int n)
{
    mat matrix, rhs, c;
    matrix = randu<mat>(3*n,3*n);
    rhs = randu<mat>(3*n,1);
   
    for (int i = 0; i < n; i++)
    {
        //double x = a + (i-1.0)*(b-a)/(double)(n-1.0);
        
        //Using Chebychev points
        double x = 0.5*(a+b) + 0.5*(b-a)*cos((double)(2*(i+1)-1)*pi/(double)(2.0*n));
        rhs(i,0) = f(x);
        rhs(i+n,0) = 0;
        rhs(i+2*n,0) = 0;
        
        for (int k = 0; k < n; k++)
        {
            matrix(i, k) = up(k,x) -2*(double)(l+1) * u(k,x)/x;
            matrix(i, k + n) = -r * u(k,x);
            matrix(i, k + 2*n) = 0;
            
            matrix(i + n, k) = 2*r * u(k,x);
            matrix(i + n, k + n) = up(k,x) - 2*u(k,x)/x;
            matrix(i + n, k + 2*n) = -2*r*u(k,x);
            
            matrix(i + 2*n, k) = 0;
            matrix(i + 2*n, k + n) = r*u(k,x);
            matrix(i + 2*n, k + 2*n) = up(k,x)+2.0*(l-1.0)*u(k,x)/x;
        }
    }
    matrix = matrix.i();
    c = randu<mat>(3*n,1);
    c = matrix * rhs;

    double sum1 = 0;
    double sum2 = 0;
    double sum3 = 0;
    double sum4 = 0;
    double sum5 = 0;
    double sum6 = 0;

    for (int k = 0; k < n; k++)
    {
        sum1+=c(k,0) * u(k,b);
        sum2+=c(k,0) * u(k,a);
        sum3+=c(k+n,0) * u(k,b);
        sum4+=c(k+n,0) * u(k,a);
        sum5+=c(k+2*n,0) * u(k,b);
        sum6+=c(k+2*n,0) * u(k,a);
    }
    double j1 = sph_bessel(l,r*b);
    double j2 = sph_bessel(l,r*a);
    double j3 = sph_bessel(l-1,r*b);
    double j4 = sph_bessel(l-1,r*a);

    double result = sum1 * j1 * j1 - sum2 * j2 * j2 +\
                    sum3 * j1 * j3 - sum4 * j2 * j4 +\
                    sum5 * j3 * j3 - sum6 * j4 * j4;

    return result;
}

double Levin::integrate_2sphj_2r(double (*f)(double), double r1, double r2, int l, int n)
{
    mat matrix, rhs, c;
    matrix = randu<mat>(4*n,4*n);
    rhs = randu<mat>(4*n,1);
   
    for (int i = 0; i < n; i++)
    {
        //double x = a + (i-1.0)*(b-a)/(double)(n-1.0);
        
        //Using Chebychev points
        double x = 0.5*(a+b) + 0.5*(b-a)*cos((double)(2*(i+1)-1)*pi/(double)(2.0*n));
        rhs(i,0) = f(x);
        rhs(i+n,0) = 0;
        rhs(i+2*n,0) = 0;
        rhs(i+3*n,0) = 0;

        for (int k = 0; k < n; k++)
        {
            matrix(i, k) = up(k,x) -2*(double)(l+1) * u(k,x)/x;
            matrix(i, k + n) = -r1 * u(k,x);
            matrix(i, k + 2*n) = -r2 * u(k,x);
            matrix(i, k + 3*n) = 0;
            
            matrix(i + n, k) = r1 * u(k,x);
            matrix(i + n, k + n) = up(k,x) - 2*u(k,x)/x;
            matrix(i + n, k + 2*n) = 0;
            matrix(i + n, k + 3*n) = -r2*u(k,x);
            
            matrix(i + 2*n, k) = r2 * u(k,x);
            matrix(i + 2*n, k + n) = 0;
            matrix(i + 2*n, k + 2*n) = up(k,x)-2*u(k,x)/x;
            matrix(i + 2*n, k + 3*n) = -r1*u(k,x);
        
            matrix(i + 3*n, k) = 0;
            matrix(i + 3*n, k + n) = r2 * u(k,x);
            matrix(i + 3*n, k + 2*n) = r1 * u(k,x);
            matrix(i + 3*n, k + 3*n) = up(k,x) + 2 *(double)(l-1)*u(k,x)/x;
        }
    }
    matrix = matrix.i();
    c = randu<mat>(4*n,1);
    c = matrix * rhs;

    double sum1 = 0;
    double sum2 = 0;
    double sum3 = 0;
    double sum4 = 0;
    double sum5 = 0;
    double sum6 = 0;
    double sum7 = 0;
    double sum8 = 0;

    for (int k = 0; k < n; k++)
    {
        sum1+=c(k,0) * u(k,b);
        sum2+=c(k,0) * u(k,a);
        sum3+=c(k+n,0) * u(k,b);
        sum4+=c(k+n,0) * u(k,a);
        sum5+=c(k+2*n,0) * u(k,b);
        sum6+=c(k+2*n,0) * u(k,a);
        sum7+=c(k+3*n,0) * u(k,b);
        sum8+=c(k+3*n,0) * u(k,a);

    }
    double j1 = sph_bessel(l,r1*b);
    double j2 = sph_bessel(l,r1*a);
    double j3 = sph_bessel(l-1,r1*b);
    double j4 = sph_bessel(l-1,r1*a);
    
    double j5 = sph_bessel(l,r2*b);
    double j6 = sph_bessel(l,r2*a);
    double j7 = sph_bessel(l-1,r2*b);
    double j8 = sph_bessel(l-1,r2*a);
    

    double result = sum1 * j1 * j5 - sum2 * j2 * j6 +\
                    sum3 * j3 * j5 - sum4 * j4 * j6 +\
                    sum5 * j1 * j7 - sum6 * j2 * j8 +\
                    sum7 * j3 * j7 - sum8 * j4 * j8;

    return result;
}

double Levin::bessel_J(double nu, double x)
{
    return boost::math::cyl_bessel_j(nu,x);
}

double Levin::sph_bessel(int l, double x)
{
    return Basis->sph_bessel_camb(l,x);
}

double Levin::u(int k, double x)
{
    return pow(x-d, k-1);
}

double Levin::up(int k, double x)
{
    return (k-1) * pow(x-d, k-2);
}
