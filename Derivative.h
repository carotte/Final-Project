#ifndef DERIVATIVE_H
#define DERIVATIVE_H
#include "PROJECTHEADER.h"
/*-------------------------------------------*
 *---                                     ---*
 *---              Derivative.h           ---*   
 *---                                     ---*
 *---  This file include the header file  ---*
 *---  used to calculate derivative.      ---*
 *---  ---   ---   ---   ---   ---   ---  ---*
 *---              Tianyi Hu              ---*
 *---                                     ---*
 *-------------------------------------------*/

/*-------------------------------------------*
 *---  PURPOSE: To create a base class    ---*
 *-------------------------------------------*/
class derivative_base
{
/* I. Define empty constructor and a pure   *
      virtual function.                     */
public:
	derivative_base(){}
	virtual double f(double x)=0;
};

/*-------------------------------------------*
 *---  PURPOSE: To create a function      ---*
 *---  class that implement single para-  ---*
 )---  meter function                     ---*
 *-------------------------------------------*/
class function:public derivative_base
{
public:
/* I. Define empty constructor and implement*
      the virtual function                  */
	function(){}
	virtual double f(double x)
	{
		return (x-2)*(x-5);
	}

/* II. Define g(x) function that can be used*
       in fixed point solve function        */
	double g(double x)
	{
		return x+f(x);
	}

/* III. Define derivative function and second*
        derivative function that will be used*
		in solve&optimizer function          */
	double dx_single(double x,double h=0.00001)
	{
		return ((f(x+h)-f(x-h))/(2.0*h));
	}
	double gx_single(double x, double h=0.00001)
	{
		return ((g(x+h)-g(x-h))/(2.0*h));
	}
	double ddx_single(double x, double h=0.00001)
	{
		return (f(x+h)-f(x)*2+f(x-h))/(h*h);
	}
	double polynomial(double x,int n)
	{
		// x represent the point, n represent the function type;
		// eg: n=1 is linear, n=2 is qudratic.... 
		double sum=1.0;
		for(int i=0;i<n;i++)
			sum=sum+pow(x,n+1);
		return sum;
	}
};

/*-------------------------------------------*
 *---  PURPOSE: To create a multiple      ---*
 *---  function class                     ---*
 *-------------------------------------------*/
class function_multiple
{
/* I. Define empty constructor and implement*
      the virtual function                  */
	function_multiple(){}
	double f(vector<double>&x)
	{
		return x[0]*2;
	}

/* III. Define partial derivative function   *
        and multi-partial derivative function*
		that will be used in jacobian, hessian
		function                             */
	double partial(vector<double>&x,int i, double h=exp(-4.0))
	{
		vector<double>x_plus=x;
		x_plus[i]+=h;
		vector<double>x_minus=x;
		x_minus[i]-=h;
		return (f(x_plus)-f(x_minus))/(2*(2.0*h));
	}
	double ppartial(vector<double>&x,int i ,int k,double h= exp(-4.0))
	{
		vector<double>x_plus=x;
		x_plus[i]+=h;
		vector<double>x_minus=x;
		x_minus[i]-=h;
		vector<double>y_plus=x;
		y_plus[k]+=h;
		vector<double>y_minus=x;
		y_minus[k]-=h;	
		return 0;
	}
};


/*-------------------------------------------*
 *--- PURPOSE: To create function solver  ---*
 *---          contain all 5 method       ---*
 *-------------------------------------------*/
double solve_newton(function,double,double,double,int);
double solve_bisection(function,double,double,double Ap=exp(-6.0),double Rp=exp(-4.0),int Ns=100);
double solve_fixed_point(function , double, double Ap=exp(-6.0),double Rp=exp(-4.0),int Ns=100);
double solve_secant(function, double, double, double, int);
double solve_newton_stabilized(function,double,double,double Ap=exp(-6.0),double Rp= exp(-4.0), int Ns=20);

/*-------------------------------------------*
 *--- PURPOSE: To create function optimizer--*
 *---          contain all 5 method       ---*
 *-------------------------------------------*/
double optimize_bisection(function,double,double,double Ap=exp(-6.0),double Rp=exp(-4.0),int Ns=100);
double optimize_newton(function,double, double Ap=exp(-6.0),double Rp=exp(-4.0),int Ns=100);
double optimize_secant(function,double,double Ap=exp(-6.0),double Rp=exp(-4.0),int Ns=100);
double optimize_newton_stabilized(function,double,double,double Ap=exp(-6.0),double Rp=exp(-4.0),int Ns=100);
double optimize_golden_search(function,double,double,double Ap=exp(-6.0),double Rp=exp(-4.0),int Ns=100);


#endif