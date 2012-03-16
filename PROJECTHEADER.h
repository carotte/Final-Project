/*-------------------------------------------*
 *---                                     ---*
 *---           PROJECTHEADER.h           ---*   
 *---                                     ---*
 *---  This file include the header file  ---*
 *---  used by Numeric final project.     ---*
 *---  ---   ---   ---   ---   ---   ---  ---*
 *---              Tianyi Hu              ---*
 *---                                     ---*
 *-------------------------------------------*/

/*---     Inclusion of C++ standard       ---*
 *---          library headers:           ---*/

#ifndef PROJECTHEADER_H
#define PROJECTHEADER_H
#include <typeinfo>
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

#include "Derivative.h"
/*------------------------------------------*
 *--- PURPOSE: To create a Matrix class  ---*
 *------------------------------------------*/

class Matrix
{
/* I. Define a out stream operator overload *
      to print Matrix.                      */
	friend ostream &operator<<(ostream &,const Matrix &);   

public:
/* II. Class Constructor and Copy Constructor*/
	Matrix(int r=1,int c=1); 
	Matrix (const Matrix &); 
	
/* III.Operator Overloading including: 
       SUBSCRIPT, PLUS, MINUS, MULTIPLE, 
	   DIVIDEN, TRANSPOSE and IDENTITY.     */
	float operator()(int,int) const;
	float &operator()(int,int);
	Matrix operator+(const Matrix &);
	Matrix operator-(const Matrix &);
	Matrix operator+(float);
	Matrix operator-(float);
	Matrix operator*(float );
	Matrix operator*(const Matrix &);
	Matrix t();
	Matrix identity(int,double,double);

/* IV.Swap function used in inverse matrix   */
	void swap(float &a, float&b){float c=a;a=b;b=c;};

/* V. Access functions to read rows and columns.
	  And utility functons to set rows 
	  and columns.                          */
	int getRows() const {return rows;};
	int getCols() const {return cols;};
	void setRows(int a) { rows=a;};
	void setCols(int a) { cols=a;};
	vector<float> getData() const {return data;};

/* VI. Private parameter.                   */
private:
	int rows;
	int cols;
	vector<float> data;
};


/*-------------------------------------------*
 *--- PURPOSE: To create numeric function ---*
 *---          and implement test.        ---*
 *-------------------------------------------*/


/* I. Matrix checking function, check symmetric
      zero and positive definite.            */
Matrix inv(Matrix );
bool is_almost_symmetric(Matrix, double Ap=exp(-6.0),double Rp=exp(-4.0));
bool is_almost_zero(Matrix, double Aap=exp(-6.0),double Rp=exp(-4.0));
bool is_positive_definite(Matrix&);

/* II. Norm and condition number function,
       norm return the norm of both Matrix
	   and single numer.                    */
double norm(double,int p=1);
double norm(Matrix, double p=1.0);
double cond_number(Matrix);

/* III. Cholesky, Exponential Taylor expension
        Markovitz and fit least square. 
		struct result is used to return multiple
	    value for Markovitz and fit least square
		function.                            */
Matrix exp_taylor(Matrix &,double Ap=exp(-6.0),double Rp=exp(-4.0),int Ns=40);
Matrix cholesky(Matrix & );
struct result
{
	result(){};
	Matrix weight;
	Matrix ret;
	float rsk;
};
result markovitz(Matrix &, Matrix &, double Rf=0.05);
result fit_least_squares(Matrix &, int n);


Matrix gradient(function_multiple, vector<double> , double H=exp(-4.0));
Matrix hessian(function_multiple, vector<double>, double H=exp(-4.0));
Matrix jacobian(function_multiple, vector<double>,double H=exp(-4.0));
void Run();

#endif