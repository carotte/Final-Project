#include "PROJECTHEADER.h"
#include "Derivative.h"
#include <stdexcept>
using std::runtime_error;
#include<iomanip>
#include <cmath>
#include <iostream>
using namespace std;
/*-------------------------------------------*
 *---                                     ---*
 *---              Function.cpp           ---*   
 *---                                     ---*
 *---  This file define all the function  ---*
 *---  claimed in projectheader and       ---*
 *---  derivative header file. And it     ---*
 *---  finally provides a test function   ---*
 *---  named Run to test all these        ---*
 *---  functions.                         ---*
 *---                                     ---*
 *---  The order functions get called and ---*
 *---  the data used test are identical   ---*
 *---  to the lecture note. Therefore,    ---*
 *---  test result should also be the same---*
 *---  as that in lecture note.           ---*
 *---  ---   ---   ---   ---   ---   ---  ---*
 *---              Tianyi Hu              ---*
 *---                                     ---*
 *-------------------------------------------*/

/*--------------------------------------------------------*
 *---           PURPOSE: print the Matrix              ---*
 *--------------------------------------------------------*/
ostream &operator<<(ostream &output, const Matrix &matrix)
{
	output<<"[";
	for (int r=0;r<matrix.rows;r++)
	{
		output<<"[";
		for(int c=0;c<matrix.cols;c++)
		{
			if (c>0) 
				output<<",";
			output<<matrix(r,c);
		}
		output<<"],\n";
	}
	output<<"]\n";
	return output;
}

/*--------------------------------------------------------*
 *---      PURPOSE: initiate class and copy the value  ---*
 *--------------------------------------------------------*/
Matrix::Matrix(int rows,int cols)
{
	this->rows=rows;
	this->cols=cols;
	this->data.resize(rows*cols);
	for(int r=0;r<rows;r++)
		for(int c=0;c<cols;c++)
			this->data[r*cols+c]=0;
}
/*--------------------------------------------------------*
 *---      PURPOSE: copy constructor                   ---*
 *--------------------------------------------------------*/
Matrix::Matrix(const Matrix & A)
	:rows(A.getRows ()),cols(A.getCols ()),data(A.getData ())
{
}

/*--------------------------------------------------------*
 *---         PURPOSE: Matrix calculation              ---*
 *--------------------------------------------------------*/
float Matrix::operator()(int i, int j) const
{
	return data[i*cols+j];
}
float &Matrix::operator()(int i,int j)
{
	return data[i*cols+j];
}
Matrix Matrix::operator+(const Matrix &A)
{
	if (this->rows!=A.rows||this->cols !=A.cols )
		cout<<"Unable to ADD\n";
	Matrix C(this->rows ,this->cols );
	for(int r=0;r<this->rows;r++)
		for(int c=0;c<this->cols;c++)
			C(r,c)=A(r,c)+this->data[r*cols+c];
	return C;
}
Matrix Matrix::operator+(float a)
{
	Matrix C(this->rows ,this->cols );
	for(int r=0;r<this->getRows ();r++)
		for(int c=0;c<this->getCols ();c++)
			C(r,c)=this->data[r*cols+c]+a;
	return C;
}
Matrix Matrix::operator-(const Matrix &A)
{
	if (this->rows!=A.rows||this->cols !=A.cols )
		cout<<"Unable to ADD\n";
	Matrix C(this->rows ,this->cols );
	for(int r=0;r<this->rows;r++)
		for(int c=0;c<this->cols;c++)
			C(r,c)=this->data[r*cols+c]-A(r,c);
	return C;
}
Matrix Matrix::operator-(float a)
{
	return *this+(-1*a);
}
Matrix Matrix::operator*(float a)
{
	Matrix C(this->rows ,this->cols );
	for(int r=0;r<this->rows;r++)
		for(int c=0;c<this->cols;c++)
			C(r,c)=a*this->data[r*cols+c];
	return C;
}
Matrix Matrix::operator*(const Matrix  &B)
{
	
	if(this->getCols ()!=B.getRows ())
		cout<<"Unable to multiple\n";
	Matrix C(this->getRows (),B.getCols () );
	for(int r=0;r<this->getRows ();r++)
		for(int c=0;c<B.getCols ();c++)
			for(int k=0;k<this->getCols ();k++)
				C(r,c)+=this->data[r*cols+k]*B(k,c);
	return C;

}
Matrix inv(Matrix A)
{
	if(A.getCols ()!=A.getRows ())
		cout<<"BAD\n"<<endl;
	Matrix B(A.getRows (),A.getCols ());
	float p; 
	float q; 
	for(int r=0;r<B.getCols ();r++) B(r,r)=1;
	for(int c=0;c<A.getCols ();c++)
	{
		p=A(c,c);
		for(int i=0;i<A.getCols ();i++){
			A(c,i)/=p;
			B(c,i)/=p;
		}
		for(int r=0;r<A.getRows ();r++)
			if(r!=c){
				q=A(r,c);
				for(int i=0;i<A.getCols ();i++){
					A(r,i)-=q*A(c,i);
					B(r,i)-=q*B(c,i);
				}
			}
	}
	return B;
}

/*--------------------------------------------------------*
 *---   PURPOSE: Matrix identity and transpose matrix  ---*
 *--------------------------------------------------------*/
Matrix Matrix::identity(int r=1,double one=1.0,double fill=0.0)
{
	Matrix M;
	M=Matrix(r,r);
	for(int i=0;i<M.getCols ();i++)
		for(int j=0;j<M.getRows ();j++)
			M(j,i)=static_cast<float>(fill);
	for(int x=0;x<M.getRows ();x++)
		M(x,x)=static_cast<float>(one);
	return M;
}
Matrix Matrix::t()
{
	Matrix C(this->cols ,this->rows );
	for(int r=0;r<this->rows;r++)
		for(int c=0;c<this->cols;c++)
			C(c,r)=this->data[r*cols+c];
	return C;
}

/*--------------------------------------------------------*
 *---   PURPOSE: Matrix symmetric and zero check.      ---*
 *---            return 0 for faluse and 1 for true    ---*
 *--------------------------------------------------------*/
bool is_almost_symmetric(Matrix a, double ap,double rp)
{
	double delta;
	if (a.getRows()!=a.getCols() )
		return false;
	for(int r=0;r<(a.getRows()-1);r++)
		for(int c=0;c<r;c++)
		{
			delta=abs(a(r,c)-a(c,r));
			if (delta>ap && delta>rp*(max(abs(a(r,c)),abs(a(c,r)))))
				return false;
		}
	return true;
}
bool is_almost_zero(Matrix a , double ap,double rp)
{
	double delta;
	for(int r=0;r<a.getRows();r++)
		for(int c=0;c<a.getCols();c++)
		{
			delta=abs(a(r,c)-a(c,r));
			if (delta>ap && delta>rp*(max(abs(a(r,c)),abs(a(c,r)))))
				return false;
		}
	return true;
}
bool is_positive_definite(Matrix & a)
{
	if(is_almost_symmetric(a)==false)
		return false;
	cholesky(a);
	return true;
}

/*--------------------------------------------------------*
 *---   PURPOSE: Return norm for both Matrix and       ---*
 *---            single number                         ---*
 *--------------------------------------------------------*/
double norm(Matrix a,double p)
{
	double sum=0.0;
	if(a.getRows()==1)
	{
		for(int i=0;i<a.getCols();i++)
			sum=sum+pow(static_cast<double>(a(0,i)),p);
		return abs(pow(sum,1.0/p));
	}
	if (a.getCols ()==1)
	{
		for(int i=0;i<a.getRows ();i++)
			sum=sum+pow(static_cast<double>(a(i,0)),p);
		return abs(pow(sum,1.0/p));
	}
	if(p==1.0)
	{
		double temp=0.0;
		for(int c=0;c<a.getCols ();c++)
		{	double col_temp=0.0;
			for(int r=0;r<a.getRows ();r++)
				col_temp+=abs(a(r,c));
			if (col_temp>=temp)
				temp=abs(col_temp);
		}
		return temp;
	}
	else 
		throw runtime_error("unable to implement request.");
}
double norm(double a,int p)
{
	if(typeid(a)==typeid(double)|| typeid(a)==typeid(int))
		return abs(static_cast<double>(a));
	else 
		throw runtime_error(" unable to implement request.");
}

/*--------------------------------------------------------*
 *---   PURPOSE: Return the condition number of        ---*
 *---            a Matrix                              ---*
 *--------------------------------------------------------*/
double cond_number(Matrix a)
{
	return norm(a)*norm(inv(a));
}


Matrix exp_taylor(Matrix & a ,double ap,double rp,int ns)
{	
	Matrix x;
	Matrix t=x.identity(a.getCols());
	Matrix s=x.identity(a.getCols ());
	for (int k=1;k<ns;k++)
	{
		t=t*a*(1/static_cast<double>(k));	
		s=s+t;
		if (norm(t)<max(ap,rp*norm(s)))
			return s;
	} 
	throw runtime_error("No Convergence! ");
}
Matrix cholesky(Matrix & a)
{
	if (is_almost_symmetric(a)==false)
		throw runtime_error("Not symmetric");
	Matrix L(a);
	for(int k=0;k<L.getCols ();k++)
	{
		if(L(k,k)<=0)
			throw runtime_error("Not Postitive definitive! ");
		L(k,k)=sqrt(L(k,k));
		double p=L(k,k);
		for(int i=k+1;i<L.getRows ();i++)
			L(i,k)/=p;
		for(int j=k+1;j<L.getRows ();j++)
		{
			p=L(j,k);
			for(int i=k+1;i<L.getRows ();i++)
				L(i,j)-=L(i,k)*p;
		}
	}  
	for(int i =0;i<L.getRows ();i++)
		for(int j=i+1;j<L.getCols ();j++)
			L(i,j)=0;
	return L;
}


result markovitz(Matrix & mu, Matrix & cov, double rf)
{
	Matrix x(cov.getRows (),1);	
	x=(inv(cov))*(mu-static_cast<float>(rf));
	
	double sum=0.0;
	for(int r=0;r<x.getRows ();r++)
		sum+=x(r,0);
	x=x*(1/sum);
	vector<float> risk;
	risk=(x.t()*(cov*x)).getData();
	result ans;
	ans.weight=x.t();
	ans.ret=(x.t()*mu);
	ans.rsk=sqrt(risk[0]);
	return ans;
}
result fit_least_squares(Matrix & points, int n)
{
	Matrix A(points.getRows (),n+1);
	Matrix B(points.getRows (),1);
	result ans;
	for(int i=0;i<A.getRows ();i++)
		{
			float weight=0.0;
			if (points.getRows ()>2)
				weight=1.0/points(i,2);
			else
				weight=1.0;		
			B(i,0)=points(i,1)*weight;
			
			for(int j=0;j<A.getCols ();j++)
				A(i,j)=pow(points(i,0)*weight,j);
	}
	
	Matrix c=inv((A.t()*A))*(A.t()*B);
	Matrix chi=A*c-B;
	double chi2=norm(chi,2)*norm(chi,2);
	ans.ret=c;
	ans.rsk=chi2;
	return ans;
}



double solve_newton(function F,double x ,double ap=exp(-6.0),double rp=exp(-4.0),int ns=100)
{
	double fx; double dfx;double x_old;
	for(int k=0;k<ns;k++)
	{
		fx=F.f(x);
		dfx=F.dx_single(x);

		if(norm(dfx,1)<ap)
			throw runtime_error("unstable solution! ");
		x_old=x;
		x=x-fx/dfx;
		if(k>2 && norm(x-x_old,1)<max(ap,norm(x,1)*rp))
			return x;
	}
	throw runtime_error("no convergence! ");
	
}
double solve_bisection(function F ,double a,double b,double ap,double rp,int ns)
{	

		double fa=F.f(a);
		double fb=F.f(b);
		double x;double fx;
		if(fa==0)
			return a;
		if(fb==0)
			return b;
		if(fa*fb>0)
			throw runtime_error("f(a) and f(b) must have opposite sign! ");
		for(int k=0;k<ns;k++)
		{
			x=(a+b)/2;
			fx=F.f(x);
			if(fx==0 || norm(b-a,1)<max(ap,norm(x,1)*rp))
				return x;
			else if ((fx*fa)<0)
			{
				b=x;
				fb=fx;
			}
			else
			{
				a=x;
				fa=fx;
			}
		}
		throw runtime_error("no convergence! ");
}
double solve_fixed_point(function F, double x, double ap,double rp,int ns)
{
		double g=F.g(x);
		double dg=F.dx_single(x)+1.0;
		double x_old;
		for(int k=0;k<ns;k++)
		{
			if(abs(dg)>=1)
				throw runtime_error("error Derivative of g larger than 1");
			x_old=x;
			x=F.g(x);
			if(k>2 && norm(x_old-x,1) < max(ap,norm(x,1)*rp))
				return x;
		}
		throw runtime_error("No convergence! ");
}
double solve_secant(function F,double x, double ap=exp(-6.0), double rp=exp(-4.0),int ns=20)
{

	double fx=F.f(x);
	double dfx=F.dx_single(x);double fx_old; double x_old;
	for(int k=0;k<ns;k++)
	{	
		if (norm(dfx)<ap)
			throw runtime_error("unstable solution ");
		x_old=x; fx_old=fx; x=x-fx/dfx;
		if(k>2 &&norm(x-x_old)<max(ap,norm(x)*rp))
			return x;
		fx=F.f(x);
		dfx=(fx-fx_old)/(x-x_old);
	}
	throw runtime_error("No convergence");

}

double solve_newton_stabilized(function F,double a ,double b,double ap,double rp, int ns)
{
	double fa=F.f(a);
	double fb=F.f(b);
	if(fa==0)
		return a;
	if(fb==0)
		return b;
	if(fa*fb > 0 )
		throw runtime_error("fa fb same sign");
	double x=(a+b)/2.0; double x_old;double fx_old; 
	double fx=F.f(x); double dfx=F.dx_single(x);
	for(int k=0;k<ns;k++)
	{
		x_old=x;
		fx_old=fx;
		if (norm(dfx)>ap)
			x=x-fx/dfx;
		if(x==x_old ||  x<a || x>b)
			x=(a+b)/2;
		fx=F.f(x);
		if(fx==0 || norm(x-x_old)<max(ap,norm(x)*rp))
			return x;
		dfx=(fx-fx_old)/(x-x_old);
		if(fx*fa<0)
		{
			b=x;
			fb=fx;
		}
		else
		{
			a=x;
			fa=fx;
		}
	}
		throw runtime_error("no convergence");
}
double optimize_bisection(function F,double a ,double b,double ap,double rp,int ns)
{
	double dfa=F.dx_single(a);
	double dfb=F.dx_single(b);
	double x;double dfx;
	if(dfa*dfb>0)
		runtime_error("dfa dfb opposite sign");
	for(int k=0;k<ns;k++)
	{
		x=(a+b)/2;
		dfx=F.dx_single (x);
		if(dfx==0||norm(b-a)<max(ap,norm(x)*rp))
			return x;
		else if(dfx*dfa<0)
		{
			b=x;
			dfb=dfx;
		}
		else
		{
			a=x;
			dfa=dfx;
		}
	}
	throw runtime_error("no convengence");
		
}
double optimize_newton(function F,double x, double ap,double rp,int ns)
{	
	double dfx;double ddfx;double x_old;
	for(int k =0; k<ns;k++)
	{
		dfx=F.dx_single(x);
		ddfx=F.ddx_single(x);
		if (dfx==0)
			return x;
		if(norm(ddfx)<ap)
			runtime_error("unstable solution");
		x_old=x;
		x=x-dfx/ddfx;
		if(norm(x-x_old)<max(ap,norm(x)*rp))
			return x;
	}
		throw runtime_error("no convengence");
}
double optimize_secant(function F,double x,double ap,double rp,int ns)
{
	double fx=F.f(x);
	double dfx=F.dx_single(x);
	double ddfx=F.ddx_single(x);
	double x_old;double dfx_old;
	for(int k=0;k<ns;k++)
	{
		if(dfx==0)
			return x;
		if (norm(ddfx)<ap)
			runtime_error("unstable solution");
		x_old=x;
		dfx_old=dfx;
		x=x-dfx/ddfx;
		if(norm(x-x_old)<max(ap,norm(x)*rp))
			return x;
		fx=F.f(x);
		dfx=F.dx_single(x);
		ddfx=(dfx-dfx_old)/(x-x_old);
	}
	throw runtime_error("no convengence");
}
double optimize_newton_stabilized(function F,double a,double b,double ap,double rp,int ns)
{
	double dfa=F.dx_single(a);
	double dfb=F.dx_single(b);
	if(dfa==0)
		return a;
	if(dfb==0)
		return b;
	if(dfa*dfb > 0 )
		throw runtime_error("fa fb have same sign");
	double x=(a+b)/2;
	double x_old;double fx_old;double dfx_old;
	double fx=F.f(x);double dfx=F.dx_single(x);double ddfx=F.ddx_single(x);
	for(int k=0;k<ns;k++)
	{
		if(dfx==0)
			return x;
		x_old=x; fx_old=fx; dfx_old=dfx;
		if(norm(ddfx)>ap)
			x=x-dfx/ddfx;
		if (x==x_old || x<a|| x>b)
			x=(a+b)/2;
		if(norm(x-x_old)<max(ap,norm(x)*rp))
			return x;
		fx=F.f(x);
		dfx=(fx-fx_old)/(x-x_old);
		ddfx=(dfx-dfx_old)/(x-x_old);
		if(dfx*dfa<0)
		{
			b=x;
			dfb=dfx;
		}
		else
		{
			a=x;
			dfa=dfx;
		}
	}
	throw runtime_error("no convergence");
}
double optimize_golden_search(function F,double a,double b,double ap,double rp,int ns)
{
	double tau=(sqrt(5.0)-1.0)/2.0;
	double x1=a+(1-tau)*(b-a);
	double x2=a+tau*(b-a);
	double fa=F.f(a);
	double fb=F.f(b);
	double f1=F.f(x1);
	double f2=F.f(x2);
	for(int k=0;k<ns;k++)
	{
		if(f1>f2)
		{
			a=x1; fa=f1; x1=x2; f1=f2;
			x2= a+tau*(b-a);
			f2=F.f(x2);
		}
		else
		{
			b=x2; fb=f2; x2=x1; f2=f1;
			x1=a+(1-tau)*(b-a);
			f1=F.f(x1);
		}
		if (k>2 && norm(b-a)<max(ap,norm(b)*rp))
			return b;
	}
	throw runtime_error("no convergence");
	
}





Matrix gradient( function_multiple F, vector<double> x, double h)
{
	int size=x.size();
	Matrix g(size);
	for(int i=0;i<size;i++)
	{
		
	}
	return g;
}
Matrix hessian(function_multiple F, vector<double> x, double h)
{
	int size=x.size();
	Matrix g(size,size);
	return g;
}
Matrix jacobian(function_multiple F, vector<double> x,double h)
{
	int size=x.size();
	Matrix g(size,size);
	//for(int i=0;i<size;i++)
	return g;

}


void Run()
{
	Matrix A(3,3);
	A(0,0)=1;A(0,1)=2;A(0,2)=3;
	A(1,0)=4;A(1,1)=5;A(1,2)=6;
	A(2,0)=4;A(2,1)=4;A(2,2)=4;
	
	Matrix B(2,2);
	B(0,0)=1;B(0,1)=2;
	B(1,0)=3;B(1,1)=4;
	Matrix C(3,3);
	C(0,0)=4;C(0,1)=2;C(0,2)=1;
	C(1,0)=2;C(1,1)=9;C(1,2)=3;
	C(2,0)=1;C(2,1)=3;C(2,2)=16;

	Matrix cov(3,3);
	cov(0,0)=0.04;cov(0,1)=0.006;cov(0,2)=0.02;
	cov(1,0)=0.006;cov(1,1)=0.09;cov(1,2)=0.06;
	cov(2,0)=0.02;cov(2,1)=0.06;cov(2,2)=0.16;

	double rf=0.05;
	Matrix mu(3,1);
	mu(0,0)=0.1;mu(1,0)=0.12;mu(2,0)=0.15;
	function F;

	cout<<"     "<<"Matrix and Derivative test"<<endl;
	cout<<"     "<<"Matrix A is:"<<endl<<A<<endl;
	cout<<"     "<<"Matrix symmetric test \n     If Matrix A is symmetric? (0 for No and 1 for Yes):"<<is_almost_symmetric(A)<<endl;
	cout<<"\n     Matrix zero test \n";
	cout<<"     If Matrix A is almost zero? (0 for No and 1 for Yes):"<<is_almost_zero(A)<<endl;
	cout<<"\n     Matrix norm test \n     Matrix B is:\n"<<B<<endl;
	cout<<"     Matrix B's condition number is:";
	try{cout<<cond_number(B)<<endl;}
	catch(runtime_error & error)
	{
		cout<<"     "<<"Exception Occured: "<<error.what()<<endl;
	}
	cout<<"\n     Exp function test:\n"<<"     "<<"Exponential Taylor series result is:\n";
	try{cout<<exp_taylor(B)<<endl;}
	catch(runtime_error & error)
	{
		cout<<"     "<<"Exception Occured: "<<error.what()<<endl;
	}
	cout<<"\n     Cholesky function test:\n";
	cout<<"     "<<"Matrix used in test is Matrix C:\n"<<C<<endl;
	cout<<"     "<<"Matrix L=cholesky(C), L is:\n";
	try{cout<<cholesky(C)<<endl;}
	catch(runtime_error & error)
	{
		cout<<"Exception Occured: "<<error.what()<<endl;
	}
	Matrix L=cholesky(C); 
	cout<<"     "<<"And check if A-L*L.t is almost zero. (0 for No and 1 for Yes):"<<is_almost_zero(C-L*L.t())<<endl;
	cout<<"     "<<"Markovitz Test\n";
	cout<<"     "<<"Risk free rate is: 5%."<<endl;
	cout<<"     "<<"Portfolio cov is:\n"<<cov<<endl;
	cout<<"     "<<"Portfolio mu is:\n"<<mu<<endl;
	cout<<"     "<<"Stock weight is:";
	cout<<markovitz(mu,cov,rf).weight;
	cout<<"\n     Stock return is:";
	cout<<markovitz(mu,cov,rf).ret;
	cout<<"\n     Stock risk is:";
	cout<<markovitz(mu,cov,rf).rsk<<endl;

	cout<<"\n     Fixed point method test:"<<endl;
	cout<<"     "<<"Tested function f(x) is (x-2)*(x-5)"<<endl;
	cout<<"     "<<"The fixed point method result when point x=1 is: ";
	try{cout<<solve_fixed_point(F,1.0)<<endl;}
	catch(runtime_error & error)
	{
		cout<<"     "<<"Exception Occured:\n"<<"     "<<error.what()<<endl;
	}

	cout<<"\n     Bisection method test:"<<endl;
	cout<<"     "<<"Tested function f(x) is (x-2)*(x-5)"<<endl;
	cout<<"     "<<"The bisection method result when range is 1 to 3 is: ";
	try{cout<<solve_bisection(F,1.0,3.0)<<endl;}
	catch(runtime_error & error)
	{
		cout<<"     "<<"Exception Occured: "<<error.what()<<endl;
	}
	
	cout<<"\n     Newton method test:"<<endl;
	cout<<"     "<<"Tested function f(x) is (x-2)*(x-5)"<<endl;
	cout<<"     "<<"The Newton method result when point is 1 is: ";
	try{cout<<solve_newton(F,1.0)<<endl;}
	catch(runtime_error & error)
	{
		cout<<"     "<<"Exception Occured: "<<error.what()<<endl;
	}
	
	cout<<"\n     Secant method test:"<<endl;
	cout<<"     "<<"Tested function f(x) is (x-2)*(x-5)"<<endl;
	cout<<"     "<<"The Secant method result when range is from 1 to 3 is: ";
	try{cout<<solve_newton(F,1.0, 3.0)<<endl;}
	catch(runtime_error & error)
	{
		cout<<"     "<<"Exception Occured: "<<error.what()<<endl;
	}

	cout<<"\n     Newton Stabilized method test:"<<endl;
	cout<<"     "<<"Tested function f(x) is (x-2)*(x-5)"<<endl;
	cout<<"     "<<"The Newton Stabilized method result when range is from 1 to 3 is: ";
	try{cout<<solve_newton_stabilized(F,1.0, 3.0)<<endl;}
	catch(runtime_error & error)
	{
		cout<<"     "<<"Exception Occured: "<<error.what()<<endl;
	}

	cout<<"\n     Bisection Optimize method test:"<<endl;
	cout<<"     "<<"Tested function f(x) is (x-2)*(x-5)"<<endl;
	cout<<"     "<<"The Bisection Optimize method result when range is from 2 to 5 is: ";
	try{cout<<optimize_bisection(F,2.0, 5.0)<<endl;}
	catch(runtime_error & error)
	{
		cout<<"     "<<"Exception Occured: "<<error.what()<<endl;
	}

	cout<<"\n     Newton Optimize method test:"<<endl;
	cout<<"     "<<"Tested function f(x) is (x-2)*(x-5)"<<endl;
	cout<<"     "<<"The Newton Optimize method result when point is 3 is: ";
	try{cout<<optimize_newton(F,3.0)<<endl;}
	catch(runtime_error & error)
	{
		cout<<"     "<<"Exception Occured: "<<error.what()<<endl;
	}

	cout<<"\n     Secant Optimize method test:"<<endl;
	cout<<"     "<<"Tested function f(x) is (x-2)*(x-5)"<<endl;
	cout<<"     "<<"The Secant Optimize method result when point is 3 is: ";
	try{cout<<optimize_secant(F,3.0)<<endl;}
	catch(runtime_error & error)
	{
		cout<<"     "<<"Exception Occured: "<<error.what()<<endl;
	}

	cout<<"\n     Newton Stabilized Optimize method test:"<<endl;
	cout<<"     "<<"Tested function f(x) is (x-2)*(x-5)"<<endl;
	cout<<"     "<<"The Newton Stabilized Optimize method result when range is 2 to 5 is: ";
	try{cout<<optimize_newton_stabilized(F,2.0,5.0)<<endl;}
	catch(runtime_error & error)
	{
		cout<<"     "<<"Exception Occured: "<<error.what()<<endl;
	}

	cout<<"\n     Golden Search Optimize method test:"<<endl;
	cout<<"     "<<"Tested function f(x) is (x-2)*(x-5)"<<endl;
	cout<<"     "<<"The Golden Search Optimize method result when range is 2 to 5 is: ";
	try{cout<<optimize_golden_search(F,2.0,5.0)<<endl;}
	catch(runtime_error & error)
	{
		cout<<"     "<<"Exception Occured: "<<error.what()<<endl;
	}
}