#include<bits/stdc++.h>
#include "Matrix.h"
#include "Range.h"
#include "Exception.h"
#include <memory>
using namespace usr;
using namespace std;
int main()
{
	{
	auto A = Matrix<double>(2,2,1,0),B = Matrix<double>(2,2,1,0);
	auto C = Matrix<complex<double>>(2,2,1,0);
	C[0][0][0] = complex<double>(1,1);
	C[0][0][1] = complex<double>(1,1);
	C[0][1][0] = complex<double>(1,1);
	C[0][1][1] = complex<double>(1,1);
	cin>>A;
	cout<<"A = "<<endl<<A<<endl;

	cout<<"transpose(A) = "<<endl<<A.transpose()<<endl;

	cout<<"conjugate(C) = "<<endl<<C.conjugate()<<endl;
	
	cout<<"determinant(A) = "<<A.determinant()<<endl;
	cout<<"rank(A) = "<<A.rank()<<endl;
	cout<<"trace(A) = "<<A.trace()<<endl;
	cout<<"max(A) = "<<A.max()<<endl;
	cout<<"min(A) = "<<A.min()<<endl;
	cout<<"sum(A) = "<<A.sum()<<endl;
	cout<<"avg(A) = "<<A.avg()<<endl;
	
	cout<<"toDiagnal(A) = "<<endl<<A.toDiagnal()<<endl;
	
	cout<<"inverse(A) = "<<endl<<A.inverse()<<endl;
	
	cout<<A*A.inverse()<<endl<<A.inverse()*A<<endl;

	cout<<"submatrix(A) = "<<endl<<A.subMatrix(Range(0, 2), Range(0,1))<<endl;

	cout<<"slice(A) = "<<endl<<A.slice(Range(0,2), 1)<<endl;

	cout<<"reshape(A)= "<<endl<<A.reshape(2,2)<<endl;

	cout<<"power(A)= "<<endl<<(A^2)<<endl;
	
	vector<complex<double>> eigenValue_A =A.eigenValue();
	for(int i=0;i<eigenValue_A.size();i++)
	{
		cout<<eigenValue_A[i]<<endl<<A.eigenVector(eigenValue_A[i])<<endl;
	}
	cout<<endl; 
	//cout<<"eigenValue(A) = "<<A.eigenValue()<<endl;
	//cout<<"eigenVector(A) = "<<A.eigenVector()<<endl;
	cin>>B;
	cout<<"B = "<<endl<<B<<endl;

	cout<<"A*B = "<<endl<<A*B<<endl;

	cout<<"A/B = "<<endl<<A/B<<endl;

	cout<<"dotProduct(AB) = "<<A.dotProduct(B)<<endl;

	// address1.reset(A.crossProduct(B));
	// cout<<"crossProduct(AB) = "<<*address1<<endl;
	}
	MemoryDetech::instance().show();
	return 0;
}
/*

*/