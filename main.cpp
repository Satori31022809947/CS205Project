#include<bits/stdc++.h>
#include "Matrix.h"
#include "Range.h"
#include "Exception.h"
#include <memory>
using namespace usr;
using namespace std;
typedef unique_ptr<Matrix<double>> ptr_double;
typedef unique_ptr<Matrix<complex<double>>> ptr_cd;
int main()
{
	ptr_double address1, address2;
	ptr_cd address3, address4;

	auto A = Matrix<double>(2,2,1,0),B = Matrix<double>(2,2,1,0);
	auto C = Matrix<complex<double>>(2,2,1,0);
	C[0][0][0] = complex<double>(1,1);
	C[0][0][1] = complex<double>(1,1);
	C[0][1][0] = complex<double>(1,1);
	C[0][1][1] = complex<double>(1,1);
	cin>>A;
	cout<<"A = "<<endl<<A<<endl;

	address1.reset(A.transpose());
	cout<<"transpose(A) = "<<endl<<*address1<<endl;

	address3.reset(C.conjugate());
	cout<<"conjugate(A) = "<<endl<<*address3<<endl;
	
	cout<<"determinant(A) = "<<A.determinant()<<endl;
	cout<<"rank(A) = "<<A.rank()<<endl;
	cout<<"trace(A) = "<<A.trace()<<endl;
	cout<<"max(A) = "<<A.max()<<endl;
	cout<<"min(A) = "<<A.min()<<endl;
	cout<<"sum(A) = "<<A.sum()<<endl;
	cout<<"avg(A) = "<<A.avg()<<endl;
	
	address1.reset(A.toDiagnal());
	cout<<"toDiagnal(A) = "<<endl<<*address1<<endl;
	
	address1.reset(A.inverse());
	cout<<"inverse(A) = "<<endl<<*address1<<endl;
	
	address1.reset(A*(*address1));
	address2.reset(A.inverse());
	address2.reset((*address2)*A);
	cout<<*address1<<endl<<*address2<<endl;

	address1.reset(A.subMatrix(Range(0, 2), Range(0,1)));
	cout<<"submatrix(A) = "<<endl<<*address1<<endl;

	address1.reset(A.slice(Range(0,2), 1));
	cout<<"slice(A) = "<<endl<<*address1<<endl;

	address1.reset(A.reshape(2,2));
	cout<<"reshape(A)= "<<endl<<*address1<<endl;

	address1.reset(A^5);
	cout<<"power(A)= "<<endl<<*address1<<endl;
	
	vector<complex<double>> eigenValue_A =A.eigenValue();
	for(int i=0;i<eigenValue_A.size();i++)
	{
		address3.reset(A.eigenVector(eigenValue_A[i]));
		cout<<eigenValue_A[i]<<endl<<*address3<<endl;
	}
	cout<<endl; 
	//cout<<"eigenValue(A) = "<<A.eigenValue()<<endl;
	//cout<<"eigenVector(A) = "<<A.eigenVector()<<endl;
	cin>>B;
	cout<<"B = "<<endl<<B<<endl;

	address1.reset(A*B);
	cout<<"A*B = "<<endl<<*address1<<endl;

	address1.reset(A/B);
	cout<<"A/B = "<<endl<<*address1<<endl;

	cout<<"dotProduct(AB) = "<<A.dotProduct(B)<<endl;

	// address1.reset(A.crossProduct(B));
	// cout<<"crossProduct(AB) = "<<*address1<<endl;

	A.release(), B.release(), C.release(), address1.reset(), address2.reset(), address3.reset(), address4.reset();
	MemoryDetech::instance().show();
	return 0;
}
/*

*/