#include<bits/stdc++.h>
#include "Matrix.h"
#include "Range.h"
#include "Exception.h"
using namespace usr;
int main()
{
	Matrix<double> A = Matrix<double>(2,2,0),B = Matrix<double>(2,2,0);
	std::cin>>A;
	std::cout<<"A = "<<std::endl<<A<<std::endl;
	std::cout<<"AT = "<<std::endl<<(*A.transpose())<<std::endl;
	std::cout<<"determinant(A) = "<<A.determinant()<<std::endl;
	vector<std::complex<double>> eigenValue_A =A.eigenValue();
	for(int i=0;i<eigenValue_A.size();i++)
		std::cout<<eigenValue_A[i]<<std::endl<<(*A.eigenVector(eigenValue_A[i]))<<std::endl;
	std::cout<<endl; 
	//std::cout<<"eigenValue(A) = "<<A.eigenValue()<<std::endl;
	//std::cout<<"eigenVector(A) = "<<A.eigenVector()<<std::endl;
	std::cin>>B;
	std::cout<<"B = "<<std::endl<<B<<std::endl;

	return 0;
}
/*

*/