/*
 * @Author: Satori 3102809947@qq.com
 * @Date: 2022-05-27 23:01:36
 * @LastEditors: Satori 3102809947@qq.com
 * @LastEditTime: 2022-06-03 09:57:41
 * @FilePath: \CS205Project\matrix.cpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include "interface.h"
#include <iostream>
using namespace std;
using namespace usr;
int main()
{

    uint32 row1, row2, col1, col2, channel1, channel2;
    cout << "Test1: basic operation" << endl;
    {
        cin >> row1 >> col1; 
        Matrix<double> usr1(row1, col1);
        cin >> usr1;
        cout<<"determinant = "<<usr1.determinant()<<endl;
        cout<<"rank = "<<usr1.rank()<<endl;
        cout<<"trace = "<<usr1.trace()<<endl;
        cout<<"max = "<<usr1.max()<<endl;
        cout<<"min = "<<usr1.min()<<endl;
        cout<<"sum = "<<usr1.sum()<<endl;
        cout<<"avg = "<<usr1.avg()<<endl;
        cout<<"transpose = \n"<<usr1.transpose()<<endl;
        cout<<"inverse = \n"<<usr1.inverse()<<endl;
        cout<<"gaussian elimination = \n"<<usr1.toDiagnal()<<endl;
        cout<<"submatrix = \n"<<usr1.subMatrix(Range(0, 2), Range(col1/2,col1))<<endl;
        cout<<"slice = \n"<<usr1.slice(Range(0,2), col1/2)<<endl;
        cout<<"reshape = \n"<<usr1.reshape(2,2)<<endl;
    }

    cout << "Test2: eigenvalue and eigenvector" << endl;
    {
        cin >> row1 >> col1; 
        Matrix<double> usr1(row1, col1);
        cin >> usr1;
        usr1 = usr1 + usr1.transpose();
        vector<complex<double>> eigenvalue = usr1.eigenValue();
        cout<<"eigenvalue = \n";
        for (int i = 0; i < eigenvalue.size(); i++)
            cout << eigenvalue[i] << " ";
        putchar('\n');
        cout<<"eigenvector = \n";
        for (int i = 0; i < eigenvalue.size(); i++)
            cout << usr1.eigenVector(eigenvalue[i]) << endl;
    }

    cout << "Test3: cross product" << endl;
    {
        cin >> row1 >> col1; 
        Matrix<double> usr1(row1, col1);
        cin >> usr1;
        cin >> row2 >> col2;
        Matrix<double> usr2(row2, col2);
        cin >> usr2;
        cout << "crossproduct = \n";
        cout << usr1.crossProduct(usr2) << endl;
    }

    cout << "Test4: dot product" << endl;
    {
        cin >> row1 >> col1; 
        Matrix<double> usr1(row1, col1);
        cin >> usr1;
        cin >> row2 >> col2;
        Matrix<double> usr2(row2, col2);
        cin >> usr2;
        cout << "dotproduct = \n";
        cout << usr1.dotProduct(usr2) << endl;
    }

    cout << "Test5: conjugation" << endl;
    {
        cin >> row1 >> col1; 
        Matrix<complex<double>> usr1(row1, col1);
        cin >> usr1;
        cout << "conjugate = \n";
        cout << usr1.conjugate() << endl;
    }

    cout << "Test6: convolution" << endl;
    {
        cin >> row1 >> col1 >> channel1; 
        Matrix<double> src(row1, col1, channel1);
        cin >> src;
        cin >> row2 >> col2 >> channel2;
        Matrix<double> kernel(row2, col2, channel2);
        cin >> kernel;
        cout << "convolution = \n";
        cout << src.convolute(kernel) << endl;
    }

    
    MemoryDetech::instance().show();
}




