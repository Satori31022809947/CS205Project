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
#include <sys/time.h>
#include <stdlib.h>
using namespace std;
using namespace usr;
int main()
{
    freopen("data.in", "r", stdin);
    srand(time(NULL));
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
        cout<<"reshape = \n"<<usr1.reshape(0,row1/2)<<endl;
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
            cout << usr1.eigenVector(eigenvalue[i]);
        putchar('\n');
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
        putchar('\n');
    }

    cout << "Test5: convolution" << endl;
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

    cout << "Test6: conjugation" << endl;
    {
        cin >> row1 >> col1; 
        Matrix<complex<double>> usr1(row1, col1);
        cin >> usr1;
        cout << "conjugate = \n";
        cout << usr1.conjugate() << endl;
    }

    cout << "Test7: sparse matrix" << endl;
    {
        cin >> row1 >> col1;
        SparseMatrix<double> usr1(row1, col1);
        double t;
        for (int i = max(row1, col1); i > max(row1, col1)-10; i--)
        {
            cin >> t;
            usr1[0][make_pair(rand()%row1,rand()%col1)] = t;
        }
        cin >> row2 >> col2;
        SparseMatrix<double> usr2(row2, col2);
        for (int i = max(row1, col1); i > max(row1, col1)-10; i--)
        {
            cin >> t;
            usr2[0][make_pair(rand()%row1,rand()%col1)] = t;
        }    
        Matrix<double> usr3 = usr1;
        Matrix<double> usr4 = usr2;
        struct timeval t1,t2;
        double timeuse;
        gettimeofday(&t1,NULL);
        cout << "sparse matrix addition = \n";
        cout << usr1 + usr2 << endl;
        gettimeofday(&t2,NULL);
        timeuse = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
        cout<<"time = "<<timeuse << "s" <<endl;
        gettimeofday(&t1,NULL);
        cout << "general matrix addition = \n";
        usr3 + usr4;
        gettimeofday(&t2,NULL);
        timeuse = (t2.tv_sec - t1.tv_sec) + (double)(t2.tv_usec - t1.tv_usec)/1000000.0;
        cout<<"time = "<<timeuse << "s" <<endl;    
        putchar('\n');    
    }

    cout << "Test8: conversion between usr and cv" << endl;
    {
        cv::Mat cv1(4,4,CV_32SC3, cv::Scalar(1,2,3));
        cout << "cv::Mat(CV_32SC3) \n" << cv1 << endl << endl;
        Matrix<double> usr1(1,1,1);
        cout << "usr::Matrix<double> \n" << usr1 << endl;
        cvTousr(cv1, usr1);
        cout << "cvTousr \n" << "usr::Matrix<double> \n" << usr1 << endl;
        cv::Matx<char, 4, 4> cv2;
        cout << "cv::Matx \n" << cv2 << endl;
        usrTocv(usr1, cv2);
        cout << "usrTocv \n" << "cv::Matx \n" << cv2 << endl << endl;
    }

    cout << "Test9: exception" << endl;
    {
        Matrix<double> usr1(2,3,1);
        Matrix<double> usr2(1,3,1);
        usr1 + usr2;
        usr1 * usr2;
        usr2 = Matrix<double>(1,3,2);
        usr1.inverse();
        usr1.crossProduct(usr2);
        usr1.convolute(usr2);
    }

    MemoryDetech::instance().show();
}




