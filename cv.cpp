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
using namespace cv;
int main()
{
    freopen("data.in", "r", stdin);
    uint32 row1, row2, col1, col2, channel1, channel2;
    cout << "Test1: basic operation" << endl;
    {
        cin >> row1 >> col1; 
        usr::Matrix<double> temp(row1, col1);
        cin >> temp;
        Mat cv1(row1, col1, CV_64FC1);
        usr::usrTocv(temp, cv1);
        cout<<"determinant = "<< determinant(cv1)<<endl;
        // cout<<"rank = "<<cv1.rank()<<endl;
        cout<<"trace = "<<trace(cv1)<<endl;
        double min, max;
        minMaxLoc(cv1, &min, &max);
        cout<<"max = "<<max<<endl;
        cout<<"min = "<<min<<endl;
        cout<<"sum = "<<sum(cv1)<<endl;
        cout<<"avg = "<<mean(cv1)<<endl;
        transpose(cv1, cv1);
        cout<<"transpose = \n"<<cv1<<endl<<endl;
        cout<<"inverse = \n"<<cv1.inv()<<endl<<endl;
        // cout<<"gaussian elimination = \n"<<cv1.toDiagnal()<<endl;
        cout<<"submatrix = \n"<<cv1(Range(0, 2), Range(col1/2,col1))<<endl<<endl;
        cout<<"slice = \n"<<cv1(Rect(0,0,2, col1/2))<<endl<<endl;
        cout<<"reshape = \n"<<cv1.reshape(0,row1/2)<<endl;
        putchar('\n');
    }

    cout << "Test2: eigenvalue and eigenvector" << endl;
    {
        cin >> row1 >> col1; 
        usr::Matrix<double> temp(row1, col1);
        cin >> temp;
        Mat cv1(row1, col1, CV_64FC1);
        usr::usrTocv(temp, cv1);
        Mat cv2;
        transpose(cv1, cv2);
        cv1 = cv1 + cv2;
        Mat eValuesMat;
        Mat eVectorsMat;
        eigen(cv1, eValuesMat, eVectorsMat);
        cout<<"eigenvalue = \n";
        for(auto i=0; i<eValuesMat.rows; i++)
            for(auto j=0; j<eValuesMat.cols; j++)
                cout << eValuesMat.at<double>(i,j) << " ";
        putchar('\n');
        cout<<"eigenvector = \n";
        for(auto i=0; i<eVectorsMat.rows; i++)
            for(auto j=0; j<eVectorsMat.cols; j++)
                cout << eVectorsMat.at<double>(i,j) << " ";
        putchar('\n');
    }

    cout << "\nTest3: cross product" << endl;
    {
        cin >> row1 >> col1; 
        usr::Matrix<double> temp(row1, col1);
        cin >> temp;
        Mat cv1(row1, col1, CV_64FC1);
        usr::usrTocv(temp, cv1);
        cin >> row2 >> col2;
        temp = usr::Matrix<double>(row2, col2);
        cin >> temp;
        Mat cv2(row2, col2, CV_64FC1);
        usr::usrTocv(temp, cv2);
        cout << "crossproduct = \n";
        cout << cv1.cross(cv2) << endl;
        putchar('\n');
    }

    cout << "Test4: dot product" << endl;
    {
        cin >> row1 >> col1; 
        usr::Matrix<double> temp(row1, col1);
        cin >> temp;
        Mat cv1(row1, col1, CV_64FC1);
        usr::usrTocv(temp, cv1);
        cin >> row2 >> col2;
        temp = usr::Matrix<double>(row2, col2);
        cin >> temp;
        Mat cv2(row2, col2, CV_64FC1);
        usr::usrTocv(temp, cv2);
        cout << "dotproduct = \n";
        cout << cv1.dot(cv2) << endl;
        putchar('\n');
    }

    cout << "Test5: convolution" << endl;
    {
        cin >> row1 >> col1 >> channel1; 
        usr::Matrix<double> temp(row1, col1, channel1);
        cin >> temp;
        Mat src(row1, col1, CV_64FC(channel1));
        usr::usrTocv(temp, src);
        cin >> row2 >> col2 >> channel2;
        temp = usr::Matrix<double>(row2, col2, channel2);
        cin >> temp;
        Mat kernel(row2, col2, CV_64FC(channel2));
        usr::usrTocv(temp, kernel);
        cout << "convolution = \n";
        Mat dst;
        filter2D(src, dst, src.depth(), kernel, Point(-1,-1));
        cout << dst << endl;
        putchar('\n');
    }

    MemoryDetech::instance().show();
}





