#ifndef INTERFACE
#define INTERFACE
#include <opencv2/opencv.hpp>
#include "Matrix.h"

namespace usr
{

template <typename _Tp> static inline
void dataCopy(const cv::Mat& src, Matrix<_Tp>& dst)
{
    for(int i = 0; i < src.rows; i++)
    {
        const _Tp * ptr = src.ptr<_Tp>(i);
        for(int j = 0; j < src.cols * src.channels(); j++)
            dst[j%src.channels()][i][j/src.channels()] = ptr[j];
    } 
}

template <typename _Tp> static inline
void cvTousr(const cv::Mat& src, usr::Matrix<_Tp>& dst)
{
    int r = src.rows, c = src.cols, channel = src.channels();
    int type = src.type()&7;
    if (type==0)
    {
        usr::Matrix<uchar> temp1(r,c,channel);
        dataCopy<uchar>(src, temp1);
        dst = (usr::Matrix<_Tp>)temp1;
    }
    else if (type==1)
    {
        usr::Matrix<schar> temp2(r,c,channel);
        dataCopy(src, temp2);
        dst = (usr::Matrix<_Tp>)temp2;
    }
    else if (type==2)
    {
        usr::Matrix<ushort> temp3(r,c,channel);
        dataCopy(src, temp3);
        dst = (usr::Matrix<_Tp>)temp3;
    }
    else if (type==3)
    {
        usr::Matrix<short> temp4(r,c,channel);
        dataCopy(src, temp4);
        dst = (usr::Matrix<_Tp>)temp4;
    }
    else if (type==4)
    {
        usr::Matrix<int> temp5(r,c,channel);
        dataCopy(src, temp5);
        dst = (usr::Matrix<_Tp>)temp5;
    }
    else if (type==5)
    {
        usr::Matrix<float> temp6(r,c,channel);
        dataCopy(src, temp6);
        dst = (usr::Matrix<_Tp>)temp6;
    }                     
    else if (type==6)
    {
        usr::Matrix<double> temp7(r,c,channel);
        dataCopy(src, temp7);
        dst = (usr::Matrix<_Tp>)temp7;
    }
}

template<typename _Tp1, typename _Tp2, int _rows, int _cols> static inline
void usrTocv(const usr::Matrix<_Tp1>& src, cv::Matx<_Tp2, _rows, _cols>& dst)
{
    if (_rows != src.getRow() || _cols != src.getCol())
        throw(SizeMismatchException("Two matrices must have the same size", HERE));
    _Tp2* p = new _Tp2[_rows*_cols];
    for (int i = 0; i < _rows; i++)
        for (int j = 0; j < _cols; j++)
            p[i*_cols+j] = static_cast<_Tp2>(const_cast<usr::Matrix<_Tp1>&>(src)[0][i][j]);
    dst = cv::Matx<_Tp2, _rows, _cols>(p);
    delete p;
}

template<typename _Tp1, typename _Tp2> inline
void convert(_Tp1* src, const _Tp2* dst, int size)
{
    for (int i = 0; i < size; i++)
        src[i] = static_cast<_Tp1>(dst[i]);
}

template<typename _Tp> static inline
void usrTocv(const usr::Matrix<_Tp>& src, cv::Mat& dst)
{
    if (src.getChannel() < dst.channels())
        throw(ChannelMismatchException("usrMatrix's channel must greater than or equal to cvMatrix's", HERE));
    uint32 size = src.getRow()*src.getCol()*src.getChannel();
    _Tp array[size];
    for (int i = 0; i < src.getRow(); i++)
    {
        for (int j = 0; j < src.getCol(); j++)
        {
            for (int k = 0; k < dst.channels(); k++)
            {
                array[i*src.getCol()*dst.channels()+j*dst.channels()+k] = const_cast<Matrix<_Tp>&>(src)[k][i][j];
            }
        }
    }
    int r = dst.rows, c = dst.cols, channel = dst.channels();
    int type = dst.type()&7;
    if (type==0)
    {
        uchar arr[size];
        convert(arr, array, size);
        cv::Mat temp = cv::Mat(src.getRow(), src.getCol(), dst.type(), arr);
        temp.copyTo(dst);
    }
    else if (type==1)
    {
        schar arr[size];
        convert(arr, array, size);
        cv::Mat temp = cv::Mat(src.getRow(), src.getCol(), dst.type(), arr);
        temp.copyTo(dst);
    }
    else if (type==2)
    {
        ushort arr[size];
        convert(arr, array, size);
        cv::Mat temp = cv::Mat(src.getRow(), src.getCol(), dst.type(), arr);
        temp.copyTo(dst);
    }
    else if (type==3)
    {
        short arr[size];
        convert(arr, array, size);
        cv::Mat temp = cv::Mat(src.getRow(), src.getCol(), dst.type(), arr);
        temp.copyTo(dst);
    }
    else if (type==4)
    {
        int arr[size];
        convert(arr, array, size);
        cv::Mat temp = cv::Mat(src.getRow(), src.getCol(), dst.type(), arr);
        temp.copyTo(dst);
    }
    else if (type==5)
    {
        float arr[size];
        convert(arr, array, size);
        cv::Mat temp = cv::Mat(src.getRow(), src.getCol(), dst.type(), arr);
        temp.copyTo(dst);
    }                     
    else if (type==6)
    {
        double arr[size];
        convert(arr, array, size);
        cv::Mat temp = cv::Mat(src.getRow(), src.getCol(), dst.type(), arr);
        temp.copyTo(dst);
    }
}

} // namespace usr
#endif
