/*
 * @Author: Satori 3102809947@qq.com
 * @Date: 2022-05-26 11:51:58
 * @LastEditors: Satori 3102809947@qq.com
 * @LastEditTime: 2022-06-04 17:53:54
 * @FilePath: \CS205Project\SparseMatrix.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#ifndef SparseMatrix_H_
#define SparseMatrix_H_
#include <complex>
#include <vector>
#include <iostream>
#include "Matrix.h"
// #include <opencv2/opencv.hpp>
#include "Range.h"
#include "Exception.h"
#define DEBUG
#include "MemoryDetect.h"
#define new DEBUG_NEW
#define is_same(P,Q) std::is_same<P,Q>::value
namespace usr
{
#define Mp std::map<std::pair<uint32,uint32>,T>

template <class T>
class SparseMatrix {
    private:
        uint32 col;
        uint32 row;
        uint32 mod;
        std::vector<Mp>data;
        inline void allocate(const uint32 _row, const uint32 _col, const uint32 channel, const uint32 _mod);
        inline void release();
        SparseMatrix Addition(const SparseMatrix& src1, const SparseMatrix& src2) const;
        SparseMatrix Subtraction(const SparseMatrix& src1, const SparseMatrix& src2) const;
        SparseMatrix Multiplication(const SparseMatrix& src1, const SparseMatrix& src2) const;
        SparseMatrix Multiplication(const SparseMatrix& src1, const T& src2) const;
        SparseMatrix Strassen(const SparseMatrix& src1, const SparseMatrix& src2) const;
        SparseMatrix Division(const SparseMatrix& src1, const SparseMatrix& src2) const;
        SparseMatrix Division(const SparseMatrix& src1, const T& src2) const;

    public:
        /** @brief constructor
         */
        SparseMatrix(uint32 row=1,uint32 col=1,uint32 channel=1, uint32 mod=0);
        /** 
         * @description: deep copy 
         */
        SparseMatrix(const SparseMatrix&);
        /** @brief destructor
         *  release if is not empty
         */
        ~SparseMatrix();

        inline uint32 getRow()const {return row;}
        inline uint32 getCol()const {return col;}
        inline uint32 getChannel()const {return data.size();}
        inline void setMod(uint32 m) {mod = m;}
        inline void setVal(uint32 id,uint32 x,uint32 y,T v){
            try{
                if (id>=getChannel()){
                    throw(MultiChannelException("Channel greater than one is not supported", HERE));
                }
                if (x >= getRow() || y >= getCol())
                    throw(RangeOutOfBoundException("range is out of bound", HERE));
                data[id][make_pair(x,y)]=v;
            }
            catch(const RangeOutOfBoundException& e)
            {
                std::cerr << "RangeOutOfBoundException: " << e.what() << '\n';
            }
            catch(const MultiChannelException& e)
            {
                std::cerr << "MultiChannelException: " << e.what() << '\n';
            }
        }
        
        /** 
         * return true if the matrix is empty (data is empty)
         */
        inline bool isEmpty()const {return data.empty();};

        /**
         * return true if the matrix is a vector
         */
        inline bool isVector()const {return row==1||col==1;};

        /**
         * return true if the matrix is square
         */
        inline bool isSquare()const {return row==col;}; 

        /**
         * return true if the matrix is invertible
         */
        inline bool isInvertible()const {return determinant()!=0;};

        /**
         * @description: deep copy
         */
        void clone(const SparseMatrix& m);

        SparseMatrix operator+(const SparseMatrix& b);
        SparseMatrix operator-(const SparseMatrix& b);
        SparseMatrix operator*(const SparseMatrix& b);
        SparseMatrix operator*(const T& b);
        template <class _T>
        friend SparseMatrix<_T>* operator*(const _T& a, const SparseMatrix<_T>& b);    
        SparseMatrix operator/(const SparseMatrix& b);
        SparseMatrix operator/(const T& b);
        SparseMatrix operator^(int b);
        void operator=(const SparseMatrix& m);
		void operator+=(const SparseMatrix& m);
		void operator-=(const SparseMatrix& m);
		void operator*=(const SparseMatrix& m);
        void operator*=(const T t);
        template <typename _T>
        operator SparseMatrix<_T>();
        template <typename _T>
        operator usr::Matrix<_T>();
        bool operator==(const SparseMatrix& m)const;
        Mp operator[](uint32 i){ return data.at(i); }
        template <class _T>
		friend std::ostream& operator<<(std::ostream&, const SparseMatrix<_T>&);

        T determinant()const;
        uint32 rank()const;
        T trace()const;
        T max(const Range &row=Range::all(), const Range &col=Range::all())const;
        T min(const Range &row=Range::all(), const Range &col=Range::all())const;
        T avg(const Range &row=Range::all(), const Range &col=Range::all())const;
        T sum(const Range &row=Range::all(), const Range &col=Range::all())const;

        /** @brief Computes a dot-product of two vectors.
        The method computes a dot-product of two matrices. If the matrices are not single-column or
        single-row vectors, the top-to-bottom left-to-right scan ordering is used to treat them as 1D
        vectors. The vectors must have the same size and type.
        */        
        T dotProduct(const SparseMatrix&)const;

        /** @brief Computes a cross-product of two 3-element vectors.
        The method computes a cross-product of two 3-element vectors. The vectors must be
        3-element vectors of the same shape and size. The result is another 3-element vector
        of the same shape and type as operands.
        */
        SparseMatrix crossProduct(const SparseMatrix&)const;

        std::vector<std::complex<double>> eigenValue()const;
        Matrix<std::complex<double>> eigenVector(std::complex<double>)const;
        SparseMatrix subMatrix(const Range& row=Range::all(), const Range& col=Range::all())const;
        SparseMatrix inverse()const;
        SparseMatrix transpose()const;
        SparseMatrix conjugate()const;
        SparseMatrix toDiagnal()const;
        SparseMatrix reshape(const uint32 _row, const uint32 _col)const;
        SparseMatrix slice(const Range& row, const uint32 col)const;
        SparseMatrix convolute(const SparseMatrix& kernal,uint32 anchor_x,uint32 anchor_y)const;

};

template <class T>
inline void SparseMatrix<T>::allocate(const uint32 _row, const uint32 _col, const uint32 channel, const uint32 _mod)
{
    row = _row, col = _col, mod = _mod;
    for (int i=0;i<channel;i++){
        data.push_back(Mp());
    }
}

template <class T>
SparseMatrix<T>::SparseMatrix(uint32 _row, uint32 _col, uint32 _channel, uint32 _mod)
{
    if (_row==0 || _col==0 || _channel==0)
        throw(InvalidArgsException("row or col or channel can not be zero", HERE));
    allocate(_row, _col, _channel, _mod);
}

template <class T>
SparseMatrix<T>::SparseMatrix(const SparseMatrix<T>& m)
{
    allocate(m.row, m.col, m.getChannel(), m.mod);
    clone(m);
}

template <class T>
SparseMatrix<T>::~SparseMatrix()
{
    release();
}

template <class T>
inline void SparseMatrix<T>::release()
{
    if (isEmpty())
        return;
    data.clear();
    row = 0, col = 0, mod = 0;
}

template <class T>
void SparseMatrix<T>::clone(const SparseMatrix<T>& m)
{
    for (int k = 0; k < m.getChannel(); k++)
        data[k] = const_cast<SparseMatrix<T>&>(m)[k];
}

template <class T>
SparseMatrix<T> SparseMatrix<T>::Addition(const SparseMatrix<T>& src1, const SparseMatrix<T>& src2) const
{
    try
    {
        if (src1.isEmpty() || src2.isEmpty())
            throw(EmptyMatrixException("SparseMatrix is empty", HERE)); 
        else if (src1.getChannel()!=src2.getChannel())
            throw(ChannelMismatchException("Both of the Matrices must have the same Channels", HERE));
        else if (src1.getRow()!=src2.getRow() || src1.getCol()!= src2.getCol())
            throw(SizeMismatchException("Matrices do not match in size", HERE)); 
        else
        {
            uint32 r = src1.getRow(), c = src1.getCol(), channel = src1.getChannel();
            SparseMatrix<T> rst(r, c, channel);
            for (int k = 0; k < channel; k++){
                for (auto it=const_cast<SparseMatrix<T>&>(src1)[k].begin();it!=const_cast<SparseMatrix<T>&>(src1)[k].end();){
                    rst[k][it->first]+=it->second;
                    it++;
                }
                for (auto it=const_cast<SparseMatrix<T>&>(src2)[k].begin();it!=const_cast<SparseMatrix<T>&>(src2)[k].end();it++){
                    rst[k][it->first]+=it->second;
                }
            }
            return rst;
        }
    }
    catch(const EmptyMatrixException& e)
    {
        std::cerr << "EmptyMatrixException: " << e.what() << '\n';
    }
    catch(const ChannelMismatchException& e)
    {
        std::cerr << "ChannelMismatchException: " << e.what() << '\n';
    }
    catch(const SizeMismatchException& e)
    {
        std::cerr << "SizeMismatchException: " << e.what() << '\n';
    }
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return SparseMatrix<T>();
}

template <class T>
SparseMatrix<T> SparseMatrix<T>::Subtraction(const SparseMatrix<T>& src1, const SparseMatrix<T>& src2) const
{
    try
    {
        if (src1.isEmpty() || src2.isEmpty())
            throw(EmptyMatrixException("SparseMatrix is empty", HERE)); 
        else if (src1.getChannel()!=src2.getChannel())
            throw(ChannelMismatchException("Both of the Matrices must have the same Channels", HERE));
        else if (src1.getRow()!=src2.getRow() || src1.getCol()!= src2.getCol())
            throw(SizeMismatchException("Matrices do not match in size", HERE)); 
        else
        {
            uint32 r = src1.getRow(), c = src1.getCol(), channel = src1.getChannel();
            SparseMatrix<T> rst(r, c, channel);
            for (int k = 0; k < channel; k++){
                for (auto it=const_cast<SparseMatrix<T>&>(src1)[k].begin();it!=const_cast<SparseMatrix<T>&>(src1)[k].end();it++){
                    rst[k][it->first]+=it->second;
                }
                for (auto it=const_cast<SparseMatrix<T>&>(src2)[k].begin();it!=const_cast<SparseMatrix<T>&>(src2)[k].end();it++){
                    rst[k][it->first]-=it->second;
                }
            }
            return rst;
        }
    }
    catch(const EmptyMatrixException& e)
    {
        std::cerr << "EmptyMatrixException: " << e.what() << '\n';
    }
    catch(const ChannelMismatchException& e)
    {
        std::cerr << "ChannelMismatchException: " << e.what() << '\n';
    }
    catch(const SizeMismatchException& e)
    {
        std::cerr << "SizeMismatchException: " << e.what() << '\n';
    }
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return SparseMatrix<T>();
}

template <class T>
SparseMatrix<T> SparseMatrix<T>::Multiplication(const SparseMatrix<T>& src1, const SparseMatrix<T>& src2) const
{
    try
    {
        if (src1.isEmpty() || src2.isEmpty())
            throw(EmptyMatrixException("SparseMatrix is empty", HERE)); 
        else if (src1.getChannel()!=src2.getChannel())
            throw(ChannelMismatchException("Both of the Matrices must have the same Channels", HERE));
        else if (src1.getChannel()>1)
            throw(MultiChannelException("Channels greater than one is not supported", HERE));
        else if (src1.getCol() != src2.getRow())
            throw(SizeMismatchException("Matrices do not match in size", HERE));
        else
        {   
            Matrix<T> a=src1;
            Matrix<T> b=src2;
            SparseMatrix<T> res=Multiplication(a,b);
            return res;
        }
    }
    catch(const EmptyMatrixException& e)
    {
        std::cerr << "EmptyMatrixException: " << e.what() << '\n';
    }
    catch(const SizeMismatchException& e)
    {
        std::cerr << "SizeMismatchException: " << e.what() << '\n';
    }
    catch(const ChannelMismatchException& e)
    {
        std::cerr << "ChannelMismatchException: " << e.what() << '\n';
    }
    catch(const MultiChannelException& e)
    {
        std::cerr << "MultiChannelException: " << e.what() << '\n';
    }
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return SparseMatrix<T>();
}

template <class T>
SparseMatrix<T> SparseMatrix<T>::Multiplication(const SparseMatrix<T>& src1, const T& src2) const
{
    try
    {
        if (src1.isEmpty())
            throw(EmptyMatrixException("SparseMatrix is empty", HERE)); 
        else
        {   
            Matrix<T> a=src1;
            SparseMatrix<T> res=Multiplication(a,src2);
            return res;
        }
    }
    catch(const EmptyMatrixException& e)
    {
        std::cerr << "EmptyMatrixException: " << e.what() << '\n';
    }
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return SparseMatrix<T>();
}

template <class T>
SparseMatrix<T> SparseMatrix<T>::Division(const SparseMatrix<T>& src1, const SparseMatrix<T>& src2) const
{
    try
    {
        if (src1.isEmpty() || src2.isEmpty())
            throw(EmptyMatrixException("SparseMatrix is empty", HERE)); 
        else if (src1.getChannel()!=src2.getChannel())
            throw(ChannelMismatchException("Both of the Matrices must have the same Channels", HERE));
        else if (src1.getCol()!=src1.getCol() || src1.getRow()!=src2.getRow())
            throw(SizeMismatchException("Matrices do not match in size", HERE));
        else
        {
            Matrix<T> a=src1;
            Matrix<T> b=src2;
            SparseMatrix<T>res=Division(a,b);
            return res;
        }
    }
    catch(const EmptyMatrixException& e)
    {
        std::cerr << "EmptyMatrixException: " << e.what() << '\n';
    }
    catch(const SizeMismatchException& e)
    {
        std::cerr << "SizeMismatchException: " << e.what() << '\n';
    }
    catch(const ChannelMismatchException& e)
    {
        std::cerr << "ChannelMismatchException: " << e.what() << '\n';
    }
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return SparseMatrix<T>();
}

template <class T>
SparseMatrix<T> SparseMatrix<T>::Division(const SparseMatrix<T>& src1, const T& src2) const
{
    try
    {
        if (src1.isEmpty())
            throw(EmptyMatrixException("SparseMatrix is empty", HERE)); 
        else if (src2 == 0)
            throw(ArithmeticException("Divide by Zero", HERE));
        else
        {
            Matrix<T> a=src1;
            SparseMatrix<T>res=Division(a,src2);
            return res;
        }
    }
    catch(const EmptyMatrixException& e)
    {
        std::cerr << "EmptyMatrixException: " << e.what() << '\n';
    }
    catch(const ArithmeticException& e)
    {
        std::cerr << "ArithmeticException: " << e.what() << '\n';
    }
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return SparseMatrix<T>();
}

template <class T>
SparseMatrix<T> SparseMatrix<T>::operator+(const SparseMatrix<T>& m)
{
    return Addition(*this, m);
}

template <class T>
SparseMatrix<T> SparseMatrix<T>::operator-(const SparseMatrix<T>& m)
{
    return Subtraction(*this, m);
}

template <class T>
SparseMatrix<T> SparseMatrix<T>::operator*(const SparseMatrix<T>& m)
{
    return Multiplication(*this, m);
}

template <class T>
SparseMatrix<T> SparseMatrix<T>::operator*(const T& a)
{
    return Multiplication(*this, a);
}

template <class T>
SparseMatrix<T> operator*(const T& a, const SparseMatrix<T>& b)
{
    return b*a;
}    

template <class T>
SparseMatrix<T> SparseMatrix<T>::operator/(const SparseMatrix<T>& m)
{
    return Division(*this, m);
}

template <class T>
SparseMatrix<T> SparseMatrix<T>::operator/(const T& a)
{
    return Division(*this, a);
}

template <class T>
SparseMatrix<T> SparseMatrix<T>::operator^(int b)
{
    try
    {
        if (isEmpty())
            throw(EmptyMatrixException("SparseMatrix is empty, can't find determinant", HERE)); 
        else if (!isSquare())
            throw(MatrixNotSquareException("SparseMatrix must be square", HERE));
        else if (getChannel()>2)
            throw(MultiChannelException("Channels greater than two is not supported", HERE));
        else
        {
            Matrix<T> res=*this; 
            SparseMatrix<T> ans=res^b;
            return ans;
        }
    }
    catch(const EmptyMatrixException& e)
    {
        std::cerr << "EmptyMatrixException: " << e.what() << '\n';
    }
    catch(const MatrixNotSquareException& e)
    {
        std::cerr << "MatrixNotSquareException: " << e.what() << '\n';
    }
    catch(const MultiChannelException& e)
    {
        std::cerr << "MultiChannelException: " << e.what() << '\n';
    }
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return SparseMatrix<T>();
}

template <class T>
void SparseMatrix<T>::operator=(const SparseMatrix<T>& m)
{
    if (this != &m)
    {
        if (!(row==m.row && col==m.col && getChannel()==m.getChannel()))
        {
            release();
            allocate(m.row, m.col, m.getChannel() ,m.mod);
        }
        clone(m);
    }
}

template <class T>
void SparseMatrix<T>::operator+=(const SparseMatrix<T>& m)
{
    (*this) = (*this)+m;
}

template <class T>
void SparseMatrix<T>::operator-=(const SparseMatrix<T>& m)
{
    (*this) = (*this)-m;
}

template <class T>
void SparseMatrix<T>::operator*=(const SparseMatrix<T>& m)
{
    (*this) = (*this)*m;
}

template <class T>
void SparseMatrix<T>::operator*=(const T t)
{
    (*this) = (*this)*t;
}

template <class T>
bool SparseMatrix<T>::operator==(const SparseMatrix<T>& m)const
{
    if(row!=m.getRow()||col!=m.getCol()||getChannel()!=m.getChannel()) return 0;
    Matrix<T> x=this;
    Matrix<T> y=m;
    return x==y;
}

template <class T>
template <typename _T>
SparseMatrix<T>::operator SparseMatrix<_T>()
{
    SparseMatrix<_T> t(row, col, getChannel());
    for (int i = 0; i < getChannel(); i++){
        for (auto it=data[i].begin();it!=data[i].end();it++){
            t[i][it->first]=static_cast<_T>(it->second);
        }
    }
    return t;
}

template <class T>
template <typename _T>
SparseMatrix<T>::operator Matrix<_T>()
{
    Matrix<_T> t(row, col, getChannel());
    for (int i = 0; i < getChannel(); i++){
        for (auto it=data[i].begin();it!=data[i].end();it++){
            t[i][it->first.first][it->first.second]=static_cast<_T>(it->second);
        }
    }
    return t;
}

template <class _T>
std::ostream& operator<<(std::ostream& os, const SparseMatrix<_T>& m)
{
    for (int i = 0; i < m.getChannel(); i++)
    {
        for (auto j = const_cast<SparseMatrix<_T>&>(m)[i].begin(); j != const_cast<SparseMatrix<_T>&>(m)[i].end(); j++)
        {
            os << j->first.first << " " << j->first.second << " " << j->second << " ";
        }
    }
    return os;   
}


template<class T>
SparseMatrix<T> SparseMatrix<T>::toDiagnal()const
{
    Matrix<T>a=*this;
    SparseMatrix<T>res=a.toDiagnal();
    return res;
}
template<class T>
T SparseMatrix<T>::determinant()const
{
    try
    {
        if (isEmpty())
            throw(EmptyMatrixException("SparseMatrix is empty, can't find determinant", HERE)); 
        else if (getChannel()>1)
            throw(MultiChannelException("Channel greater than one is not supported", HERE));
        else if (!isSquare())
            throw(MatrixNotSquareException("SparseMatrix must be square", HERE));
        else 
        {   
            Matrix<T>a=*this;
            return a.determinant();
        }
    }
    catch(const EmptyMatrixException& e)
    {
        std::cerr << "EmptyMatrixException: " << e.what() << '\n';
    }
    catch(const MultiChannelException& e)
    {
        std::cerr << "MultiChannelException: " << e.what() << '\n';
    }
    catch(const MatrixNotSquareException& e)
    {
        std::cerr << "MatrixNotSquareException: " << e.what() << '\n';
    }
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return 0;
}

template <class T>
uint32 SparseMatrix<T>::rank()const{
    try 
    {
        if (isEmpty())
            throw(EmptyMatrixException("SparseMatrix is empty, can't find rank", HERE)); 
        else if (getChannel()>1)
            throw(MultiChannelException("Channel greater than one is not supported", HERE));
        
        Matrix<T>a=*this;
        return a.rank();
    }
    catch(const EmptyMatrixException& e)
    {
        std::cerr << "EmptyMatrixException: " << e.what() << '\n';
    }
    catch(const MultiChannelException& e)
    {
        std::cerr << "MultiChannelException: " << e.what() << '\n';
    }
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return 0;
}

template<class T>
T SparseMatrix<T>::trace()const
{
    try
    {
        if (isEmpty())
            throw(EmptyMatrixException("SparseMatrix is empty, can't find trace", HERE)); 
        else if (getChannel()>1)
            throw(MultiChannelException("Channel greater than one is not supported", HERE));
        else if (!isSquare())
            throw(MatrixNotSquareException("SparseMatrix must be square", HERE));
        else 
        {   
            Matrix<T> a=*this;
            return a.trace();
        }
    }
    catch(const EmptyMatrixException& e)
    {
        std::cerr << "EmptyMatrixException: " << e.what() << '\n';
    }
    catch(const MultiChannelException& e)
    {
        std::cerr << "MultiChannelException: " << e.what() << '\n';
    }
    catch(const MatrixNotSquareException& e)
    {
        std::cerr << "MatrixNotSquareException: " << e.what() << '\n';
    }
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return 0;
}

template <class T>
T SparseMatrix<T>::max(const Range &row, const Range &col)const
{
    try
    {
        if (row.empty() || col.empty())
            throw(InvalidArgsException("range can not be zero", HERE));
        else if (row.start >= getRow() || col.start >= getCol())
            throw(RangeOutOfBoundException("range is out of bound", HERE));
        else
        {
            Matrix<T> a=*this;
            return a.max(row,col);
        }
    }
    catch(const InvalidArgsException& e)
    {
        std::cerr << "InvalidArgsException: " << e.what() << '\n';
    }
    catch(const RangeOutOfBoundException& e)
    {
        std::cerr << "RangeOutOfBoundException: " << e.what() << '\n';
    }
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return 0;
}

template <class T>
T SparseMatrix<T>::min(const Range &row, const Range &col)const
{
    try
    {
        if (row.empty() || col.empty())
            throw(InvalidArgsException("range can not be zero", HERE));
        else if (row.start >= getRow() || col.start >= getCol())
            throw(RangeOutOfBoundException("range is out of bound", HERE));
        else
        {
            Matrix<T> a=*this;
            return a.min(row,col);
        }
    }
    catch(const InvalidArgsException& e)
    {
        std::cerr << "InvalidArgsException: " << e.what() << '\n';
    }
    catch(const RangeOutOfBoundException& e)
    {
        std::cerr << "RangeOutOfBoundException: " << e.what() << '\n';
    }
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return 0;
}

template <class T>
T SparseMatrix<T>::avg(const Range &row, const Range &col)const
{
    try
    {
        if (row.empty() || col.empty())
            throw(InvalidArgsException("range can not be zero", HERE));
        else if (row.start >= getRow() || col.start >= getCol())
            throw(RangeOutOfBoundException("range is out of bound", HERE));
        else
        {
            Matrix<T> a=*this;
            return a.avg(row,col);
        }
    }
    catch(const InvalidArgsException& e)
    {
        std::cerr << "InvalidArgsException: " << e.what() << '\n';
    }
    catch(const RangeOutOfBoundException& e)
    {
        std::cerr << "RangeOutOfBoundException: " << e.what() << '\n';
    }
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return 0;
}

template <class T>
T SparseMatrix<T>::sum(const Range &row, const Range &col)const
{
    try
    {
        if (row.empty() || col.empty())
            throw(InvalidArgsException("range can not be zero", HERE));
        else if (row.start >= getRow() || col.start >= getCol())
            throw(RangeOutOfBoundException("range is out of bound", HERE));
        else
        {
            Matrix<T> a=*this;
            return a.sum(row,col);
        }
    }
    catch(const InvalidArgsException& e)
    {
        std::cerr << "InvalidArgsException: " << e.what() << '\n';
    }
    catch(const RangeOutOfBoundException& e)
    {
        std::cerr << "RangeOutOfBoundException: " << e.what() << '\n';
    }
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return 0;
}
template <class T>
T SparseMatrix<T>::dotProduct(const SparseMatrix<T>& m)const
{
    try{
        if(row*col!=m.getRow()*m.getCol())
            throw(SizeMismatchException("the total size must be equal", HERE));
        else if (getChannel() != m.getChannel())
            throw(ChannelMismatchException("Both of the Matrices must have the same Channels", HERE));
        else
        {
            Matrix<T> a=*this;
            Matrix<T> b=m;
            return a.dotProduct(b);
        }
    }catch(const SizeMismatchException& e)
    {
        std::cerr << "SizeMismatchException: "<<e.what()<<'\n';
    }
    catch(const ChannelMismatchException& e)
    {
        std::cerr << "ChannelMismatchException: " << e.what() << '\n';
    }
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return 0;
}
template <class T>
SparseMatrix<T> SparseMatrix<T>::crossProduct(const SparseMatrix<T>& m)const
{
    try{
        if(row!=1||col!=3||m.getRow()!=1||m.getCol()!=3) 
            throw(SizeMismatchException("row of both matrix must be 1 and col must be 3", HERE));
        else if (getChannel()!=m.getChannel())
            throw(ChannelMismatchException("Both of the Matrices must have the same Channels", HERE));
        else if (getChannel()>1)
            throw(MultiChannelException("Channel greater than one is not supported", HERE));
        else
        {
            Matrix<T> a=*this;
            Matrix<T> b=m;
            return a.crossProduct(b);
        }
    }catch(SizeMismatchException& e)
    {
        std::cerr << "SizeMismatchException: "<<e.what()<<'\n';
    }
    catch(const ChannelMismatchException& e)
    {
        std::cerr << "ChannelMismatchException: " << e.what() << '\n';
    }
    catch(const MultiChannelException& e)
    {
        std::cerr << "MultiChannelException: " << e.what() << '\n';
    }
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return SparseMatrix<T>();
}
template <class T>
SparseMatrix<T> SparseMatrix<T>::transpose()const
{   
    try{
        if (isEmpty())
            throw(EmptyMatrixException("SparseMatrix is empty, can't transpose", HERE)); 
        else
        {
            Matrix<T> a=*this;
            SparseMatrix<T> res=a.transpose();
            return res;
        }
    }
    catch(const EmptyMatrixException& e)
    {
        std::cerr << "EmptyMatrixException: " << e.what() << '\n';
    }
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return SparseMatrix<T>();
}
template <class T>
SparseMatrix<T> SparseMatrix<T>::conjugate()const
{
    try{
        if (isEmpty())
            throw(EmptyMatrixException("SparseMatrix is empty, can't conjugate", HERE)); 
        else
        {
            Matrix<T> a=*this;
            SparseMatrix<T> res=a.conjugate();
            return res;
        }
    }
    catch(const EmptyMatrixException& e)
    {
        std::cerr << "EmptyMatrixException: " << e.what() << '\n';
    }
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return SparseMatrix<T>();
}
template <class T>
SparseMatrix<T> SparseMatrix<T>::inverse()const
{
    try{
        if (isEmpty())
            throw(EmptyMatrixException("SparseMatrix is empty, can't find inverse", HERE)); 
        else if (!isSquare())
            throw(MatrixNotSquareException("non-square matrix is not invertible", HERE));
        else if (!isInvertible())
            throw(InverseNotExistException("singular matrix is not invertible", HERE));
        else
        {
            Matrix<T> a=*this;
            SparseMatrix<T> res=a.inverse();
            return res;
        }
    }
    catch(const EmptyMatrixException& e)
    {
        std::cerr << "EmptyMatrixException: " << e.what() << '\n';
    }
    catch(const MatrixNotSquareException& e)
    {
        std::cerr << "MatrixNotSquareException: " << e.what() << '\n';
    }
    catch(const InverseNotExistException& e)
    {
        std::cerr << "InverseNotExistException: " << e.what() <<'\n';
    }
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return SparseMatrix<T>();
}
template <class T>
SparseMatrix<T> SparseMatrix<T>::subMatrix(const Range& row, const Range& col)const
{
    try
    {
        if (row.empty() || col.empty())
            throw(InvalidArgsException("range can not be zero", HERE));
        else if (row.start >= getRow() || col.start >= getCol())
            throw(RangeOutOfBoundException("range is out of bound", HERE));
        else
        {
            
            Matrix<T> a=*this;
            SparseMatrix<T> res=a.subMatrix(row,col);
            return res;
        }
    }
    catch(const InvalidArgsException& e)
    {
        std::cerr << "InvalidArgsException: " << e.what() << '\n';
    }
    catch(const RangeOutOfBoundException& e)
    {
        std::cerr << "RangeOutOfBoundException: " << e.what() << '\n';
    }
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return SparseMatrix<T>();
}
template <class T>
SparseMatrix<T> SparseMatrix<T>::reshape(const uint32 _row, const uint32 _col)const
{
    try{
        if(1ll*_row*_col!=1ll*row*col)
            throw(SizeMismatchException("cannot fit into this shape", HERE));
        else
        {
            Matrix<T> a=*this;
            SparseMatrix<T> res=a.reshape(row,col);
            return res;
        }
    }
    catch(const SizeMismatchException& e)
    {
        std::cerr << "SizeMismatchException: " << e.what() <<'\n';
    }
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return SparseMatrix<T>();
}
template <class T>
SparseMatrix<T> SparseMatrix<T>::slice(const Range& row, const uint32 col)const
{
    try{
        if (row.empty())
            throw(InvalidArgsException("range can not be zero", HERE));
        if(col>=this->col||row.start>=this->row)
            throw(RangeOutOfBoundException("range out of bound", HERE));
        else
        {
            Matrix<T> a=*this;
            SparseMatrix<T> res=a.slice(row,col);
            return res;
        }
    }
    catch(const InvalidArgsException& e)
    {
        std::cerr << "InvalidArgsException: " << e.what() << '\n';
    }
    catch(const RangeOutOfBoundException& e)
    {
        std::cerr << "RangeOutOfBoundException: " << e.what() <<'\n';
    }
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return SparseMatrix<T>();
}

template <class T>
SparseMatrix<T> SparseMatrix<T>::convolute(const SparseMatrix& kernal,uint32 anchor_x,uint32 anchor_y)const
{
    try{
        if (isEmpty()||kernal.isEmpty())
            throw(EmptyMatrixException("SparseMatrix is empty", HERE)); 
        else if (getChannel()!=kernal.getChannel())
            throw(ChannelMismatchException("Both of the Matrices must have the same Channels", HERE));
        else if(anchor_x>=kernal.getRow()||anchor_y>=kernal.getCol())
            throw(RangeOutOfBoundException("range out of bound", HERE));
        else{
            Matrix<T> a=*this;
            Matrix<T> b=kernal;
            SparseMatrix<T> res=a.convolute(b,row,col);
            return res;
        }
    }
    catch(const EmptyMatrixException& e)
    {
        std::cerr << "EmptyMatrixException: " << e.what() << '\n';
    }
    catch(const RangeOutOfBoundException& e)
    {
        std::cerr << "RangeOutOfBoundException: " << e.what() <<'\n';
    }
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return SparseMatrix<T>();
}

template <class T>
std::vector<std::complex<double>> SparseMatrix<T>::eigenValue() const
{
    try{
        if (isEmpty())
            throw(EmptyMatrixException("SparseMatrix is empty, can't find eigenvalue", HERE)); 
        else if (!isSquare())
            throw(MatrixNotSquareException("Non-square matrix does not have eigenValue", HERE));
        else if (getChannel()>1)
            throw(MultiChannelException("Channels greater than one is not supported", HERE));
        else
        {
            
            Matrix<T> a=*this;
            return a.eigenValue();
        }
    }
    catch(const EmptyMatrixException& e)
    {
        std::cerr << "EmptyMatrixException: " << e.what() << '\n';
    }
    catch(const MatrixNotSquareException& e)
    {
        std::cerr << "MatrixNotSquareException: " << e.what() << '\n';
    }
    catch(const MultiChannelException& e)
    {
        std::cerr << "MultiChannelException: " << e.what() << '\n';
    }
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return std::vector<std::complex<double>>();
}

template <class T>
Matrix<std::complex<double>> SparseMatrix<T>::eigenVector(std::complex<double> val)const{
    
    try{
        Matrix<T> a=*this;
        return a.eigenVector(val);
    }
    catch(const WrongEigenValueException& e)
    {
        std::cerr << "WrongEigenValueException: " << e.what() << '\n';
    }
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return Matrix<double>();
}
}

#undef Mp
#endif