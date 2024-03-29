/*
 * @Author: Satori 3102809947@qq.com
 * @Date: 2022-05-26 11:51:58
 * @LastEditors: Satori 3102809947@qq.com
 * @LastEditTime: 2022-06-04 17:53:54
 * @FilePath: \CS205Project\Matrix.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#ifndef Matrix_H_
#define Matrix_H_
#include <complex>
#include <vector>
#include <iostream>
#include "Range.h"
#include "Exception.h"
#define DEBUG
#include "MemoryDetect.h"
#define new DEBUG_NEW
#define is_same(P,Q) std::is_same<P,Q>::value
namespace usr
{

template <class T>
class Matrix {
    private:
        uint32 col;
        uint32 row;
        uint32 mod;
        std::vector<T**> data;
        inline void allocate(const uint32 _row, const uint32 _col, const uint32 channel, const uint32 _mod);
        inline void release();
        Matrix Addition(const Matrix& src1, const Matrix& src2) const;
        Matrix Subtraction(const Matrix& src1, const Matrix& src2) const;
        Matrix Multiplication(const Matrix& src1, const Matrix& src2) const;
        Matrix Multiplication(const Matrix& src1, const T& src2) const;
        Matrix Strassen(const Matrix& src1, const Matrix& src2) const;
        Matrix Division(const Matrix& src1, const Matrix& src2) const;
        Matrix Division(const Matrix& src1, const T& src2) const;

    public:
        /** @brief constructor
         */
        Matrix(uint32 row=1,uint32 col=1,uint32 channel=1, uint32 mod=0);
        /** 
         * @description: deep copy 
         */
        Matrix(const Matrix&);
        /** @brief destructor
         *  release if is not empty
         */
        ~Matrix();

        inline uint32 getRow()const {return row;}
        inline uint32 getCol()const {return col;}
        inline uint32 getChannel()const {return data.size();}
        inline void setMod(uint32 m) {mod = m;}
        
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
        void clone(const Matrix& m)const;

        Matrix operator+(const Matrix& b);
        Matrix operator-(const Matrix& b);
        Matrix operator*(const Matrix& b);
        Matrix operator*(const T& b);
        template <class _T>
        friend Matrix<_T> operator*(const _T& a, const Matrix<_T>& b);    
        Matrix operator/(const Matrix& b);
        Matrix operator/(const T& b);
        Matrix operator^(int b);
        void operator=(const Matrix& m);
		void operator+=(const Matrix& m);
		void operator-=(const Matrix& m);
		void operator*=(const Matrix& m);
        void operator*=(const T t);
        template <typename _T>
        operator Matrix<_T>();
        bool operator==(const Matrix& m)const;
        T** operator[](uint32 i){ return data.at(i); }
        template <class _T>
        friend std::istream& operator>>(std::istream&, Matrix<_T>&);
        template <class _T>
        friend std::istream& operator>>(std::istream&, Matrix<std::complex<_T>>&);
        template <class _T>
		friend std::ostream& operator<<(std::ostream&, const Matrix<_T>&);

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
        T dotProduct(const Matrix&)const;

        /** @brief Computes a cross-product of two 3-element vectors.
        The method computes a cross-product of two 3-element vectors. The vectors must be
        3-element vectors of the same shape and size. The result is another 3-element vector
        of the same shape and type as operands.
        */
        Matrix crossProduct(const Matrix&)const;

        std::vector<std::complex<double>> eigenValue()const;
        Matrix<std::complex<double>> eigenVector(std::complex<double>)const;
        Matrix subMatrix(const Range& row=Range::all(), const Range& col=Range::all())const;
        Matrix inverse()const;
        Matrix transpose()const;
        Matrix conjugate()const;
        Matrix toDiagnal()const;
        Matrix reshape(const uint32 n, const uint32 _row)const;
        Matrix slice(const Range& row, const uint32 col)const;
        Matrix convolute(const Matrix& kernel,uint32 anchor_x=-1,uint32 anchor_y=-1)const;

};

template <class T>
inline void Matrix<T>::allocate(const uint32 _row, const uint32 _col, const uint32 channel, const uint32 _mod)
{
    row = _row, col = _col, mod = _mod;
    for (int i = 0; i < channel; i++)
    {
        T** tem = new T* [row];
        for (int i = 0; i < row; i++)
            tem[i] = new T[col];
        data.push_back(tem);
    }
}

template <class T>
Matrix<T>::Matrix(uint32 _row, uint32 _col, uint32 _channel, uint32 _mod)
{
    if (_row==0 || _col==0 || _channel==0)
        throw(InvalidArgsException("row or col or channel can not be zero", HERE));
    allocate(_row, _col, _channel, _mod);
}

template <class T>
Matrix<T>::Matrix(const Matrix<T>& m)
{
    allocate(m.row, m.col, m.getChannel(), m.mod);
    clone(m);
}

template <class T>
Matrix<T>::~Matrix()
{
    release();
}

template <class T>
inline void Matrix<T>::release()
{
    if (isEmpty())
        return;
    for (int i = 0; i < getChannel(); i++)
    {
        for (int j = 0; j < row; j++)
            delete[] data[i][j];
        delete[] data[i];
    }
    data.clear();
    row = 0, col = 0, mod = 0;
}

template <class T>
void Matrix<T>::clone(const Matrix<T>& m) const
{
    for (int k = 0; k < m.getChannel(); k++)
        for (int i = 0; i < row; i++)
            for (int j = 0; j < col; j++)
                data[k][i][j] = const_cast<Matrix<T>&>(m)[k][i][j];
}

template <class T>
Matrix<T> Matrix<T>::Addition(const Matrix<T>& src1, const Matrix<T>& src2) const
{
    try
    {
        if (src1.isEmpty() || src2.isEmpty())
            throw(EmptyMatrixException("Matrix is empty", HERE)); 
        else if (src1.getChannel()!=src2.getChannel())
            throw(ChannelMismatchException("Both of the Matrices must have the same Channels", HERE));
        else if (src1.getRow()!=src2.getRow() || src1.getCol()!= src2.getCol())
            throw(SizeMismatchException("Matrices do not match in size", HERE)); 
        else
        {
            uint32 r = src1.getRow(), c = src1.getCol(), channel = src1.getChannel();
            Matrix<T> rst(r, c, channel);
            for (int k = 0; k < channel; k++)
                for (int i = 0; i < r; i++)
                    for (int j = 0; j < c; j++)
                        rst[k][i][j] = const_cast<Matrix<T>&>(src1)[k][i][j] + const_cast<Matrix<T>&>(src2)[k][i][j];
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
    return Matrix<T>();
}

template <class T>
Matrix<T> Matrix<T>::Subtraction(const Matrix<T>& src1, const Matrix<T>& src2) const
{
    try
    {
        if (src1.isEmpty() || src2.isEmpty())
            throw(EmptyMatrixException("Matrix is empty", HERE)); 
        else if (src1.getChannel()!=src2.getChannel())
            throw(ChannelMismatchException("Both of the Matrices must have the same Channels", HERE));
        else if (src1.getRow()!=src2.getRow() || src1.getCol()!= src2.getCol())
            throw(SizeMismatchException("Matrices do not match in size", HERE)); 
        else
        {
            uint32 r = src1.getRow(), c = src1.getCol(), channel = src1.getChannel();
            Matrix<T> rst(r, c, channel);
            for (int k = 0; k < channel; k++)
                for (int i = 0; i < r; i++)
                    for (int j = 0; j < c; j++)
                        rst[k][i][j] = const_cast<Matrix<T>&>(src1)[k][i][j] - const_cast<Matrix<T>&>(src2)[k][i][j];
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
    return Matrix<T>();
}

template <class T>
Matrix<T> Matrix<T>::Multiplication(const Matrix<T>& src1, const Matrix<T>& src2) const
{
    try
    {
        if (src1.isEmpty() || src2.isEmpty())
            throw(EmptyMatrixException("Matrix is empty", HERE)); 
        else if (src1.getChannel()!=src2.getChannel())
            throw(ChannelMismatchException("Both of the Matrices must have the same Channels", HERE));
        else if (src1.getChannel()>1)
            throw(MultiChannelException("Channels greater than one is not supported", HERE));
        else if (src1.getCol() != src2.getRow())
            throw(SizeMismatchException("Matrices do not match in size", HERE));
        else
        {
            uint32 r = src1.getRow(), c = src2.getCol(), m = src1.getCol();
            if (src1.isSquare() && src2.isSquare() && r&1==0 && r >= 256)
                return Strassen(src1, src2);
            else 
            {
                Matrix<T> rst(r, c);
                for(int k=0;k<m;k++)
                    for(int i=0;i<r;i++)
                    {
                        T t = const_cast<Matrix<T>&>(src1)[0][i][k];
                        for(int j=0;j<c;j++)
                        {
                            #if mod>0
                                rst[0][i][j]=(rst[0][i][j]+1ll*t*const_cast<Matrix<T>&>(src2)[0][k][j])%mod;
                            #else
                                rst[0][i][j]+=t*const_cast<Matrix<T>&>(src2)[0][k][j];
                            #endif
                        }
                    }
                return rst;
            }
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
    return Matrix<T>();
}

template <class T>
Matrix<T> Matrix<T>::Multiplication(const Matrix<T>& src1, const T& src2) const
{
    try
    {
        if (src1.isEmpty())
            throw(EmptyMatrixException("Matrix is empty", HERE)); 
        else
        {
            uint32 r = src1.getRow(), c = src1.getCol(), channel = src1.getChannel();
            Matrix<T> rst(r, c, channel);
            for (int k = 0; k < channel; k++)
                for (int i = 0; i < r; i++)
                    for (int j = 0; j < c; j++)
                        rst[k][i][j] = const_cast<Matrix<T>&>(src1)[k][i][j] * src2;
            return rst;
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
    return Matrix<T>();
}

template <class T>
Matrix<T> Matrix<T>::Strassen(const Matrix<T>& src1, const Matrix<T>& src2) const
{
    uint32 size = src1.getRow()/2;
    Matrix<T> matrices[21] = {
        Matrix<T>(size, size), Matrix<T>(size, size), Matrix<T>(size, size), 
        Matrix<T>(size, size), Matrix<T>(size, size), Matrix<T>(size, size), 
        Matrix<T>(size, size), Matrix<T>(size, size), Matrix<T>(size, size),   
        Matrix<T>(size, size), Matrix<T>(size, size), Matrix<T>(size, size), 
        Matrix<T>(size, size), Matrix<T>(size, size), Matrix<T>(size, size), 
        Matrix<T>(size, size), Matrix<T>(size, size), Matrix<T>(size, size), 
        Matrix<T>(size, size), Matrix<T>(size, size), Matrix<T>(size, size),                                               
    };

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            matrices[0][0][i][j] = const_cast<Matrix<T>&>(src1)[0][i][j];
            matrices[1][0][i][j] = const_cast<Matrix<T>&>(src1)[0][i][j+size];
            matrices[2][0][i][j] = const_cast<Matrix<T>&>(src1)[0][i+size][j];
            matrices[3][0][i][j] = const_cast<Matrix<T>&>(src1)[0][i+size][j+size];
            matrices[4][0][i][j] = const_cast<Matrix<T>&>(src2)[0][i][j];
            matrices[5][0][i][j] = const_cast<Matrix<T>&>(src2)[0][i][j+size]; 
            matrices[6][0][i][j] = const_cast<Matrix<T>&>(src2)[0][i+size][j];
            matrices[7][0][i][j] = const_cast<Matrix<T>&>(src2)[0][i+size][j+size];                       
        }
    }

    matrices[19] = matrices[0] + matrices[3];
    matrices[20] = matrices[4] + matrices[7];
    matrices[12] = matrices[19] * matrices[20];

    matrices[19] = matrices[2] + matrices[3];
    matrices[13] = matrices[19] * matrices[4];

    matrices[20] = matrices[5] - matrices[7];
    matrices[14] = matrices[0] * matrices[20];

    matrices[20] = matrices[6] - matrices[4];
    matrices[15] = matrices[3] * matrices[20];
    
    matrices[19] = matrices[0] + matrices[1];
    matrices[16] = matrices[19] * matrices[7];

    matrices[20] = matrices[2] - matrices[0];
    matrices[19] = matrices[4] + matrices[5];
    matrices[17] = matrices[20] * matrices[19];
    
    matrices[20] = matrices[1] - matrices[3];
    matrices[19] = matrices[6] + matrices[7];
    matrices[18] = matrices[20] * matrices[19];
    
    matrices[20] = matrices[18] - matrices[16];
    matrices[19] = matrices[12] + matrices[15];
    matrices[8] = matrices[19] + matrices[20];
    
    matrices[9] = matrices[14] + matrices[16];
 
    matrices[10] = matrices[13] + matrices[15];

    matrices[20] = matrices[12] - matrices[13];
    matrices[19] = matrices[14] + matrices[17];
    matrices[17] = matrices[20] + matrices[19];

    Matrix<T> ret(2*size, 2*size);
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            ret[0][i][j] = matrices[8][0][i][j];
            ret[0][i][j+size] = matrices[9][0][i][j];
            ret[0][i+size][j] = matrices[10][0][i][j];
            ret[0][i+size][j+size] = matrices[11][0][i][j];                                    
        }
    }
    return ret;
}

template <class T>
Matrix<T> Matrix<T>::Division(const Matrix<T>& src1, const Matrix<T>& src2) const
{
    try
    {
        if (src1.isEmpty() || src2.isEmpty())
            throw(EmptyMatrixException("Matrix is empty", HERE)); 
        else if (src1.getChannel()!=src2.getChannel())
            throw(ChannelMismatchException("Both of the Matrices must have the same Channels", HERE));
        else if (src1.getCol()!=src1.getCol() || src1.getRow()!=src2.getRow())
            throw(SizeMismatchException("Matrices do not match in size", HERE));
        else
        {
            uint32 r = src1.getRow(), c = src1.getCol(), channel = src1.getChannel();
            Matrix<T> rst(r, c, channel);
            for (int k = 0; k < channel; k++)
                for (int i = 0; i < row; i++)
                    for (int j = 0; j < col; j++)
                        rst[k][i][j] = const_cast<Matrix<T>&>(src1)[k][i][j] / const_cast<Matrix<T>&>(src2)[k][i][j];
            return rst;
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
    return Matrix<T>();
}

template <class T>
Matrix<T> Matrix<T>::Division(const Matrix<T>& src1, const T& src2) const
{
    try
    {
        if (src1.isEmpty())
            throw(EmptyMatrixException("Matrix is empty", HERE)); 
        else if (src2 == 0)
            throw(ArithmeticException("Divide by Zero", HERE));
        else
        {
            uint32 r = src1.getRow(), c = src1.getCol(), channel = src1.getChannel();
            Matrix<T> rst(r, c, channel);
            for (int k = 0; k < channel; k++)
                for (int i = 0; i < r; i++)
                    for (int j = 0; j < c; j++)
                        rst[k][i][j] = const_cast<Matrix<T>&>(src1)[k][i][j] / src2;
            return rst;
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
    return Matrix<T>();
}

template <class T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& m)
{
    return Addition(*this, m);
}

template <class T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& m)
{
    return Subtraction(*this, m);
}

template <class T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& m)
{
    return Multiplication(*this, m);
}

template <class T>
Matrix<T> Matrix<T>::operator*(const T& a)
{
    return Multiplication(*this, a);
}

template <class _T>
Matrix<_T> operator*(const _T& a, const Matrix<_T>& b)
{
    return b*a;
}    

template <class T>
Matrix<T> Matrix<T>::operator/(const Matrix<T>& m)
{
    return Division(*this, m);
}

template <class T>
Matrix<T> Matrix<T>::operator/(const T& a)
{
    return Division(*this, a);
}

template <class T>
Matrix<T> Matrix<T>::operator^(int b)
{
    try
    {
        if (isEmpty())
            throw(EmptyMatrixException("Matrix is empty, can't find determinant", HERE)); 
        else if (!isSquare())
            throw(MatrixNotSquareException("Matrix must be square", HERE));
        else if (getChannel()>2)
            throw(MultiChannelException("Channels greater than two is not supported", HERE));
        else
        {
            Matrix<T> ret(row,col);
            for(int i=0;i<row;i++)
                for(int j=0;j<col;j++) ret[0][i][j]=(i==j);
            Matrix<T> cur=*this;
            for(;b;b>>=1,cur*=cur)
                if(b&1) ret*=cur;
            return ret;
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
    return Matrix<T>();
}

template <class T>
void Matrix<T>::operator=(const Matrix<T>& m)
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
void Matrix<T>::operator+=(const Matrix<T>& m)
{
    (*this) = (*this)+m;
}

template <class T>
void Matrix<T>::operator-=(const Matrix<T>& m)
{
    (*this) = (*this)-m;
}

template <class T>
void Matrix<T>::operator*=(const Matrix<T>& m)
{
    (*this) = (*this)*m;
}

template <class T>
void Matrix<T>::operator*=(const T t)
{
    (*this) = (*this)*t;
}

template <class T>
bool Matrix<T>::operator==(const Matrix<T>& m)const
{
    if(row!=m.getRow()||col!=m.getCol()||getChannel()!=m.getChannel()) return 0;
    for(int k=0;k<getChannel();k++)
        for(int i=0;i<m.row;i++)
            for(int j=0;j<m.col;j++)
                if(data[k][i][j]!=const_cast<Matrix<T>&>(m)[k][i][j]) return 0;
    return 1;
}

template <class T>
template <typename _T>
Matrix<T>::operator Matrix<_T>()
{
    Matrix<_T> t(row, col, getChannel());
    for (int i = 0; i < getChannel(); i++)
        for (int j = 0; j < row; j++)
            for (int k = 0; k < col; k++)
                t[i][j][k] = static_cast<_T>(data[i][j][k]);
    return t;
}

template <class _T>
std::istream& operator>>(std::istream& is, Matrix<_T>& m)
{
    for (int k = 0; k < m.getChannel(); k++)
        for (int i = 0; i < m.row; i++)
            for (int j = 0; j < m.col; j++)
                is >> m[k][i][j];
    return is;
}

template <class _T>
std::istream& operator>>(std::istream& is, Matrix<std::complex<_T>>& m)
{
    for (int k = 0; k < m.getChannel(); k++)
        for (int i = 0; i < m.row; i++)
            for (int j = 0; j < m.col; j++)
            {
                _T rea, ima;
                is >> rea >> ima;
                m[k][i][j].real(rea);
                m[k][i][j].imag(ima);
            }
    return is;
}

template <class _T>
std::ostream& operator<<(std::ostream& os, const Matrix<_T>& m)
{
    os << "[";
    for (int i = 0; i < m.row; i++)
    {
        if(i>0) os<<" ";
        for (int j = 0; j < m.col; j++)
            for (int k = 0; k < m.getChannel(); k++)
                if (i == m.row-1 && j == m.col-1 && k == m.getChannel()-1)
                    os << const_cast<Matrix<_T>&>(m)[k][i][j] << "]";
                else if(j == m.col-1 && k == m.getChannel()-1) os<<const_cast<Matrix<_T>&>(m)[k][i][j]<<";";
                else 
                    os << const_cast<Matrix<_T>&>(m)[k][i][j] << ", ";
        os << "\n";
    }
    return os;
}

uint32 ksm(int a,uint32 x,uint32 mod)
{
    uint32 ret = 1;
    for(;x;x>>=1,a=1ll*a*a%mod) if(x&1) ret = 1ll*ret*a%mod;
    return ret;
}

template<class T>
Matrix<T> Matrix<T>::toDiagnal()const
{
    #if (mod==0)
    {
        uint32 r = row, c = col;
        Matrix<T> rst = *this;
        for (int i = 0; i < r-1; i++)
        {
            int id=-1;
            for (int j = i; j < r; j++)
                if (rst[0][j][i])
                {
                    id=j;
                    break;
                }
            if (id!=i)
            {
                for (int j = 0; j < c; j++)
                {
                    std::swap(rst[0][id][j],rst[0][i][j]);
                }
            }
            for (int j = i + 1; j < r; j++)
            {
                for (int k = c - 1; k >= i; k--)
                {
                    rst[0][j][k] -= rst[0][i][k] * rst[0][j][i] / rst[0][i][i];
                }
            }
        }
        for(int i=r-1;i>=0;i--)
        {
            if(rst[0][i][i]==0) continue;
            for(int j=i-1;j>=0;j--)
            {
                for (int k = c - 1; k >= i; k--)
                {
                    rst[0][j][k] -= rst[0][i][k] * rst[0][j][i] / rst[0][i][i];
                }
            }
        }        
        // for(int i=r-1;i>=0;i--)
        // {
        //     for(int j=r-1;j>i;j--)
        //     {
        //         if(rst[0][j][j]==0) continue;
        //         T t = rst[0][i][j]/rst[0][j][j];
        //         for(int k=i+1;k<col;k++) rst[0][i][k]-=t*rst[0][j][j];
        //     }
        // }
        return rst;
    }
    #else 
    {
        uint32 r = row, c = col, mo = mod;
        Matrix<T> rst(r, c, 1, mo);
        for (int i = 0; i < r; i++)
        {
            int id=-1;
            for (int j = i; j < r; j++)
                if (rst[0][j][i])
                {
                    id=j;
                    break;
                }
            if (id!=i)
            {
                for (int j = 0; j < c; j++)
                {
                    std::swap(rst[0][id][j],rst[0][i][j]);
                }
            }
            uint32 p=1,k=mo-2,a=rst[0][i][i];
            for (;k;k>>=1){
                if (k&1)p=1ull* p * a %mo;
                a=1ull * a *a % mo;
            }
            for (int j= i + 1; j < c; j++){
                rst[0][i][j] = 1ull * rst[0][i][j] * p %mo;   
            }
            for (int j = i + 1; j < r; j++)
            {
                for (int k = c - 1; k > i; k--)
                {
                    rst[0][j][k] = (rst[0][j][k] - rst[0][i][k] * rst[0][j][i]) % mo;
                }
            }
        }
        for(int i=r-1;i>=0;i--)
        {
            for(int j=r-1;j>i;j--)
            {
                if(rst[0][j][j]==0) continue;
                T t = rst[0][i][j]*ksm(rst[0][j][j],mod-2,mod);
                for(int k=i+1;k<col;k++) rst[0][i][k]=((rst[0][i][k]-1ll*t*rst[0][j][j])%mod+mod)%mod;
            }
        }
        return rst;
    }
    #endif
}
template<class T>
T Matrix<T>::determinant()const
{
    try
    {
        if (isEmpty())
            throw(EmptyMatrixException("Matrix is empty, can't find determinant", HERE)); 
        else if (getChannel()>1)
            throw(MultiChannelException("Channel greater than one is not supported", HERE));
        else if (!isSquare())
            throw(MatrixNotSquareException("Matrix must be square", HERE));
        else 
        {   
            #if mod==0
            {
                uint32 r = row, c = col;
                Matrix<T> rst(r,c);
                rst = (*this);
                T res=1,res_=1;
                for (int i = 0; i < r; i++)
                {
                    int id=-1;
                    for (int j = i; j < r; j++)
                        if (rst[0][j][i]!=0.0)
                        {
                            id=j;
                            break;
                        }
                    if (id==-1)
                    {
                        res=0;
                        return res;
                    }
                    if (id!=i)
                    {
                        res = -res;
                        for (int j = 0; j < c; j++)
                        {
                            std::swap(rst[0][id][j],rst[0][i][j]);
                        }
                    }
                    res = res * rst[0][i][i];
                    for (int j = i + 1; j < r; j++)
                    {
                        res_*=rst[0][i][i];
                        for (int k = c - 1; k >= i; k--){
                            rst[0][j][k] *= rst[0][i][i];
                        }
                        for (int k = c - 1; k >= i; k--)
                        {
                            rst[0][j][k] -= rst[0][i][k] * rst[0][j][i] / rst[0][i][i];
                        }
                    }
                }
                //std::cerr<<"-------------------"<<cerr<<std::endl<<rst<<std::endl;
                return res/res_;
            }
            #else 
            {
                uint32 r = row, c = col, mo = mod;
                Matrix<T> rst(r, c, mo);
                uint32 res = 1;
                for (int i = 0; i < r; i++)
                {
                    int id=-1;
                    for (int j = i; j < r; j++)
                        if (rst[0][j][i])
                        {
                            id=j;
                            break;
                        }
                    if (id==-1)
                    {
                        res=0;
                        return res;
                    }
                    if (id!=i)
                    {
                        res = -res;
                        for (int j = 0; j < c; j++)
                        {
                            std::swap(rst[0][id][j],rst[0][i][j]);
                        }
                    }
                    res = 1ull* res * rst[0][i][i] % mo;
                    uint32 p=1,k=mo-2,a=rst[0][i][i];
                    for (;k;k>>=1){
                        if (k&1)p=1ull* p * a %mo;
                        a=1ull * a *a % mo;
                    }
                    for (int j= i + 1; j < c; j++){
                        rst[0][i][j] = 1ull * rst[0][i][j] * p %mo;   
                    }

                    for (int j = i + 1; j < r; j++)
                    {
                        for (int k = c - 1; k >= i; k--)
                        {
                            rst[0][j][k] = (rst[0][j][k] - rst[0][i][k] * rst[0][j][i]) % mo;
                        }
                    }
                }
                res %= mo;
                if (res<0)res += mo;
                return res;
            }
            #endif
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
uint32 Matrix<T>::rank()const{
    try 
    {
        if (isEmpty())
            throw(EmptyMatrixException("Matrix is empty, can't find rank", HERE)); 
        else if (getChannel()>1)
            throw(MultiChannelException("Channel greater than one is not supported", HERE));
        uint32 res=0;
        uint32 r = row, c = col;
        Matrix<T> rst(r, c);
        rst = (*this);
        for (int i = 0; i < r && i < c; i++)
        {
            int id=-1;
            for (int j = i; j < r; j++)
                if (rst[0][j][i])
                {
                    id=j;
                    break;
                }
            if (id==-1)
            {
                continue;
            }
            res++;
            if (id!=i)
            {
                for (int j = 0; j < c; j++)
                {
                    std::swap(rst[0][id][j],rst[0][i][j]);
                }
            }
            for (int j = i + 1; j < r; j++)
            {
                for (int k = c - 1; k >= i; k--)
                    rst[0][j][k] = rst[0][j][k] * rst[0][i][i];
                for (int k = c - 1; k >= i; k--)
                {
                    rst[0][j][k] -= rst[0][i][k] * rst[0][j][i] / rst[0][i][i];
                }
            }
        }
        return res;
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
T Matrix<T>::trace()const
{
    try
    {
        if (isEmpty())
            throw(EmptyMatrixException("Matrix is empty, can't find trace", HERE)); 
        else if (getChannel()>1)
            throw(MultiChannelException("Channel greater than one is not supported", HERE));
        else if (!isSquare())
            throw(MatrixNotSquareException("Matrix must be square", HERE));
        else 
        {   
            T res=0;
            for (int i=0;i<row;i++){
                res=res+data[0][i][i];
            }
            return res;
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
T Matrix<T>::max(const Range &row, const Range &col)const
{
    try
    {
        if (row.empty() || col.empty())
            throw(InvalidArgsException("range can not be zero", HERE));
        else if (row.start >= getRow() || col.start >= getCol())
            throw(RangeOutOfBoundException("range is out of bound", HERE));
        else
        {
            T m = data[0][row.start][col.start];
            for (int k = 0; k < getChannel(); k++)
                for (int i = row.start; i < std::min(getRow(), row.end); i++)
                    for (int j = col.start; j < std::min(getCol(), col.end); j++)
                        if (m < data[k][i][j])
                            m = data[k][i][j];
            return m;
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
T Matrix<T>::min(const Range &row, const Range &col)const
{
    try
    {
        if (row.empty() || col.empty())
            throw(InvalidArgsException("range can not be zero", HERE));
        else if (row.start >= getRow() || col.start >= getCol())
            throw(RangeOutOfBoundException("range is out of bound", HERE));
        else
        {
            T m = data[0][row.start][col.start];
            for (int k = 0; k < getChannel(); k++)
                for (int i = row.start; i < std::min(getRow(), row.end); i++)
                    for (int j = col.start; j < std::min(getCol(), col.end); j++)
                        if (m > data[k][i][j])
                            m = data[k][i][j];
            return m;
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
T Matrix<T>::avg(const Range &row, const Range &col)const
{
    try
    {
        if (row.empty() || col.empty())
            throw(InvalidArgsException("range can not be zero", HERE));
        else if (row.start >= getRow() || col.start >= getCol())
            throw(RangeOutOfBoundException("range is out of bound", HERE));
        else
        {
            T tot = 0;
            uint32 cnt = 0;
            for (int k = 0; k < getChannel(); k++)
                for (int i = row.start; i < std::min(getRow(), row.end); i++)
                    for (int j = col.start; j < std::min(getCol(), col.end); j++)
                        tot += data[k][i][j],cnt++;                 
            return tot/cnt;
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
T Matrix<T>::sum(const Range &row, const Range &col)const
{
    try
    {
        if (row.empty() || col.empty())
            throw(InvalidArgsException("range can not be zero", HERE));
        else if (row.start >= getRow() || col.start >= getCol())
            throw(RangeOutOfBoundException("range is out of bound", HERE));
        else
        {
            T tot = 0;
            for (int k = 0; k < getChannel(); k++)
                for (int i = row.start; i < std::min(getRow(), row.end); i++)
                    for (int j = col.start; j < std::min(getCol(), col.end); j++)
                        tot += data[0][i][j];
            return tot;
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
T Matrix<T>::dotProduct(const Matrix<T>& m)const
{
    try{
        if(row*col!=m.getRow()*m.getCol())
            throw(SizeMismatchException("the total size must be equal", HERE));
        else if (getChannel() != m.getChannel())
            throw(ChannelMismatchException("Both of the Matrices must have the same Channels", HERE));
        else
        {
            T tot = 0;
            for (int k=0;k<getChannel();k++)
                for(int i=0;i<row*col;i++)
                {
                    int rowx = i/col, colx = i%col, rowy = i/m.getCol(), coly = i%m.getCol();
                    tot += data[k][rowx][colx]*const_cast<Matrix<T>&>(m)[k][rowy][coly];
                }
            return tot;
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
Matrix<T> Matrix<T>::crossProduct(const Matrix<T>& m)const
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
            Matrix<T> ret(1,3);
            ret[0][0][0] = data[0][0][1]*const_cast<Matrix<T>&>(m)[0][0][2]-data[0][0][2]*const_cast<Matrix<T>&>(m)[0][0][1];
            ret[0][0][1] = data[0][0][2]*const_cast<Matrix<T>&>(m)[0][0][0]-data[0][0][0]*const_cast<Matrix<T>&>(m)[0][0][2];
            ret[0][0][2] = data[0][0][0]*const_cast<Matrix<T>&>(m)[0][0][1]-data[0][0][1]*const_cast<Matrix<T>&>(m)[0][0][0];
            return ret;
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
    return Matrix<T>();
}
template <class T>
Matrix<T> Matrix<T>::transpose()const
{   
    try{
        if (isEmpty())
            throw(EmptyMatrixException("Matrix is empty, can't transpose", HERE)); 
        else
        {
            Matrix<T> ret(col,row,getChannel(),mod);
            for (int k=0;k<getChannel();k++)
                for(int i=0;i<col;i++)
                    for(int j=0;j<row;j++) 
                        ret[k][i][j]=data[k][j][i];
            return ret;
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
    return Matrix<T>();
}
template <class T>
Matrix<T> Matrix<T>::conjugate()const
{
    try{
        if (isEmpty())
            throw(EmptyMatrixException("Matrix is empty, can't conjugate", HERE)); 
        else
        {
            Matrix<T> ret(row,col,getChannel(),mod);
            for (int k=0;k<getChannel();k++)
                for(int i=0;i<row;i++)
                    for(int j=0;j<col;j++) 
                        if constexpr 
                        (is_same(T, std::complex<double>) || 
                        is_same(T, std::complex<int>) || 
                        is_same(T, std::complex<float>))
                            ret[k][i][j]=std::conj(data[k][i][j]);
                        else 
                            ret[k][i][j]=std::conj(data[k][i][j]).real();
            return ret;
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
    return Matrix<T>();
}
template <class T>
Matrix<T> Matrix<T>::inverse()const
{
    try{
        if (isEmpty())
            throw(EmptyMatrixException("Matrix is empty, can't find inverse", HERE)); 
        else if (!isSquare())
            throw(MatrixNotSquareException("non-square matrix is not invertible", HERE));
        else if (!isInvertible())
            throw(InverseNotExistException("singular matrix is not invertible", HERE));
        else
        {
            Matrix<T> t(row,row+col,1,mod);
            for(int i=0;i<row;i++)
                for(int j=0;j<row+col;j++)
                {
                    if(j<col) t[0][i][j]=data[0][i][j];
                    else t[0][i][j]=(j-col==i); 
                }
            Matrix<T> m = t.toDiagnal();
            Matrix<T> ret(row,col,1,mod);
            for(int i=0;i<row;i++)
                for(int j=0;j<col;j++)
                {
                    #if mod>0
                        ret[0][i][j] = 1ll*m[0][i][row+j]*ksm(m[0][i][i],mod-2)%mod;
                    #else 
                        ret[0][i][j] = m[0][i][row+j]/m[0][i][i];
                    #endif
                }
            return ret;
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
    return Matrix<T>();
}
template <class T>
Matrix<T> Matrix<T>::subMatrix(const Range& row, const Range& col)const
{
    try
    {
        if (row.empty() || col.empty())
            throw(InvalidArgsException("range can not be zero", HERE));
        else if (row.start >= getRow() || col.start >= getCol())
            throw(RangeOutOfBoundException("range is out of bound", HERE));
        else
        {
            Matrix<T> ret(std::min(getRow(), row.end)-row.start,std::min(getCol(), col.end)-col.start,getChannel(),mod);
            for (int k=0;k<getChannel();k++)
                for(int i=0;i<ret.getRow();i++)
                    for(int j=0;j<ret.getCol();j++)
                        ret[k][i][j]=data[k][i+row.start][j+col.start];
            return ret;
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
    return Matrix<T>();
}
template <class T>
Matrix<T> Matrix<T>::reshape(const uint32 n, const uint32 _row)const
{
    try{
        if(n != 0 && 1ll*row*col*getChannel() % 1ll*_row*n != 0 || row*col % _row != 0)
            throw(SizeMismatchException("cannot fit into this shape", HERE));
        else
        {
            uint32 channel = (n!=0)?n:getChannel();
            Matrix<T> ret(_row, 1ll*row*col*getChannel()/1ll*_row*channel, channel,mod);
            T element[row*col*getChannel()];
            int t = 0;
            for (int i = 0; i < row; i++)
                for (int j = 0; j < col; j++)
                    for (int k = 0; k < getChannel(); k++)
                        element[t++] = data[k][i][j];
            t = 0;
            for (int i = 0; i < ret.getRow(); i++)
                for (int j = 0; j < ret.getCol(); j++)
                    for (int k = 0; k < ret.getChannel(); k++)
                        ret[k][i][j] = element[t++];
            return ret;
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
    return Matrix<T>();
}
template <class T>
Matrix<T> Matrix<T>::slice(const Range& row, const uint32 col)const
{
    try{
        if (row.empty())
            throw(InvalidArgsException("range can not be zero", HERE));
        if(col>=this->col||row.start>=this->row)
            throw(RangeOutOfBoundException("range out of bound", HERE));
        else
        {
            Matrix<T> ret(std::min(this->getRow(), row.end) - row.start,1,getChannel(),mod);
            for (int k=0;k<getChannel();k++)
                for(int i=row.start;i< ret.getRow();i++)
                    ret[k][i-row.start][0] = data[k][i][col];
            return ret;
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
    return Matrix<T>();
}

template <class T>
Matrix<T> Matrix<T>::convolute(const Matrix& kernel,uint32 anchor_x,uint32 anchor_y)const
{
    if (anchor_x==-1)
        anchor_x = kernel.getRow()/2;
    if (anchor_y==-1)
        anchor_y = kernel.getCol()/2;
    try{
        if (isEmpty()||kernel.isEmpty())
            throw(EmptyMatrixException("Matrix is empty", HERE)); 
        else if (getChannel()!=kernel.getChannel())
            throw(ChannelMismatchException("Both of the Matrices must have the same Channels", HERE));
        else if(anchor_x>=kernel.getRow()||anchor_y>=kernel.getCol())
            throw(RangeOutOfBoundException("range out of bound", HERE));
        else{
            //cout<<"anchor_x = "<<anchor_x<<" anchor_y = "<<anchor_y<<" row = "<<row<<" col = "<<col<<endl;
            Matrix<T> ret(row,col);
            for(int i=0;i<row;i++)
                for(int j=0;j<col;j++) ret[0][i][j] = 0;
            int R = kernel.getRow(),C = kernel.getCol();
            for (int m=0;m<getChannel();m++)
                for(int i=0;i<row;i++)
                    for(int j=0;j<col;j++)
                    {
                        for(int k=0;k<R;k++)
                        {
                            for(int l=0;l<C;l++)
                            {
                                int x = i-anchor_x+k,y=j-anchor_y+l;
                                if(x<0) x=-x;
                                else if(x>=row) x = 2*(row-1)-x;
                                if(y<0) y=-y;
                                else if(y>=col) y = 2*(col-1)-y;
                                //int x = i-anchor_x+k,y = j-anchor_y+l;
                                #if mod>0
                                    ret[0][i][j]=(ret[0][i][j]+1ll*const_cast<Matrix<T>&>(kernel)[m][k][l]*data[m][x][y])%mod;
                                #else
                                    ret[0][i][j]+=const_cast<Matrix<T>&>(kernel)[m][k][l]*data[m][x][y];
                                #endif
                            }
                        }
                    }
            return ret;
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
    return Matrix<T>();
}

void QR_Decomposition(Matrix<std::complex<double> >& A,Matrix<std::complex<double> >& Q,Matrix<std::complex<double> >& R){
    uint32 n = A.getRow();
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            Q[0][i][j]=R[0][i][j]=0;
        }    
    }
    for (int k = 0; k < n; k++){
		double MOD = 0;
		for (int i = 0; i < n; i++)
		{
			MOD += abs(A[0][i][k] * A[0][i][k]); 
		}
		R[0][k][k] = sqrt(MOD);
		for (int i = 0; i < n; i++)
		{
			Q[0][i][k] = A[0][i][k] / R[0][k][k]; 
		}

		for (int i = k + 1; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				R[0][k][i] += A[0][j][i] * Q[0][j][k];
			}
			for (int j = 0; j < n; j++)
			{
				A[0][j][i] -= R[0][k][i] * Q[0][j][k]; 
			}
		}
	}
}

template <class T>
std::vector<std::complex<double>> Matrix<T>::eigenValue() const
{
    try{
        if (isEmpty())
            throw(EmptyMatrixException("Matrix is empty, can't find eigenvalue", HERE)); 
        else if (!isSquare())
            throw(MatrixNotSquareException("Non-square matrix does not have eigenValue", HERE));
        else if (getChannel()>1)
            throw(MultiChannelException("Channels greater than one is not supported", HERE));
        else
        {
            int n = row;
            Matrix<std::complex<double> > A(n,n);
            for (int i=0;i<n;i++){
                for (int j=0;j<n;j++){
                    A[0][i][j]=data[0][i][j];
                }
            }
            Matrix<std::complex<double> > U(n,n);
            for (int i=0;i<n;i++){
                for (int j=0;j<n;j++){
                    if (i==j){
                        U[0][i][j]=1;
                    }
                    else{
                        U[0][i][j]=0;
                    }
                }
            }
            Matrix<std::complex<double> > Q(n,n);
            Matrix<std::complex<double> > R(n,n);
            Matrix<std::complex<double> > temp(n,n);
            for (int tim=0;tim<50;tim++){
                QR_Decomposition(A,Q,R);
                temp = R*Q;
                A = temp;
                temp = U*Q;
                U = temp;
            }
            std::vector<std::complex<double>> res;
            for (int i=0;i<n;i++){
                res.push_back(A[0][i][i]);
            }
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
Matrix<std::complex<double>> Matrix<T>::eigenVector(std::complex<double> val)const{
    
    try{
        int n=col;
        Matrix<std::complex<double>> a (this->row,this->col);
        for (int i=0;i<n;i++){
            for (int j=0;j<n;j++){
                a[0][i][j]=-this->data[0][i][j];
                if (i==j)a[0][i][j]=-this->data[0][i][j]+val;
            }
        }
        if (std::abs(a.determinant())>1e-9)
            throw(WrongEigenValueException("eigenvalue is wrong", HERE));
        for (int i=0;i<n;i++){
            int id=-1;
            for (int j=i;j<n;j++){
                if (a[0][j][i] != 0.0){
                    id=j;
                    break;
                }
            }
            if (id==-1){
                continue;
            }
            if (id!=i){
                for (int j=0;j<n;j++){
                    swap(a[0][i][j],a[0][id][j]);
                }
            }
            for (int j=i+1;j<n;j++){
                for (int k=n-1;k>=i;k--){
                    a[0][j][k]-=a[0][i][k]*a[0][j][i]/a[0][i][i];
                }
            }
        }
        Matrix<std::complex<double>> res(1,n);
        res[0][0][n-1]=1;
        for (int i=n-2;i>=0;i--){
            std::complex<double>tmp=0;
            for (int j=i+1;j<n;j++){
                tmp-=a[0][i][j]*res[0][0][j];
            }
            if (std::abs(a[0][i][i])>1e-9){
                res[0][0][i]=tmp/(a[0][i][i]);
            }
        }
        return res;
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


} // namespace usr


#endif