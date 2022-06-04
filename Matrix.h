/*
 * @Author: Satori 3102809947@qq.com
 * @Date: 2022-05-26 11:51:58
 * @LastEditors: Satori 3102809947@qq.com
 * @LastEditTime: 2022-06-03 17:48:32
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
// #define DEBUG
#include "MemoryDetect.h"
#define new DEBUG_NEW
namespace usr
{

template <class T>
class Matrix {
    private:
        uint32 col;
        uint32 row;
        uint32 mod;
        T** data;
        inline void allocate(const uint32 _row, const uint32 _col, const uint32 _mod);
        Matrix* Addition(const Matrix& src1, const Matrix& src2) const;
        Matrix* Subtraction(const Matrix& src1, const Matrix& src2) const;
        Matrix* Multiplication(const Matrix& src1, const Matrix& src2) const;
        Matrix* Multiplication(const Matrix& src1, const T& src2) const;
        Matrix* Strassen(const Matrix& src1, const Matrix& src2) const;
        Matrix* Division(const Matrix& src1, const Matrix& src2) const;
        Matrix* Division(const Matrix& src1, const T& src2) const;

    public:
        /**
         * @description: construct a Matrix
         * @param {int} col 
         * @param {int} row
         * @return {*}
         */
        Matrix(uint32 row=1,uint32 col=1,uint32 mod=0);
        /** 
         * @description: Shallow copy 
         */
        Matrix(const Matrix&);
        /** @brief destructor
         *  release if is not empty
         */
        ~Matrix();

        inline uint32 getRow()const {return row;}
        inline uint32 getCol()const {return col;}
        
        /** 
         * return true if the matrix is empty (data is equal to nullptr)
         */
        inline bool isEmpty()const {return data==nullptr;};

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
         * @description: Deep copy
         */
        void clone(const Matrix& m)const;

        /** @brief Free data
         *  Set data pointer as nullptr after deleting 
         */
        inline void release();

        Matrix* operator+(const Matrix& b);
        Matrix* operator-(const Matrix& b);
        Matrix* operator*(const Matrix& b);
        Matrix* operator*(const T& b);
        template <class _T>
        friend Matrix<_T>* operator*(const _T& a, const Matrix<_T>& b);    
        Matrix* operator/(const Matrix& b);
        Matrix* operator/(const T& b);
        Matrix* operator^(int b);
        void operator=(const Matrix& m);
		void operator+=(const Matrix& m);
		void operator-=(const Matrix& m);
		void operator*=(const Matrix& m);
        void operator*=(const T t);
        bool operator==(const Matrix& m)const;
        T* operator[](uint32 i){ return *(data+i); }
        template <class _T>
        friend std::istream& operator>>(std::istream&, Matrix<_T>&);
        template <class _T>
		friend std::ostream& operator<<(std::ostream&, Matrix<_T>&);

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
        Matrix* crossProduct(const Matrix&)const;

        std::vector<std::complex<double>> eigenValue()const;
        std::vector<Matrix*> eigenVector()const;
        Matrix<std::complex<double>>* eigenVector(std::complex<double>)const;
        Matrix* subMatrix(const Range& row=Range::all(), const Range& col=Range::all())const;
        Matrix* inverse()const;
        Matrix* transpose()const;
        Matrix* hermite()const;
        Matrix* toDiagnal()const;
        Matrix* reshape(const uint32 _row, const uint32 _col)const;
        Matrix* slice(const Range& row, const uint32 col)const;
        Matrix* convolute(const Matrix& kernal,uint32 anchor_x,uint32 anchor_y)const;

};
template <class T>
inline void Matrix<T>::allocate(const uint32 _row, const uint32 _col, const uint32 _mod)
{
    row = _row, col = _col, mod = _mod;
    data = new T* [row];
    for (int i = 0; i < row; i++)
        data[i] = new T[col];
}

template <class T>
Matrix<T>::Matrix(uint32 _row, uint32 _col, uint32 _mod)
{
    if (_row==0 || _col==0)
        throw(InvalidArgsException("row or col can not be zero", __func__, __FILE__, __LINE__));
    allocate(_row, _col, _mod);
}

template <class T>
Matrix<T>::Matrix(const Matrix<T>& m)
{
    allocate(m.row, m.col, m.mod);
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
    for (int i = 0; i < row; i++)
        delete[] data[i];
    delete[] data;
    data = nullptr;
    row = 0, col = 0, mod = 0;
}

template <class T>
void Matrix<T>::clone(const Matrix<T>& m) const
{
    for (int i = 0; i < row; i++)
        for (int j = 0; j < col; j++)
            data[i][j] = const_cast<Matrix<T>&>(m)[i][j];
}

template <class T>
Matrix<T>* Matrix<T>::Addition(const Matrix<T>& src1, const Matrix<T>& src2) const
{
    try
    {
        if (src1.isEmpty() || src2.isEmpty())
            throw(EmptyMatrixException("Matrix is empty", __func__, __FILE__, __LINE__)); 
        else if (src1.getRow()!=src2.getRow() || src1.getCol()!= src2.getCol())
            throw(SizeMismatchException("Matrices do not match in size", __func__, __FILE__, __LINE__)); 
        else
        {
            uint32 r = src1.getRow(), c = src1.getCol();
            Matrix<T>* rst = new Matrix<T>(r, c);
            for (int i = 0; i < r; i++)
                for (int j = 0; j < c; j++)
                    (*rst)[i][j] = const_cast<Matrix<T>&>(src1)[i][j] + const_cast<Matrix<T>&>(src2)[i][j];
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
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return nullptr;
}

template <class T>
Matrix<T>* Matrix<T>::Subtraction(const Matrix<T>& src1, const Matrix<T>& src2) const
{
    try
    {
        if (src1.isEmpty() || src2.isEmpty())
            throw(EmptyMatrixException("Matrix is empty", __func__, __FILE__, __LINE__)); 
        else if (src1.getRow()!=src2.getRow() || src1.getCol()!= src2.getCol())
            throw(SizeMismatchException("Matrices do not match in size", __func__, __FILE__, __LINE__)); 
        else
        {
            uint32 r = src1.getRow(), c = src1.getCol();
            Matrix<T>* rst = new Matrix<T>(r, c);
            for (int i = 0; i < r; i++)
                for (int j = 0; j < c; j++)
                    (*rst)[i][j] = const_cast<Matrix<T>&>(src1)[i][j] - const_cast<Matrix<T>&>(src2)[i][j];
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
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return nullptr;
}

template <class T>
Matrix<T>* Matrix<T>::Multiplication(const Matrix<T>& src1, const Matrix<T>& src2) const
{
    try
    {
        if (src1.isEmpty() || src2.isEmpty())
            throw(EmptyMatrixException("Matrix is empty", __func__, __FILE__, __LINE__)); 
        else if (src1.getCol() != src2.getRow())
            throw(SizeMismatchException("Matrices do not match in size", __func__, __FILE__, __LINE__));
        else
        {
            uint32 r = src1.getRow(), c = src2.getCol(), m = src1.getCol();
            Matrix<T>* rst = nullptr;
            if (src1.isSquare() && src2.isSquare() && r&1==0 && r >= 256)
                rst = Strassen(src1, src2);
            else 
            {
                rst = new Matrix<T>(r, c);
                for(int k=0;k<m;k++)
                    for(int i=0;i<r;i++)
                    {
                        T t = const_cast<Matrix<T>&>(src1)[i][k];
                        for(int j=0;j<c;j++)
                        {
                            #if mod>0
                                (*rst)[i][j]=((*rst)[i][j]+1ll*t*const_cast<Matrix<T>&>(src2)[k][j])%mod;
                            #else
                                (*rst)[i][j]+=t*const_cast<Matrix<T>&>(src2)[k][j];
                            #endif
                        }
                    }
            }
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
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return nullptr;
}

template <class T>
Matrix<T>* Matrix<T>::Multiplication(const Matrix<T>& src1, const T& src2) const
{
    try
    {
        if (src1.isEmpty())
            throw(EmptyMatrixException("Matrix is empty", __func__, __FILE__, __LINE__)); 
        else
        {
            uint32 r = src1.getRow(), c = src1.getCol();
            Matrix<T>* rst = new Matrix<T>(r, c);
            for (int i = 0; i < r; i++)
                for (int j = 0; j < c; j++)
                    (*rst)[i][j] = const_cast<Matrix<T>&>(src1)[i][j] * src2;
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
    return nullptr;
}

template <class _T>
inline void add(Matrix<_T>& a, Matrix<_T>& b, Matrix<_T>& c)
{
    Matrix<_T>* temp = a + b;
    c = *temp;
    delete temp;
}

template <class _T>
inline void sub(Matrix<_T>& a, Matrix<_T>& b, Matrix<_T>& c)
{
    Matrix<_T>* temp = a - b;
    c = *temp;
    delete temp;
}

template <class _T>
inline void mul(Matrix<_T>& a, Matrix<_T>& b, Matrix<_T>& c)
{
    Matrix<_T>* temp = a * b;
    c = *temp;
    delete temp;
}

template <class T>
Matrix<T>* Matrix<T>::Strassen(const Matrix<T>& src1, const Matrix<T>& src2) const
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
            matrices[0][i][j] = const_cast<Matrix<T>&>(src1)[i][j];
            matrices[1][i][j] = const_cast<Matrix<T>&>(src1)[i][j+size];
            matrices[2][i][j] = const_cast<Matrix<T>&>(src1)[i+size][j];
            matrices[3][i][j] = const_cast<Matrix<T>&>(src1)[i+size][j+size];
            matrices[4][i][j] = const_cast<Matrix<T>&>(src2)[i][j];
            matrices[5][i][j] = const_cast<Matrix<T>&>(src2)[i][j+size]; 
            matrices[6][i][j] = const_cast<Matrix<T>&>(src2)[i+size][j];
            matrices[7][i][j] = const_cast<Matrix<T>&>(src2)[i+size][j+size];                       
        }
    }

    add(matrices[0], matrices[3], matrices[19]); 
    add(matrices[4], matrices[7], matrices[20]); 
    mul(matrices[19], matrices[20], matrices[12]); 

    add(matrices[2], matrices[3], matrices[19]); 
    mul(matrices[19], matrices[4], matrices[13]); 
    
    sub(matrices[5], matrices[7], matrices[20]);
    mul(matrices[0], matrices[20], matrices[14]); 
    
    sub(matrices[6], matrices[4], matrices[20]);
    mul(matrices[3], matrices[20], matrices[15]); 
    
    add(matrices[0], matrices[1], matrices[19]); 
    mul(matrices[19], matrices[7], matrices[16]); 
    
    sub(matrices[2], matrices[0], matrices[20]);
    add(matrices[4], matrices[5], matrices[19]); 
    mul(matrices[20], matrices[19], matrices[17]); 
    
    sub(matrices[1], matrices[3], matrices[20]);
    add(matrices[6], matrices[7], matrices[19]); 
    mul(matrices[20], matrices[19], matrices[18]);
    
    sub(matrices[18], matrices[16], matrices[20]);
    add(matrices[12], matrices[15], matrices[19]); 
    add(matrices[19], matrices[20], matrices[8]); 

    add(matrices[14], matrices[16], matrices[9]); 

    add(matrices[13], matrices[15], matrices[10]); 

    sub(matrices[12], matrices[13], matrices[20]);
    add(matrices[14], matrices[17], matrices[19]); 
    add(matrices[20], matrices[19], matrices[11]);

    Matrix<T>* ret = new Matrix<T>(2*size, 2*size);
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            (*ret)[i][j] = matrices[8][i][j];
            (*ret)[i][j+size] = matrices[9][i][j];
            (*ret)[i+size][j] = matrices[10][i][j];
            (*ret)[i+size][j+size] = matrices[11][i][j];                                    
        }
    }
    return ret;
}

template <class T>
Matrix<T>* Matrix<T>::Division(const Matrix<T>& src1, const Matrix<T>& src2) const
{
    try
    {
        if (src1.isEmpty() || src2.isEmpty())
            throw(EmptyMatrixException("Matrix is empty", __func__, __FILE__, __LINE__)); 
        else if (!src1.isSquare() || !src2.isSquare())
            throw(MatrixNotSquareException("Not square matrix", __func__, __FILE__, __LINE__));
        else if (src1.getRow() == src2.getRow())
            throw(SizeMismatchException("Matrices do not match in size", __func__, __FILE__, __LINE__));
        else if (!src2.isInvertible())
            throw(InverseNotExistException("Not invertible matrix", __func__, __FILE__, __LINE__));
        else
        {
            Matrix<T>* inv = src2.inverse();
            Matrix<T>* rst = src1 * inv;
            delete inv;
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
    catch(const MatrixNotSquareException& e)
    {
        std::cerr << "MatrixNotSquareException: " << e.what() << '\n';
    }
    catch(const InverseNotExistException& e)
    {
        std::cerr << "InverseNotExistException: " << e.what() << '\n';
    }
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return nullptr;
}

template <class T>
Matrix<T>* Matrix<T>::Division(const Matrix<T>& src1, const T& src2) const
{
    try
    {
        if (src1.isEmpty())
            throw(EmptyMatrixException("Matrix is empty", __func__, __FILE__, __LINE__)); 
        else if (src2 == 0)
            throw(ArithmeticException("Divide by Zero", __func__, __FILE__, __LINE__));
        else
        {
            uint32 r = src1.getRow(), c = src1.getCol();
            Matrix<T>* rst = new Matrix<T>(r, c);
            for (int i = 0; i < r; i++)
                for (int j = 0; j < c; j++)
                    (*rst)[i][j] = const_cast<Matrix<T>&>(src1)[i][j] / src2;
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
    return nullptr;
}

template <class T>
Matrix<T>* Matrix<T>::operator+(const Matrix<T>& m)
{
    return Addition(*this, m);
}

template <class T>
Matrix<T>* Matrix<T>::operator-(const Matrix<T>& m)
{
    return Subtraction(*this, m);
}

template <class T>
Matrix<T>* Matrix<T>::operator*(const Matrix<T>& m)
{
    return Multiplication(*this, m);
}

template <class T>
Matrix<T>* Matrix<T>::operator*(const T& a)
{
    return Multiplication(*this, a);
}

template <class T>
Matrix<T>* operator*(const T& a, const Matrix<T>& b)
{
    return b*a;
}    

template <class T>
Matrix<T>* Matrix<T>::operator/(const Matrix<T>& m)
{
    return Division(*this, m);
}

template <class T>
Matrix<T>* Matrix<T>::operator/(const T& a)
{
    return Division(*this, a);
}

template <class T>
Matrix<T>* Matrix<T>::operator^(int b)
{
    try
    {
        if (isEmpty())
            throw(EmptyMatrixException("Matrix is empty, can't find determinant", __func__, __FILE__, __LINE__)); 
        else if (!isSquare())
            throw(MatrixNotSquareException("Matrix must be square", __func__, __FILE__, __LINE__));
        else
        {
            Matrix<T>* ret = new Matrix<T>(row,col);
            for(int i=0;i<row;i++)
                for(int j=0;j<col;j++) (*ret)[i][j]=(i==j);
            Matrix<T> cur(row,col);
            for(int i=0;i<row;i++)
                for(int j=0;j<col;j++) cur[i][j]=(i==j);
            for(;b;b>>=1,cur*=(*this))
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
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
}

template <class T>
void Matrix<T>::operator=(const Matrix<T>& m)
{
    if (this != &m)
    {
        if (!(row==m.row && col==m.col))
        {
            release();
            allocate(m.row, m.col, m.mod);
        }
        clone(m);
    }
}

template <class T>
void Matrix<T>::operator+=(const Matrix<T>& m)
{
    Matrix<T>* ret = (*this)+m;
    (*this) = *ret;
    delete ret;
}

template <class T>
void Matrix<T>::operator-=(const Matrix<T>& m)
{
    Matrix<T>* ret = (*this)-m;
    (*this) = *ret;
    delete ret;
}

template <class T>
void Matrix<T>::operator*=(const Matrix<T>& m)
{
    Matrix<T>* ret = (*this)*m;
    (*this) = *ret;
    delete ret;
}

template <class T>
void Matrix<T>::operator*=(const T t)
{
    Matrix<T>* ret = (*this)*t;
    (*this) = *ret;
    delete ret;
}

template <class T>
bool Matrix<T>::operator==(const Matrix<T>& m)const
{
    if(row!=m.getRow()||col!=m.getCol()) return 0;
    for(int i=0;i<m.row;i++)
        for(int j=0;j<m.col;j++)
            if(data[i][j]!=const_cast<Matrix<T>&>(m)[i][j]) return 0;
    return 1;
}

template <class _T>
std::istream& operator>>(std::istream& is, Matrix<_T>& m)
{
    for (int i = 0; i < m.row; i++)
        for (int j = 0; j < m.col; j++)
            is >> m[i][j];
    return is;
}

template <class _T>
std::ostream& operator<<(std::ostream& os, Matrix<_T>& m)
{
    os << "[";
    for (int i = 0; i < m.row; i++)
    {
        if(i>0) os<<" ";
        for (int j = 0; j < m.col; j++)
            if (i == m.row-1 && j == m.col-1)
                os << m[i][j] << "]";
            else if(j == m.col-1) os<<m[i][j]<<";";
            else 
                os << m[i][j] << ",";
        os << "\n";
    }
    return os;
}

template<class T>
Matrix<T>* Matrix<T>::toDiagnal()const
{
    if (mod==0)
    {
        uint32 r = row, c = col;
        Matrix<T>* rst = new Matrix<T>(r,c);
        (*rst) = *this;
        for (int i = 0; i < r; i++)
        {
            int id=-1;
            for (int j = i; j < r; j++)
                if ((*rst)[j][i])
                {
                    id=j;
                    break;
                }
            if (id!=i)
            {
                for (int j = 0; j < c; j++)
                {
                    std::swap((*rst)[id][j],(*rst)[i][j]);
                }
            }
            for (int j = i + 1; j < r; j++)
            {
                for (int k = c - 1; k > i; k--)
                {
                    (*rst)[j][k] -= (*rst)[i][k] * (*rst)[j][i] / (*rst)[i][i];
                }
            }
        }
        for(int i=r-1;i>=0;i--)
        {
            for(int j=r-1;j>i;j--)
            {
                if((*rst)[j][j]==0) continue;
                T t = (*rst)[i][j]/(*rst)[j][j];
                for(int k=i;k<=j;k++) (*rst)[i][k]-=t*(*rst)[j][j];
            }
        }
        return rst;
    }
    else {
        uint32 r = row, c = col, mo = mod;
        Matrix<T>* rst = new Matrix<T>(r, c, mo);
        for (int i = 0; i < r; i++)
        {
            int id=-1;
            for (int j = i; j < r; j++)
                if ((*rst)[j][i])
                {
                    id=j;
                    break;
                }
            if (id!=i)
            {
                for (int j = 0; j < c; j++)
                {
                    std::swap((*rst)[id][j],(*rst)[i][j]);
                }
            }
            uint32 p=1,k=mo-2,a=(*rst)[i][i];
            for (;k;k>>=1){
                if (k&1)p=1ull* p * a %mo;
                a=1ull * a *a % mo;
            }
            for (int j= i + 1; j < c; j++){
                (*rst)[i][j] = 1ull * (*rst)[i][j] * p %mo;   
            }
            for (int j = i + 1; j < r; j++)
            {
                for (int k = c - 1; k > i; k--)
                {
                    (*rst)[j][k] = ((*rst)[j][k] - (*rst)[i][k] * (*rst)[j][i]) % mo;
                }
            }
        }
        for(int i=r-1;i>=0;i--)
        {
            for(int j=r-1;j>i;j--)
            {
                if((*rst)[j][j]==0) continue;
                T t = (*rst)[i][j]*ksm((*rst)[j][j],mod-2);
                for(int k=i;k<=j;k++) (*rst)[i][k]=(((*rst)[i][k]-1ll*t*(*rst)[j][j])%mod+mod)%mod;
            }
        }
        return rst;
    }
}
template<class T>
T Matrix<T>::determinant()const
{
    try
    {
        if (isEmpty())
            throw(EmptyMatrixException("Matrix is empty, can't find determinant", __func__, __FILE__, __LINE__)); 
        else if (!isSquare())
            throw(MatrixNotSquareException("Matrix must be square", __func__, __FILE__, __LINE__));
        else 
        {   
            #if mod==0
            {
                uint32 r = row, c = col;
                Matrix<T> rst(r,c);
                rst = (*this);
                T res=1;
                for (int i = 0; i < r; i++)
                {
                    int id=-1;
                    for (int j = i; j < r; j++)
                        if (rst[j][i])
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
                            std::swap(rst[id][j],rst[i][j]);
                        }
                    }
                    res = res * rst[i][i];
                    for (int j = i + 1; j < r; j++)
                    {
                        for (int k = c - 1; k >= i; k--){
                            rst[j][k] *= rst[i][i];
                        }
                        for (int k = c - 1; k >= i; k--)
                        {
                            rst[j][k] -= rst[i][k] * rst[j][i] / rst[i][i];
                        }
                    }
                }
                return res;
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
                        if (rst[j][i])
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
                            std::swap(rst[id][j],rst[i][j]);
                        }
                    }
                    res = 1ull* res * rst[i][i] % mo;
                    uint32 p=1,k=mo-2,a=rst[i][i];
                    for (;k;k>>=1){
                        if (k&1)p=1ull* p * a %mo;
                        a=1ull * a *a % mo;
                    }
                    for (int j= i + 1; j < c; j++){
                        rst[i][j] = 1ull * rst[i][j] * p %mo;   
                    }

                    for (int j = i + 1; j < r; j++)
                    {
                        for (int k = c - 1; k >= i; k--)
                        {
                            rst[j][k] = (rst[j][k] - rst[i][k] * rst[j][i]) % mo;
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
            throw(EmptyMatrixException("Matrix is empty, can't find rank", __func__, __FILE__, __LINE__)); 
        uint32 res=0;
        uint32 r = row, c = col;
        Matrix<T> rst(r, c);
        rst = (*this);
        for (int i = 0; i < r && i < c; i++)
        {
            int id=-1;
            for (int j = i; j < r; j++)
                if (rst[j][i])
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
                    std::swap(rst[id][j],rst[i][j]);
                }
            }
            for (int j = i + 1; j < r; j++)
            {
                for (int k = c - 1; k >= i; k--)
                    rst[j][k] = rst[j][k] * rst[i][i];
                for (int k = c - 1; k >= i; k--)
                {
                    rst[j][k] -= rst[i][k] * rst[j][i] / rst[i][i];
                }
            }
        }
        return res;
    }

    catch(const EmptyMatrixException& e)
    {
        std::cerr << "EmptyMatrixException: " << e.what() << '\n';
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
            throw(EmptyMatrixException("Matrix is empty, can't find trace", __func__, __FILE__, __LINE__)); 
        else if (!isSquare())
            throw(MatrixNotSquareException("Matrix must be square", __func__, __FILE__, __LINE__));
        else 
        {   
            T res=0;
            for (int i=0;i<row;i++){
                res=res+data[i][i];
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
            throw(InvalidArgsException("range can not be zero", __func__, __FILE__, __LINE__));
        else if (row.start >= getRow() || col.start >= getCol())
            throw(RangeOutOfBoundException("range is out of bound", __func__, __FILE__, __LINE__));
        else
        {
            T m = data[row.start][col.start];
            for (int i = row.start; i < std::min(getRow(), row.end); i++)
                for (int j = col.start; j < std::min(getCol(), col.end); j++)
                    if (m < data[i][j])
                        m = data[i][j];
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
            throw(InvalidArgsException("range can not be zero", __func__, __FILE__, __LINE__));
        else if (row.start >= getRow() || col.start >= getCol())
            throw(RangeOutOfBoundException("range is out of bound", __func__, __FILE__, __LINE__));
        else
        {
            T m = data[row.start][col.start];
            for (int i = row.start; i < std::min(getRow(), row.end); i++)
                for (int j = col.start; j < std::min(getCol(), col.end); j++)
                    if (m > data[i][j])
                        m = data[i][j];
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
            throw(InvalidArgsException("range can not be zero", __func__, __FILE__, __LINE__));
        else if (row.start >= getRow() || col.start >= getCol())
            throw(RangeOutOfBoundException("range is out of bound", __func__, __FILE__, __LINE__));
        else
        {
            T tot = 0;
            uint32 cnt = 0;
            for (int i = row.start; i < std::min(getRow(), row.end); i++)
                for (int j = col.start; j < std::min(getCol(), col.end); j++)
                    tot += data[i][j],cnt++;
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
            throw(InvalidArgsException("range can not be zero", __func__, __FILE__, __LINE__));
        else if (row.start >= getRow() || col.start >= getCol())
            throw(RangeOutOfBoundException("range is out of bound", __func__, __FILE__, __LINE__));
        else
        {
            T tot = 0;
            for (int i = row.start; i < std::min(getRow(), row.end); i++)
                for (int j = col.start; j < std::min(getCol(), col.end); j++)
                    tot += data[i][j];
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
            throw(SizeMismatchException("the total size must be equal", __func__, __FILE__, __LINE__));
        else
        {
            T tot = 0;
            for(int i=0;i<row*col;i++)
            {
                int rowx = i/col, colx = i%col, rowy = i/m.getCol(), coly = i%m.getCol();
                tot += data[rowx][colx]*const_cast<Matrix<T>&>(m)[rowy][coly];
            }
            return tot;
        }
    }catch(SizeMismatchException& e)
    {
        std::cerr << "SizeMismatchException: "<<e.what()<<'\n';
    }
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return 0;
}
template <class T>
Matrix<T>* Matrix<T>::crossProduct(const Matrix<T>& m)const
{
    try{
        if(row!=1||col!=3||m.getRow()!=1||m.getCol()!=3) 
            throw(SizeMismatchException("row of both matrix must be 1 and col must be 3", __func__, __FILE__, __LINE__));
        else if(typeid(data[0][0])!=typeid(m[0][0]))
            throw(TypeMismatchException("data type in the matrix does not match",__func__, __FILE__, __LINE__));
        else
        {
            Matrix<T>* ret = new Matrix<T>(1,3);
            (*ret)[0][0] = data[0][1]*const_cast<Matrix<T>&>(m)[0][2]-data[0][2]*const_cast<Matrix<T>&>(m)[0][1];
            (*ret)[0][1] = data[0][2]*const_cast<Matrix<T>&>(m)[0][0]-data[0][0]*const_cast<Matrix<T>&>(m)[0][2];
            (*ret)[0][2] = data[0][0]*const_cast<Matrix<T>&>(m)[0][1]-data[0][1]*const_cast<Matrix<T>&>(m)[0][0];
            return ret;
        }
    }catch(SizeMismatchException& e)
    {
        std::cerr << "SizeMismatchException: "<<e.what()<<'\n';
    }
    catch(TypeMismatchException& e)
    {
        std::cerr << "TypeMismatchException: "<<e.what()<<'\n';
    }
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return nullptr;
}
template <class T>
Matrix<T>* Matrix<T>::transpose()const
{   
    try{
        if (isEmpty())
            throw(EmptyMatrixException("Matrix is empty, can't transpose", __func__, __FILE__, __LINE__)); 
        else
        {
            Matrix<T>* ret = new Matrix<T>(col,row,mod);
            for(int i=0;i<col;i++)
                for(int j=0;j<row;j++) 
                    (*ret)[i][j]=data[j][i];
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
    return nullptr;
}
template <class T>
Matrix<T>* Matrix<T>::inverse()const
{
    try{
        if (isEmpty())
            throw(EmptyMatrixException("Matrix is empty, can't find inverse", __func__, __FILE__, __LINE__)); 
        else if (!isSquare())
            throw(MatrixNotSquareException("non-square matrix is not invertible", __func__, __FILE__, __LINE__));
        else if (!isInvertible())
                throw(InverseNotExistException("singular matrix is not invertible", __func__, __FILE__, __LINE__));
        else
        {
            Matrix<T> t(row,row+col,mod);
            for(int i=0;i<row;i++)
                for(int j=0;j<row+col;j++)
                {
                    if(j<col) t[i][j]=data[i][j];
                    else t[i][j]=(j-col==i); 
                }
            Matrix<T>* m = t.toDiagnal();
            Matrix<T>* ret = new Matrix<T>(row,col,mod);
            for(int i=0;i<row;i++)
                for(int j=0;j<col;j++)
                {
                    if(mod) (*ret)[i][j] = 1ll*(*m)[i][row+j]*ksm((*m)[i][i],mod-2)%mod;
                    else (*ret)[i][j] = (*m)[i][row+j]/(*m)[i][i];
                }
            delete m;
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
    return nullptr;
}
template <class T>
Matrix<T>* Matrix<T>::subMatrix(const Range& row, const Range& col)const
{
    try
    {
        if (row.empty() || col.empty())
            throw(InvalidArgsException("range can not be zero", __func__, __FILE__, __LINE__));
        else if (row.start >= getRow() || col.start >= getCol())
            throw(RangeOutOfBoundException("range is out of bound", __func__, __FILE__, __LINE__));
        else
        {
            Matrix<T>* ret = new Matrix<T>(std::min(getRow(), row.end)-row.start,std::min(getCol(), col.end)-col.start,mod);
            for(int i=0;i<ret->getRow();i++)
                for(int j=0;j<ret->getCol();j++)
                    (*ret)[i][j]=data[i+row.start][j+col.start];
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
    return nullptr;
}
template <class T>
Matrix<T>* Matrix<T>::reshape(const uint32 _row, const uint32 _col)const
{
    try{
        if(1ll*_row*_col!=1ll*row*col)
            throw(SizeMismatchException("cannot fit into this shape", __func__, __FILE__, __LINE__));
        else
        {
            Matrix<T>* ret = new Matrix<T>(_row,_col,mod);
            for(int i=0;i<_row;i++)
                for(int j=0;j<_col;j++)
                {
                    long long id = i*_col+j;
                    int x = id/col,y=id%col;
                    (*ret)[i][j] = data[x][y];
                }
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
    return nullptr;
}
template <class T>
Matrix<T>* Matrix<T>::slice(const Range& row, const uint32 col)const
{
    try{
        if (row.empty())
            throw(InvalidArgsException("range can not be zero", __func__, __FILE__, __LINE__));
        if(col>=this->col||row.start>=this->row)
            throw(RangeOutOfBoundException("range out of bound", __func__, __FILE__, __LINE__));
        else
        {
            Matrix<T>* ret = new Matrix<T>(std::min(this->getRow(), row.end) - row.start,1,mod);
            for(int i=row.start;i< ret->getRow();i++)
                (*ret)[i-row.start][0] = data[i][col];
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
    return nullptr;
}

template <class T>
Matrix<T>* Matrix<T>::convolute(const Matrix& kernal,uint32 anchor_x,uint32 anchor_y)const
{
    try{
        if (kernal.isEmpty())
            throw(EmptyMatrixException("Matrix is empty", __func__, __FILE__, __LINE__)); 
        if(anchor_x>=kernal.getRow()||anchor_y>=kernal.getCol())
            throw(RangeOutOfBoundException("range out of bound", __func__, __FILE__, __LINE__));
        else{
            Matrix<T>* ret = new Matrix<T>(row,col);
            for(int i=0;i<row;i++)
                for(int j=0;j<col;j++)
                    for(int k=anchor_x-i;k<std::min(kernal.getRow(),row+anchor_x-i);k++)
                    {
                        //i-anc_x+k>=0&&i-anc_x+k<row
                        //k>=anc_x-i && k< row +anc_x-i
                        for(int l=anchor_y-j;l<std::min(kernal.getCol(),col+anchor_y-j);l++)
                        {
                            int x = i-anchor_x+k,y = j-anchor_y+l;
                            #if mod>0
                                (*ret)[i][j]=((*ret)[i][j]+1ll*const_cast<Matrix<T>&>(kernal)[x][y]*data[i][j])%mod;
                            #else
                                (*ret)[i][j]+=const_cast<Matrix<T>&>(kernal)[x][y]*data[i][j];
                            #endif
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
    return nullptr;
}

template <class T>
std::vector<std::complex<double>> Matrix<T>::eigenValue() const
{
    try{
        if (isEmpty())
            throw(EmptyMatrixException("Matrix is empty, can't find eigenvalue", __func__, __FILE__, __LINE__)); 
        else if (!isSquare())
            throw(MatrixNotSquareException("non-square matrix does not have eigenValue", __func__, __FILE__, __LINE__));
        else
        {
            int n = row;
            Matrix<T> A(n,n);
            A = (*this);
            Matrix<T> U(n,n);
            for (int i=0;i<n;i++){
                for (int j=0;j<n;j++){
                    if (i==j){
                        U[i][j]=1;
                    }
                    else{
                        U[i][j]=0;
                    }
                }
            }
            Matrix<T> Q(n,n);
            Matrix<T> R(n,n);
            Matrix<T>* temp = nullptr;
            for (int tim=0;tim<50;tim++){
                QR_Decomposition(A,Q,R);
                temp = R*Q;
                A = *temp;
                delete temp;
                temp = U*Q;
                U = *temp;
                delete temp;
            }
            std::vector<std::complex<double>> res;
            for (int i=0;i<n;i++){
                res.push_back(A[i][i]);
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
    catch(const Exception& e)
    {
        std::cerr << "Fatal: " << e.what() << '\n';
    }
    return std::vector<std::complex<double>>();
}

template <class T>
Matrix<std::complex<double>>* Matrix<T>::eigenVector(std::complex<double> val)const{
    
    try{
        int n=col;
        Matrix<std::complex<double>> a (this->row,this->col);
        for (int i=0;i<n;i++){
            for (int j=0;j<n;j++){
                a[i][j]=-this->data[i][j];
                if (i==j)a[i][j]=-this->data[i][j]+val;
            }
        }
        for (int i=0;i<n;i++){
            int id=-1;
            for (int j=i;j<n;j++){
                if (a[j][i] != 0.0){
                    id=j;
                    break;
                }
            }
            if (id==-1){
                throw(WrongEigenValueException("eigenvalue is wrong", __func__, __FILE__, __LINE__));
            }
            if (id!=i){
                for (int j=0;j<n;j++){
                    swap(a[i][j],a[id][j]);
                }
            }
            for (int j=i+1;j<n;j++){
                for (int k=n;k>=i;k--){
                    a[j][k]-=a[i][k]*a[j][i]/a[i][i];
                }
            }
        }
        Matrix<std::complex<double>>* res = new Matrix<std::complex<double>>(1,n);
        (*res)[0][n-1]=1;
        for (int i=n-2;i>=0;i--){
            std::complex<double>tmp=0;
            for (int j=i+1;j<n;j++){
                tmp-=a[i][j]*(*res)[0][j];
            }
            (*res)[0][i]=tmp;
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
    return nullptr;
}
template <class T>
void QR_Decomposition(Matrix<T>& A,Matrix<T>& Q,Matrix<T>& R){
    uint32 n = A.getRow();
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            Q[i][j]=R[i][j]=0;
        }    
    }
    for (int k = 0; k < n; k++){
		double MOD = 0;
		for (int i = 0; i < n; i++)
		{
			MOD += A[i][k] * A[i][k]; 
		}
		R[k][k] = sqrt(MOD);
		for (int i = 0; i < n; i++)
		{
			Q[i][k] = A[i][k] / R[k][k]; 
		}

		for (int i = k + 1; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				R[k][i] += A[j][i] * Q[j][k];
			}
			for (int j = 0; j < n; j++)
			{
				A[j][i] -= R[k][i] * Q[j][k]; 
			}
		}
	}
}
} // namespace usr


#endif