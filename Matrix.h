/*
 * @Author: Satori 3102809947@qq.com
 * @Date: 2022-05-26 11:51:58
 * @LastEditors: Satori 3102809947@qq.com
 * @LastEditTime: 2022-05-27 22:54:28
 * @FilePath: \CS205Project\Matrix.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#ifndef Matrix_H_
#define Matrix_H_
#include <iostream>
#include "Range.h"
#include "Exception.h"
namespace usr
{

template <class T>

class Matrix {
    private:
        uint32 col;
        uint32 row;
        T** data;

        Matrix Addition(const Matrix& src1, const Matrix& src2) const;
        Matrix Subtraction(const Matrix& src1, const Matrix& src2) const;
        Matrix Multiplication(const Matrix& src1, const Matrix& src2) const;
        Matrix Multiplication(const Matrix& src1, const T& src2) const;
        Matrix Strassen(const Matrix& src1, const Matrix& src2) const;
        Matrix Division(const Matrix& src1, const Matrix& src2) const;
        Matrix Division(const Matrix& src1, const T& src2) const;

    public:
        /**
         * @description: construct a Matrix
         * @param {int} col 
         * @param {int} row
         * @return {*}
         */
        Matrix(uint32 row=1,uint32 col=1);
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
        inline bool isEmpty()const {return data==nullptr};

        /**
         * return true if the matrix is a vector
         */
        inline bool isVector()const {return row==1||col==1};

        /**
         * return true if the matrix is square
         */
        inline bool isSquare()const {return row==col;}; 

        /**
         * return true if the matrix is invertible
         */
        inline bool canInverse()const {};

        /**
         * @description: Deep copy
         */
        Matrix clone()const;

        /** @brief Free data
         *  Set data pointer as nullptr after deleting 
         */
        inline void release();

        Matrix operator+(const Matrix& b);
        Matrix operator-(const Matrix& b);
        Matrix operator*(const Matrix& b);
        Matrix operator*(const T& b);
        friend Matrix operator*(const T& a, const Matrix& b);    
        Matrix operator/(const Matrix& b);
        Matrix operator/(const T& b);
        Matrix operator^(const int b);
        Matrix operator=(const Matrix& m);
		Matrix& operator+=(const Matrix& m);
		Matrix& operator-=(const Matrix& m);
		Matrix& operator*=(const Matrix& m);
        Matrix& operator*=(const T t);
        bool operator==(const Matrix& m)const;
        T* operator[](uint32 i){ return data + i*col; }
        void* operator new(size_t sz, char* filename, int line);
        void* operator new[](size_t sz, char* filename, int line);
        void operator delete(void* p);
        void operator delete[](void* p);
        friend istream& operator>>(istream&, Matrix&)const;
		friend ostream& operator<<(ostream&, Matrix&)const;

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

        T eigenValue()const;
        Matrix eigenVector()const;
        Matrix subMatrix(const Range& row=Range::all(), const Range& col=Range::all())const;
        Matrix inverse()const;
        Matrix transpoe()const;
        Martix hermite()const;
        Matrix reshape(const uint32 _row, const uint32 _col)const;
        Matrix slice(const Range& row, const uint32 col)const;
        Matrix convolute(const Matrix& kernel)const;
};

template <class T>
Matrix<T>::Matrix(uint32 _row, uint32 _col)
{
    if (_row==0 || _col==0)
        throw(InvalidArgsException("row or col can not be zero", __func__, __FILE__, __LINE__));
    row = _row, col = _col;
    data = new T* [row];
    for (int i = 0; i < col; i++)
        data[i] = new T[col];
}

template <class T>
Matrix<T>::Matrix(const Matrix<T>& m)
{
    row = m.row, col = m.col;
    data = m.data;
}

template <class T>
Matrix<T>::~Matrix()
{
    if (!this->isEmpty())
        release();
}

template <class T>
inline void Matrix<T>::release()
{
    for (int i = 0; i < row; i++)
        delete[] data[i];
    delete[] data;
    data = nullptr;
}

template <class T>
Matrix<T> Matrix<T>::clone() const
{
    Matrix<T> m(row, col);
    for (int i = 0; i < row; i++)
        for (int j = 0; j < col; j++)
            m[i][j] = this->data[i][j];
    return m;
}

template <class T>
Matrix<T> Matrix<T>::Addition(const Matrix<T>& src1, const Matrix<T>& src2) const
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
            Matrix<T> rst(r, c);
            for (int i = 0; i < r; i++)
                for (int j = 0; j < c; j++)
                    rst[i][j] = src1[i][j] + src2[i][j];
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
    return Matrix<T>();
}

template <class T>
Matrix<T> Matrix<T>::Subtraction(const Matrix<T>& src1, const Matrix<T>& src2) const
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
            Matrix<T> rst(r, c);
            for (int i = 0; i < r; i++)
                for (int j = 0; j < c; j++)
                    rst[i][j] = src1[i][j] - src2[i][j];
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
    return Matrix<T>();
}

template <class T>
Matrix<T> Matrix<T>::Multiplication(const Matrix<T>& src1, const Matrix<T>& src2) const
{
    try
    {
        if (src1.isEmpty() || src2.isEmpty())
            throw(EmptyMatrixException("Matrix is empty", __func__, __FILE__, __LINE__)); 
        else if (src1.getCol() != src2.getRow())
            throw(SizeMismatchException("Matrices do not match in size", __func__, __FILE__, __LINE__));
        else
        {
            uint32 r = src1.getRow(), c = src2.getCol();
            Matrix<T> rst(r, c);
            /* TODO */
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
Matrix<T> Matrix<T>::Multiplication(const Matrix<T>& src1, const T& src2) const
{
    try
    {
        if (src1.isEmpty())
            throw(EmptyMatrixException("Matrix is empty", __func__, __FILE__, __LINE__)); 
        else
        {
            uint32 r = src1.getRow(), c = src1.getCol();
            Matrix<T> rst(r, c);
            for (int i = 0; i < r; i++)
                for (int j = 0; j < c; j++)
                    rst[i][j] = src1[i][j] * src2;
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

}

template <class T>
Matrix<T> Matrix<T>::Division(const Matrix<T>& src1, const Matrix<T>& src2) const
{
    try
    {
        if (src1.isEmpty() || src2.isEmpty())
            throw(EmptyMatrixException("Matrix is empty", __func__, __FILE__, __LINE__)); 
        else if (!src1.isSquare() || !src2.isSquare())
            throw(ArithmeticException("Not square matrix", __func__, __FILE__, __LINE__));
        else if (src1.getRow() == src2.getRow())
            throw(SizeMismatchException("Matrices do not match in size", __func__, __FILE__, __LINE__));
        else if (!src2.canInverse())
            throw(ArithmeticException("Not invertible matrix", __func__, __FILE__, __LINE__));
        else
        {
            Matrix<T> inv = src2.inverse();
            Matrix<T> rst = src1 * inv;
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
Matrix<T> Matrix<T>::Division(const Matrix<T>& src1, const T& src2) const
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
            Matrix<T> rst(r, c);
            for (int i = 0; i < r; i++)
                for (int j = 0; j < c; j++)
                    rst[i][j] = src1[i][j] / src2;
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

template <class T>
Matrix<T> operator*(const T& a, const Matrix<T>& b)
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
Matrix<T> Matrix<T>::operator^(const int b)
{

}

template <class T>
Matrix<T> Matrix<T>::operator=(const Matrix<T>& m)
{

}

template <class T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& m)
{

}

template <class T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& m)
{

}

template <class T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& m)
{

}

template <class T>
Matrix<T>& Matrix<T>::operator*=(const T t)
{

}

template <class T>
bool Matrix<T>::operator==(const Matrix<T>& m)const
{

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

}

template <class T>
T Matrix<T>::avg(const Range &row, const Range &col)const
{

}

template <class T>
T Matrix<T>::sum(const Range &row, const Range &col)const
{

}

} // namespace usr
#endif