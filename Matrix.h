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
#include "Range.h"
namespace usr
{

template <class T>

class Matrix {
    private:
        int col;
        int row;
        T** data;

        Matrix Addition(const Matrix& src1, const Matrix& src2) const;
        Matrix Subtraction(const Matrix& src1, const Matrix& src2) const;
        Matrix Multiplication(const Matrix& src1, const Matrix& src2) const;
        Matrix Multiplication(const Matrix& src1, const T& src2) const;
        Matrix Strassen(const Matrix& src1, const Matrix& src2) const;
        Matrix Division(const Matrix& src1, const Matrix& src2) const;
        Matrix Division(const Matrix& src1, const T& src2) const;
        Matrix Division(const T& src1, const Matrix& src2) const;

    public:
        /**
         * @description: construct a Matrix
         * @param {int} col 
         * @param {int} row
         * @return {*}
         */
        Matrix(int row=0,int col=0);
        /** 
         * @description: Shallow copy 
         */
        Matrix(const Matrix&);
        /** @brief destructor
         *  release if is not empty
         */
        ~Matrix();

        inline int getRow()const {return row;}
        inline int getCol()const {return col;}
        
        /** 
         * return true if the matrix is empty (data is equal to nullptr)
         */
        inline bool isEmpty()const;

        /**
         * return true if the matrix is a vector
         */
        inline bool isVector()const;

        /**
         * return true if the matrix is square
         */
        inline bool isSquare()const; 

        /**
         * return true if the matrix is invertible
         */
        bool canInverse()const;

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
        friend Matrix operator/(const T& a, const Matrix& b);
        Matrix operator^(const int b);
        Matrix operator=(const Matrix& m);
		Matrix& operator+=(const Matrix& m);
		Matrix& operator-=(const Matrix& m);
		Matrix& operator*=(const Matrix& m);
        Matrix& operator*=(const T t);
        bool operator==(const Matrix& m)const;
        T* operator[](int i){ return data + i*col; }
        void* operator new(size_t sz, char* filename, int line);
        void* operator new[](size_t sz, char* filename, int line);
        void operator delete(void* p);
        void operator delete[](void* p);
        friend istream& operator>>(istream&, Matrix&)const;
		friend ostream& operator<<(ostream&, Matrix&)const;

        T determinant()const;
        int rank()const;
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
        Matrix reshape(const int _row, const int _col)const;
        Matrix slice(const Range& row, const int col)const;
        Matrix convolute(const Matrix& kernel)const;
};

} // namespace usr
#endif