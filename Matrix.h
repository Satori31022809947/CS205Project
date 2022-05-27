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

template <class T>

class Matrix {
    private:
        int col;
        int row;

        
        T** a;

    public:
        /**
         * @description: construct a Matrix
         * @param {int} col 
         * @param {int} row
         * @return {*}
         */
        Matrix(int row=0,int col=0);
        ~Matrix();
        Matrix operator+(const Matrix& b);
        Matrix operator-(const Matrix& b);
        Matrix operator*(const Matrix& b);
        Matrix operator/(const Matrix& b);
        
        /**
         * @description: return true if is vector
         * @return {*}
         */
        bool isVector()const;
        /**
         * @description: return true if is square
         * @return {*}
         */
        bool isSquare()const; 

        T determinant()const;
        int rank()const;
        T max()const;
        T min()const;
        T avg()const;
        T sum()const;

        Matrix subMatrix(int stx,int sty,int lenx=row-stx,int leny=col-sty)const;
        Matrix transposition()const;
        Martix Hermite()const;
        friend istream &operator>>(istream&, Matrix&);
		friend ostream &operator<<(ostream&, Matrix&);
        Matrix operator=(T *);
		Matrix& operator+=(const Matrix &m);
		Matrix& operator-=(const Matrix &m);
		Matrix& operator*=(const Matrix &m);
        T*operator[](int i){ return matrix + i*col; }

};
#endif Matrix_H_