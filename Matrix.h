/*
 * @Author: Satori 3102809947@qq.com
 * @Date: 2022-05-26 11:51:58
 * @LastEditors: Satori 3102809947@qq.com
 * @LastEditTime: 2022-05-27 22:29:01
 * @FilePath: \CS205Project\Matrix.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#ifndef Matrix_H_
#define Matrix_H_

template <class T>

class Matrix {
    private:
        int row;
        int col;

        
        T** a;

    public:
        /**
         * @description: construct a Matrix
         * @param {int} col 
         * @param {int} row
         * @return {*}
         */
        Matrix(int col=0,int row=0);
        Matrix operator+(const Matrix& b);
        Matrix operator-(const Matrix& b);
        Matrix operator*(const Matrix& b);
        Matrix operator/(const Matrix& b);

        bool isVector()const;
        bool isSquare()const; 


        
        
        


        
        

};
#endif Martix_H_