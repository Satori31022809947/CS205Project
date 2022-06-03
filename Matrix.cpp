/*
 * @Author: Satori 3102809947@qq.com
 * @Date: 2022-05-27 23:01:36
 * @LastEditors: Satori 3102809947@qq.com
 * @LastEditTime: 2022-06-03 09:57:41
 * @FilePath: \CS205Project\matrix.cpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include "Matrix.h"
#include "Exception.h"


Matrix::Matrix Addition(const Matrix& src1, const Matrix& src2){
    if (src1.row!=src2.row||src1.col!=src2.col){
        throw SizeMismatchException();
    }

}
Martix::Matrix Subtraction(const Matrix& src1, const Matrix& src2) const;
        Matrix Multiplication(const Matrix& src1, const Matrix& src2) const;
        Matrix Multiplication(const Matrix& src1, const T& src2) const;
        Matrix Strassen(const Matrix& src1, const Matrix& src2) const;
        Matrix Division(const Matrix& src1, const Matrix& src2) const;
        Matrix Division(const Matrix& src1, const T& src2) const;


