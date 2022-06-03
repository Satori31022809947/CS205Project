/*
 * @Author: Satori 3102809947@qq.com
 * @Date: 2022-06-03 09:41:11
 * @LastEditors: Satori 3102809947@qq.com
 * @LastEditTime: 2022-06-03 10:12:01
 * @FilePath: \CS205Project\Exception.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#ifndef EXCEPTION_H_
#define EXCEPTION_H_
#include <iostream>
#include <sstream>
namespace usr
{
class Exception: public std::exception
{
    private:
        const char* err;
        std::string func;
        std::string file;
        int line;

    public:
        std::string msg;
        Exception();
        Exception(const char* _err, const std::string& _func, const std::string& _file, const int _line)
        : err(_err), func(_func), file(_file), line(_line)
        {
            std::stringstream ss;
            ss << _err << "\n" << "\tat " << _func << "(" << _file << ":" << _line << ")";
            msg = ss.str();
        }
        virtual inline const char* what() const throw() {return msg.c_str();}
};

class InvalidArgsException: public Exception 
{
    public:
        InvalidArgsException(const char* _err, const std::string& _func, const std::string& _file, const int _line)
        : Exception(_err, _func, _file, _line) {};
};

class SizeMismatchException: public Exception
{
    public:
        SizeMismatchException(const char* _err, const std::string& _func, const std::string& _file, const int _line)
        : Exception(_err, _func, _file, _line) {};
};

class TypeMismatchException: public Exception
{
    public:
        TypeMismatchException(const char* _err, const std::string& _func, const std::string& _file, const int _line)
        : Exception(_err, _func, _file, _line) {};
};

class ArithmeticException: public Exception
{
    public:
        ArithmeticException(const char* _err, const std::string& _func, const std::string& _file, const int _line)
        : Exception(_err, _func, _file, _line) {};
};

class EmptyMatrixException: public Exception
{
    public:
        EmptyMatrixException(const char* _err, const std::string& _func, const std::string& _file, const int _line)
        : Exception(_err, _func, _file, _line) {};
};

class RangeOutOfBoundException: public Exception
{
    public:
        RangeOutOfBoundException(const char* _err, const std::string& _func, const std::string& _file, const int _line)
        : Exception(_err, _func, _file, _line) {};
};

class MatrixNotSquareException: public Exception
{
    public:
        MatrixNotSquareException(const char* _err, const std::string& _func, const std::string& _file, const int _line)
        : Exception(_err, _func, _file, _line) {};
};

} // namespace usr
#endif