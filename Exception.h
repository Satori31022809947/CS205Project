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
        : err(_err), func(_func), file(_file), line(_line) {};
        virtual ~Exception() throw();
        virtual const char* what() const throw();
};

Exception::Exception(const char* _err, const std::string& _func, const std::string& _file, const int _line)
{
    std::stringstream ss;
    ss << _err << "\n" << "\tat " << _func << "(" << _file << ":" << _line << ")";
    msg = ss.str();
}

const char* Exception::what() const throw()
{
    return msg.c_str();
}

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

class NullPointerException: public Exception
{
    public:
        NullPointerException(const char* _err, const std::string& _func, const std::string& _file, const int _line)
        : Exception(_err, _func, _file, _line) {};
};

} // namespace usr
#endif