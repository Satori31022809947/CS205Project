#ifndef DETECT_H_
#define DETECT_H_
#include <map>
#include <sstream>
#include <iostream>
using namespace std;
typedef map<uint64_t, string> Map;
#ifdef DEBUG
#define DEBUG_NEW new(__FILE__, __LINE__)

class MemoryDetech
{
    private:
        MemoryDetech() {}
        Map mp;
    public:
        static MemoryDetech& instance()
        {
            static MemoryDetech md;
            return md;
        }
        void add(void* ptr, const char* file, int line)
        {
            ostringstream ss;
            ss << file << " " <<line;
            mp.insert(make_pair(reinterpret_cast<uint64_t>(ptr), ss.str()));
        }
        void remove(void* ptr)
        {
            mp.erase(reinterpret_cast<uint64_t>(ptr));
        }
        void show()
        {
            if (mp.empty())
            {
                cout << "Well done, no memory leak!" << endl;
                return;
            }
            for(auto pair : mp){
                cout<<pair.second<<" memory leakage"<<endl;
            }
        }
};

void * operator new(std::size_t size, const char *file, int line)
{
    void *ptr = malloc(size);
    MemoryDetech::instance().add(ptr, file, line);
    return ptr;
}

void* operator new[](std::size_t size, const char* file, int line)
{
    return operator new(size, file, line);
}

void operator delete(void* ptr)
{
    MemoryDetech::instance().remove(ptr);
    free(ptr);
    ptr = nullptr;
}
void operator delete[](void* ptr) 
{
    return operator delete(ptr);
}
#else
#define DEBUG_NEW new
#endif

#endif