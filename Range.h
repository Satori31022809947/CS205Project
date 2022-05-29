#ifndef RANGE_H_
#define RANGE_H_
typedef unsigned int uint32;
namespace usr
{
/**
 * start is inclusive while end is exclusive
 * interval: [start, end)
 */
class Range
{
public:
    uint32 start, end;
    Range();
    Range(uint32 start, uint32 end);
    uint32 size() const;
    bool empty() const;
    static Range all();
};

inline
Range::Range(): start(0), end(1) {}

inline
Range::Range(uint32 _start, uint32 _end): start(_start), end(_end) {}

inline
uint32 Range::size() const
{
    return end - start;
}

inline
bool Range::empty() const
{
    return start == end;
}

inline
Range Range::all()
{
    return Range(0, __UINT32_MAX__);
}

} // namespace usr
#endif
