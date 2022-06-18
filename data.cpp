#include <bits/stdc++.h>
#include <cstdlib>
#include <time.h>
using namespace std;
inline int read();
int main()
{
    freopen("datainput.txt", "r", stdin);
    freopen("data.in", "w", stdout);
    int max = read(), min = read();
    random_device rd;
    default_random_engine eng(rd());
    uniform_real_distribution<double> distr(min, max);

    uint32_t row, col, channel;

    // test1
    row = read(), col = read();
    printf("%d %d\n", row, col);
    for (uint32_t i = 0; i < row*col; i++)
    {
        printf("%f ", distr(eng));
    }
    putchar('\n');

    // test2
    row = read(), col = read();
    printf("%d %d\n", row, col);
    for (uint32_t i = 0; i < row; i++)
    {
        for (uint32_t j = 0; j < col; j++)
        {
            if (j >= i)
                printf("%f ", distr(eng));
            else 
                printf("%f ", 0.0);
        }
    }
    putchar('\n');

    // test3
    row = read(), col = read();
    printf("%d %d\n", row, col);
    for (uint32_t i = 0; i < row*col; i++)
    {
        printf("%f ", distr(eng));
    }
    putchar('\n');
    printf("%d %d\n", row, col);
    for (uint32_t i = 0; i < row*col; i++)
    {
        printf("%f ", distr(eng));
    }
    putchar('\n');

    // test4
    row = read(), col = read();
    printf("%d %d\n", row, col);
    for (uint32_t i = 0; i < row*col; i++)
    {
        printf("%f ", distr(eng));
    }
    putchar('\n');
    printf("%d %d\n", row, col);
    for (uint32_t i = 0; i < row*col; i++)
    {
        printf("%f ", distr(eng));
    }
    putchar('\n');

    // test5
    row = read(), col = read(), channel = read();
    printf("%d %d %d\n", row, col, channel);
    for (int i = 0; i < row*col*channel; i++)
    {
        printf("%f ", distr(eng));
    }
    putchar('\n');
    row = read(), col = read(), channel = read();
    printf("%d %d %d\n", row, col, channel);
    for (int i = 0; i < row*col*channel; i++)
    {
        printf("%f ", distr(eng));
    }
    putchar('\n');

    // test6
    row = read(), col = read();
    printf("%d %d\n", row, col);
    for (uint32_t i = 0; i < 2*row*col; i++)
    {
        printf("%f ", distr(eng));
    }
    putchar('\n');

}

inline int read(){
    register int x=0, f=1;
    register char ch = getchar();
    while(ch<'0'||ch>'9')
    {
        if(ch=='-') f=-1;
        ch=getchar();
    }
    while(ch>='0'&&ch<='9')
    {
        x=(x<<1)+(x<<3)+(ch^48);
        ch=getchar();
    }
    return x*f;
}