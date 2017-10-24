#include "neuron.h"
#include <iostream>
#include <vector>
#include <cmath>
#define INPUT       2
#define SIZE        5
#define SL_LENGTH   50

double fun(double x, double y);

int main()
{
    using std::vector;
    using std::cout;
    int * structure = new int[SIZE];
    double * u;
    double * y;
    structure[0]=5;
    structure[1]=5;
    structure[2]=4;
    structure[3]=6;
    structure[4]=1;

    Network net(SIZE,structure,INPUT);

    int ut              =           INPUT*SL_LENGTH;
    int yt              =    structure[SIZE-1]*SL_LENGTH;

    u                   = new double[ut];
    y                   = new double[yt];

    double x            = 0.01;
    double a            = 0.01;

    int j        =    0;
    bool cond    = true;
    for(int i=0;i<ut;i+=2)
    {
        u[i]     =    x;
        u[i+1]   =    a;
        cout<<"\n"<<i;
        if(cond)
        {
            x   +=  0.1;
            cond = false;
        }else{
            a   +=  0.1;
            cond = true;
        }
        y[j]=fun(x,a);
        cout<<"\n "<<x<<" "<<a<<" "<<y[j];
        j++;
    }
    std::cin.get();
    Learning_stream LS(u,ut,y,yt,SL_LENGTH,0.01,0.0001);

    grad_desc(net,LS);

    delete [] u;
    delete [] y;
    return 0;
}

double fun(double x, double y)
{
    return (x*x+y*y);
}

/*double fun(double x, double y)
{
    return (1+(1/(x*x)+(1/(sqrt(y)*sqrt(y)*sqrt(y)))))*(1+(1/(x*x)+(1/(sqrt(y)*sqrt(y)*sqrt(y)))));
}*/
