#include "neuron.h"
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>

using std::vector;
using std::cout;
using std::cin;

double fun(double x);

int main()
{
    std::srand(std::time(0));

    /*
    cout<<"\nLAYER_TEST_BEGGINING";
    Layer layer(10,20,'l');

    vector <double> x;
    vector <double> N;

    for(int i=0;i<20;i++)
    {
        x.push_back(i);
    }

    for(int i=0;i<10;i++)
    {
        N.push_back(i);
    }
    layer.process(x);
    layer.compute_delta(N);
    layer.der(0.1);
    vector <double> temp;
    layer>>temp;
    cout<<"\nLAYER_TEST_PASS";

                                */

    /*
    int length = 10;
    string typos;
    vector <int> structure;
    for(int i=5;i<=(length+5);i++)
        structure.push_back(i);
    for(int i=0;i<length;i++)
        typos[i]='l';

    vector <double> x;
    vector <double> y;

    for(int i=0;i<5;i++)
        x.push_back(i);

    Network net(structure, typos, length);
    cout<<"\nINPUT="<<net.input_length;
    cout<<"\nOUTPUT="<<net.output_length;

    y=net.process(x);

    for(int i=0;i<net.output_length;i++)
    cout<<"\nY="<<y[i];
                                            */
    string typos;
    vector <int> structure(3);
    int length = 2;
    typos = "bb";
    structure[0]=1;
    structure[1]=3;
    structure[2]=1;

    Network net(structure, typos, length);

    vector < vector <double> > x;
    vector < vector <double> > d;

    vector < vector <double> > x_testowy;
    vector < vector <double> > d_testowy;

    int stream_length = 4;

    double A = 0.0;

    int M=net.input_length;
    int N=net.output_length;

    /*
    for(int i=0;i<stream_length;i++)
    {

        vector <double> temp_x;
        vector <double> temp_d;
        for(int j=0;j<M;j++)
        {
            temp_x.push_back(A);
        }
        for(int j=0;j<N;j++)
        {
            temp_d.push_back(A);
        }
        x.push_back(temp_x);
        d.push_back(temp_d);
        A+=1.0;
    }*/

    for(int i=0;i<stream_length;i++)
    {
        vector <double> temp_x;
        vector <double> temp_d;
        for(int j=0;j<M;j++)
        {
            temp_x.push_back(A);
        }
        for(int j=0;j<N;j++)
        {
            temp_d.push_back(sin(A));
        }
        x.push_back(temp_x);
        d.push_back(temp_d);
        A+=2*3.14/stream_length;
    }

    for(int i=0;i<stream_length;i++)
    {
        vector <double> temp_x;
        vector <double> temp_d;

        temp_x=x[i];
        temp_d=d[i];

        cout<<"\n SIN("<<temp_x[0]<<")="<<temp_d[0];
    }
    A = 0.0;
    stream_length = 20;
    for(int i=0;i<stream_length;i++)
    {
        vector <double> temp_x;
        vector <double> temp_d;
        for(int j=0;j<M;j++)
        {
            temp_x.push_back(A);
        }
        for(int j=0;j<N;j++)
        {
            temp_d.push_back(fun(A));
        }
        x_testowy.push_back(temp_x);
        d_testowy.push_back(temp_d);
        A+=2*3.14/stream_length;
    }
    stream_length = 4;

    //cout<<"\nQ="<<gradient_descent(net, x, d, stream_length, 0.000000000000000005,0.8);
    cout<<"\nQ="<<gradient_descent_momentum(net, x, d, stream_length, 0.0000000005,0.2);

    stream_length = 20;
    std::cout.setf(std::ios_base::fixed,std::ios_base::floatfield);
    for(int i=0;i<stream_length;i++)
    {
            vector <double> test = x_testowy[i];
            vector <double> ETA  = d_testowy[i];
            net.process(test);
            cout<<net<<" dla PI/"<<i<<" oczekiwano "<< ETA[0];

    }

    cin.get();
    return 0;
}

double fun(double x)
{
    return sin(x);
}
