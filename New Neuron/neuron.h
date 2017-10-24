#ifndef _NEURON
#define _NEURON
#include <vector>
#include <iostream>
using std::vector;
using std::string;


/* Warstwa neuronów która dzi³a na zasadzie *
 * y_i=fun(s_i) gdzie s_i=sum(w_i_j*x_j)    *
 * a fun zale¿y do char type                */
class Layer
{
    public:
    vector <vector <double> > w; // macierz M * N gdzie M - iloœæ wejœæ N - iloœæ wyjœæ  j=0....M i=1....N
    vector <vector <double> > sigma;
    vector <double> s;           // Zsumowane wejœcia   i = 1...N
    vector <double> delta;       // Pochodna cz¹stkowa b³êdu po s i = 1...N
    vector <double> y;           // Wyjœcie i=1...N
    vector <double> x;           // Wejœcie j=1...M
    char            type;        // Typ neuronu (l - liniowy, u - sigmoidalny unipolarny, b - liniowy bipolarny
    int             N;
    int             M;

    // Wska¿nik na funkcjê kiedyœ
    // WskaŸnik na pochodn¹ kiedyœ

    Layer(int N, int M, char m_type);
    Layer(vector <vector <double> > & m_w, char m_type);
    ~Layer();

    void            sum();
    void            fun();            // funkcja  fun
    double          der(double m_x);            // pochodna fun
    void            change_wage(double eta); // w_ij(t+1)=w_ij(t)+2*eta*delta_i(t)*x_j(t)
    void            change_wage_momentum(double eta, Layer l);
    void            show_wage();
    void            show_delta();
    void            reset();
    vector <double> process(const vector <double> & m_x);
    vector <double> compute_delta(const vector <double> & epsilon);

    void operator>>(vector <double> & output);
};

class Network
{
    public:
    vector <Layer>  layer;
    vector <double> x;
    vector <double> y;

    int             input_length;
    int             output_length;
    int             L;

    Network(const vector <int> & structure,const string & structure_of_types, int length);
    ~Network();

    vector <double> process(const vector <double> & m_x);
    void            reset();

    friend std::ostream & operator<<(std::ostream & os,Network net);

};

std::ostream & operator<<(std::ostream & os,Network net);

double gradient_descent(Network & net, vector < vector <double> > x, const vector < vector <double> > d,int T, double error, double eta);

double gradient_descent_momentum(Network & net, vector < vector <double> > x, const vector < vector <double> > d,int T, double error, double eta);

double sigm(double x);

double bipo(double x);

double line(double x);

double sigm_prim(double temp);

double bipo_prim(double temp);

double line_prim(double x);

double lambda(Network & net, vector < vector <double> > x, const vector < vector <double> > d, int T);

double rand_normal(double mean, double stddev);

struct pop
{
    double *        x;
    double *        s;
    int static      l;
    int static      m;
    double static   a;
    double static   tau;
    double static   tau_prim;
} ;

#endif
