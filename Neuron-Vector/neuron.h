#ifndef _NEURON
#define _NEURON
#include <iostream>
#include <cstdlib>      // Dla rand()
#include <ctime>        // Dla time()
#include <cmath>        // Dla exp()
#include <vector>       // Dla Vector

using std::vector;

class Learning_stream
{
    public:
    double *  u;
    double *  y;
    double   ni; //Wspó³czynnik uczenia
    double  eps; //Warunek zatrzymania
    int       L; //Dlugosc ciagu
    int      lu;
    int      ly;

    Learning_stream ();
    Learning_stream (double * m_u, int length_m_u, double * m_y, int length_m_y, int m_L, double m_ni=0.5, double m_eps=0.5);
    Learning_stream (const Learning_stream & m_ls);
    ~Learning_stream();

    operator=(const Learning_stream & m_ls);
};

class Neuron
{
    public: // --- enkapsulacja robocza
    double  * wage; // Do modyfikacji dynamicznej
    double  wage_0;
    int          M;
    double       S;
    static int num;

    Neuron  ();
    Neuron  (int N); // --- wersja dynamiczna
    Neuron  (int N, double * m_wage); // --- jw.
    Neuron  (const Neuron & m_n);
    ~Neuron ();

    double process(double * u);

    Neuron & operator=(const Neuron & m_n);

    //friend double teach(l_stream s) --- funkcja uczaca
};

class Summer
{
    public: // --- enkapsulacja robocza
    double  *   wage;
    double      wage_0;
    double      S;
    int         M;

    Summer();
    Summer(int N);
    Summer(int N, double * m_wage);
    Summer(const Summer & m_s);
    ~Summer();

    double process(double * u);

    void operator=(const Summer & m_s);

    //friend double teah(l_stream s) --- funkcja uczaca
};

class Layer
{
    public: // --- enkapsulacja robocza
    vector<Neuron> neuron; // Do modyfjikacji dynamicznej
    int N;
    int M;
    double * output; // Robocze --- do usuniecia przy wersji dynamicznej

    Layer();
    Layer(int m_N,int m_M);
    Layer(const Layer & m_l);
    ~Layer();

    void process(double * u); // Do modyfikacji dynamicznej


    void operator=(const Layer & m_l);

    //friend double teach(l_stream s) --- funkcja uczaca
};

class S_Layer
{
    public:
    vector<Summer> neuron;
    int N;
    int M;
    double * output;

    S_Layer ();
    S_Layer (int m_N,int M);
    S_Layer (const S_Layer & m_l);
    ~S_Layer();

    void process(double * u); // Do modyfikacji dynamicznej


    void operator=(const S_Layer & m_l);

    //friend double teach(l_stream s) --- funkcja uczaca
};

class Network
{
    public: // --- enkapsulacja robocza
    vector<Layer>       layer;
    vector<S_Layer>     s_layer;
    int              L;
    int         y_size;
    double  *   output; // Robocze --- do usuniecia przy wersji dynamicznej

    Network ();
    Network (int m_L,int * N,int M);
    Network (const Network & m_n);
    ~Network();

    void show_result();
    void process(double * u); // Do modyfikacji dynamicznej double * proce  ss(double * u)

    operator=(const Network & m_n);

    //friend double teach(l_stream s) --- funkcja uczaca
};

double grad_desc(Network & net, const Learning_stream & ls);

double sigm_prim(double x);

double bipo_prim(double x);

double line_prim(double x);

struct pop
{
    vector<Network> x;
    vector<double>  s;
    int static      l;
    int static      m;
    double static   a;
    double static   tau;
    double static   tau_prim;
} ;

    int pop::l = 20;
    int pop::m = 10;
    double pop::a = 0.5;
    double pop::tau = 0.0;
    double pop::tau_prim = 0.0;

double rand_normal(double mean, double stddev);

double * lambda(int dimensions, int population);
#endif // _NEURON

