#include "neuron.h"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <conio.h>
#define LIMIT 10
#define BETA  0.5
#define ALFA  1
#define PARAM 0.9

using std::vector;

int pop::l = 20;
int pop::m = 10;
double pop::a = 0.5;
double pop::tau = 0.0;
double pop::tau_prim = 0.0;


Layer::Layer(int m_N, int m_M, char m_type)
{
    N=m_N;
    M=m_M;
    type = m_type;
    for(int i=0;i<N;i++)
    {
        vector <double> temp_w_row;
        vector <double> temp_sigma_row;
        for(int j=0;j<=M;j++) // j={1...M+1} poniewa¿ 1 to waga zerowa
        {
            //temp_w_row.push_back(((rand()%LIMIT)/10.0)*(-1*rand()%2));
            temp_w_row.push_back(0.0);
            temp_sigma_row.push_back(1.0);
        }
        w.push_back(temp_w_row);
        sigma.push_back(temp_sigma_row);
        s.push_back(0);
        delta.push_back(0);
        y.push_back(0);
    }
}

Layer::Layer(vector <vector <double> > & m_w, char m_type)
{
    w=m_w;
    type=m_type;
}

Layer::~Layer()
{

}

void Layer::sum()
{
    for(int i=0;i<N;i++)
    {
        s[i]=w[i][0];
        //std::cout<<"\nw0=="<<w[i][0];
        for(int j=1;j<=M;j++)
        {
            s[i]+=w[i][j]*x[j-1];
            //std::cout<<"\nw"<<j<<"=="<<w[i][j];
            //std::cout<<"\nx"<<j<<"=="<<x[j-1];
            //std::cout<<"\ns"<<j<<"=="<<s[i];
        }
        //std::cin.get();
    }
}

vector <double> Layer::process(const vector <double> & m_x)
{
    x=m_x;
    this -> sum();
    this -> fun();
    return y;
}

vector <double> Layer::compute_delta(const vector <double> & epsilon)
{
    for(int i=0;i<N;i++)
    {
        delta[i]=der(y[i])*epsilon[i];
        //std::cout<<"\nDelta("<<i<<")="<<der(y[i])<<"*"<<epsilon[i]<<"="<<delta[i];
    }
    return delta;
}

void Layer::show_delta()
{
    for(int i=0;i<N;i++)
    {
        std::cout<<"\nDelta Neuronu "<<i<<" = "<<delta[i];
    }
}

void Layer::fun()
{
        switch(type)
        {
            case 'l':
                for(int i=0;i<N;i++)
                {
                    y[i]=line(s[i]);
                }
                break;
            case 'u':
                for(int i=0;i<N;i++)
                {
                    y[i]=sigm(s[i]);
                }
                break;
            case 'b':
                for(int i=0;i<N;i++)
                {
                    y[i]=bipo(s[i]);
                }
                break;
            default:
                exit(EXIT_FAILURE);
        }
}

double Layer::der(double m_x)
{
        switch(type)
        {
            case 'l':
                return line_prim(m_x);
                break;
            case 'u':
                return sigm_prim(m_x);
                break;
            case 'b':
                return bipo_prim(m_x);
                break;
            default:
                exit(EXIT_FAILURE);
                return 0;
        }
}

void Layer::change_wage(double eta)
{
    for(int i=0;i<N;i++)
    {
        //std::cout<<"\nw0 przed zmianą"<<w[i][0];
        w[i][0]     += 2*eta*delta[i];
        //std::cout<<"\nw0 po zmianie 2*"<<eta<<"*"<<delta[i]<<"="<<w[i][0];
        for(int j=1;j<=M;j++)
        {
            //std::cout<<"\nw("<<j<<") przed zmianą "<<w[i][j];
            w[i][j] += 2*eta*delta[i]*x[j-1];
            //std::cout<<"\nw("<<j<<") po zmianie 2*"<<eta<<"*"<<delta[i]<<"*"<<x[j-1]<<"="<<w[i][j];
        }
    }
}

void Layer::change_wage_momentum(double eta, Layer l)
{
    for(int i=0;i<N;i++)
    {
        //std::cout<<"\nw0 przed zmianą"<<w[i][0];
        w[i][0]     += 2*eta*delta[i] + PARAM*(w[i][0]-l.w[i][0]);
        //std::cout<<"\nw0 po zmianie 2*"<<eta<<"*"<<delta[i]<<"="<<w[i][0];
        for(int j=1;j<=M;j++)
        {
            //std::cout<<"\nw("<<j<<") przed zmianą "<<w[i][j];
            w[i][j] += 2*eta*delta[i]*x[j-1] + PARAM*(w[i][j]-l.w[i][j]);
            //std::cout<<"\nw("<<j<<") po zmianie 2*"<<eta<<"*"<<delta[i]<<"*"<<x[j-1]<<"="<<w[i][j];
        }
    }
}
void Layer::show_wage()
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<=M;j++)
        {
            std::cout<<"\nWaga "<<j<<",Neuron "<<i<<" = "<<w[i][j];
        }
    }
}

void Layer::reset()
{
    //N=m_N;
    //M=m_M;
    //type = m_type;

    vector <vector <double> > n;
    w     = n;
    sigma = n;
    vector <double> q;
    s     = q;
    delta = q;
    y     = q;

    for(int i=0;i<N;i++)
    {
        vector <double> temp_w_row;
        vector <double> temp_sigma_row;
        for(int j=0;j<=M;j++) // j={1...M+1} poniewa¿ 1 to waga zerowa
        {
            //temp_w_row.push_back(((rand()%LIMIT)/10.0)*(-1*rand()%2));
            temp_w_row.push_back(0.0);
            temp_sigma_row.push_back(1.0);
        }
        w.push_back(temp_w_row);
        sigma.push_back(temp_sigma_row);
        s.push_back(0);
        delta.push_back(0);
        y.push_back(0);
    }
}

/* sigma(S) */
double sigm(double x)
{
    double temp = (1/(1+exp(-BETA*x)));
    if(temp==INFINITY)
        temp = 1.0;
    if(temp == NAN || temp == -NAN || temp== -INFINITY)
        temp = 0.0;
    return temp;
}

/* tgh(S)   */
double bipo(double x)
{
    //double temp = ((1-exp(BETA*x))/(1+exp(-BETA*x)));
    //return temp;
    return tanh(x);
}

/* (a*s)    */
double line(double x)
{
    return ALFA*x;
}

/* sigma'(S) */
double sigm_prim(double temp)
{
    temp        = (BETA*temp*(1-temp));
    return temp;
}

/* tgh'(S)   */
double bipo_prim(double temp)
{
    return (BETA*(1-(temp*temp)));
}

/* (a*s)'    */
double line_prim(double x)
{
    return ALFA;
}

void Layer::operator>>(vector <double> & output)
{
    output=y;
}

/*********************************************************/
/*********************NETWORK*****************************/
/*********************************************************/
/*********************************************************
    const vector <int> & structure:
        - Pierwszy element to d³ugoœæ wejœcia
        - Drugi element to iloœc neuronów w warstwie
        - Ostatni element L+1 to iloœc neuronów w ostatniej
          warstwie == d³ugoœæ wyjœcia
        Struktura powinna mieæ L+1 elementów
    const string & structure_of_types
        - Ci¹g o d³ugoœci L
        - i-ty znak oznacza typ warstwy
    int length
        - równoznaczny z int L tj. iloœæ warstw
 *********************************************************/
Network::Network(const vector <int> & structure, const string & structure_of_types, int length)
{
    L=length;

    for(int i=0;i<L;i++)
    {
        layer.push_back(Layer(structure[i+1],structure[i],structure_of_types[i]));
    }

    input_length=layer[0].M;
    output_length=layer[L-1].N;
}

Network::~Network()
{

}

void Network::reset()
{
    for(int i=0;i<L;i++)
    {
        layer[i].reset();
    }
}

vector <double> Network::process(const vector <double> & m_x)
{
    x=m_x;
    vector <double> temp;
    temp=x;
    for(int i=0;i<L;i++)
    {
        layer[i].process(temp);
        layer[i]>>temp;
    }
    y=temp;
    return y;
}

std::ostream & operator<<(std::ostream & os,Network net)
{
    int N = net.output_length;
    for(int i=0;i<N;i++)
        os<<"\nWyjscie["<<i<<"]="<<net.y[i];
    return os;
}

double gradient_descent(Network & net, vector < vector <double> > x, const vector < vector <double> > d,int T, double error, double eta)
{
    std::cout.setf(std::ios_base::scientific,std::ios_base::floatfield);
    int                 N = net.output_length;
    int                 L = net.L;
    double              Q = 0;
    double              q = 0;
    double              pre_Q;
    int                 M;
    unsigned long int   iter = 0;
    Network             temp_net = net;

    for(int t=0;t<T;t++)
    {
        net.process(x[t]);
        for(int i=0;i<N;i++)
        {
            //std::cout << "\ny("<<i<<")="<<net.y[i];
            //std::cout << "\nd("<<i<<")="<<d[t][i];
            //std::cout << "\ny-d="<<(net.y[i]-d[t][i]);
            q  = (d[t][i]-net.y[i]);
            //std::cout << "\n(y-d)^2="<<q*q;
            Q += q*q;
            //std::cin.get();
        }
    }



    std::cout<<"\nQ="<<Q;
    double err=error;
    error*=eta;
    pre_Q=2*(Q+10000000000000000*error);
    while(error<sqrt((Q-pre_Q)*(Q-pre_Q)))
    {
        temp_net = net;
        pre_Q=Q;
        Q=0;
        for(int t=0;t<T;t++)
        {
            vector <double> epsilon;
            net.process(x[t]);

            for(int i=0;i<N;i++)
            {
                //std::cout << "\ny("<<i<<")="<<net.y[i];
                //std::cout << "\nd("<<i<<")="<<d[t][i];
                //std::cout << "\ny-d="<<(d[t][i]-net.y[i]);
                q  = (d[t][i]-net.y[i]);
                epsilon.push_back(q);
            }

            net.layer[L-1].compute_delta(epsilon);

            for(int k=L-2;k>=0;k--)
            {
                vector <double> epsilon;
                N=net.layer[k].N;
                M=net.layer[k+1].N;

                double eps;

                for(int i=0;i<N;i++)
                {
                    eps = 0;
                    for(int m=0;m<M;m++)
                    {
                        //std::cout<<"\nk="<<k+1<<" delta("<<m<<")="<<net.layer[k+1].delta[m]<<" w("<<m<<")("<<i<<")="<<net.layer[k+1].w[m][i];
                        eps+= net.layer[k+1].delta[m] * net.layer[k+1].w[m][i];
                    }
                    epsilon.push_back(eps);
                }

                net.layer[k].compute_delta(epsilon);

            }

            for(int k=0;k<L;k++)
            {
                /*std::cout<<"\nWarstwa "<<k<<" przed zmiana";
                net.layer[k].show_wage();
                net.layer[k].show_delta();*/
                net.layer[k].change_wage(eta);
                /*std::cout<<"\nWarstwa "<<k<<" po zmianie";
                net.layer[k].show_wage();
                net.layer[k].show_delta();
                std::cin.get();*/
            }



            for(int i=0;i<N;i++)
            {
                q  = (net.y[i]-d[t][i]);
                Q += q*q;
            }
        }

       /* if(Q>pre_Q)
        {
            //net=temp_net;
            eta/=2.0;
            error=err*eta;
            std::cout<<"\neta="<<eta;
        }else{

            //eta+=eta*0.1;
        }*/
        if(iter%100000==0)
            std::cout<<"\nQ="<<Q;
        if(iter%1000000==0)
            std::cout<<"\nSRE="<<sqrt((Q-pre_Q)*(Q-pre_Q));
        iter++;
        if(kbhit())
            break;
    }

    std::cout<<"\nZakończono uczenie warunkiem "<<error<<">"<<sqrt((Q-pre_Q)*(Q-pre_Q))<<" w "<<iter<<" iteracji.";
    return Q;
}

double gradient_descent_momentum(Network & net, vector < vector <double> > x, const vector < vector <double> > d,int T, double error, double eta)
{
    std::cout.setf(std::ios_base::scientific,std::ios_base::floatfield);
    int                 N = net.output_length;
    int                 L = net.L;
    double              Q = 0;
    double              q = 0;
    double              pre_Q;
    int                 M;
    unsigned long int   iter = 0;
    Network             temp_net = net;
    double              small_Q =0;

    for(int t=0;t<T;t++)
    {
        net.process(x[t]);
        for(int i=0;i<N;i++)
        {
            //std::cout << "\ny("<<i<<")="<<net.y[i];
            //std::cout << "\nd("<<i<<")="<<d[t][i];
            //std::cout << "\ny-d="<<(net.y[i]-d[t][i]);
            q  = (d[t][i]-net.y[i]);
            //std::cout << "\n(y-d)^2="<<q*q;
            Q += q*q;
            //std::cin.get();
        }
    }



    std::cout<<"\nQ="<<Q;
    double err      = error;
    error          *= eta;
    pre_Q           = 2*(Q+10000000000000000*error);
    while(error<sqrt((Q-pre_Q)*(Q-pre_Q)))
    {
        pre_Q   = Q;
        Q       = 0;


        for(int t=0;t<T;t++)
        {

                for(int i=0;i<N;i++)
                {
                    q  = (net.y[i]-d[t][i]);
                    small_Q += q*q;
                }

                temp_net = net;
                small_Q     = 0;
                vector <double> epsilon;
                net.process(x[t]);
                N = net.output_length;
                for(int i=0;i<N;i++)
                {
                    /*std::cout << "\nN="<<N;
                    std::cout << "\nT="<<t;
                    std::cout << "\ny("<<i<<")="<<net.y[i];
                    std::cout << "\nd("<<i<<")="<<d[t][i];
                    std::cout << "\ny-d="<<(d[t][i]-net.y[i]);
                    q  = (d[t][i]-net.y[i]);*/
                    epsilon.push_back(q);
                }

                net.layer[L-1].compute_delta(epsilon);

                for(int k=L-2;k>=0;k--)
                {
                    vector <double> epsilon;
                    N=net.layer[k].N;
                    M=net.layer[k+1].N;

                    double eps;

                    for(int i=0;i<N;i++)
                    {
                        eps = 0;
                        for(int m=0;m<M;m++)
                        {
                            //std::cout<<"\nk="<<k+1<<" delta("<<m<<")="<<net.layer[k+1].delta[m]<<" w("<<m<<")("<<i<<")="<<net.layer[k+1].w[m][i];
                            eps+= net.layer[k+1].delta[m] * net.layer[k+1].w[m][i];
                        }
                        epsilon.push_back(eps);
                    }

                    net.layer[k].compute_delta(epsilon);

                }

                for(int k=0;k<L;k++)
                {
                    /*std::cout<<"\nWarstwa "<<k<<" przed zmiana";
                    net.layer[k].show_wage();
                    net.layer[k].show_delta();*/
                    //net.layer[k].change_wage_momentum(eta,net.layer[k]);
                    net.layer[k].change_wage(eta);
                    /*std::cout<<"\nWarstwa "<<k<<" po zmianie";
                    net.layer[k].show_wage();
                    net.layer[k].show_delta();
                    std::cin.get();*/
                }
                //std::cin.get();




            for(int i=0;i<N;i++)
            {
                q  = (net.y[i]-d[t][i]);
                Q += q*q;
            }
            //temp_net = net;
        }

        /*if(Q>pre_Q)
        {
            //net=temp_net;
            eta/=2.0;
            error=err*eta;
            std::cout<<"\neta="<<eta;
        }else{

            //eta+=eta*0.1;
            //error=err*eta;
        }*/
        if(iter%10==0)
            std::cout<<"\nQ="<<Q;
        if(iter%100==0)
            std::cout<<"\nSRE="<<sqrt((Q-pre_Q)*(Q-pre_Q));
        iter++;
        if(kbhit())
            break;
    }

    std::cout<<"\nZakończono uczenie warunkiem "<<error<<">"<<sqrt((Q-pre_Q)*(Q-pre_Q))<<" w "<<iter<<" iteracji.";
    return Q;
}


double rand_normal(double mean, double stddev)
{//Box muller method
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        double x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;

            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            double d        = sqrt(-2.0*log(r)/r);
            double n1       = x*d;
            n2              = y*d;
            double result   = n1*stddev + mean;
            n2_cached       = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}

double lambda(Network & net, vector < vector <double> > x, const vector < vector <double> > d, int T)
{

    int      iter           = 0;

    vector<Network> population;

    std::srand(std::time(0));

    for(int i=0;i<pop::m;i++)
    {
        population.push_back(net);
        population[i].reset();
    }

    vector <double>  score;
    vector <Network> temp_population;
    vector <double>  temp_score;
    double       change             = 1.0;
    double       best_score         = 0.0;

    for(int i=0;i<pop::m;i++)
    {
        best_score += (net.process(x[i]) - d[i]) * (net.process(x[i]) - d[i]);
    }

    for(int i=0;i<pop::l;i++)
    {
        temp_population.push_back(net);
    }

    pop::tau        = 1 / sqrt(2.0*sqrt(dimensions));
    pop::tau_prim   = 1 / sqrt(2.0*dimensions);

    /* Tutaj zaczyna siê WHILE  */
    while(change>0.0)
    {
        /*  Ocena populacji P za pomoc¹ funkcji fun(double * ) */
        pt_pop     = population;
        pt_score   = score;

        for(int i=0;i<pop::m;i++)
        {
            *pt_score = fun(pt_pop -> x);
            pt_pop++;
            pt_score++;
        }

        /*  Tworzenie nowej populacji O za pomoc¹ krzy¿owania i mutacji populacji P */
        /*  Najpierw krzyżowanie   */
        pt_temp_population  = temp_population;
        pt_x                = pt_temp_population -> x;
        pt_s                = pt_temp_population -> s;

        for(int i=0;i<pop::l;i++)
        {
            for(int j=0;j<dimensions;j++)
            {
                *pt_x = pop::a * *((population+rand()%pop::m)->x) + (1 - pop::a) * *((population+rand()%pop::m)->x);
                *pt_s = pop::a * *((population+rand()%pop::m)->s) + (1 - pop::a) * *((population+rand()%pop::m)->s);
                pt_x++;
                pt_s++;
            }
            pt_temp_population++;
            pt_x    = pt_pop -> x;
            pt_s    = pt_pop -> s;
        }

        /* Potem mutacja */
        pt_temp_population  = temp_population;
        pt_x                = pt_temp_population -> x;
        pt_s                = pt_temp_population -> s;
        for(int i=0;i<pop::l;i++)
        {
            /*Najpierw mutacja sigma*/
            /*Potem mutacja x*/
            for(int j=0;j<dimensions;j++)
            {
                *pt_s *= exp(pop::tau_prim * rand_normal(0,1) + pop::tau * rand_normal(0,1));
                *pt_x += *pt_s * rand_normal(0,1);
                pt_x++;
                pt_s++;
            }
            pt_temp_population++;
            pt_x = pt_temp_population -> x;
            pt_s = pt_temp_population -> s;
        }

        /*  Ocena populacji O za pomoc¹ funkcji fun(double * )  */
        pt_temp_score       = temp_score;
        pt_temp_population  = temp_population;
        for(int i=0;i<pop::l;i++)
        {
            *pt_temp_score = fun(pt_temp_population->x);
            if(*pt_temp_score<best_score)
            {
                std::cout<<"\nNowy najlepszy wynik "<<*pt_temp_score;
                best_score = *pt_temp_score;
            }

            pt_temp_population++;
            pt_temp_score++;
        }

        /*  Wybór najlepszych z populacji O                 */

        bool temp = true;

        while(temp)
        {
            temp               = false;
            pt_temp_population = temp_population + 1;
            pt_temp_score      = temp_score + 1;

            for(int i=1;i<pop::l;i++)
            {
                if(*(pt_temp_score-1)<*pt_temp_score)
                {
                    temp                = true;
                    double t_score      = *pt_temp_score;
                    *pt_temp_score      = *(pt_temp_score-1);
                    *(pt_temp_score-1)  = t_score;

                    double t_x;
                    double t_s;

                    pt_x = pt_temp_population -> x;
                    pt_s = pt_temp_population -> s;

                    double * temp_pt_x = (pt_temp_population-1) -> x;
                    double * temp_pt_s = (pt_temp_population-1) -> s;

                    for(int j=0;j<dimensions;j++)
                    {
                        t_x         = *pt_x;
                        *pt_x       = *temp_pt_x;
                        *temp_pt_x  = t_x;


                        t_s         = *pt_s;
                        *pt_s       = *temp_pt_s;
                        *temp_pt_s  = t_s;

                        pt_x++;
                        pt_s++;
                        temp_pt_s++;
                        temp_pt_x++;
                    }

                }
                pt_temp_score++;
                pt_temp_population++;
            }
        }

        pt_pop              = population;
        pt_temp_population  = temp_population;

        for(int i=0;i<pop::m;i++) // Do modyfikacji
        {
            pt_x                = pt_pop -> x;
            pt_s                = pt_pop -> s;

            double * temp_pt_x  = pt_temp_population -> x;
            double * temp_pt_s  = pt_temp_population -> s;

            for(int j=0;j<dimensions;j++)
            {
                *pt_x = *temp_pt_x;
                *pt_s = *temp_pt_s;
                pt_x++;
                pt_s++;
                temp_pt_x++;
                temp_pt_s++;
            }

            pt_pop++;
            pt_temp_population++;
        }

        iter++;

        if(kbhit())
            break;
    }
    /* Koniec While     */
    /*  Ocena populacji P za pomoc¹ funkcji fun(double * ) */
    pt_pop     = population;
    pt_score   = score;

    for(int i=0;i<pop::m;i++)
    {
        *pt_score = fun(pt_pop->x);
        pt_pop++;
        pt_score++;
    }

    /* Wybór najlepszego wektora
       Zwrócenie go              */
    best_score = *score;
    pt_score   = score;
    int k = 0;
    pt_score++;
    for(int i=1;i<pop::m;i++)
    {
        if(*pt_score<best_score)
        {
            best_score = *pt_score;
            k=i;
        }
        pt_score++;
    }

    pt_pop = population + k;
    pt_x   = pt_pop->x;
    double * x = new double[dimensions];
    std::cout<<"\nZwrócono wektor (";
    for(int i=0;i<dimensions;i++)
    {
        *x = *pt_x;
        std::cout<<*x<<",";
        pt_x++;
    }
    std::cout<<") o wyniku "<<best_score<<" w "<<iter<<" iteracjach.";
    return x;
}
