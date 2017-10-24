#include "neuron.h"
#include <conio.h>
#define LIMIT 1
#define BETA  2.0
#define ALFA  1.0
/********************************************************
    Implementacja klasy Neuron
    double * wage - M wag neuronów
    int M         - ilosc wejsc
    double wage_0 - waga zerowa
 ********************************************************/

/********************************************************
    Konstruktor domyœlny
    Brak praktycznego zastosowania,istnieje jedynie ¿eby
    byæ kompatybilnym z destruktorem.
 ********************************************************/

int Neuron::num=0;

Neuron::Neuron()
{
    wage    = new double[1];
    wage_0  = 0;
    M       = 1;
    num++;
}

/********************************************************
    Konstruktor wag losowych
    Tworzy Neuron gdzie argument int N jest iloœci¹ wejœæ
    danego neuronu
 ********************************************************/
Neuron::Neuron(int N)
{
    wage          = new double[N];
    double * temp = wage;
    num++;

    for(int i=0;i<N;i++)
    {
        *temp   =  std::rand()%LIMIT;
        //*temp   /= (std::rand()%LIMIT)+1.0;
        //*temp   *= (-(std::rand()%LIMIT));
        temp++;
    }
    wage_0      =  std::rand()%LIMIT;
    //wage_0      /= (std::rand()%LIMIT)+1.0;
    //wage_0      *= (-(std::rand()%LIMIT));
    M           =  N;
}

/*********************************************************
    Konstruktor wag dobranych
    Tworzy Neuron gdzie argument int N jest iloœci¹ wejœæ
    danego Neuronu, natomiast double * m_wage jest wskaŸ-
    nikiem wskazuj¹cym na adres pierwszego elementu adresu
    tablicy dobranych wag których musi byæ (N+1)
 *********************************************************/
Neuron::Neuron(int N, double * m_wage)
{
    wage            =new double[N];
    double * temp   = wage;
    double * m_temp = m_wage;

    for(int i=0;i<N;i++)
    {
        *temp = *m_temp;
        temp++;
        m_temp++;
    }

    wage_0  = *m_temp;
    M       = N;
    num++;
}

/*********************************************************
    Konstruktor kopiuj¹cy
    Istnieje z powodu dynamicznego przydzielania pamiêci
    zmiennej wage
 *********************************************************/
Neuron::Neuron(const Neuron & m_n)
{
    wage_0=m_n.wage_0;
    M=m_n.M;
    wage            = new double[M];

    double * temp   = wage;
    double * m_temp = m_n.wage;

    for(int i=0;i<m_n.M;i++)
    {
        *temp = *m_temp;
        temp++;
        m_temp++;
    }
    num++;
}
/********************************************************
    Destruktor
    Usuwa pamiêæ przydzielon¹ zmiennej wage
 ********************************************************/
Neuron::~Neuron()
{
    delete [] wage;
    num--;
}

/********************************************************
    Metoda process
    Pobiera jako arguemnt N wejœæ które s¹ mno¿one przez
    wagi i sumowane wraz z waga zerow¹, nastêpnie ten
    wynik jest przetwarzany przez funkcjê sigmoidaln¹ i
    zwracany
 ********************************************************/
double Neuron::process(double * u)
{
    double s                = wage_0;
    double * temp           = wage;
    const double * m_temp   = u;

    using std::cout;
    //cout<<"\nNEURON w0  ="<<wage_0;

    for(int i=0;i<M;i++)
    {
        //cout <<"\nw="<<*temp<<"*u="<<*m_temp;
        s   += *temp*(*m_temp);
        temp++;
        m_temp++;
    }

    S   = s;

    double x = (1/(1+exp(-BETA*s)));

    if(x==INFINITY)
        x==1;
    if(x==NAN || x==-NAN || x==-INFINITY)
        x==0;

    return x;
}

/********************************************************
    operator= (operator przypisania)
    Istnieje z powodu dynamicznego przydzielania pamiêci
    zmiennej wage
 ********************************************************/
Neuron & Neuron::operator=(const Neuron & m_n)
{
    if(this == & m_n)
        return *this;
    delete [] wage;
    wage            = new double[m_n.M];

    double * temp   = wage;
    double * m_temp = m_n.wage;

    for(int i=0;i<m_n.M;i++)
    {
        *temp       = *m_temp;
        temp++;
        m_temp++;
    }

    wage_0          = m_n.wage_0;
    M               = m_n.M;
    return *this;
}

/********************************************************
    Implementacja klasy Summer
    double * wage - M wag neuronów
    int M         - ilosc wejsc
    double wage_0 - waga zerowa
    Wszystkie metody poni¿ej s¹ analogiczne do metod
    klasy Neuron za wyj¹tkiem metody process
  ********************************************************/
/********************************************************
    Konstruktor domyœlny
    Brak praktycznego zastosowania,istnieje jedynie ¿eby
    byæ kompatybilnym z destruktorem.
 ********************************************************/
Summer::Summer()
{
    wage    = new double[1];
    wage_0  = 0;
    M       = 1;
}

/********************************************************
    Konstruktor wag losowych
    Tworzy Neuron gdzie argument int N jest iloœci¹ wejœæ
    danego neuronu
 ********************************************************/
Summer::Summer(int N)
{
    wage            = new double[N];
    double * temp   = wage;

    for(int i=0;i<N;i++)
    {
        *temp   =  std::rand()%LIMIT;
        //*temp   /= (std::rand()%LIMIT)+1.0;
        //*temp   *= (-(std::rand()%LIMIT));
        temp++;
    }

    wage_0      =  std::rand()%LIMIT;
    //wage_0      /= (std::rand()%LIMIT)+1.0;
    //wage_0      *= (-(std::rand()%LIMIT));
    M           =  N;
}

/*********************************************************
    Konstruktor wag dobranych
    Tworzy Neuron gdzie argument int N jest iloœci¹ wejœæ
    danego Neuronu, natomiast double * m_wage jest wskaŸ-
    nikiem wskazuj¹cym na adres pierwszego elementu adresu
    tablicy dobranych wag których musi byæ (N+1) gdzie
    ostatni element tablicy jest wag¹ zerow¹
 *********************************************************/
Summer::Summer(int N, double * m_wage)
{
    wage=new double[N];
    double * temp   = wage;
    double * m_temp = m_wage;

    for(int i=0;i<N;i++)
    {
        *temp = *m_temp;
        temp++;
        m_temp++;
    }

    wage_0  = *m_temp;
    M       = N;
}

/*********************************************************
    Konstruktor kopiuj¹cy
    Istnieje z powodu dynamicznego przydzielania pamiêci
    zmiennej wage
 *********************************************************/
Summer::Summer(const Summer & m_n)
{
    wage            = new double[m_n.M];

    double * temp   = wage;
    double * m_temp = m_n.wage;

    for(int i=0;i<m_n.M;i++)
    {
        *temp = *m_temp;
        temp++;
        m_temp++;
    }

    wage_0  = m_n.wage_0;
    M       = m_n.M;
}

/********************************************************
    Destruktor
    Usuwa pamiêæ przydzielon¹ zmiennej wage
 ********************************************************/
Summer::~Summer()
{
    delete [] wage;
}

/********************************************************
    Metoda process
    Pobiera jako arguemnt N wejœæ które s¹ mno¿one przez
    wagi i sumowane wraz z waga zerow¹, nastêpnie ten
    wynik jest przetwarzany przez funkcjê liniow¹ i
    zwracany
 ********************************************************/
double Summer::process(double * u)
{
    double         s        = wage_0;
    double       * temp     = wage;
    const double * m_temp   = u;

    using std::cout;
    //cout<<"\nSUMMER w_0 ="<<s;
    for(int i=0;i<M;i++)
    {
        //cout<<"\nW="<<*temp<<"*x="<<*m_temp;
        s   += *temp*(*m_temp);
        temp++;
        m_temp++;
    }

    S = s;

    return (ALFA*s);
}

/********************************************************
    operator= (operator przypisania)
    Istnieje z powodu dynamicznego przydzielania pamiêci
    zmiennej wage
 ********************************************************/
void Summer::operator=(const Summer & m_n)
{
    delete [] wage;
    wage            = new double[m_n.M];

    double * temp   = wage;
    double * m_temp = m_n.wage;

    for(int i=0;i<m_n.M;i++)
    {
        *temp = *m_temp;
        temp++;
        m_temp++;
    }

    wage_0  = m_n.wage_0;
    M       = m_n.M;
}

/*********************************************************
    Implementacja klasy Layer
    Neuron * neuron - WskaŸnik na pierwszy element tablicy
                      neuronów
    int N           - Iloœæ neuronów w wartswie
    int M           - Iloœæ wejœæ
    double * output - Wska¿nik na pierwszy element tablicy
                      wyjœæ neuronów
 *********************************************************/

/*********************************************************
    Konstruktor domyœlny
    Brak praktycznego zastosowania,istnieje jedynie ¿eby
    byæ kompatybilnym z destruktorem.
 *********************************************************/
Layer::Layer()
{
    output = new double[1];
}

/*********************************************************
    Konstruktor standardowy
    Tworzy warstwê m_N neuronów które posiadaj¹ po m_M
    synaps
 *********************************************************/
Layer::Layer(int m_N, int m_M)
{
    for(int i=0;i<m_N;i++)
    {
        neuron.push_back(Neuron(m_M));
    }
    N       = m_N;
    M       = m_M;
    output  = new double[m_N];
}

/********************************************************
    Konstruktor kopiuj¹cy
    Istnieje z powodu dynamicznego przydzielania pamiêci
 ********************************************************/
Layer::Layer(const Layer & m_l)
{
    N       = m_l.N;
    M       = m_l.M;
    for(int i=0;i<N;i++)
    {
        neuron.push_back(m_l.neuron[i]);
    }

    output  = new double[N];
}

/********************************************************
    Destruktor
    Usuwa pamiêæ tablic neuron oraz output
 ********************************************************/
Layer::~Layer()
{
    delete [] output;
}

/********************************************************
    Metoda process
    Przesy³a do wszystkich neuronów warstwy wskaŸnik
    pierwszezgo elementu tablic u a nastêpnie z ich wyjœæ
    tworzy tablicê output
 ********************************************************/
void Layer::process(double * u)
{
    double * m_temp = u;
    double * temp   = output;

    for(int i=0;i<N;i++)
    {
        *temp       = neuron[i].process(m_temp);
        temp++;
    }
}

/********************************************************
    operator= (operator przypisania)
    Istnieje z powodu dynamicznego przydzielania pamiêci
    tablicy output oraz tablicy obiektów klasy neuron
 ********************************************************/
void Layer::operator=(const Layer & m_l)
{
    delete [] output;
    for(int i=0;i<m_l.N;i++)
        neuron[i]=m_l.neuron[i];

    N      = m_l.N;
    M      = m_l.M;

    output = new double[N];
}

/*********************************************************
    Implementacja klasy S_Layer
    Analogiczna do klasy Layer jedna wykorzystuje
    wektro obiektów klasy Summer miast Neuron
    Neuron * neuron - WskaŸnik na pierwszy element tablicy
                      neuronów
    int N           - Iloœæ neuronów w wartswie
    int M           - Iloœæ wejœæ
    double * output - Wska¿nik na pierwszy element tablicy
                      wyjœæ neuronów
 *********************************************************/

/*********************************************************
    Konstruktor domyœlny
    Brak praktycznego zastosowania,istnieje jedynie ¿eby
    byæ kompatybilnym z destruktorem.
 *********************************************************/
S_Layer::S_Layer()
{
    output = new double[1];
}

/*********************************************************
    Konstruktor standardowy
    Tworzy warstwê m_N neuronów które posiadaj¹ po m_M
    synaps
 *********************************************************/
S_Layer::S_Layer(int m_N, int m_M)
{
    for(int i=0;i<m_N;i++)
    {
        neuron.push_back(Summer(m_M));
    }
    N      = m_N;
    M      = m_M;
    output = new double[m_N];
}

/********************************************************
    Konstruktor kopiuj¹cy
    Istnieje z powodu dynamicznego przydzielania pamiêci
 ********************************************************/
S_Layer::S_Layer(const S_Layer & m_l)
{
    for(int i=0;i<m_l.N;i++)
    {
        neuron.push_back(m_l.neuron[i]);
    }

    N = m_l.N;
    M = m_l.M;

    output = new double[N];
}

/********************************************************
    Destruktor
    Usuwa pamiêæ tablic neuron oraz output
 ********************************************************/
S_Layer::~S_Layer()
{
    delete [] output;
}

/********************************************************
    Metoda process
    Przesy³a do wszystkich neuronów warstwy wskaŸnik
    pierwszezgo elementu tablic u a nastêpnie z ich wyjœæ
    tworzy tablicê output
 ********************************************************/
void S_Layer::process(double * u)
{
    double * m_temp = u;
    double * temp   = output;

    for(int i=0;i<N;i++)
    {
        *temp = neuron[i].process(m_temp);
        temp++;
    }
}

/********************************************************
    operator= (operator przypisania)
    Istnieje z powodu dynamicznego przydzielania pamiêci
    tablicy output oraz tablicy obiektów klasy neuron
 ********************************************************/
void S_Layer::operator=(const S_Layer & m_l)
{
    delete [] output;
    N = m_l.N;
    M = m_l.M;

    for(int i=0;i<N;i++)
        neuron.push_back(m_l.neuron[i]);

    output = new double[N];
}

/***********************************************************
    Klasa Network
    Składa się z L warstw neuronów zawartych w tablicy layer
    warstwy sumującej s_layer oraz tablicy wyjściowej output
 ***********************************************************/

/********************************************************
    Konstruktor domyślny
    Brak praktycznego zastosowania,istnieje jedynie żeby
    być kompatybilnym z destruktorem.
 ********************************************************/
Network::Network()
{
    output = new double[1];
}

/********************************************************
    Konstruktor standardowy
    Tworzy m_L warstw neuronów
    Jest M wejść dla sieci, ilość neuronów dla każdej
    warstwy jest opisana w tablicy N, przekazywanej jako
    wskaźnik na adres pierwszego elementu tablicy
 ********************************************************/
Network::Network(int m_L, int * N, int M)
{
    std::srand(std::time(0));
    int * temp = N;

    int temp_M = M;
    for(int i=1;i<m_L;i++)
    {
        layer.push_back(Layer(*temp,temp_M));
        temp_M = *temp;
        temp++;
    }
    s_layer.push_back(S_Layer(*temp,temp_M));
    L      = m_L;
    y_size = *temp;
    output = new double[y_size];
}

/********************************************************
    Destruktor
    Usuwa tablicy layer oraz output, a także zmienną
    s_layer
 ********************************************************/
Network::~Network()
{
    delete [] output;
}

/********************************************************
    Metoda process
    przesyła wejście u otrzymywane argumentem jako
    wskaźnik na adres pierwszego elementu tablicy u o
    wielkośi M dla warstw
 ********************************************************/
void Network::process(double * u)
{
    double * temp       = u;
    L--;
    for(int i=0;i<L;i++)
    {
        layer[i].process(temp);
        temp            = layer[i].output;
    }

    s_layer[0].process(temp);
    int j               = s_layer[0].N;
    temp                = s_layer[0].output;
    double * temp_out   = output;

    for(int i=0;i<j;i++)
    {
        *temp_out = *temp;
        temp_out++;
        temp++;
    }
    L++;
}

void Network::show_result()
{
    double * temp = output;

    for(int i=0;i<y_size;i++)
    {
        std::cout<<"\nY["<<i<<"]="<<*temp;
        temp++;
    }
}
/***************************************************************
    Implementacja klasy Learning_stream
    - klasa ciągu uczącego
    - wskaźnik double * u wskazuje na tablicę wejść
    - wskażnik double * y wskazuje na tablicę oczekiwanych wyjść
    - double ni to współczynnik uczenia
    - double eps to maksymalny błąd
    - int L to długość ciagu
 ***************************************************************/

/***************************************************************
    Konstruktor domyślny Learning_stream()
 ***************************************************************/
 Learning_stream::Learning_stream()
 {
     u   = new double[1];
     y   = new double[1];
     ni  = 1.5;
     eps = 0.5;
     L   =   1;
 }

/***************************************************************
    Konstruktor standardowy
 ***************************************************************/
 Learning_stream::Learning_stream(double * m_u, int length_m_u, double * m_y, int length_m_y, int m_L, double m_ni, double m_eps)
 {
    u               = new double[length_m_u];
    double * temp   = m_u;
    double * pt     =   u;

    for(int i=0;i<length_m_u;i++)
    {
        *pt=*temp;
        temp++;
        pt++;
    }
    y               = new double[length_m_y];
    temp            = m_y;
    pt              = y;
    for(int i=0;i<length_m_y;i++)
    {
        *pt=*temp;
        temp++;
        pt++;
    }
    L=m_L;
    ni=m_ni;
    eps=m_eps;
    lu=length_m_u;
    ly=length_m_y;
 }
/***************************************************************
    Konstruktor kopiujący
 ***************************************************************/
 Learning_stream::Learning_stream(const Learning_stream & m_ls)
 {
    u               = new double[m_ls.lu];
    lu=m_ls.lu;
    ly=m_ls.ly;
    double * temp   = m_ls.u;
    double * pt     =      u;

    for(int i=0;i<lu;i++)
    {
        *pt=*temp;
        temp++;
        pt++;
    }
    y               = new double[ly];
    temp            = m_ls.y;
    pt              = y;
    for(int i=0;i<ly;i++)
    {
        *pt=*temp;
        temp++;
        pt++;
    }
    L=m_ls.L;
    ni=m_ls.ni;
    eps=m_ls.eps;
 }

 /***************************************************************
    Destruktor
 ***************************************************************/
 Learning_stream::~Learning_stream()
 {
    delete [] u;
    delete [] y;
 }

/***************************************************************
    Funkcja teach
    - uczy sieć neuronową przekazywaną przez referencję jako
    argument net za pomocą ciągu uczącego przekazywanego
    przez referencję jako argument ls algorytmem gradientu
    prostego (ang. Gradient Descent/Ascent)
        Network & net        - uczona sieć neuronowa
        Learning_stream & ls - ciąg uczący
 ***************************************************************/
double grad_desc(Network & net, const Learning_stream & ls)
{
    using std::vector;
    double *    temp_u = ls.u;                 // Wskażnik na tablicę wejść length_u X ls.L
    double *    temp_y = ls.y;                 // Wskażnik na tablicę oczekiwanych wyjść length_y X ls.L
    double *    temp_w;                        // Wskażnik na wagę neuronu
    double *      temp;                        // Wskażnik na wyjście warstwy
    double *        pt;                        // Dodatkowy wskaźnik przy sumowaniu iloczynów delt i wag
    double     epsilon;                        // Zmienna na sumę iloczynów delt i wag
    double         eps  = ls.eps;              // Warunek STOPU
    double          ni  = ls.ni;               // Współczynnik uczenia
    double *     delta;                        // Dynamiczna tablica Delt ze wzoru
                                               /* Tj. Q[i][L](t)*f'(s[i][L](t)) dla ostatniej warstwy
                                                    f'(s[i][k](t))*SUM{m=1:N_k+1}(delta[m][k+1](t)*w[m][i][k+1](t)) */
    double *    temp_d;                        // Wskażnik na deltę
    double *   n_delta;                        // Dodatkowa dynamiczna tablica delt dla nieostatnich warstw
    double *   temp_nd;                        // Wskaźnik na nową deltę
    double           Q;                        // Błąd ze wzoru na Qi = (y-d)^2 gdzie y - wyjście d - oczekiwane wyjście
    double           q  = eps   + 1.0;         // Średni błąd całościowy tj. dla całego ciągu uczącego tj. średnia suma wszystkich błędów Qi
    double           S;                        // Zmienna tymczasowa dla Sumy z neuronu
    int              T  = net.L - 1;           // Ilość warstw bez sumatorów
    int              t;                        // Ilość wyjść z warstwy
    int              n;                        // Ilość wejść do warstwy
    double       pre_q;

    using std::cout;
    using std::ios_base;
    //cout.setf(ios_base::scientific,ios_base::floatfield);


    t       = net.s_layer[0].N;
    temp    = net.s_layer[0].output;
    for(int l=0;l<t;l++)
    {
        temp_w  = net.s_layer[0].neuron[l].wage;
        S       = net.s_layer[0].neuron[l].S;

        Q       = (*temp_y-*temp);
        Q      *= Q;
        q      += Q;

        temp++;
        temp_y++;
    }
    cout<<"\nq="<<q;
    pre_q=q;

    int iter=0;

    while(q>eps)
    {
        temp_u  = ls.u;                // Wskażnik na tablicę wejść length_u X ls.L
        temp_y  = ls.y;                // Wskażnik na tablicę oczekiwanych wyjść length_y X ls.L
        q       =    0;                // Wyzerowanie błędu

        for(int i=0;i<ls.L;i++)
        {
            net.process(temp_u);

            // pt = Y
            pt = net.s_layer[0].output;
            // Warstwa L ma zawsze inny wzór -> Q[i][L](t)*f'(s[i][L](t))
            n       = net.s_layer[0].M;
            t       = net.s_layer[0].N;
            temp    = net.layer[T-1].output; // TEMP = X
            delta   = new  double[t];
            temp_d  = delta;

            vector  <vector <double> > neuron_holder;
            for(int l=0;l<t;l++)
            {
                vector <double> temp_nh;
                temp_w  = net.s_layer[0].neuron[l].wage;
                S       = net.s_layer[0].neuron[l].S;

                Q       = 0;
                Q       = (*temp_y-*pt);
                *temp_d = Q*line_prim(S); // DELTA DLA L==K
                for(int m=0;m<n;m++)
                {
                    temp_nh.push_back(*temp_w);
                    *temp_w += 2*ni*(*temp_d)*(*temp); // W_ijk(t+1)=w_ijk(t)+2*ni*delta_ik(t)*x_jk(t)
                    temp_w++;
                }
                neuron_holder.push_back(temp_nh);

                net.s_layer[0].neuron[l].wage_0 += 2*ni*(*temp_d);
                Q      *= Q;
                q      += Q;
                temp_d++;
                temp++;
                temp_y++;

            }

            // Warstwy inne niż L -> f'(s[i][k](t))*SUM{m=1:N_k+1}(delta[m][k+1](t)*w[m][i][k+1](t))
            for(int j=T-1;j>=0;j--)
            {
                n            = net.layer[j].M;
                t            = net.layer[j].N;
                temp_d       = delta;
                n_delta      = new double[t];
                temp_nd      = n_delta;

                vector  <vector <double> > new_neuron_holder;
                for(int l=0;l<t;l++)
                {
                    vector <double> temp_nnh;
                    temp_w      = net.layer[j].neuron[l].wage;
                    S           = net.layer[j].neuron[l].S;
                    *temp_nd    = sigm_prim(S);
                    //cout<<"\nS   ="<<S;
                    //cout<<"\nS'  ="<<sigm_prim(S);
                    epsilon     = 0;

                    if(j==0)
                        temp    = temp_u;
                    else
                        temp    = net.layer[j].output;

                    if(j==T-1)
                    {
                        for(int o=0;o<net.s_layer[0].N;o++)
                        {
                            //pt       = net.s_layer[0].neuron[o].wage;
                            //pt      += l;
                            epsilon += *temp_d*neuron_holder[o][l];
                            temp_d++;
                        }
                    }else{
                        j++;
                        for(int o=0;o<net.layer[j].N;o++)
                        {
                            //pt       = net.layer[j].neuron[o].wage;
                            //pt      += l;
                            epsilon += *temp_d*(neuron_holder[o][i]);
                            temp_d++;
                        }
                        j--;
                    }

                    *temp_nd *= epsilon;
                    for(int m=0;m<n;m++)
                    {
                        temp_nnh.push_back(*temp_w);
                        *temp_w += 2*ni*(*temp_nd)*(*temp);
                        temp_w++;
                    }
                        new_neuron_holder.push_back(temp_nnh);

                    net.layer[j].neuron[l].wage_0 += 2*ni*(*temp_nd);

                    temp_nd++;
                    temp++;
                }
                neuron_holder = new_neuron_holder;
                delete []      delta;
                delta        = n_delta;
            }

            temp_u += net.layer[0].M;
        }

        delete [] delta;
        cout <<"\nq = "<<q;
        //std::cin.get();
        //std::cin.get();
        //if(q>pre_q&&iter>3)
        //   break;
        //else
        //   pre_q=q;
        iter++;
    }
    return q;
}

/* sigma'(S) */
double sigm_prim(double x)
{
    double temp = (1/(1+exp(-BETA*x)));
    temp        = (BETA*temp*(1-temp));
    return temp;
}

/* tgh'(S)   */
double bipo_prim(double x)
{
    double temp = ((1-exp(BETA*x))/(1+exp(-BETA*x)));
    return (BETA*(1-(temp*temp)));
}

/* (a*s)'    */
double line_prim(double x)
{
    return ALFA;
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

/**********************************************************/
double * lambda(int SIZE,int * structure,int INPUT,Learning_stream learning_stream)
{

    int      iter           = 0;

    vector<net> population  = new pop[pop::m];

    std::srand(std::time(0));
    pop.x    = vector[dimensions];
    pt_pop -> s    = new double[dimensions];
    double * pt_x           = pt_pop -> x;
    double * pt_s           = pt_pop -> s;

    for(int i=0;i<population_size;i++)
    {

        for(int j=0;j<dimensions;j++)
        {
            *pt_x   = (rand()%RAND_MAX)/100.0;
            *pt_s   = 1;
            pt_x++;
            pt_s++;
        }

        pt_pop++;
        pt_pop -> x    = new double[dimensions];
        pt_pop -> s    = new double[dimensions];
        pt_x    = pt_pop -> x;
        pt_s    = pt_pop -> s;
    }

    double     * score              = new double[pop::m];
    double     * pt_score           = score;
    pop        * temp_population    = new pop[pop::l];
    pop        * pt_temp_population = temp_population;
    double     * temp_score         = new double[pop::l];
    double     * pt_temp_score      = temp_score;
    double       change             = 1.0;
    double       best_score         = fun(population->x);

    for(int i=0;i<pop::l;i++)
    {
        pt_temp_population  -> x    = new double[dimensions];
        pt_temp_population  -> s    = new double[dimensions];
        pt_temp_population++;
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
/**********************************************************/
