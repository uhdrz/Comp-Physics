#ifndef hattree_H
#define hattree_H
#include <stdlib.h>
#include <vector>

using namespace std;

class hattree
{
public:
    hattree(int n_x, int n_y, int s, double w, int c);
    hattree();
    vector<int> state;
    int cut, m_nx, m_ny, m_s;
    int pos;
    double m_w;
    double energy();
    void numtopos();
    void postonum();
    double hermite (int n, double x);
    double wavefunc (double x, double y);
    void genstate();
    void get (int x, int y, int ms);
    void print();
    long fac(int n);

    /*Integration*/
    void GaussHermiteQuadrature(double *x, double *w, int n);
    double potential (double x1, double x2, double y1, double y2);
    double Integration(int n,class hattree p,class hattree q,class hattree r, class hattree s); /*<pq|V|rs>*/
    double ASMatrixElement(int n,class hattree p,class hattree q,class hattree r, class hattree s);/*<pq|V|rs>-<pq|V|sr>*/




};


/*destruktor*, genstate*/


#endif // hattree_H
