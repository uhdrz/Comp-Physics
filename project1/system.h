#ifndef SYSTEM_H
#define SYSTEM_H
#include <stdlib.h>
#include <vector>

using namespace std;

class system
{
public:
    system(int c);
    vector<int> state;
    int cut, m_nx, m_ny, m_s;
    int pos;
    void numtopos();
    void postonum();
    double hermite (int n, double x);
    double wavefunc (double x, double y, double w, double A);
    void get (int x, int y, int ms);
    void print();

};




#endif // SYSTEM_H
