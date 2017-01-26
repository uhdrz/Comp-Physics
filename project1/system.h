#ifndef SYSTEM_H
#define SYSTEM_H


class system
{
public:
    system(int c);
    int state[];
    int cut, m_nx, m_ny, m_s;
    int numtopos();
    void postonum(int pos);
    double hermite (int n, double x);
    double wavefunc (double x, double w,double y, double A);
    void genstate();
    int cutoff(int n1, int n2, int cut_test);
    void get (int x, int y, int ms);

};





#endif // SYSTEM_H
