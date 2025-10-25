class mxmat_real
{
private:
    unsigned int M, N;
    double *data;
public:
    mxmat_real(double *data_input, unsigned long Minput, unsigned long Ninput)
    {
        *data = *data_input;
        M = Minput;
        N = Ninput;
    }
    double &term(unsigned long m, unsigned long n)
    {
        return data[m + n*M];
    }
    
};

class pcomplex
{
public:
    double *real, *imag;
    pcomplex()
    pcomplex(double *real_in, double *imag_in)
    {
        real = real_in;
        imag = imag_in;
    }
    friend const pcomplex operator+ (const pcomplex &left, const pcomplex &right);
    friend const pcomplex operator+ (double left, const pcomplex &right);
    friend const pcomplex operator+ (const pcomplex &left, double right);    
};

const pcomplex operator+ (const pcomplex &left, const pcomplex &right)
{
    *this.real = *left.real + *right.real;
    *this.imag = *left.imag + *right.imag;
    return *this;
}