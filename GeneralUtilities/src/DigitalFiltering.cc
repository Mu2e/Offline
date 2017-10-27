#include "GeneralUtilities/inc/DigitalFiltering.hh"
#include <stddef.h>
#include <math.h>

namespace mu2e {
  namespace DigitalFiltering {

    void zpk2tf(std::vector<double> &b, std::vector<double> &a, std::vector<double> &za, std::vector<double> &pa)
    {
      a[0] = 1.0;
      if (za.size() > 0){
        a[1] = 0.0;
        for (size_t i=0;i<za.size();i++){
          a[1] += -1*za[i];
        }
      }
      if (za.size() > 1){
        a[2] = 0.0;
        for (size_t i=0;i<za.size();i++){
          for (size_t j=i+1;j<za.size();j++){
            a[2] += za[i]*za[j];
          }
        }
      }
      if (za.size() > 2){
        a[3] = 0.0;
        for (size_t i=0;i<za.size();i++){
          for (size_t j=i+1;j<za.size();j++){
            for (size_t k=j+1;k<za.size();k++){
              a[3] += -1*za[i]*za[j]*za[k];
            }
          }
        }
      }
      b[0] = 1.0;
      if (pa.size() > 0){
        b[1] = 0.0;
        for (size_t i=0;i<pa.size();i++){
          b[1] += -1*pa[i];
        }
      }
      if (pa.size() > 1){
        b[2] = 0.0;
        for (size_t i=0;i<pa.size();i++){
          for (size_t j=i+1;j<pa.size();j++){
            b[2] += pa[i]*pa[j];
          }
        }
      }
      if (pa.size() > 2){
        b[3] = 0.0;
        for (size_t i=0;i<pa.size();i++){
          for (size_t j=i+1;j<pa.size();j++){
            for (size_t k=j+1;k<pa.size();k++){
              b[3] += -1*pa[i]*pa[j]*pa[k];
            }
          }
        }
      }
      if (pa.size() > 3){
        b[4] = 0.0;
        for (size_t i=0;i<pa.size();i++){
          for (size_t j=i+1;j<pa.size();j++){
            for (size_t k=j+1;k<pa.size();k++){
              for (size_t l=k+1;l<pa.size();l++){
                b[4] += pa[i]*pa[j]*pa[k]*pa[l];
              }
            }
          }
        }
      }

    }

    unsigned int iter_factorial(unsigned int n)
    {
      unsigned int ret = 1;
      for(unsigned int i = 1; i <= n; ++i)
        ret *= i;
      return ret;
    }

    double comb(double n, double k)
    {
      return iter_factorial(n)/(iter_factorial(k)*iter_factorial(n-k));
    }

    void bilinear(std::vector<double> &bprime, std::vector<double> &aprime, std::vector<double> &b, std::vector<double> &a, double fs)
    { 
      int D = a.size()-1;
      int N = b.size()-1;

      int M = D;
      if (N > D)
        M = N;

      for (int j=0;j<M+1;j++){
        double val = 0;
        for (int i=0;i<N+1;i++){
          for (int k=0;k<i+1;k++){
            for (int l=0;l<M-i+1;l++){
              if (k+l==j)
                val += comb(i,k)*comb(M-i,l)*b[N-i]*pow(2*fs,i)*pow(-1,k);
            }
          }
        }
        bprime[j] = val;
      }
      for (int j=0;j<M+1;j++){
        double val = 0;
        for (int i=0;i<D+1;i++){
          for (int k=0;k<i+1;k++){
            for (int l=0;l<M-i+1;l++){
              if (k+l==j)
                val += comb(i,k)*comb(M-i,l)*a[D-i]*pow(2*fs,i)*pow(-1,k);
            }
          }
        }
        aprime[j] = val;
      }

      for (int j=M;j>=0;j--){
        bprime[j] /= aprime[0];
        aprime[j] /= aprime[0];
      }
    }

  }
}
