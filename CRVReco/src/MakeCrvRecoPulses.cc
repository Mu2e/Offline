#include "Offline/CRVReco/inc/MakeCrvRecoPulses.hh"

#include <algorithm>
#include <cmath>
#include <limits>

namespace
{
  typedef std::array<double,3> Vec3;
  typedef std::array<std::array<double,3>,3> Mat3;

  constexpr double kE = 2.718281828459045;
  constexpr int    kLMIterations = 40;        // Levenberg-Marquardt iterations
  constexpr int    kPolishIterations = 8;     // Gauss-Newton iterations with step halving after LM
  constexpr double kLimitTolerance = 0.01;    // a parameter within this fraction of its range from a limit fails the fit
  constexpr double kEdmTolerance = 1e-3;      // converged if the expected chi2 decrease is below this x max(chi2,1)

  double Gumbel(double t, const Vec3 &p)
  {
    double const z = (t-p[1])/p[2];
    return p[0]*std::exp(-z-std::exp(-z));
  }

  // Minuit's sine transformation between bounded external and unbounded internal parameters
  double ToExternal(double u, double lo, double hi) {return lo+(hi-lo)*(std::sin(u)+1.0)*0.5;}
  double ToInternal(double p, double lo, double hi) {return std::asin(std::clamp(2.0*(p-lo)/(hi-lo)-1.0,-1.0,1.0));}

  // solves M x = b by Gaussian elimination with partial pivoting; false if M is singular
  bool Solve3(Mat3 M, Vec3 b, Vec3 &x)
  {
    for(int c=0; c<3; ++c)
    {
      int piv=c;
      for(int r=c+1; r<3; ++r) if(std::abs(M[r][c])>std::abs(M[piv][c])) piv=r;
      if(M[piv][c]==0.0 || !std::isfinite(M[piv][c])) return false;
      std::swap(M[c],M[piv]);
      std::swap(b[c],b[piv]);
      for(int r=c+1; r<3; ++r)
      {
        double const f=M[r][c]/M[c][c];
        for(int k=c; k<3; ++k) M[r][k]-=f*M[c][k];
        b[r]-=f*b[c];
      }
    }
    for(int r=2; r>=0; --r)
    {
      double s=b[r];
      for(int k=r+1; k<3; ++k) s-=M[r][k]*x[k];
      x[r]=s/M[r][r];
    }
    return std::isfinite(x[0]) && std::isfinite(x[1]) && std::isfinite(x[2]);
  }

  // smallest and largest eigenvalue of a symmetric 3x3 matrix (closed form)
  void EigenRange3(const Mat3 &A, double &emin, double &emax)
  {
    double const p1=A[0][1]*A[0][1]+A[0][2]*A[0][2]+A[1][2]*A[1][2];
    double const q=(A[0][0]+A[1][1]+A[2][2])/3.0;
    if(p1==0.0)
    {
      emin=std::min({A[0][0],A[1][1],A[2][2]});
      emax=std::max({A[0][0],A[1][1],A[2][2]});
      return;
    }
    double const p2=(A[0][0]-q)*(A[0][0]-q)+(A[1][1]-q)*(A[1][1]-q)+(A[2][2]-q)*(A[2][2]-q)+2.0*p1;
    double const p=std::sqrt(p2/6.0);
    Mat3 B;
    for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) B[i][j]=(A[i][j]-(i==j?q:0.0))/p;
    double const detB=B[0][0]*(B[1][1]*B[2][2]-B[1][2]*B[2][1])-B[0][1]*(B[1][0]*B[2][2]-B[1][2]*B[2][0])
                     +B[0][2]*(B[1][0]*B[2][1]-B[1][1]*B[2][0]);
    double const phi=std::acos(std::clamp(detB/2.0,-1.0,1.0))/3.0;
    emax=q+2.0*p*std::cos(phi);
    emin=q+2.0*p*std::cos(phi+2.0*M_PI/3.0);
  }
}

namespace mu2eCrv
{

MakeCrvRecoPulses::MakeCrvRecoPulses(float minADCdifference, float defaultBeta, float minBeta, float maxBeta,
                                     float maxTimeDifference, float minPulseHeightRatio, float maxPulseHeightRatio,
                                     float LEtimeFactor, int samplesBefore, int samplesAfter) :
                                     _minADCdifference(minADCdifference),
                                     _defaultBeta(defaultBeta), _minBeta(minBeta), _maxBeta(maxBeta),
                                     _maxTimeDifference(maxTimeDifference),
                                     _minPulseHeightRatio(minPulseHeightRatio),
                                     _maxPulseHeightRatio(maxPulseHeightRatio),
                                     _LEtimeFactor(LEtimeFactor),
                                     _samplesBefore(samplesBefore), _samplesAfter(samplesAfter)
{
}

bool MakeCrvRecoPulses::FindNextPeak(const std::vector<double> &y, int start, int &peakStart, int &peakEnd, int &fitStart, int &fitEnd) const
{
  int const n=y.size();
  bool foundRisingEdge=false;
  bool foundPeak=false;
  for(int bin=start+1; bin<n; ++bin)
  {
    if(y[bin-1]<y[bin]) //rising edge
    {
      foundRisingEdge=true;
      peakStart=bin;
      peakEnd=bin;
    }
    if(y[bin-1]==y[bin] && foundRisingEdge) //potentially a peak where consecutive ADC values are equal
    {
      peakEnd=bin;
    }
    if(y[bin-1]>y[bin] && foundRisingEdge) //falling edge
    {
      if(y[peakEnd]>_minADCdifference)
      {
        foundPeak=true;
        break; //found the peak. stop searching.
      }
      foundRisingEdge=false;  //so that the loop has to wait for the next rising edge
    }
  }

  if(!foundPeak) return false;

  //select a range of up to 3 points before and 2 after the peak (to avoid hidden reflection pulses)
  //-find up to 4/3 points before and after the peak for which the waveform is stricly decreasing
  //-remove 1 point on each side. this removes potentially "bad points" belonging to a second pulse (i.e. in double pulses)

  fitStart=peakStart-1;
  fitEnd=peakEnd+1;
  for(int i=peakStart-1; i>=peakStart-_samplesBefore; --i)
  {
    if(y[i]<=y[i+1]) fitStart=i;
    else break;
    if(i==0) break;
  }
  for(int i=peakEnd+1; i<n && i<=peakEnd+_samplesAfter; ++i)
  {
    if(y[i-1]>=y[i]) fitEnd=i;
    else break;
  }
  if(peakStart-fitStart>1) fitStart++;
  if(fitEnd-peakEnd>1) fitEnd--;

  return true;
}

//bounded least-squares Gumbel fit: Levenberg-Marquardt in sine-mapped parameters, then a Gauss-Newton polish
MakeCrvRecoPulses::FitResult MakeCrvRecoPulses::FitGumbel(const std::vector<double> &t, const std::vector<double> &y,
                                                          const Vec3 &p0, const Vec3 &lower, const Vec3 &upper) const
{
  size_t const n=t.size();
  std::vector<double> r(n), rn(n);
  std::vector<Vec3> J(n), Jn(n);

  auto residualsAndJacobian=[&](const Vec3 &u, std::vector<double> &res, std::vector<Vec3> &jac)
  {
    Vec3 p, dpdu;
    for(int k=0; k<3; ++k)
    {
      p[k]=ToExternal(u[k],lower[k],upper[k]);
      dpdu[k]=(upper[k]-lower[k])*0.5*std::cos(u[k]);
    }
    double chi2=0;
    for(size_t i=0; i<n; ++i)
    {
      double const z=(t[i]-p[1])/p[2];
      double const ez=std::exp(-z);
      double const g=std::exp(-z-ez);
      double const f=p[0]*g;
      double const dfdz=f*(ez-1.0);
      res[i]=y[i]-f;
      jac[i]={g*dpdu[0], -dfdz/p[2]*dpdu[1], -dfdz*z/p[2]*dpdu[2]};
      chi2+=res[i]*res[i];
    }
    return chi2;
  };
  auto normalEquations=[&](Mat3 &JTJ, Vec3 &JTr)
  {
    for(int a=0; a<3; ++a)
    {
      JTr[a]=0;
      for(int b=0; b<3; ++b) JTJ[a][b]=0;
    }
    for(size_t i=0; i<n; ++i)
    {
      for(int a=0; a<3; ++a)
      {
        JTr[a]+=J[i][a]*r[i];
        for(int b=0; b<3; ++b) JTJ[a][b]+=J[i][a]*J[i][b];
      }
    }
  };

  Vec3 u;
  for(int k=0; k<3; ++k) u[k]=ToInternal(p0[k],lower[k],upper[k]);
  double chi2=residualsAndJacobian(u,r,J);
  double lambda=1e-3;
  Mat3 JTJ;
  Vec3 JTr, du;

  for(int iter=0; iter<kLMIterations; ++iter)
  {
    normalEquations(JTJ,JTr);
    Mat3 M=JTJ;
    for(int k=0; k<3; ++k) M[k][k]+=lambda*(JTJ[k][k]+1e-12);
    if(!Solve3(M,JTr,du)) du={0,0,0};
    Vec3 un={u[0]+du[0],u[1]+du[1],u[2]+du[2]};
    double const chi2n=residualsAndJacobian(un,rn,Jn);
    bool const better=chi2n<chi2;
    if(better)
    {
      u=un;
      r.swap(rn);
      J.swap(Jn);
      chi2=chi2n;
    }
    lambda=std::clamp(better?lambda*0.3:lambda*10.0,1e-9,1e9);
  }

  for(int iter=0; iter<kPolishIterations; ++iter)
  {
    normalEquations(JTJ,JTr);
    Mat3 M=JTJ;
    for(int k=0; k<3; ++k) M[k][k]+=1e-12;
    if(!Solve3(M,JTr,du)) break;
    for(double scale : {1.0,0.5,0.25,0.125})
    {
      Vec3 un={u[0]+scale*du[0],u[1]+scale*du[1],u[2]+scale*du[2]};
      double const chi2n=residualsAndJacobian(un,rn,Jn);
      if(chi2n<chi2)
      {
        u=un;
        r.swap(rn);
        J.swap(Jn);
        chi2=chi2n;
        break;
      }
    }
  }

  FitResult result;
  for(int k=0; k<3; ++k) result.par[k]=ToExternal(u[k],lower[k],upper[k]);
  result.chi2=chi2;
  normalEquations(JTJ,JTr);
  Mat3 M=JTJ;
  for(int k=0; k<3; ++k) M[k][k]+=1e-12;
  Vec3 step;
  if(Solve3(M,JTr,step))
  {
    double const edm=0.5*(JTr[0]*step[0]+JTr[1]*step[1]+JTr[2]*step[2]);
    double emin, emax;
    EigenRange3(JTJ,emin,emax);
    result.converged=std::isfinite(edm) && edm<=kEdmTolerance*std::max(chi2,1.0);
    result.hessianOk=emin>1e-12*std::max(emax,1e-300);
  }
  else
  {
    result.converged=false;
    result.hessianOk=false;
  }
  return result;
}

void MakeCrvRecoPulses::SetWaveform(const std::vector<int16_t> &waveform,
                                    uint16_t startTDC, float digitizationPeriod, float pedestal,
                                    float calibrationFactor, float calibrationFactorPulseHeight)
{
  _pulseTimes.clear();
  _pulseHeights.clear();
  _pulseBetas.clear();
  _pulseFitChi2s.clear();
  _PEs.clear();
  _PEsPulseHeight.clear();
  _LEtimes.clear();
  _zeroNdf.clear();
  _failedFits.clear();

  size_t nBins=waveform.size();
  std::vector<double> y(nBins);
  for(size_t bin=0; bin<nBins; ++bin) y[bin]=waveform[bin]-pedestal;
  std::vector<double> tFit, yFit;

  //loop through all peaks
  int start=0;
  int peakStart=0;
  int peakEnd=0;
  int fitStart=0;
  int fitEnd=0;
  while(FindNextPeak(y, start, peakStart, peakEnd, fitStart, fitEnd))
  {
    double peakStartTime=(startTDC+peakStart)*digitizationPeriod;
    double peakEndTime=(startTDC+peakEnd)*digitizationPeriod;
    double peakTime=0.5*(peakStartTime+peakEndTime);

    double const peakValue=y[peakStart];
    Vec3 const p0={peakValue*kE, peakTime, _defaultBeta};
    Vec3 const lower={peakValue*kE*_minPulseHeightRatio, peakStartTime-_maxTimeDifference, _minBeta};
    Vec3 const upper={peakValue*kE*_maxPulseHeightRatio, peakEndTime+_maxTimeDifference, _maxBeta};

    tFit.clear();
    yFit.clear();
    for(int bin=fitStart; bin<=fitEnd; ++bin)
    {
      tFit.push_back((startTDC+bin)*digitizationPeriod);
      yFit.push_back(y[bin]);
    }

    //do the fit
    FitResult const fr=FitGumbel(tFit, yFit, p0, lower, upper);
    int const ndf=static_cast<int>(tFit.size())-3;
    bool atLimit=false;
    for(int k=0; k<3; ++k)
    {
      double const range=upper[k]-lower[k];
      if((fr.par[k]-lower[k])/range<kLimitTolerance || (upper[k]-fr.par[k])/range<kLimitTolerance) atLimit=true;
    }

    //collect fit information for the first peak
    float  PEs          = fr.par[0]*fr.par[2] / calibrationFactor;
    double pulseTime    = fr.par[1];
    float  pulseHeight  = fr.par[0]/kE;
    float  pulseBeta    = fr.par[2];
    float  pulseFitChi2 = (ndf>0?fr.chi2/ndf:-1);
    bool   zeroNdf      = (ndf>0?false:true);
    bool   failedFit    = !fr.converged || !fr.hessianOk || atLimit;

    if(failedFit)
    {
      PEs          = peakValue * kE * _defaultBeta / calibrationFactor;
      pulseTime    = peakTime;
      pulseHeight  = peakValue;
      pulseBeta    = _defaultBeta;
      pulseFitChi2 = -1;
    }

    double LEtime         = pulseTime-_LEtimeFactor*pulseBeta;  //50% pulse height is reached at -0.985*beta before the peak
    float  PEsPulseHeight = pulseHeight / calibrationFactorPulseHeight;

    _pulseTimes.push_back(pulseTime);
    _pulseHeights.push_back(pulseHeight);
    _pulseBetas.push_back(pulseBeta);
    _pulseFitChi2s.push_back(pulseFitChi2);
    _PEs.push_back(PEs);
    _PEsPulseHeight.push_back(PEsPulseHeight);
    _LEtimes.push_back(LEtime);
    _zeroNdf.push_back(zeroNdf);
    //if one pulse reconstruction fails, all subsequent pulses should also get the error flag,
    //because a wrong pulse subtraction will create wrong secondary pulses
    if(_failedFits.size()>0)
    {
      if(_failedFits.back()) _failedFits.push_back(true);
      else _failedFits.push_back(failedFit);
    }
    else _failedFits.push_back(failedFit);

    //subtract the fitted pulse (also for a failed fit) before searching for the next peak
    for(size_t bin=0; bin<nBins; ++bin) y[bin]-=std::round(Gumbel((startTDC+bin)*digitizationPeriod, fr.par));
    start=fitEnd;
  }

}

}
