#ifndef TRKEXCHANGEPAR_HH
#define TRKEXCHANGEPAR_HH
#include "BTrk/TrkBase/TrkParams.hh"

// Class interface //
class HelixParams : public TrkParams {
public:
  enum ParIndex {d0Index=0, phi0Index, omegaIndex, z0Index, tanDipIndex, NHLXPRM};

  HelixParams(const CLHEP::HepVector&, const CLHEP::HepSymMatrix&);
  ~HelixParams();
  
  double d0() const                              {return parvec[d0Index];}
  double phi0() const                            {return parvec[phi0Index];}
  double omega() const                           {return parvec[omegaIndex];}
  double z0() const                              {return parvec[z0Index];}
  double tanDip() const                          {return parvec[tanDipIndex];}

  const CLHEP::HepVector& params() const                {return parvec;}
  CLHEP::HepVector& params()                            {return parvec;}
  const CLHEP::HepSymMatrix& covariance() const            {return parcov;}
  CLHEP::HepSymMatrix& covariance()                        {return parcov;}

  void setD0(double in)                          {parvec[d0Index] = in;} 
  void setPhi0(double in)                        {parvec[phi0Index] = in;} 
  void setOmega(double in)                       {parvec[omegaIndex] = in;} 
  void setZ0(double in)                          {parvec[z0Index] = in;} 
  void setTanDip(double in)                      {parvec[tanDipIndex] = in;} 
  void setError(const CLHEP::HepSymMatrix& in)          {parcov = in;}

  void print(std::ostream& o) const;		// Print parameters on one line
  void printAll(std::ostream& o) const;	// Print parameters and error matrix

private:	
};

// Output operator, useful for debugging
std::ostream& operator<<(std::ostream& o, const HelixParams& helix);

#endif
