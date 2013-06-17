//
//  $Id: 
//  $Author: 
//  $Date: 
//
//  Original author Vadim Rusu

#ifndef PIDProduct_HH
#define PIDProduct_HH

#include <utility>

namespace mu2e {



  class PIDProduct {

  public:
    PIDProduct(); 
    PIDProduct (const PIDProduct & p) ; 
    ~PIDProduct() {;} 
    PIDProduct & operator = (const PIDProduct & p) ;

    void clear() ;

    void SetTrkID(int id){_trkid=id;}
    void SetResidualsSlope(double s){_residualsSlope=s;}
    void SetResidualsSlopeError(double s){_residualsSlopeError=s;}

    void SetLogEProb(double d) { _logeprob = d;}
    void SetLogMProb(double d) { _logmprob = d;}
    
    double GetResidualsSlope() const { return _residualsSlope;}
    double GetResidualsSlopeError() const { return _residualsSlopeError;}
    double GetTrkID() const { return _trkid;}
    double GetLogEProb() const { return _logeprob;}
    double GetLogMProb() const { return _logmprob;}


  private:
    double _residualsSlope;
    double _residualsSlopeError;
    double _logeprob;
    double _logmprob;
    int _trkid;

  };



} // end namespace mu2e


#endif
