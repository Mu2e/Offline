#ifndef FOFM_FOFM_H
#define FOFM_FOFM_H

// ----------------------------------------------------------------------
//
// FofM.h
//
// Figure of Merit computation class
//
// ----------------------------------------------------------------------

#include <vector>
#include <iostream>
#include <string>

#include "splines/Grid.h"
#include "splines/Spline.h"
#include "FigureOfMerit/inc/FeldmanCousinsSensitivity.h"
// include of root THF1 header and maybe others



namespace mu2e {

enum MeritFunctionChoice {
    SmoothedPunziMeritFunction
  , FCsensitivityMeritFunction
};

class FofM {

// Utility sub-classes for spectra
// -------------------------------

public:

  class Spectrum;
  class Smearing;

public:

  static double defaultSpectrumLowerLimit() {return  98.0;}
  static double defaultSpectrumUpperLimit() {return 108.0;}
  static double defaultSmearingLowerLimit() {return  -8.0;}
  static double defaultSmearingUpperLimit() {return   8.0;}
  static int    defaultSpectrumNodes()      {return 1000;}
  static int    defaultSmearingNodes()      {return  800;}

class Spectrum {
  public:
  // ctors
  // From a function, or a function object
  explicit Spectrum ( std::function<double (double)> const & shape );
           Spectrum ( std::function<double (double)> const & shape, 
                      double range_low, double range_high );
  // explicit Spectrum ( root::THF1 const & shape );  // TODO -- not yet implemented
  //          Spectrum ( root::THF1 const & shape,
  //                     double range_low, double range_high );  
           Spectrum ( std::vector<double> const &  strengths, 
                      double range_low, double range_high );
  // Convolution
  FofM::Spectrum convolve( FofM::Smearing const & resolution ) const;
  // Usage by FofM class
  splines::Grid<1> grid() const;
  splines::Spline<1> representation() const;
  double a() const {return a_;}
  double b() const {return b_;}
  private:
  double a_;  // range lower limit
  double b_;  // range upper limit
  unsigned int nodes_;
  double gridSpacing_;
  splines::Grid<1> grid_;
  splines::Spline<1> spline_;
}; // Spectrum

class Smearing {
  public:
  // ctors
  // From a function, or a function object
  explicit Smearing ( std::function<double (double)> const & shape );
           Smearing ( std::function<double (double)> const & shape, 
                      double range_low, double range_high );
  // explicit Smearing ( root::THF1 const & shape );  // TODO -- not yet implemented
  //          Smearing ( root::THF1 const & shape,
  //                     double range_low, double range_high );  
           Smearing ( std::vector<double> const &  shape, 
                      double range_low, double range_high );
  // Convolution
  FofM::Smearing convolve( FofM::Smearing const & s2 ) const;
  // Usage by FofM class
  splines::Grid<1> grid() const;
  splines::Spline<1> representation() const;
  double a() const {return a_;}
  double b() const {return b_;}
  private:
  double a_;  // range lower limit
  double b_;  // range upper limit
  unsigned int nodes_;
  double gridSpacing_;
  splines::Grid<1> grid_;
  splines::Spline<1> spline_;
}; // Smearing 

// Utility sub-class for summary
// -----------------------------

public:

  struct Summary;

struct Summary {
  double singleEventSensitivity;
  double CL90sensitivity;
  double smoothedPunziSensitivity;
  double pCutLo;
  double pCutHi;
  double figureOfMerit;
};  


// The main FofM class
// -------------------

public:
  // ctors
  
  FofM ( FofM::Spectrum signalEfficiency, 
         double protonsOnProductionTarget,
         MeritFunctionChoice mfunction,
         std::ostream & ost);
  FofM ( FofM::Spectrum signalEfficiency, 
         FofM::Spectrum backgroundStrength,
         double protonsOnProductionTarget,
         MeritFunctionChoice mfunction,
         std::ostream & ost); 
  FofM ( FofM::Spectrum signalEfficiency, 
         FofM::Spectrum backgroundStrength,
         FofM::Smearing resolution,
         double protonsOnProductionTarget,
         MeritFunctionChoice mfunction,
         std::ostream & ost);

  // Preparing signal and background spectra

  void addSignal       ( FofM::Spectrum signalEfficiency );
  void addSignal       ( FofM::Spectrum signalEfficiency, 
                         FofM::Smearing resolution );
  void addBackground   ( FofM::Spectrum bkg ); 
  void addBackground   ( FofM::Spectrum bkg, 
                         FofM::Smearing resolution );
  void smearSignal     ( FofM::Smearing resolution );  
  void smearBackground ( FofM::Smearing resolution );  
  void smearAll        ( FofM::Smearing resolution );  
  
  // Obtain Figure of Merit and hypothetical results table
  double figureOfMerit () const;
  double worstCaseBranchingRatio () const;
  double lowMomentumCut () const;
  double highMomentumCut () const;
  std::string tables (int maximumSignalCount, double & lowCut, double & highCut,
                      std::ostream & os, FofM::Summary & summary) const; 
  void results (int maximumSignalCount,
                double &lowCut,  double &highCut, 
                double &signalEff,
                double &bkgd,
                double &figOfM, 
                int    &discoveryThresholdN, 
                double &discoveryThreshold,
                double &smoothedWorstCaseSensitivity,
                double &punziSensitivity,
                double &CL90Sensitivity,
                double &singleEventSensitivity,
                std::vector<double> & upperLimitBR,
                std::vector<double> & CLlowerBR,
                std::vector<double> & CLupperBR,
                std::ostream & os ) const; 

  // Obtain Figure of Merit and hypothetical results table,pinning momentum cuts

  // fix high cut only
  double figureOfMerit_fixed_highCut (double momcutHigh) const;
  double worstCaseBranchingRatio_fixed_highCut (double momcutHigh) const;
  double lowMomentumCut_fixed_highCut (double momcutHigh) const;
  std::string tables_fixed_highCut 
        (int maximumSignalCount, double& lowCut, double momcutHigh, 
                     std::ostream & os, FofM::Summary & summary) const; 
  void results_fixed_highCut (int maximumSignalCount, 
                double &lowCut,  double highCut, 
                double &signalEff,
                double &bkgd,
                double &figOfM, 
                int    &discoveryThresholdN, 
                double &discoveryThreshold,
                double &smoothedWorstCaseSensitivity,
                double &punziSensitivity,
                double &CL90Sensitivity,
                double &singleEventSensitivity,
                std::vector<double> & upperLimitBR,
                std::vector<double> & CLlowerBR,
                std::vector<double> & CLupperBR ,
                std::ostream & os ) const;  

  // fix low cut only
  double figureOfMerit_fixed_lowCut (double momcutLow) const;
  double worstCaseBranchingRatio_fixed_lowCut (double momcutLow) const;
  double highMomentumCut_fixed_highCut (double momcutLow) const;
  std::string tables_fixed_lowCut 
        (int maximumSignalCount, double momcutLow, double& highCut,
                     std::ostream & os, FofM::Summary & summary) const; 
  void results_fixed_lowCut (int maximumSignalCount, 
                double lowCut,  double &highCut, 
                double &signalEff,
                double &bkgd,
                double &figOfM, 
                int    &discoveryThresholdN, 
                double &discoveryThreshold,
                double &smoothedWorstCaseSensitivity,
                double &punziSensitivity,
                double &CL90Sensitivity,
                double &singleEventSensitivity,
                std::vector<double> & upperLimitBR,
                std::vector<double> & CLlowerBR,
                std::vector<double> & CLupperBR,
                std::ostream & os ) const; 

  // fix both cuts
  double figureOfMerit_fixed_cuts 
        (double momcutLow, double momcutHigh) const;
  double worstCaseBranchingRatio_fixed_cuts 
        (double momcutLow, double momcutHigh) const;
  std::string tables_fixed_cuts (int maximumSignalEventCount, 
                double momcutLow, double momcutHigh, 
                     std::ostream & os, FofM::Summary & summary) const; 
  void results_fixed_cuts  (int maximumSignalCount, 
                double lowCut, double highCut, 
                double &signalEff,
                double &bkgd,
                double &figOfM,  
                int    &discoveryThresholdN,
                double &discoveryThreshold,
                double &smoothedWorstCaseSensitivity,
                double &punziSensitivity,
                double &CL90Sensitivity,
                double &singleEventSensitivity,
                std::vector<double> & upperLimitBR,
                std::vector<double> & CLlowerBR,
                std::vector<double> & CLupperBR,
                std::ostream & os ) const;  

  // The following may be useful for understanding calculations or debugging
  static double punziMeritDenominator ( double B );
  static double discoveryCount   ( double B );   
  void displayBackground(double pmin, double pmax, int n);
  splines::Grid<1>   getGrid      () const {return grid;}
  splines::Spline<1> getBackground() const {return background;}
  splines::Spline<1> getSignal    () const {return signal;}
 
  void computeSensitivities( 
                double eps, double B,  
                int    &discoveryThresholdN,
                double &discoveryThreshold,
                double &smoothedWorstCaseSensitivity,
                double &punziSensitivity,
                double &CL90Sensitivity,
                double &singleEventSensitivity ) const;

  void fillHypotheticalVectors (int maximumSignalCount,
                double eps, double B, int discoveryThresholdN, 
                std::vector<double> & upperLimitBR,
                std::vector<double> & CLlowerBR,
                std::vector<double> & CLupperBR ) const;
  
private:
  std::ostream & os;
  double L;
  double a; // momentum range bottom
  double b; // momentum range top
  splines::Spline<1> signal;
  splines::Spline<1> background;
  splines::Grid<1>   grid;
  MeritFunctionChoice mfunc;
  mutable bool computationIsDone;
  // TODO - other mutables

// helper methods
private:
  double optimize_low_cut  (double fixed_high_cut) const;
  double optimize_high_cut (double fixed_low_cut)  const;
  std::string tablesString (int maximumSignalCount,
                double lowCut, double highCut, 
                double epsilon,  // integratedSignalEfficiency
                double B,        // integratedBackground
                double figOfM, 
                int    discoveryThresholdN, 
                double discoveryThreshold,
                double smoothedWorstCaseSensitivity,
                double punziSensitivity,
                double CL90Sensitivity,
                double singleEventSensitivity,
                std::vector<double> const & upperLimitBR,
                std::vector<double> const & CLlowerBR,
                std::vector<double> const & CLupperBR,
                bool lowCutFixed, bool highCutFixed) const;
  FeldmanCousins90pctSensitivity fc90MeritDenominator;
   
};  // FofM

// --------------------
// Pre-packaged Spectra
// --------------------

// The spectra we provide are:
//
//    For use as signals
//      CEdeltaFunction
//
//    For use as backgrounds
//      DIOmomentumSpectrum (= CzarneckiDIOSpectrum)
//      CzarneckiDIOSpectrum
//      FlatBackground
//      LinearBackground
//
//    For use in smearing
//      InnerBremsstrahlungCorrection
//      GaussianSmearing
//      TwoGaussianAsymmetricSmearing
//       

// Functions used to help create splines

double zeroFunctionForNullBackgroundSpectrum (double x) {return 0.0;}

double branchingFractionForCapturedMuonsToDIO();

struct CzarneckiDIOspectrumFunction  {
public:
  explicit CzarneckiDIOspectrumFunction(double strength); 
  double operator() (double p) const;  
  double integrated (double p) const;  // integral from p to endpoint   
  double operator() (double *p, double*) const {return operator()(*p);}  
private:
  double multiplier;
  double alpha_5;
  double alpha_6;
  double alpha_7;
  double alpha_8;
  double Emax;
  double E_mu;
  double M_Al;
  double E_mue;
  double beta_6;
  double beta_7;
  double beta_8;
  double beta_9;
  double beta_10;
};

// Spectra -- each is a class inheriting from Spectrum

class NullBackgroundSpectrum : public FofM::Spectrum {
  public:
  // ctors
  NullBackgroundSpectrum();

  // All functions are done by the base class
  // Note - convolve can be done more efficiently by returning 
  // NullBackgroundSpectrum().  But convolving with the null spectrum is not 
  // going to be a usual mode of usage anyway.
    
}; // NullBackgroundSpectrum

class CzarneckiDIOSpectrum : public FofM::Spectrum {
  public:
  // ctors
  explicit CzarneckiDIOSpectrum( 
        double strength, double a = 98.0, double b = 108.0 ); 
  // strength is the number of muons stopped in stopping target, 
  // per proton on production target.

  // All functions are done by the base class

}; //  CzarneckiDIOSpectrum

typedef CzarneckiDIOSpectrum DIOmomentumSpectrum;





// TODO --
//      prepackaged spectra
//      prepackaged smearings


} // end of namespace mu2e

// Output results table  
std::ostream & operator<< (std::ostream & os, mu2e::FofM const & fm) {
  mu2e::FofM::Summary s;
  double lowCut;
  double highCut;
  std::string t = fm.tables(25, lowCut, highCut, os, s);
  return os << t << "\nFigure of merit = " << s.figureOfMerit << "\n";
}

// -----------------------------
// Helper Functions and Functors
// -----------------------------

namespace mu2e {

class FofMwithLowerLimitVarying {
public:
  FofMwithLowerLimitVarying( splines::Spline<1> const & s, 
                             splines::Spline<1> const & b,
                             double luminosity, 
                             mu2e::MeritFunctionChoice mfunction )
      : negIntegratedSignal(s)
      , negIntegratedBackground(b)
      , L(luminosity)
      , mfunc (mfunction) {}
  double operator()(double a) const;
private:
  splines::Spline<1> negIntegratedSignal;
  splines::Spline<1> negIntegratedBackground;
  double L;
  MeritFunctionChoice mfunc;
  FeldmanCousins90pctSensitivity fc90MeritDenominator;
};

class FofMwithUpperLimitVarying {
public:
  FofMwithUpperLimitVarying( splines::Spline<1> const & s, 
                             splines::Spline<1> const & b,
                             double luminosity, 
                             MeritFunctionChoice mfunction )
      : integratedSignal(s)
      , integratedBackground(b)
      , L(luminosity) 
      , mfunc (mfunction) {}
  double operator()(double a) const;
private:
  splines::Spline<1> integratedSignal;
  splines::Spline<1> integratedBackground;
  double L;
  MeritFunctionChoice mfunc;
  FeldmanCousins90pctSensitivity fc90MeritDenominator;
};

} // end namespace mu2e

#endif  /* FOFM_FOFM_H */
