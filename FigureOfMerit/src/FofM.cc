// FofM.cc
//
// Implementation for methods of FofM>
//

#include "splines/Spline.h"
#include "FigureOfMerit/inc/FofM.h"
#include "FigureOfMerit/inc/PoissonCDFsolver.h"
#include "FigureOfMerit/inc/FeldmanCousinsSensitivity.h"

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>

using std::setw;

namespace mu2e {

// --------
// Spectrum
// --------

FofM::Spectrum::Spectrum ( std::function<double (double)> const & shape )
  : a_           (defaultSpectrumLowerLimit())
  , b_           (defaultSpectrumUpperLimit())
  , nodes_       (defaultSpectrumNodes())
  , gridSpacing_ ( (b_-a_)/nodes_ )
  , grid_        ( a_ , b_+gridSpacing_ , nodes_ )
  , spline_      ( shape, grid_ )
{
} // ctor from function

FofM::Spectrum::Spectrum ( std::function<double (double)> const & shape, 
                           double range_low, double range_high )
  : a_           (range_low)
  , b_           (range_high)
  , nodes_       (defaultSpectrumNodes())
  , gridSpacing_ ( (b_-a_)/nodes_ )
  , grid_        ( a_ , b_+gridSpacing_ , nodes_ )
  , spline_      ( shape, grid_ )
{
} // ctor from function with specified range

FofM::Spectrum::Spectrum ( std::vector<double> const &  strengths, 
                           double range_low, double range_high )
  : a_           (range_low)
  , b_           (range_high)
  , nodes_       (strengths.size())
  , gridSpacing_ ( (b_-a_)/nodes_ )
  , grid_        ( a_ , b_+gridSpacing_ , nodes_ )
  , spline_      ( strengths, grid_ )
{
} // ctor from vector with specified range

splines::Grid<1> FofM::Spectrum::grid() const {
  return grid_;
}
  
splines::Spline<1> FofM::Spectrum::representation() const {
  return spline_;
}

FofM::Spectrum FofM::Spectrum::convolve( FofM::Smearing const & resolution ) const
{
  return FofM::Spectrum( spline_.convolution(resolution.representation()), a_, b_ );
  // OPTIMIZATION -- This makes a spline from a spline.
  //         It is almost twice as efficient copy the spline but I do not want to 
  //         deal with range issues at this moment.   
} // convolve

// --------
// Smearing
// --------

FofM::Smearing::Smearing ( std::function<double (double)> const & shape )
  : a_           (defaultSmearingLowerLimit())
  , b_           (defaultSmearingUpperLimit())
  , nodes_       (defaultSmearingNodes())
  , gridSpacing_ ( (b_-a_)/nodes_ )
  , grid_        ( a_ , b_+gridSpacing_ , nodes_ )
  , spline_      ( shape, grid_ )
{
  // TODO -- normalize spline_
} // ctor from function

FofM::Smearing::Smearing ( std::function<double (double)> const & shape, 
                           double range_low, double range_high )
  : a_           (range_low)
  , b_           (range_high)
  , nodes_       (defaultSmearingNodes())
  , gridSpacing_ ( (b_-a_)/nodes_ )
  , grid_        ( a_ , b_+gridSpacing_ , nodes_ )
  , spline_      ( shape, grid_ )
{
  // TODO -- normalize spline_
} // ctor from function with specified range

FofM::Smearing::Smearing ( std::vector<double> const &  shape, 
                           double range_low, double range_high )
  : a_           (range_low)
  , b_           (range_high)
  , nodes_       (shape.size())
  , gridSpacing_ ( (b_-a_)/nodes_ )
  , grid_        ( a_ , b_+gridSpacing_ , nodes_ )
  , spline_      ( shape, grid_ )
{

} // ctor from vector with specified range

FofM::Smearing FofM::Smearing::convolve( FofM::Smearing const & resolution ) const
{
  return FofM::Smearing( spline_.convolution(resolution.representation()), a_, b_ );
  // OPTIMIZATION -- This makes a spline from a spline.
  //         It is almost twice as efficient copy the spline but I do not want to 
  //         deal with range issues at this moment.   
} // convolve

splines::Grid<1> FofM::Smearing::grid() const {
  return grid_;
}
  
splines::Spline<1> FofM::Smearing::representation() const {
  return spline_;
}

// ----
// FofM
// ----

FofM::FofM ( FofM::Spectrum signalEfficiency, 
             double protonsOnProductionTarget,
             MeritFunctionChoice mfunction )
  : L (protonsOnProductionTarget)
  , a (NullBackgroundSpectrum().a())
  , b (NullBackgroundSpectrum().b())
  , signal (signalEfficiency.representation())
  , background (NullBackgroundSpectrum().representation())
  , grid (NullBackgroundSpectrum().grid())
  , mfunc(mfunction)
  , computationIsDone(false)
{
} // ctor from signal and L; null background implied
  
FofM::FofM ( FofM::Spectrum signalEfficiency, 
             FofM::Spectrum backgroundStrength,
             double protonsOnProductionTarget,
             MeritFunctionChoice mfunction )
  : L (protonsOnProductionTarget)
  , a (backgroundStrength.a())
  , b (backgroundStrength.b())
  , signal (signalEfficiency.representation())
  , background (backgroundStrength.representation())
  , grid (backgroundStrength.grid())
  , mfunc(mfunction)
  , computationIsDone(false)
{
  // TODO -- We are taking the grid from the background only.
  //         That is probably right but we should think it through.
  // std::cout << "FofM ctor from signal and background\n";
} // ctor from signal, background  and L
         
FofM::FofM ( FofM::Spectrum signalEfficiency, 
             FofM::Spectrum backgroundStrength,
             FofM::Smearing resolution,
             double protonsOnProductionTarget,
             MeritFunctionChoice mfunction )
  : L (protonsOnProductionTarget)
  , a (backgroundStrength.a())
  , b (backgroundStrength.b())
  , signal 
     (signalEfficiency.representation().convolution(resolution.representation()))
  , background 
     (backgroundStrength.representation().convolution(resolution.representation()))
  , grid (backgroundStrength.grid())
  , mfunc(mfunction)
  , computationIsDone(false)
{
  // TODO -- We are taking the grid from the background only.
  //         That is probably right but we should think it through.
} // ctor from signal, bcakground, resolution and L

void FofM::addSignal       ( FofM::Spectrum signalEfficiency )
{
  // TODO -- implement
}
void FofM::addSignal  ( FofM::Spectrum signalEfficiency, FofM::Smearing resolution )
{
  // TODO -- implement
}

void FofM::addBackground   ( FofM::Spectrum bkg ) 
{
  std::vector<double> newBackground(grid.nPoints());
  for (unsigned int i = 0; i < grid.nPoints(); ++i) {
    double x = grid[i];
    newBackground[i] = background(x)+bkg.representation()(x);
  }
  splines::Spline<1> b(newBackground, grid);
  background = b;
}

void FofM::addBackground  ( FofM::Spectrum background, FofM::Smearing resolution )
{
  // TODO -- implement
}

void FofM::smearSignal     ( FofM::Smearing resolution )  
{
  // TODO -- implement
}
void FofM::smearBackground ( FofM::Smearing resolution )  
{
  // TODO -- implement
}
void FofM::smearAll        ( FofM::Smearing resolution )  
{
  // TODO -- implement
}

void FofM::results (int maximumSignalCount,
                double &lowCut, double &highCut, 
                double &epsilon,  // integratedSignalEfficiency
                double &B,        // integratedBackground
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
                std::ostream & os ) const
{
  // starting guesses for momentum cuts
  double old_pmin = 103.6;
  double pmin = old_pmin;
  double old_pmax = 105.0;
  double pmax = old_pmax;
  
  int cycleCount = 0;
  int maxCycles = 3;
  while (cycleCount < maxCycles) 
  {
    pmin = optimize_low_cut ( old_pmax );
#ifdef TRACE_PRANGE_OPTIMIZATION
    std::cout << "Low momentum cut = " << pmin << ";  ";
#endif
    pmax = optimize_high_cut( pmin );
#ifdef TRACE_PRANGE_OPTIMIZATION
    std::cout << "High momentum cut = " << pmax << "\n";
#endif
    if ( (std::fabs(pmin-old_pmin) < 0.01) && 
         (std::fabs(pmax-old_pmax) < 0.01) ) {
      break;
    }
    old_pmin = pmin;
    old_pmax = pmax;
  }
  lowCut  = pmin;
  highCut = pmax;
  double integratedSignalEfficiency = signal.integrate(pmin,pmax);
  double integratedBackground = L * background.integrate(pmin,pmax);
#ifdef TRACE_PRANGE_OPTIMIZATION
  std::cout << "Optimal pmin = "   << pmin 
            << "  Optimal pmax = " << pmax << "\n";
#endif
#ifdef TRACE_INTEGRATED_STRENGTHS
  std::cout << "Integrated signal efficiency = " << integratedSignalEfficiency
            << "    Background counts = " << integratedBackground << "\n";
#endif
  computeSensitivities (integratedSignalEfficiency, integratedBackground,
       discoveryThresholdN,
       discoveryThreshold, smoothedWorstCaseSensitivity, punziSensitivity,
       CL90Sensitivity, singleEventSensitivity);
   
  epsilon = integratedSignalEfficiency;
  B =  integratedBackground;
  
      os    << "(optimized) integratedSignalEfficiency = " << epsilon << "\n"
            << "(optimized) integratedBackground = " << B << "\n"
            << "Smoothed Punzi merit denominator(B) = " 
            << punziMeritDenominator(B) << "\n"
            << "fc90 merit denominator(B) = " << fc90MeritDenominator(B) 
            << "\n";

  if (mfunc == SmoothedPunziMeritFunction) {
    figOfM = epsilon / FofM::punziMeritDenominator (B);  // equation (2)
  } else {
    figOfM = epsilon / fc90MeritDenominator (B);
  }

  // figOfM = epsilon / meritDenominator(B); // equation (2)        

  fillHypotheticalVectors( maximumSignalCount,
                epsilon, B, discoveryThresholdN, 
                upperLimitBR, CLlowerBR, CLupperBR );
} // results (free-floating range)

#define TRACE_MERIT_DENOMINATOR

void FofM::results_fixed_highCut (int maximumSignalCount, 
                double &lowCut,  double highCut, 
                double &epsilon,
                double &B,
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
                std::ostream & os ) const
{
  double pmax = highCut;
  double pmin = optimize_low_cut ( pmax );
  lowCut = pmin;
  double integratedSignalEfficiency = signal.integrate(pmin,pmax);
  double integratedBackground = L * background.integrate(pmin,pmax);
  #ifdef TRACE_INTEGRATED_STRENGTHS
  std::cout << "Integrated signal efficiency = " << integratedSignalEfficiency
            << "    Background counts = " << integratedBackground << "\n";
  #endif
  computeSensitivities (integratedSignalEfficiency, integratedBackground,
       discoveryThresholdN,
       discoveryThreshold, smoothedWorstCaseSensitivity, punziSensitivity,
       CL90Sensitivity, singleEventSensitivity);
  epsilon = integratedSignalEfficiency;
  B =  integratedBackground;
  
  #ifdef TRACE_MERIT_DENOMINATOR
       os   << "(optimized) integratedSignalEfficiency = " << epsilon << "\n"
            << "(optimized) integratedBackground = " << B << "\n"
            << "Smoothed Punzi merit denominator(B) = " 
            << punziMeritDenominator(B) << "\n"
            << "fc90 merit denominator(B) = " << fc90MeritDenominator(B) 
            << "\n";
  #endif
  
  if (mfunc == SmoothedPunziMeritFunction) {
    figOfM = epsilon / FofM::punziMeritDenominator (B);  // equation (2)
  } else {
    figOfM = epsilon / fc90MeritDenominator (B);
  }
  
  fillHypotheticalVectors( maximumSignalCount,
                epsilon, B, discoveryThresholdN, 
                upperLimitBR, CLlowerBR, CLupperBR );  
} // results (fixed high cut)

void FofM::results_fixed_lowCut (int maximumSignalCount, 
                double lowCut,  double &highCut, 
                double &epsilon,
                double &B,
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
                std::ostream & os ) const
{
  double pmin = lowCut;
  double pmax = optimize_high_cut ( pmin );
  highCut = pmax;
  double integratedSignalEfficiency = signal.integrate(pmin,pmax);
  double integratedBackground = L * background.integrate(pmin,pmax);
  #ifdef TRACE_INTEGRATED_STRENGTHS
  std::cout << "Integrated signal efficiency = " << integratedSignalEfficiency
            << "    Background counts = " << integratedBackground << "\n";
  #endif
  computeSensitivities (integratedSignalEfficiency, integratedBackground,
       discoveryThresholdN,
       discoveryThreshold, smoothedWorstCaseSensitivity, punziSensitivity,
       CL90Sensitivity, singleEventSensitivity);
  epsilon = integratedSignalEfficiency;
  B =  integratedBackground;
  
  #ifdef TRACE_MERIT_DENOMINATOR
       os   << "(optimized) integratedSignalEfficiency = " << epsilon << "\n"
            << "(optimized) integratedBackground = " << B << "\n"
            << "Smoothed Punzi merit denominator(B) = " 
            << punziMeritDenominator(B) << "\n"
            << "fc90 merit denominator(B) = " << fc90MeritDenominator(B) 
            << "\n";
  #endif
  
  if (mfunc == SmoothedPunziMeritFunction) {
    figOfM = epsilon / FofM::punziMeritDenominator (B);  // equation (2)
  } else {
    figOfM = epsilon / fc90MeritDenominator (B);
  }
  
  // figOfM = epsilon / meritDenominator(B); // equation (2)        
  fillHypotheticalVectors( maximumSignalCount,
                epsilon, B, discoveryThresholdN, 
                upperLimitBR, CLlowerBR, CLupperBR );  
} // results (fixed low cut)                

void FofM::results_fixed_cuts  (int maximumSignalCount, 
                double lowCut, double highCut, 
                double &epsilon,
                double &B,
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
                std::ostream & os ) const
{
  double pmin = lowCut;
  double pmax = highCut;
  double integratedSignalEfficiency = signal.integrate(pmin,pmax);
  double integratedBackground = L * background.integrate(pmin,pmax);
  #ifdef TRACE_INTEGRATED_STRENGTHS
  std::cout << "Integrated signal efficiency = " << integratedSignalEfficiency
            << "    Background counts = " << integratedBackground << "\n";
  #endif
  computeSensitivities (integratedSignalEfficiency, integratedBackground,
       discoveryThresholdN,
       discoveryThreshold, smoothedWorstCaseSensitivity, punziSensitivity,
       CL90Sensitivity, singleEventSensitivity);
   
  epsilon = integratedSignalEfficiency;
  B =  integratedBackground;
  
  #ifdef TRACE_MERIT_DENOMINATOR
       os   << "(optimized) integratedSignalEfficiency = " << epsilon << "\n"
            << "(optimized) integratedBackground = " << B << "\n"
            << "Smoothed Punzi merit denominator(B) = " 
            << punziMeritDenominator(B) << "\n"
            << "fc90 merit denominator(B) = " << fc90MeritDenominator(B) 
            << "\n";
  #endif
  
  if (mfunc == SmoothedPunziMeritFunction) {
    figOfM = epsilon / FofM::punziMeritDenominator (B);  // equation (2)
  } else {
    figOfM = epsilon / fc90MeritDenominator (B);
  }
  
  // figOfM = epsilon / meritDenominator(B); // equation (2)        

  fillHypotheticalVectors( maximumSignalCount,
                epsilon, B, discoveryThresholdN, 
                upperLimitBR, CLlowerBR, CLupperBR );
} // results (fixed momentum cuts)

double FofM::punziMeritDenominator ( double B )  {
  if (B < 0) return 3.148;
  double s = std::sqrt(B);
  return 3.148 + 3.538 * std::sqrt(s) + 5.390 * s;
}

void FofM::computeSensitivities( 
                double eps, double B,  
                int    & discoveryThresholdN,
                double &discoveryThreshold,
                double &smoothedWorstCaseSensitivity,
                double &punziSensitivity,
                double &CL90Sensitivity,
                double &singleEventSensitivity ) const
{
  double p5sigma = 1.0 - 2.86652e-7;
  discoveryThresholdN =  poissonInverseCDF(B, p5sigma);
  if (discoveryThresholdN < 5) discoveryThresholdN = 5;
  discoveryThreshold = (discoveryThresholdN-B)/(eps*L);
  smoothedWorstCaseSensitivity = punziMeritDenominator(B)/(eps*L);
  double p90pct = 0.90;
  double punziSensitivityMu = 
        poissonMeanSolver(discoveryThresholdN-1,1-p90pct);
  punziSensitivity = (punziSensitivityMu-B)/(eps*L);
//  FeldmanCousins90pctSensitivity fc90;
  CL90Sensitivity = fc90MeritDenominator(B)/(eps*L);
  singleEventSensitivity = 1.0/(eps*L);
  
  // To check - since the quoted BR limit is based on total counts, and that
  // includes B, it seems that we needed to subtract B from the count level
  // to translate to a BR.  The two that are not background-based, of course,
  // do not subtract B.  So we need to subtract B when computing the Punzi
  // sensitivity.  But the merit denominator formulat already subtracts B,
  // so the Punzi sensitivity is the only place where explicit subtraction is
  // needed.  Is this reasoning OK?
  
} // computeSensitivities

void FofM::fillHypotheticalVectors (int maximumSignalCount,
                double eps, double B, int discoveryThresholdN, 
                std::vector<double> & upperLimitBR,
                std::vector<double> & CLlowerBR,
                std::vector<double> & CLupperBR ) const
{
  double p90pct = 0.90;
  double p90pct_twoSided = .95;
  // Note - for hypothetical counts below the expected background,
  // this computation is NOT doing Cousins-Feldman.  Should be no problem:
  // background is probably going to be less than one count, so the only
  // issue is the upper limit we announce if we see 0 events, and that is
  // pretty cleanly the usual 90% CL number.
  upperLimitBR = std::vector<double>(maximumSignalCount+1);
  CLlowerBR    = std::vector<double>(maximumSignalCount+1);
  CLupperBR    = std::vector<double>(maximumSignalCount+1);
  // Up to n = discoveryThresholdN, we would state an upper limit on 
  // Branching Ratio.
  //#define TRACE_TABLE_CREATION
  for (int i=0;  i < discoveryThresholdN; ++i) {
    if (i > maximumSignalCount) break;
    upperLimitBR[i] = (poissonMeanSolver(i, 1-p90pct)-B)/(eps*L);
    #ifdef TRACE_TABLE_CREATION
    std::cout << "n = " << i 
              << "  upperLimitBR:  poissonMeanSolver(" << i << ", " 
              << 1-p90pct << ") - B";
    std::cout <<  " = "  << L*eps*upperLimitBR[i] 
              << " --> " << upperLimitBR[i] << "\n";
    #endif
  }
  // Past n = discoveryThresholdN, we sould state a 90% CL range for 
  // an announced Branching Ratio.
  if ( discoveryThresholdN > maximumSignalCount ) return;
  for (int i=discoveryThresholdN; i <= maximumSignalCount; ++i) {
    CLlowerBR[i] = std::max(0.0,
                (poissonMeanSolver(i, p90pct_twoSided)-B)/(eps*L));
    #ifdef TRACE_TABLE_CREATION
    std::cout << "n = " << i 
              << "  CLlowerBR:  poissonMeanSolver(" << i << ", " 
              << p90pct_twoSided << ") - B";
    std::cout <<  " = "  << L*eps*CLlowerBR[i] 
              << " --> " << CLlowerBR[i] << "\n";
    #endif
    if ( CLlowerBR[i] == 0 ) {
      CLupperBR[i] = (poissonMeanSolver(i, 1-p90pct)-B)/(eps*L);
      #ifdef TRACE_TABLE_CREATION
      std::cout << "n = " << i 
              << "  CLupperBR:  poissonMeanSolver(" << i << ", " 
              << 1-p90pct << ") - B";
      std::cout <<  " = "  << L*eps*CLupperBR[i] 
                << " --> " << CLupperBR[i] << "\n";
      #endif
    } else {
      CLupperBR[i] = (poissonMeanSolver(i, 1-p90pct_twoSided)-B)/(eps*L);
      #ifdef TRACE_TABLE_CREATION
      std::cout << "n = " << i 
              << "  CLupperBR:  poissonMeanSolver(" << i << ", " 
              << 1-p90pct_twoSided << ") - B";
      std::cout <<  " = "  << L*eps*CLupperBR[i] 
                << " --> " << CLupperBR[i] << "\n";
      #endif
    }
  }
}  // fillHypotheticalVectors             

std::string FofM::tables (int maximumSignalCount, 
             std::ostream & os, FofM::Summary & summary) const
{
  double lowCut;
  double highCut; 
  double epsilon;
  double B;
  double figOfM;
  int    discoveryThresholdN;
  double discoveryThreshold;
  double smoothedWorstCaseSensitivity;
  double punziSensitivity;
  double CL90sensitivity;
  double singleEventSensitivity;
  std::vector<double>  upperLimitBR;
  std::vector<double>  CLlowerBR;
  std::vector<double>  CLupperBR;
  results (maximumSignalCount,
           lowCut, highCut, epsilon, B, figOfM, 
           discoveryThresholdN, discoveryThreshold,
           smoothedWorstCaseSensitivity,punziSensitivity,
           CL90sensitivity, singleEventSensitivity,
           upperLimitBR, CLlowerBR, CLupperBR, os );
  summary.singleEventSensitivity = singleEventSensitivity;
  summary.CL90sensitivity = CL90sensitivity;
  summary.smoothedPunziSensitivity = smoothedWorstCaseSensitivity;
  summary.pCutLo = lowCut;
  summary.pCutHi = highCut;
  summary.figureOfMerit = figOfM;
  return tablesString (maximumSignalCount,
           lowCut, highCut, epsilon, B, figOfM, 
           discoveryThresholdN, discoveryThreshold,
           smoothedWorstCaseSensitivity,punziSensitivity,
           CL90sensitivity, singleEventSensitivity,
           upperLimitBR, CLlowerBR, CLupperBR,
           false, false);
} // tables (floating cuts)

std::string FofM::tables_fixed_highCut 
        (int maximumSignalCount, double highCut, 
             std::ostream & os, FofM::Summary & summary) const
{
  double lowCut;
  double epsilon;
  double B;
  double figOfM;
  int    discoveryThresholdN;
  double discoveryThreshold;
  double smoothedWorstCaseSensitivity;
  double punziSensitivity;
  double CL90sensitivity;
  double singleEventSensitivity;
  std::vector<double>  upperLimitBR;
  std::vector<double>  CLlowerBR;
  std::vector<double>  CLupperBR;
  results_fixed_highCut (maximumSignalCount,
           lowCut, highCut, epsilon, B, figOfM, 
           discoveryThresholdN, discoveryThreshold,
           smoothedWorstCaseSensitivity,punziSensitivity,
           CL90sensitivity, singleEventSensitivity,
           upperLimitBR, CLlowerBR, CLupperBR, os );
  summary.singleEventSensitivity = singleEventSensitivity;
  summary.CL90sensitivity = CL90sensitivity;
  summary.smoothedPunziSensitivity = smoothedWorstCaseSensitivity;
  summary.pCutLo = lowCut;
  summary.pCutHi = highCut;
  summary.figureOfMerit = figOfM;
  return tablesString (maximumSignalCount,
           lowCut, highCut, epsilon, B, figOfM, 
           discoveryThresholdN, discoveryThreshold,
           smoothedWorstCaseSensitivity,punziSensitivity,
           CL90sensitivity, singleEventSensitivity,
           upperLimitBR, CLlowerBR, CLupperBR,
           false, true);
} // tables (fixed high cut)

std::string FofM::tables_fixed_lowCut 
        (int maximumSignalCount, double lowCut, 
             std::ostream & os, FofM::Summary & summary) const
{
  double highCut;
  double epsilon;
  double B;
  double figOfM;
  int    discoveryThresholdN;
  double discoveryThreshold;
  double smoothedWorstCaseSensitivity;
  double punziSensitivity;
  double CL90sensitivity;
  double singleEventSensitivity;
  std::vector<double>  upperLimitBR;
  std::vector<double>  CLlowerBR;
  std::vector<double>  CLupperBR;
  results_fixed_lowCut (maximumSignalCount,
           lowCut, highCut, epsilon, B, figOfM, 
           discoveryThresholdN, discoveryThreshold,
           smoothedWorstCaseSensitivity,punziSensitivity,
           CL90sensitivity, singleEventSensitivity,
           upperLimitBR, CLlowerBR, CLupperBR, os );
  summary.singleEventSensitivity = singleEventSensitivity;
  summary.CL90sensitivity = CL90sensitivity;
  summary.smoothedPunziSensitivity = smoothedWorstCaseSensitivity;
  summary.pCutLo = lowCut;
  summary.pCutHi = highCut;
  summary.figureOfMerit = figOfM;
  return tablesString (maximumSignalCount,
           lowCut, highCut, epsilon, B, figOfM, 
           discoveryThresholdN, discoveryThreshold,
           smoothedWorstCaseSensitivity,punziSensitivity,
           CL90sensitivity, singleEventSensitivity,
           upperLimitBR, CLlowerBR, CLupperBR,
           true, false );
} // tables (fixed low cut)

std::string FofM::tables_fixed_cuts (int maximumSignalCount, 
                  double momcutLow, double momcutHigh, 
             std::ostream & os, FofM::Summary & summary) const
{
  double lowCut  = momcutLow;
  double highCut = momcutHigh; 
  double epsilon;
  double B;
  double figOfM;
  int    discoveryThresholdN;
  double discoveryThreshold;
  double smoothedWorstCaseSensitivity;
  double punziSensitivity;
  double CL90sensitivity;
  double singleEventSensitivity;
  std::vector<double>  upperLimitBR;
  std::vector<double>  CLlowerBR;
  std::vector<double>  CLupperBR;
  results_fixed_cuts (maximumSignalCount,
           lowCut, highCut, epsilon, B, figOfM, 
           discoveryThresholdN, discoveryThreshold,
           smoothedWorstCaseSensitivity,punziSensitivity,
           CL90sensitivity, singleEventSensitivity,
           upperLimitBR, CLlowerBR, CLupperBR, os);
  summary.singleEventSensitivity = singleEventSensitivity;
  summary.CL90sensitivity = CL90sensitivity;
  summary.smoothedPunziSensitivity = smoothedWorstCaseSensitivity;
  summary.pCutLo = lowCut;
  summary.pCutHi = highCut;
  summary.figureOfMerit = figOfM;
  return tablesString (maximumSignalCount,
           momcutLow, momcutHigh, epsilon, B, figOfM, 
           discoveryThresholdN, discoveryThreshold,
           smoothedWorstCaseSensitivity,punziSensitivity,
           CL90sensitivity, singleEventSensitivity,
           upperLimitBR, CLlowerBR, CLupperBR,
           true, true);
} // tables (fixed cuts)
           
std::string FofM::tablesString (int maximumSignalCount,
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
                bool lowCutFixed, bool highCutFixed) const
{
  std::ostringstream t;
  t << "\nSensitivity tables: All branching ratios given  are * 10^{-16} \n\n";
  unsigned int prec = t.precision(6);
  double displayScale = 1.0e16;
  if ( (!lowCutFixed) && (!highCutFixed) ) {
    t << "(optimized) momentum range  " << lowCut 
      << " --- " << highCut << "\n";
  } else if (lowCutFixed && !(highCutFixed) )  {
    t << "momentum range  " << lowCut 
      << " (fixed) --- " << highCut << " (optimized)\n";
  } else if (highCutFixed && !(lowCutFixed)) {
    t << "momentum range  " << lowCut 
      << " (optimized) --- " << highCut << " (fixed)\n";
  } else {
    t << "momentum range  " << lowCut 
      << " (fixed) --- " << highCut << " (fixed)\n";
  }
  if ( (!lowCutFixed) || (!highCutFixed) ) {
    t << "optimization based on ";
    if (mfunc == SmoothedPunziMeritFunction) {
      t << "smoothed worst-case sensitivity \n";
    } else {
      t << "90% CL sensitivity \n";
    }
  }
    
  t << "Signal Efficiency with those cuts:  " << epsilon << "\n";
  t << "Background mean with those cuts:  " << B 
    << " (" << L << " protons on target)\n";
  t << "Single event sensitivity (cutting at stated momentum range): " 
    << singleEventSensitivity * displayScale << "\n";
  t << "90% CL sensitivity (cutting at stated momentum range)      : " 
    << CL90Sensitivity * displayScale << "\n";
  t << "BR sensitivity by Punzi's worst-case definition           : " 
    << punziSensitivity * displayScale << "\n";
  t << "Smoothed worst-case BR sensitivity                         : " 
    << smoothedWorstCaseSensitivity * displayScale << "\n";

  if (mfunc == SmoothedPunziMeritFunction) {
    t << "\nFigure of merit (based on 90% CL sensitivity) =              "
      << figOfM * smoothedWorstCaseSensitivity/CL90Sensitivity << "\n";     
    t << "\nFigure of merit (based on smoothed worst-case sensitivity) = "
      << figOfM << "\n\n";     
  } else {
    t << "\nFigure of merit (based on 90% CL sensitivity) =              "
      << figOfM << "\n";     
    t << "\nFigure of merit (based on smoothed worst-case sensitivity) = "
      << figOfM * CL90Sensitivity/smoothedWorstCaseSensitivity << "\n\n";     
  }  
  
  t << "(5-sigma) Discovery threshold: " << discoveryThresholdN << "\n";
  t << "BR quoted mean if minimum counts for a discovery are seen  : " 
    << discoveryThreshold * displayScale << "\n";

  t << "90% CL upper limits for BR if below discovery threshold: \n"
    << " N    BR (* 10^{16}) \n";
  t.precision(4);
  for (int i=0;  i < discoveryThresholdN; ++i) {
    if (i > maximumSignalCount) break;
    if (i < 10) {
      t << " " << i;
    } else {
      t << i;
    } 
    t << ":  " <<  upperLimitBR[i]  * displayScale << "\n";
  }
  t << "90% CL (two-side) range for BR (* 10^{16}) "
    << "if above discovery threshold: \n"
    << " N    BR Lower Limit    BR Upper Limit (* 10^{16}) \n";
  for (int i=discoveryThresholdN; i <= maximumSignalCount; ++i) {
    if (i < 10) {
      t << " " << i;
    } else {
      t << i;
    } 
    t << ":    " <<  std::setw(5) << CLlowerBR[i]  * displayScale 
      << "              " << CLupperBR[i]  * displayScale << "\n";
  }
  t.precision (prec);
  return t.str();

} // tables (both cuts fixed)                                


double FofM::discoveryCount ( double B )  {
  
  if ( B > 20 ) return B + 5.0 * std::sqrt(B);
  // The above is the Gaussian approximation for a Poisson distribution with 
  // "large" mean, which gives sllightly worse p-values than the true
  // Gaussian 5-sigma.  Below 20, we use a fit (below) to the actual needed 
  // count to get the actual correct p-value; this fit is more accurate.
  // There is a bit of discontinuity here which we can fix later.
  // But this is not urgent; the experiment lives in a region where B << 20.
  return 0.759*B + 6.054* std::sqrt(B) + 2.764;
}

double FofM::optimize_low_cut  (double fixed_high_cut) const 
{
  splines::Spline<1> negIntegratedSignal = 
        signal.indefiniteIntegral(fixed_high_cut);

#ifdef BACKGROUND_TRACE
  std::cout << "\nBackground spline about to be used at high cut of " 
            << fixed_high_cut << ": \n\n";
  for (unsigned int i=0; i< grid.nPoints(); ++i) {
    if (i%1 == 0) {
      double x = grid[i];
      std::cout << x << "    " << background(x) << "\n";
    }
  }
#endif
  splines::Spline<1> negIntegratedBackground = 
        background.indefiniteIntegral(fixed_high_cut);
  // We are about to form a spline representing the figure of merit.
  // However, this is pathological if a >= b.  In fact, we would never
  // consider a range of less than, say 0.1.  Since Spline<1> does not
  // have a function maximum(range),we will need to prepare a different
  // grid for the new spline.
  double minimumAcceptanceRange = 0.1;
  double gridMin     = grid[0];
  double gridSpacing = grid[1] - gridMin;
  double aMax = fixed_high_cut - minimumAcceptanceRange;
  int nPointsInOptimizingGrid = (aMax - gridMin)/gridSpacing;
  double gridMax = gridMin + nPointsInOptimizingGrid * gridSpacing;
#ifdef OPTIMIZATION_OUTPUT
  std::cout <<  "aMax = " << aMax
            <<  "  nPointsInOptimizingGrid = " << nPointsInOptimizingGrid
            <<  "  gridMax = " << gridMax << "\n";
#endif
  
  splines::Grid<1> optimizingGrid( gridMin, gridMax, nPointsInOptimizingGrid );
  FofMwithLowerLimitVarying f(negIntegratedSignal, negIntegratedBackground, 
                                L, mfunc);
  // mfunc is either SmoothedPunziMeritFunction or FCsensitivityMeritFunction
#ifdef ANALYSIS_OUTPUT
  std::cout << "\nnegIntegratedSignal at high cut of " 
            << fixed_high_cut << ": \n\n";
  for (unsigned int i=0; i< optimizingGrid.nPoints(); ++i) {
    if (i%5 == 0) {
      double x = optimizingGrid[i];
      std::cout << x << "    " << negIntegratedSignal(x) << "\n";
    }
  }
 
  std::cout << "\nnegIntegratedBackground at high cut of " 
            << fixed_high_cut << ": \n\n";
  for (unsigned int i=0; i< optimizingGrid.nPoints(); ++i) {
    if (i%5 == 0) {
      double x = optimizingGrid[i];
      std::cout << x << "    " << negIntegratedBackground(x) << "\n";
    }
  }
 
  std::cout << "\nMerit function used in optimize_low_cut: \n\n";
  for (unsigned int i=0; i< optimizingGrid.nPoints(); ++i) {
    if (i%5 == 0) {
      double x = optimizingGrid[i];
      std::cout << x << "    " << f(x) << "\n";
    }
  }
#endif
  
  splines::Spline<1> figureOfMeritSpline ( f, optimizingGrid );
#ifdef ANALYSIS_OUTPUT
  std::cout << "\nfigureOfMeritSpline for fixed high cut of "
            <<  fixed_high_cut << ": \n\n";
  for (unsigned int i=0; i< optimizingGrid.nPoints(); ++i) {
    if (i%5 == 0) {
      double x = optimizingGrid[i];
      std::cout << x << "    " << figureOfMeritSpline(x) << "\n";
    }
  }
#endif
  double maximum =  figureOfMeritSpline.maximum();  
  // TODO - There may be a much quicker way to do this, if we modify the
  //        splines package to add the right capabilities.  That is a 
  //        performance and perhaps a code cleanup optimization for later.
  return maximum;
}

double FofM::optimize_high_cut (double fixed_low_cut) const 
{
  splines::Spline<1> integratedSignal = 
        signal.indefiniteIntegral(fixed_low_cut);
  splines::Spline<1> integratedBackground = 
        background.indefiniteIntegral(fixed_low_cut);
  // We are about to form a spline representing the figure of merit.
  // However, this is pathological if b <= a.  In fact, we would never
  // consider a range of less than, say 0.1.  Since Spline<1> does not
  // have a function maximum(range),we will need to prepare a different
  // grid for the new spline.
  double minimumAcceptanceRange = 0.1;
  unsigned int nGridPoints = grid.nPoints();
  double gridMax     = grid[nGridPoints-1];
  double gridSpacing = gridMax - grid[nGridPoints-2];
  double bMin = fixed_low_cut + minimumAcceptanceRange;
  int nPointsInOptimizingGrid = (gridMax - bMin)/gridSpacing;
  double gridMin = gridMax - nPointsInOptimizingGrid * gridSpacing;
#ifdef OPTIMIZATION_OUTPUT
  std::cout <<  "bMin = " << bMin
            <<  "  nPointsInOptimizingGrid = " << nPointsInOptimizingGrid
            <<  "  gridMin = " << gridMin << "\n";
#endif
  splines::Grid<1> optimizingGrid( gridMin, gridMax, nPointsInOptimizingGrid );
  FofMwithUpperLimitVarying f(integratedSignal, integratedBackground, 
                                L, mfunc);

#ifdef ANALYSIS_OUTPUT
  std::cout << "\nintegratedSignal at low cut of " 
            << fixed_low_cut << ": \n\n";
  for (unsigned int i=0; i< optimizingGrid.nPoints(); ++i) {
    if (i%5 == 0) {
      double x = optimizingGrid[i];
      std::cout << x << "    " << integratedSignal(x) << "\n";
    }
  }
 
   std::cout << "\nintegratedBackground at low cut of " 
            << fixed_low_cut << ": \n\n";
 for (unsigned int i=0; i< optimizingGrid.nPoints(); ++i) {
    if (i%5 == 0) {
      double x = optimizingGrid[i];
      std::cout << x << "    " << integratedBackground(x) << "\n";
    }
  }

  std::cout << "\nMerit function used in optimize_high_cut: \n\n";
  for (unsigned int i=0; i< optimizingGrid.nPoints(); ++i) {
    if (i%5 == 0) {
      double x = optimizingGrid[i];
      std::cout << x << "    " << f(x) << "\n";
    }
  }
#endif
  
  splines::Spline<1> figureOfMeritSpline ( f, optimizingGrid );
#ifdef ANALYSIS_OUTPUT
  std::cout << "\nfigureOfMeritSpline for fixed high cut of "
            <<  fixed_high_cut << ": \n\n";
  for (unsigned int i=0; i< optimizingGrid.nPoints(); ++i) {
    if (i%5 == 0) {
      double x = optimizingGrid[i];
      std::cout << x << "    " << figureOfMeritSpline(x) << "\n";
    }
  }
#endif
  double maximum =  figureOfMeritSpline.maximum();  
  return maximum;
}


// -----------------------------
// Helper Functions and Functors
// -----------------------------

double FofMwithLowerLimitVarying::operator()(double a) const
{
  double eps = - negIntegratedSignal(a);
  double B   = - L * negIntegratedBackground(a);
  double f;
  if (mfunc == SmoothedPunziMeritFunction) {
    f = eps / FofM::punziMeritDenominator (B);
  } else {
    f = eps / fc90MeritDenominator (B);
  }
  return f;
}

double FofMwithUpperLimitVarying::operator()(double a) const
{
  double eps = integratedSignal(a);
  double B   = L * integratedBackground(a);
  double f;
  if (mfunc == SmoothedPunziMeritFunction) {
    f = eps / FofM::punziMeritDenominator (B);
  } else {
    f = eps / fc90MeritDenominator (B);
  }
  return f;
}

void FofM::displayBackground(double pmin, double pmax, int n)
{
  std::cout << "Background in FofM:\n\n";
  double step = (pmax-pmin)/n;
  double roughIntegral = 0.0;
  for (int i=0; i < n; ++i) {
    double bk = background(pmin + i*step);
    std::cout << "p = " << pmin + i*step << "  background = " << bk << "\n";
    roughIntegral += bk*step;
  } 
  std::cout << "Rough background integral is " <<  roughIntegral << "\n";
}

// --------------------
// Pre-packaged Spectra
// -------------------

NullBackgroundSpectrum::NullBackgroundSpectrum()  
  : FofM::Spectrum(zeroFunctionForNullBackgroundSpectrum, 
                   FofM::defaultSpectrumLowerLimit(), 
                   FofM::defaultSpectrumUpperLimit()) 
{}

double branchingFractionForCapturedMuonsToDIO()
{
  // TODO -- get this from the conditions database  
  return 0.41;
} // branchingFractionForCapturedMuonsToDIO()


CzarneckiDIOspectrumFunction::CzarneckiDIOspectrumFunction (double strength)
  : multiplier(strength)
{
  // TODO -- obtain alpha_5_ thru alpha_8 and Emax from conditions database
  //         This is why we are setting these, temproarily, in the ctor body
  //         rather than in an initializer list.

  // The following coefficients come from equation (26) in Czarnecki,
  // arXiv:1106.4756v1 [hep-ph] 23 Jun 2011
  alpha_5 =  8.6434e-17;
  alpha_6 =  1.16874e-17;
  alpha_7 = -1.87828e-19;
  alpha_8 =  9.16327e-20;
  Emax    =  104.973;
  
  // TODO -- PERHAPS obtain E_mu (= m_u - binding energy) and M_Al 
  //         from conditions DB
  
  E_mu = 105.194;
  M_Al = 25133;
}

double CzarneckiDIOspectrumFunction::operator() (double p)
{
  double delta = E_mu - p - p*p/(2.0*M_Al); 
  if (delta < 0) return 0;
  double delta5 = delta*delta*delta*delta*delta;
  double delta6 = delta*delta5;
  double delta7 = delta*delta6;
  double delta8 = delta*delta7;
  double czarnecki_equation_25 =
         alpha_5*delta5 + alpha_6*delta6 + alpha_7*delta7 + alpha_8*delta8;
  return multiplier * czarnecki_equation_25;
}

CzarneckiDIOSpectrum::CzarneckiDIOSpectrum 
                (double strength, double a, double b)
  : FofM::Spectrum( CzarneckiDIOspectrumFunction
                    (strength*branchingFractionForCapturedMuonsToDIO()) )
{}         

} // end of namespace mu2e
