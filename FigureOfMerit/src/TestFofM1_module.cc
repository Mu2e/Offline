//
//  The TestFofM1 plugin; the first example of a module.
//
//  $Id: TestFofM1_module.cc,v 1.1 2012/03/12 18:23:25 mf Exp $
//  $Author: mf $
//  $Date: 2012/03/12 18:23:25 $
//
//  Original author Rob Kutschke
//

// C++ includes.
#include <iostream>
#include <vector>
#include <string>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/TFileService.h"

// mu2e includes
#include "FigureOfMerit/inc/FofM.h"
#include "FigureOfMerit/inc/PoissonCDFsolver.h"

// Root includes.
#include "TFile.h"
#include "TH1F.h"
#include "TNtuple.h"
#include "TChain.h"

// Helper functions

namespace mu2e {

void extractFitmom 
        ( TTree * tracks
        , int    mimimum_nactive  
        , float  maximum_t0err
        , float  maximum_fitmomerr
        , float  minimum_fitcon   
        , double lowestBin 
        , double topOfLastBin              
        , unsigned int nBins
        , std::vector<double> & counts
        )
{
    // Obtain the momentum data as a vector of momenta
    std::vector<double> momenta;
    tracks->SetBranchStyle(0);
    int   nactive;   tracks->SetBranchAddress("nactive"  , &nactive);
    float t0err;     tracks->SetBranchAddress("t0err"    , &t0err);    
    float fitmomerr; tracks->SetBranchAddress("fitmomerr", &fitmomerr);    
    float fitcon;    tracks->SetBranchAddress("fitcon"   , &fitcon);    
    int   fitstatus; tracks->SetBranchAddress("fitstatus", &fitstatus);    
    float fitmom;    tracks->SetBranchAddress("fitmom"   , &fitmom);    
    int ntracks = 0;
    for (int i = 0; true; ++i) {
      int bytesRead = tracks->GetEntry(i);
      if (bytesRead == 0) break;
      ++ntracks;
      if  (    (fitstatus == 1) 
            && (nactive   >= mimimum_nactive)
            && (t0err     <= maximum_t0err)
            && (fitmomerr <= maximum_fitmomerr)
            && (fitcon    >= minimum_fitcon) 
            && (fitmom    < 130.0) // deals with the one really screwy track    
          ) {
        momenta.push_back( fitmom );           
      }     
    } 
    std::cout << "There are " << ntracks << " entries in the tree\n";
    std::cout << "There are " << momenta.size() 
              << " entries passing the cuts\n";
    double binSize = (topOfLastBin - lowestBin)/nBins;
    int nBinnedTracks = 0;
    for (unsigned int i=0; i < momenta.size(); ++i) {
      double p = momenta[i];
      if (p >= lowestBin && p < topOfLastBin) {
        unsigned int binNumber = (p-lowestBin)/binSize;
        if (binNumber >= nBins) binNumber = nBins-1;
        counts[binNumber] += 1.0; // everything has equal weight
        ++nBinnedTracks;
      }
    }
    std::cout << nBinnedTracks << " tracks between " << lowestBin << " and "
              << topOfLastBin << "\n";
} // extractFitmom

void extractFitmom 
        ( std::string const & fileName
        , int    mimimum_nactive  
        , float  maximum_t0err
        , float  maximum_fitmomerr
        , float  minimum_fitcon   
        , double lowestBin 
        , double topOfLastBin              // top of last bin
        , unsigned int nBins
        , std::vector<double> & counts
        )
{
    std::cout << "extracting fitmom from " << fileName << "\n";
    // Obtain the momentum data as a vector of momenta
    TFile* tfp    = new TFile (fileName.c_str());
    TTree* tracks = (TTree*)tfp->Get ("ReadKalFits/trkdiag");
    extractFitmom 
        (  tracks
        ,  mimimum_nactive  
        ,  maximum_t0err
        ,  maximum_fitmomerr
        ,  minimum_fitcon   
        ,  lowestBin 
        ,  topOfLastBin         
        ,  nBins
        ,  counts );
    delete tracks;   tracks = 0;  // tracks would delete when tfp does,
                                  // but its delete informs the TFile so
                                  // this is OK as well, and may be saner
    delete tfp;      tfp = 0;
}  

void extractFitmom 
        ( std::vector<std::string> const & listOfFileNames
        , int    mimimum_nactive  
        , float  maximum_t0err
        , float  maximum_fitmomerr
        , float  minimum_fitcon   
        , double lowestBin 
        , double topOfLastBin             
        , unsigned int nBins
        , std::vector<double> & counts
        )
{
    std::cout << "extracting fitmom from a chain of files \n";
    // Obtain the momentum data as a vector of momenta
   for (unsigned int i=0; i<listOfFileNames.size(); ++i) {
      std::cout << "listOfFileNames["  << i << "].c_str() = "
                << listOfFileNames[i].c_str() << "\n";
      TFile* tfp    = new TFile (listOfFileNames[i].c_str());
      TTree* tracks = (TTree*)tfp->Get ("ReadKalFits/trkdiag");

      extractFitmom 
        (  tracks
        ,  mimimum_nactive  
        ,  maximum_t0err
        ,  maximum_fitmomerr
        ,  minimum_fitcon   
        ,  lowestBin 
        ,  topOfLastBin         
        ,  nBins
        ,  counts );
      delete tracks;   tracks = 0;
      delete tfp;      tfp = 0;
    }
} 


} // End of namespace mu2e for helper functions


using namespace std;

namespace mu2e {

  class TestFofM1 : public art::EDAnalyzer {

  public:
    explicit TestFofM1(fhicl::ParameterSet const& pset){}

    virtual void beginRun(art::Run const &);
    void analyze(const art::Event& event);

  private:

  };

  void TestFofM1::analyze(const art::Event& event){
    cout << "This is TestFofM1.  From analyze: "
         << event.id()
         << endl;
  }

  void TestFofM1::beginRun(art::Run const &){
    cout << "This is TestFofM1.  From BeginRun() "
         << endl;

    // Set filter levels for the CE data:  Cut set B
    int    mimimum_nactive   = 20;
    float  maximum_t0err     = 1.5; 
    float  maximum_fitmomerr = 0.2;
    float  minimum_fitcon    = 1.0e-4;
#ifdef OTHER_POSSIBLE_CHOICES
    // Set filter levels for the CE data:  Cut set A
    int    mimimum_nactive   = 20;
    float  maximum_t0err     = 10.0; 
    float  maximum_fitmomerr = 1.0;
    float  minimum_fitcon    = 1.0e-6;
    // Set filter levels for the CE data:  Cut set C
    int    mimimum_nactive   = 25;
    float  maximum_t0err     = 1.0; 
    float  maximum_fitmomerr = 0.18;
    float  minimum_fitcon    = 1.0e-3;
    // Set filter levels for the CE data:  Cut set D
    int    mimimum_nactive   = 30;
    float  maximum_t0err     = 0.9; 
    float  maximum_fitmomerr = 0.15;
    float  minimum_fitcon    = 1.0e-2;
#endif
    // Obtain the CE data as a vector of momenta
    double lowestBin = 102.5;  
    double topOfLastBin = 106.0;  // top of last bin
    unsigned int nBins = 350;
    double canonicalRangeLo = 103.51;
    double canonicalRangeHi = 104.70;
    double binSize = (topOfLastBin - lowestBin)/nBins;
    
    std::vector<double> sigEfficiency (nBins);
        extractFitmom  
        ( "/mu2e/data/users/mf/CEFilesMixed.root"
        , mimimum_nactive  
        , maximum_t0err
        , maximum_fitmomerr
        , minimum_fitcon   
        , lowestBin 
        , topOfLastBin              // top of last bin
        , nBins
        , sigEfficiency
        );
    // Normalization of the signal efficiencies
    // Note that our convention is to normalize the signal efficiency to a
    // per POT number.
    // Note also that only muons that are captured by the nucleus are 
    // candidates for direct conversion.  
    double protonsOnTarget = 3.6e20;  
    double numberOfGeneratedMuonConversions = 50000;
    double stoppedMuonsPerPOT = 0.0021;  // I have seen .00215 
    double capturedMuonsPerStoppedMuon = 0.609; // DocDB 48 - I have seen .59
    double liveGateFraction = .47;  // or .50 -- captures in time window in 48
    double POTsRepresented = numberOfGeneratedMuonConversions /
       ( stoppedMuonsPerPOT * capturedMuonsPerStoppedMuon * liveGateFraction );
    double sigEffCountNormalization =  1.0/POTsRepresented;
    double sigEffNormalization =  sigEffCountNormalization / binSize;
    std::cout << "Signal Efficiency: \n\n";
    unsigned int canonicalRangeCEcount = 0;
    for (unsigned int i=0; i < sigEfficiency.size(); ++i) {
      double p = lowestBin + binSize*i;
      if ( (p >= canonicalRangeLo) && (p <= canonicalRangeHi) ) {
       canonicalRangeCEcount += sigEfficiency[i];
      }
    }
    std::cout << "There are " <<  canonicalRangeCEcount
              << " CE's in the files in the range of "
              << canonicalRangeLo << " - "  << canonicalRangeHi << " \n";
    std::cout << "With normalization, this amounts to "
              << canonicalRangeCEcount * sigEffCountNormalization
                                       * protonsOnTarget
               << " CE events in the range if BR = 1\n";
    for (unsigned int i=0; i < sigEfficiency.size(); ++i) {
      sigEfficiency[i] *= sigEffNormalization;
//#define OUTPUT_FREQUENCIES
#ifdef OUTPUT_FREQUENCIES
      std::cout << lowestBin + binSize*i << ":   " 
                << sigEfficiency[i]/sigEffNormalization;
      std::cout << " --> " << sigEfficiency[i] 
                << " ==> " << sigEfficiency[i] * protonsOnTarget << "\n";
#endif
    }
    double integratedSignalEfficiencyInRange = 0;
    for (unsigned int i=0; i < sigEfficiency.size(); ++i) {
      double p = lowestBin + binSize*i ;
      if ( (p >= canonicalRangeLo) && (p <= canonicalRangeHi) ) {
        integratedSignalEfficiencyInRange += 
                sigEfficiency[i];
      } 
    }
    double CEcount = integratedSignalEfficiencyInRange * protonsOnTarget;
    std::cout << "CE count for BR = 1.0e-16 between 102.5 and 106.0 is " 
              << CEcount*1.0e-16 << "\n";

    // Obtain the DIO data as a vector of momenta 
    // Since there will be multiple files for this, we use a list.
    // (Because the TTree is placed under a directory, we cannot use
    // a TChain, we must use a list of TFiles.  Oh well...)
    std::vector<std::string> DIOfileList;
    DIOfileList.push_back(
                std::string("/mu2e/data/users/mf/DIOFilesMixed1tag.root") );
    DIOfileList.push_back(
                std::string("/mu2e/data/users/mf/DIOFilesMixed2tag.root") );
    std::vector<double> DIObackground (nBins);
        extractFitmom 
        ( DIOfileList
        , mimimum_nactive  
        , maximum_t0err
        , maximum_fitmomerr
        , minimum_fitcon   
        , lowestBin 
        , topOfLastBin              // top of last bin
        , nBins
        , DIObackground
        );
    // Normalization of the DIO background
    // Note that our convention is to normalize the background to a
    // number per proton on target, so the luminosity does not appear in 
    // the background normalization.  (But the luminosity is passed to the
    // FofM calculator because *it* needs a total number of background counts.)

    // All the DIO's in the files were created from a distribution of DIO's 
    // with p >= 102.5.
    double numberOfGeneratedDIOs = 194000;
    double DIOsPerStoppedMuon =  1.0 - capturedMuonsPerStoppedMuon;
    double fractionOfDIOsAbove102_5 = 4.42e-15;
    double DIObackgroundCountNormalization =
        ( stoppedMuonsPerPOT * DIOsPerStoppedMuon 
            * fractionOfDIOsAbove102_5  * liveGateFraction ) /
        ( numberOfGeneratedDIOs );

        // TEMPORARY TEST
        double DIOrescaler = 1.0;
        if ( DIOrescaler != 1.0 ) {
          std::cout << "\n**** For this trial, DIO normalization is scaled to "
                    << DIOrescaler << " times its actual value ****\n\n";
        }
        DIObackgroundCountNormalization *= DIOrescaler;
    
    double DIObackgroundNormalization =  
                DIObackgroundCountNormalization / binSize;  
    std::cout << "\nDIObackground: (Count normalization is " 
              << DIObackgroundCountNormalization << " ==> " 
              << DIObackgroundCountNormalization*protonsOnTarget << ")\n\n";
    unsigned int canonicalRangeDIOcount = 0;
    for (unsigned int i=0; i < DIObackground.size(); ++i) {
      double p = lowestBin + binSize*i;
      if ( p >= 103.51 && p <= 104.70 )  {
        canonicalRangeDIOcount += DIObackground[i];
      }
    }
    std::cout << "There are " <<  canonicalRangeDIOcount
              << " DIO's in the files in the range of "
              << canonicalRangeLo << " - "  << canonicalRangeHi << " \n";
    std::cout << "With normalization, this amounts to "
              << canonicalRangeDIOcount * DIObackgroundCountNormalization 
                                        * protonsOnTarget 
               << " DIO events in the range \n";
    for (unsigned int i=0; i < DIObackground.size(); ++i) {
      DIObackground[i] *= DIObackgroundNormalization;
#ifdef OUTPUT_FREQUENCIES
      std::cout << lowestBin + binSize*i << ":   "
                << DIObackground[i] / DIObackgroundCountNormalization;
      std::cout << " --> " << DIObackground[i] 
                << " ==> " << DIObackground[i] * protonsOnTarget << "\n";
#endif
    }


    // Obtain the RPC data as a vector of momenta
    // Since there will be multiple files for this, we use a list.
    std::vector<std::string> RPCfileList;
    RPCfileList.push_back(
                std::string("/mu2e/data/users/mf/RPCFilesMixed1tag.root") );
    RPCfileList.push_back(
                std::string("/mu2e/data/users/mf/RPCFilesMixed2tag.root") );
    std::vector<double> RPCbackground (nBins);
        extractFitmom 
        ( RPCfileList
        , mimimum_nactive  
        , maximum_t0err
        , maximum_fitmomerr
        , minimum_fitcon   
        , lowestBin 
        , topOfLastBin              // top of last bin
        , nBins
        , RPCbackground
        );
    // Normalization of the RPCbackground
    // Note that our convention is to normalize the background to a
    // number per proton on target, so the luminosity does not appear in 
    // the background normalization.
    
//    double suppliedRPCnormalization =  9.45e-5;   
    double suppliedRPCnormalization =  2.0 * 9.45e-5;   
        // Explanation:  This is the number to multiply the counts in the 
        // files by, to get a number of RPC's per 36.e20 POTs.  It includes
        // the .021 RPC's per stopped pion, and the live gate fraction,
        // which is very small because most stopped pions are in the first
        // 200 nsec, and the lifetime is only 21 nsec.    
        //
        // Since we want a number per POT, we need to divide that back out.  
        //
        // The number 9.45e-5 is based on just the processes Bob B. had
        // simulated to get these files.  But in fact there are other processes;
        // an email has said to just multiply the normalizatoin by 2, which
        // we do above.
    double RPCbackgroundCountNormalization =
        ( suppliedRPCnormalization ) /
        ( protonsOnTarget );

        // TEMPORARY TEST
        double RPCrescaler = 1.0;
        // std::cout << "Enter RPC rescaling factor: ";
        // std::cin  >> RPCrescaler;
        if ( RPCrescaler != 1.0 ) {
          std::cout << "\n**** For this trial, RPC normalization is scaled to "
                    << RPCrescaler << " times its actual value ****\n\n";
        }
        RPCbackgroundCountNormalization *= RPCrescaler;

    double RPCbackgroundNormalization =  
                RPCbackgroundCountNormalization / binSize;  

    std::cout << "\nRPCbackground: (count normalization is " 
              << RPCbackgroundCountNormalization << " ==> " 
              << RPCbackgroundCountNormalization*protonsOnTarget << ")\n\n";
    unsigned int canonicalRangeRPCcount = 0;
    for (unsigned int i=0; i < RPCbackground.size(); ++i) {
      double p = lowestBin + binSize*i;
      if ( (p >= canonicalRangeLo) && (p <= canonicalRangeHi) ) {
        canonicalRangeRPCcount += RPCbackground[i];
      }
    }
    std::cout << "There are " <<  canonicalRangeRPCcount
              << " RPC's in the files in the range of "
              << canonicalRangeLo << " - "  << canonicalRangeHi << " \n";
    std::cout << "With normalization, this amounts to "
              << canonicalRangeRPCcount * RPCbackgroundCountNormalization 
                                        * protonsOnTarget
              << " RPC events in the range \n";
    for (unsigned int i=0; i < RPCbackground.size(); ++i) {
      RPCbackground[i] *= RPCbackgroundNormalization;
#ifdef OUTPUT_FREQUENCIES
      std::cout << lowestBin + binSize*i << ":   " 
                << RPCbackground[i]/RPCbackgroundCountNormalization;
      std::cout << " --> " << RPCbackground[i] 
                << " ==> " << RPCbackground[i] * protonsOnTarget << "\n";
#endif
    }

    
    //
    // Use the FofM class to calculate quantities of interest
    //
    
    // Since our vectors represent values at the middle of intervals,
    // we should state our sample points vectors as such when using
    // them to construct spectra
    double midBinOffset = (topOfLastBin - lowestBin)/(2.0*nBins);
    double lowestPoint = lowestBin + midBinOffset;
    double endingPoint = topOfLastBin - midBinOffset;

    std::cout << "\n------------------------------------ \n\n"
              << "Using the FofM tool \n\n";
 

    std::cout << "lowestBin = " << lowestBin
              << "  topOfLastBin = " << topOfLastBin 
              << "\nlowestPoint = " << lowestPoint
              << "  endingPoint = " << endingPoint
              << "  nBins = " << nBins << "\n\n";    
    
    FofM::Spectrum CE_spectrum (sigEfficiency, lowestPoint, endingPoint);
    FofM::Spectrum DIO_spectrum (DIObackground, lowestPoint, endingPoint);
    FofM figureOfMeritCalculator ( CE_spectrum, 
                                   DIO_spectrum, 
                                   protonsOnTarget );  
//#define EXAMINESPECTRUM
#ifdef EXAMINESPECTRUM
    figureOfMeritCalculator.displayBackground(lowestBin, topOfLastBin, 100);
    splines::Spline<1> sbkg = figureOfMeritCalculator.getBackground();
    splines::Grid<1> sgrid = figureOfMeritCalculator.getGrid();
    std::cout << " L * DIO Background (normalized) vs L * Background \n";
    double roughInt = 0.0;
    double gridstep = sgrid[1] - sgrid[0];
    for (unsigned int i = 0; i < sgrid.nPoints(); ++i) {
      double p =  sgrid[i];
      std::cout << i << ": p = " << p  << "  " 
                << DIObackground[i] / DIObackgroundNormalization << "  "
                << protonsOnTarget * DIObackground[i] << "  "
                << protonsOnTarget * sbkg(p) << "\n";
      if ( (p >= canonicalRangeLo) && (p <= canonicalRangeHi) )  {
        roughInt += sbkg(p);
      }
    } 
    std::cout << "Rough integral in canonical range based on sbkg = " 
              << protonsOnTarget*gridstep*roughInt << "\n";
#endif
 
    FofM::Spectrum RPC_spectrum (RPCbackground, lowestPoint, endingPoint);
    figureOfMeritCalculator.addBackground (RPC_spectrum);
#ifdef EXAMINESPECTRUM
    figureOfMeritCalculator.displayBackground(lowestBin, topOfLastBin, 100);
    sbkg = figureOfMeritCalculator.getBackground();
    std::cout << " L * (DIO+RPC) Background (normalized) vs L * Background \n";
    roughInt = 0.0;
    for (unsigned int i = 0; i < sgrid.nPoints(); ++i) {
      double p =  sgrid[i];
      std::cout << i << ": p = " << p  << "  " 
                << DIObackground[i] / DIObackgroundNormalization 
                +  RPCbackground[i] / RPCbackgroundNormalization << "  "
                << protonsOnTarget * 
                        (DIObackground[i] + RPCbackground[i]) << "  "
                << protonsOnTarget * sbkg(p) << "\n";
      if ( (p >= canonicalRangeLo) && (p <= canonicalRangeHi) )  {
        roughInt += sbkg(p);
      }
    } 
    std::cout << "DIO + RPC Rough integral in canonical range based on sbkg = " 
              << protonsOnTarget*gridstep*roughInt << "\n";
#endif
    
    int maximumSignalCount = 25; // print hypothetical BR ranges for up to 
                                 // 25 accepted signal counts
    double lowCut;
    double highCut;

#ifdef EXERCISE_RESULTS_METHOD
    double sigEff;
    double bkgd;
    double figOfM;
    int    discoveryThresholdN;
    double discoveryThreshold;
    double smoothedWorstCaseSensitivity;
    double punziSensitivity;
    double CL90Sensitivity;
    double singleEventSensitivity;
    std::vector<double> upperLimitBR;
    std::vector<double> CLlowerBR;
    std::vector<double> CLupperBR;
      figureOfMeritCalculator.results 
        (maximumSignalCount, lowCut, highCut, sigEff, bkgd, figOfM,  
         discoveryThresholdN, discoveryThreshold, 
         smoothedWorstCaseSensitivity,punziSensitivity,
         CL90Sensitivity, singleEventSensitivity,
         upperLimitBR, CLlowerBR, CLupperBR); 
#endif
    
    std::cout <<
      "Enter fixed low momentum cut (or 0 for automatic optimizing): ";
    std::cin >> lowCut;
    std::cout <<
      "Enter fixed high momentum cut (or 0 for automatic optimizing): ";
    std::cin >> highCut;
    
    std::string t; 
    if ((lowCut == 0) && (highCut == 0)) {    
      t = figureOfMeritCalculator.tables(maximumSignalCount);
    } else if (lowCut == 0) { 
      t = figureOfMeritCalculator.tables_fixed_highCut 
                                (maximumSignalCount, highCut);
    } else if (highCut == 0) { 
      t = figureOfMeritCalculator.tables_fixed_lowCut 
                                (maximumSignalCount, lowCut);
    } else {
      t = figureOfMeritCalculator.tables_fixed_cuts 
                                (maximumSignalCount, lowCut, highCut);
    }

    std::cout << t << "\n";
       
    std::cout << "\n--------------------------------------- \n\n";

    // TEMPORARY: Exercise the poissonCDFSolver functions

#ifdef EXERCISE_POISSON_SOLVER    
    std::cout << "\npoissonComplementaryCDF(.4,6) = "
              << poissonComplementaryCDF(.4,6) 
              << "\n";
    std::cout << "\npoissonComplementaryCDF(.6,7) = "
              << poissonComplementaryCDF(.6,7) 
              << "\n";
    std::cout << "\npoissonInverseCDF(.4,.99999) = "
              << poissonInverseCDF(.4,.99999) 
              << "\n";
    std::cout << "\npoissonInverseCDF(.6,.999998) = "
              << poissonInverseCDF(.6,.999998) 
              << "\n";
    std::cout << "\npoissonInverseCDF(.3,.9) = "
              << poissonInverseCDF(.3,.9) 
              << "\n";
    std::cout << "\npoissonInverseCDF(.8,.9) = "
              << poissonInverseCDF(.8,.9) 
              << "\n";
    std::cout << "\npoissonMeanSolver(7,1-2.86652e-7) = "
              << poissonMeanSolver(7,1-2.86652e-7) 
              << "\n";
    std::cout << "\npoissonMeanSolver(6,1-2.86652e-7) = "
              << poissonMeanSolver(6,1-2.86652e-7) 
              << "\n";
    std::cout << "\npoissonMeanSolver(6,.9) = "
              << poissonMeanSolver(6,.9) 
              << "\n";
#endif

  } // end of TestFofM1 beginRun

} // end namespace mu2e

using mu2e::TestFofM1;
DEFINE_ART_MODULE(TestFofM1);
