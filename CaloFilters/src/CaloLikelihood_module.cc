//
// A Filter module aimed to select events using a Likelihood defined with calorimeter cluster info
//
// $Id: $
// $Author: $
// $Date: $
//
// Original author E. Castiglia, G. Pezzullo
//

#include "CLHEP/Units/SystemOfUnits.h"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"

#include "CaloCluster/inc/ClusterMoments.hh"

#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "RecoDataProducts/inc/TriggerInfo.hh"

// #include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"

#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"

#include "TDirectory.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TFile.h"

#include <cmath>
#include <string>
#include <vector>


using namespace std;

namespace mu2e {


  class CaloLikelihood : public art::EDFilter {
     
  public:
    
    enum {
      kN1DVar    = 10,
      kN2DVar    = 3,
      kNCorHist  = 10
    };

    virtual ~CaloLikelihood() { }

    virtual void beginJob();
    virtual void endJob  ();
    virtual bool filter  (art::Event& event) override;
    virtual bool endRun( art::Run& run ) override;

    explicit CaloLikelihood(const fhicl::ParameterSet& PSet);

  private:
       
    typedef art::Ptr< CaloCrystalHit> CaloCrystalHitPtr;

    int                     _diagLevel;
    int                     _nProcess;
    int                     _nPass;
    ClusterMoments::cogtype _cogType;
    art::InputTag           _clTag;
    std::string             _signalTemplateFile;
    std::string             _bkgTemplateFile;
    bool                    _dropSecondDisk;
    std::string             _signalTemplates;
    std::string             _bkgTemplates;
    double                  _minClEnergy, _clEStep;
    double                  _minRDist   , _rDistStep;
    double                  _minLH;

  //Histograms need to load the templates for the signal and background hypothesis
    TH1F*       _signalHist1D[2][kN1DVar];
    TH2F*       _signalHist2D[2][kN2DVar];

    TH1F*       _bkgHist1D   [2][kN1DVar];
    TH2F*       _bkgHist2D   [2][kN2DVar];

    TH1F*       _signalCorrHist1D[kN2DVar][kNCorHist];
    TH1F*       _bkgCorrHist1D   [kN2DVar][kNCorHist];
    
    void   buildTemplateHist(TH2F*Input, TH1F**FinalHist,  TString Label, double MinX, double StepSize);

    //initialize calclate prob in private
    double calculateProb(double &variable, TH1* templates);
  };


  CaloLikelihood::CaloLikelihood(const fhicl::ParameterSet & pset) :
    _diagLevel                   (pset.get<int>("diagLevel",0)),
    _nProcess                    (0),
    _nPass                       (0),
    _cogType                     (ClusterMoments::Linear),                
    _clTag                       (pset.get<art::InputTag>("CaloClusterModuleLabel")),
    _signalTemplateFile          (pset.get<string>("SignalTemplates")),
    _bkgTemplateFile             (pset.get<string>("BackgroundTemplates")),
    _dropSecondDisk              (pset.get<bool>  ("DropSecondDisk"       , false)),
    _minClEnergy                 (pset.get<double>("MinClusterEnergy"     ,   50.)),   // MeV
    _clEStep                     (pset.get<double>("ClusterEnergyStep"    ,   10.)),   // MeV
    _minRDist                    (pset.get<double>("MinClusterRadialDist" ,  350.)),   // mm
    _rDistStep                   (pset.get<double>("ClusterRadialDistStep",   50.)),   // mm
    _minLH                       (pset.get<double>("MinLikelihoodCut"     ,   1.)){   // likelihood threshold

    ConfigFileLookupPolicy configFile;
    _signalTemplates = configFile(_signalTemplateFile);
    _bkgTemplates    = configFile(_bkgTemplateFile);
    
    TFile* signalFile = TFile::Open(_signalTemplates.c_str());
    TFile* bkgFile    = TFile::Open(_bkgTemplates.c_str());

    //get the templates histograms
    //signal on calo-disk 0
    _signalHist1D[0][0]   = (TH1F*)signalFile->Get("ceAna/cluster_electron/clDisk0_E2");	
    _signalHist1D[0][1]   = (TH1F*)signalFile->Get("ceAna/cluster_electron/clDisk0_e1eRatio2");   
    _signalHist1D[0][2]   = (TH1F*)signalFile->Get("ceAna/cluster_electron/clDisk0_nCr2");   	
    _signalHist1D[0][3]   = (TH1F*)signalFile->Get("ceAna/cluster_electron/clDisk0_rDist2");	
    _signalHist1D[0][4]   = (TH1F*)signalFile->Get("ceAna/cluster_electron/clDisk0_e2eRatio2");	
    _signalHist1D[0][5]   = (TH1F*)signalFile->Get("ceAna/cluster_electron/clDisk0_e2e1Ratio2");
    //signal on calo-disk 1
    _signalHist1D[1][0]   = (TH1F*)signalFile->Get("ceAna/cluster_electron/clDisk1_E2");
    _signalHist1D[1][1]   = (TH1F*)signalFile->Get("ceAna/cluster_electron/clDisk1_e1eRatio2");   
    _signalHist1D[1][2]   = (TH1F*)signalFile->Get("ceAna/cluster_electron/clDisk1_nCr2");   
    _signalHist1D[1][3]   = (TH1F*)signalFile->Get("ceAna/cluster_electron/clDisk1_rDist2");
    _signalHist1D[1][4]   = (TH1F*)signalFile->Get("ceAna/cluster_electron/clDisk1_e2eRatio2");
    _signalHist1D[1][5]   = (TH1F*)signalFile->Get("ceAna/cluster_electron/clDisk1_e2e1Ratio2");

    //background on calo-disk 0
    _bkgHist1D   [0][0]   = (TH1F*)bkgFile->Get("ceAna/cluster_all/clDisk0_E0");	
    _bkgHist1D   [0][1]   = (TH1F*)bkgFile->Get("ceAna/cluster_all/clDisk0_e1eRatio0");
    _bkgHist1D   [0][2]   = (TH1F*)bkgFile->Get("ceAna/cluster_all/clDisk0_nCr0");   
    _bkgHist1D   [0][3]   = (TH1F*)bkgFile->Get("ceAna/cluster_all/clDisk0_rDist0");
    _bkgHist1D   [0][4]   = (TH1F*)bkgFile->Get("ceAna/cluster_all/clDisk0_e2eRatio0");
    _bkgHist1D   [0][5]   = (TH1F*)bkgFile->Get("ceAna/cluster_all/clDisk0_e2e1Ratio0");
    //background on calo-disk 1
    _bkgHist1D   [1][0]   = (TH1F*)bkgFile->Get("ceAna/cluster_all/clDisk1_E0");	
    _bkgHist1D   [1][1]   = (TH1F*)bkgFile->Get("ceAna/cluster_all/clDisk1_e1eRatio0");
    _bkgHist1D   [1][2]   = (TH1F*)bkgFile->Get("ceAna/cluster_all/clDisk1_nCr0");   
    _bkgHist1D   [1][3]   = (TH1F*)bkgFile->Get("ceAna/cluster_all/clDisk1_rDist0");
    _bkgHist1D   [1][4]   = (TH1F*)bkgFile->Get("ceAna/cluster_all/clDisk1_e2eRatio0");
    _bkgHist1D   [1][5]   = (TH1F*)bkgFile->Get("ceAna/cluster_all/clDisk1_e2e1Ratio0");
    

    //get the 2D histograms
    // signal
    //disk 0
    _signalHist2D[0][0]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk0_eRDist2"); 
    _signalHist2D[0][1]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk0_e1RDist2");
    _signalHist2D[0][2]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk0_e2RDist2");
    //disk 1
    _signalHist2D[1][0]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk1_eRDist2");
    _signalHist2D[1][1]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk1_e1RDist2");
    _signalHist2D[1][2]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk1_e2RDist2");

    //background
    //disk 0
    _bkgHist2D   [0][0]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk0_eRDist0"); 
    _bkgHist2D   [0][1]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk0_e1RDist0");
    _bkgHist2D   [0][2]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk0_e2RDist0");
    //disk 1
    _bkgHist2D   [1][0]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk1_eRDist0"); 
    _bkgHist2D   [1][1]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk1_e1RDist0");
    _bkgHist2D   [1][2]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk1_e2RDist0");

    //make correlation templates
    
    //cluster energy vs cluster radial distance
    buildTemplateHist(_signalHist2D[0][0], _signalCorrHist1D[0], "SignalDisk0_eRDist", _minClEnergy, _clEStep);
    buildTemplateHist(_bkgHist2D   [0][0], _bkgCorrHist1D   [0], "BkgDisk0_eRDist"   , _minClEnergy, _clEStep);
    
    //E1 (seedHitEnergy/cluster_energy) vs cluster radial distance
    buildTemplateHist(_signalHist2D[0][1], _signalCorrHist1D[1], "SignalDisk0_e1RDist", _minRDist, _rDistStep);
    buildTemplateHist(_bkgHist2D   [0][1], _bkgCorrHist1D   [1], "BkgDisk0_e1RDist"   , _minRDist, _rDistStep);

    //E2 (2nd_HitEnergy/cluster_energy) vs cluster radial distance
    buildTemplateHist(_signalHist2D[0][2], _signalCorrHist1D[2], "SignalDisk0_e2RDist", _minRDist, _rDistStep);
    buildTemplateHist(_bkgHist2D   [0][2], _bkgCorrHist1D   [2], "BkgDisk0_e2RDist"   , _minRDist, _rDistStep);
  
    //E2/E1 vs cluster radial distance
    // buildTemplateHist(_signalHist2D[0][3], _signalCorrHist1D[3], "SignalDisk0_e2e1RDist", _minRDist, _rDistStep);
    // buildTemplateHist(_bkgHist2D   [0][3], _bkgCorrHist1D   [3], "BkgDisk0_e2e1RDist"   , _minRDist, _rDistStep);
  
  }

  //--------------------------------------------------------------------------------
  // routine function used to produce the template histograms 
  // from the 2D distributions
  //--------------------------------------------------------------------------------
  void   CaloLikelihood::buildTemplateHist(TH2F*Input, TH1F**FinalHist,  TString Label, double MinX, double StepSize){
    double    binsize = Input->GetYaxis()->GetBinWidth(1);
    double    binlow  = Input->GetYaxis()->GetBinLowEdge(1);
    
    for (int i=0; i<kNCorHist; ++i){
      double  start      = (MinX + i*StepSize);
      int     binstart   = (start - binlow)/binsize + 1;
      int     binend     = binstart + (StepSize/binsize);
      FinalHist[i]       = (TH1F*)Input->ProjectionY(Form("%s_%i", Label.Data(), i), binstart, binend);
      double  area       = FinalHist[i]->Integral();
      FinalHist[i]->Scale(1/area);
    }
    
  }

  void CaloLikelihood::beginJob(){ }

  void CaloLikelihood::endJob(){}

  bool CaloLikelihood::endRun( art::Run& run ) {
    if(_diagLevel > 0 && _nProcess > 0){
      cout << *currentContext()->moduleLabel() << " passed " <<  _nPass << " events out of " << _nProcess << " for a ratio of " << float(_nPass)/float(_nProcess) << endl;
    }
    return true;
  }
  
  //--------------------------------------------------------------------------------
  // Follow the body of the Filter logic
  //--------------------------------------------------------------------------------
  bool CaloLikelihood::filter(art::Event& event) {

    ++_nProcess;
    if (_nProcess%10==0 && _diagLevel > 0) std::cout<<"Processing event from CaloLikelihood =  "<<_nProcess  <<std::endl;
   
    unique_ptr<TriggerInfo> triginfo(new TriggerInfo);
    bool   retval(false);

    //Get handle to the calorimeter
    art::ServiceHandle<GeometryService> geom;
    if( ! geom->hasElement<Calorimeter>() )       return 0;
    Calorimeter const & cal = *(GeomHandle<Calorimeter>());
  
    //Get calo cluster collection
    auto  clH = event.getValidHandle<CaloClusterCollection>(_clTag);
    const CaloClusterCollection*  caloClusters = clH.product();

    const CaloCrystalHit* crystalHit(0);
    const CaloCluster::CaloCrystalHitPtrVector* caloClusterHits(0);
    
    //for loop over the clusters in the calorimeter
    for(auto icl = caloClusters->begin();icl != caloClusters->end(); ++icl){
      auto const& cluster = *icl;
      
      int clSection = cluster.diskId();
      if (_dropSecondDisk && (clSection == 1))        continue;
      
      ClusterMoments cogCalculator(cal, cluster,clSection);
      cogCalculator.calculate(_cogType);

      double                   clEnergy = cluster.energyDep();
      if ( clEnergy < _minClEnergy)                   continue;

      const CLHEP::Hep3Vector& cog      = cluster.cog3Vector();
      double                   xpos     = cog(0);
      double                   ypos     = cog(1);
      double                   rDist    = sqrt(xpos*xpos+ypos*ypos);

      caloClusterHits = &cluster.caloCrystalHitsPtrVector();

      int       nCrystalHits = caloClusterHits->size();
      double    maxECrystal(0);
      double    secondmaxECrystal(0);
      double    indexMaxECrystal(0);

      //first loop to find the most energetic crystalHit
      for (int j=0; j<nCrystalHits; ++j){
	crystalHit = &(*caloClusterHits->at(j));
	double   crystalEnergy = crystalHit->energyDep();
	if (crystalEnergy > maxECrystal) {
	  maxECrystal    = crystalEnergy;
	  indexMaxECrystal   = j;
	}
      }

      //second loop to find the second most energetic crystalHit
      for (int j=0; j<nCrystalHits; ++j){
	if (j == indexMaxECrystal) continue;
	crystalHit = &(*caloClusterHits->at(j));
	double crystalEnergy= crystalHit->energyDep();
	if (crystalEnergy > secondmaxECrystal){
	  secondmaxECrystal =crystalEnergy;
	}
      }
            
      double   e1Ratio   = maxECrystal/clEnergy;
      double   e2Ratio   = secondmaxECrystal/clEnergy;
      double   e2e1Ratio = secondmaxECrystal/maxECrystal;

      //calculate probability for likelihood 
      // double  signalEClusterProb    = calculateProb(clEnergy, _signalHist1D[0]);
      // double  bkgEClusterProb       = calculateProb(clEnergy, _bkgHist1D[0]);
      // double  logEClusterRatio      = log(signalEClusterProb/bkgEClusterProb);

      // double  signalE1Prob          = calculateProb(e1Ratio, _signalHist1D[1]);
      // double  bkgE1Prob             = calculateProb(e1Ratio, _bkgHist1D[1]);
      // double  logE1Ratio            = log(signalE1Prob/bkgE1Prob);
      
      double  clusterSize = (double) nCrystalHits;
      
      double  signalClusterSizeProb = calculateProb(clusterSize, _signalHist1D[clSection][2]);
      double  bkgClusterSizeProb    = calculateProb(clusterSize , _bkgHist1D[clSection][2]);
      double  logClusterSizeRatio   = log(signalClusterSizeProb/bkgClusterSizeProb);
      
      double  signalRDistProb       = calculateProb(rDist, _signalHist1D[clSection][3]);
      double  bkgRDistProb          = calculateProb(rDist, _bkgHist1D   [clSection][3]);
      double  logRDistRatio         = log(signalRDistProb/bkgRDistProb);

      double  signalE2E1Prob        = calculateProb(e2e1Ratio, _signalHist1D[clSection][4]);
      double  bkgE2E1Prob           = calculateProb(e2e1Ratio, _bkgHist1D   [clSection][4]);
      double  logE2E1Ratio          = log(signalE2E1Prob/bkgE2E1Prob);
      

      //Energy versus Radial correlation
      double  logEClusterRDistRatio(-1e10);
      for (int j=0; j<kNCorHist; ++j){
	double   eCl_min = (_minClEnergy + (double)j*_clEStep);
	double   eCl_max = (_minClEnergy + (double)(j+1)*_clEStep);
	if ( ( clEnergy >= eCl_min) && ( clEnergy < eCl_max) ){
	  double signalEClusterRDistProb = calculateProb(rDist,  _signalCorrHist1D[0][j]);
	  double bkgEClusterRDistProb    = calculateProb(rDist, _bkgCorrHist1D    [0][j]);
	  logEClusterRDistRatio          = log(signalEClusterRDistProb/bkgEClusterRDistProb);
	  break;
	}
      }
   
      //E1 Ratio (Energy_most_energetic_crystalHit / cluster_energy) versus Radial correlation
      double logE1RDistRatio(-1e10);
      for (int j=0; j<kNCorHist; ++j){
	double   rDist_min = (_minRDist + (double)j*_rDistStep);
	double   rDist_max = (_minRDist + (double)(j+1)*_rDistStep);
	if ( ( rDist >= rDist_min) && (rDist < rDist_max) ){
	  double signalE1RDistProb = calculateProb(e1Ratio,  _signalCorrHist1D[1][j]);
	  double bkgE1RDistProb    = calculateProb(e1Ratio, _bkgCorrHist1D    [1][j]);
	  logE1RDistRatio          = log(signalE1RDistProb/bkgE1RDistProb);
	  break;
	}
      }

      //Second Highest Crystal Ratio versus Radial
      double logE2RDistRatio(-1e10);
      for (int j=0; j<kNCorHist; ++j){
	double   rDist_min = (_minRDist + (double)j*_rDistStep);
	double   rDist_max = (_minRDist + (double)(j+1)*_rDistStep);
	if ( ( rDist >= rDist_min) && (rDist < rDist_max) ){
	  double signalE2RDistProb = calculateProb(e2Ratio,  _signalCorrHist1D[2][j]);
	  double bkgE2RDistProb    = calculateProb(e2Ratio, _bkgCorrHist1D    [2][j]);
	  logE2RDistRatio          = log(signalE2RDistProb/bkgE2RDistProb);
	  break;
	}
      }

      // double lh_1 = logE1Ratio + logEClusterRatio + logClusterSizeRatio + logRDistRatio;
      // double lh_2 = lh_1 + logE2E1Ratio;

      double lh_4 = logE1RDistRatio  + logE2RDistRatio + logEClusterRDistRatio + logClusterSizeRatio + logRDistRatio;
      double lh_5 = lh_4 + logE2E1Ratio;
   
      if ( lh_5 > _minLH) {
	retval = true;
	++_nPass;
        // Fill the trigger info object
        triginfo->_triggerBits.merge(TriggerFlag::caloCluster);
        // associate to the caloCluster which triggers.  Note there may be other caloClusters which also pass the filter
        // but filtering is by event!
        size_t index = std::distance(caloClusters->begin(), icl);
        triginfo->_caloCluster = art::Ptr<CaloCluster>(clH, index);
        
	if(_diagLevel > 1){
          cout << *currentContext()->moduleLabel() << " passed event " << event.id() << endl;
        }
	
        break;
      }
      
      
    }

    event.put(std::move(triginfo));
    return retval;
  }
 
  double CaloLikelihood::calculateProb(double&Variable, TH1* Template){

    double     thisprob      = 1;
    //not sure how to loop over variable value   i should be event number
    double     binSize       = Template->GetBinWidth(1);
    int        binIndex      = (Variable - Template->GetBinLowEdge(1))/binSize;
    //what is the probability for this variable
    int        templateNBins = Template->GetNbinsX();
    if (binIndex > templateNBins) {
      thisprob = 1e-17;
      std::cout << "whoops out of range" << std::endl;
    } else{
      thisprob = Template->GetBinContent(binIndex);
    }
      
    if (thisprob < 1e-10) thisprob = 1e-17;

    return thisprob;

  }


}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::CaloLikelihood);


