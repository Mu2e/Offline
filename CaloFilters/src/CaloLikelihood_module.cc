//
// A Filter module aimed to select events using a Likelihood defined with calorimeter cluster info
//
//
// Original author E. Castiglia, G. Pezzullo
//

#include "CLHEP/Units/SystemOfUnits.h"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"

#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/CaloCluster/inc/ClusterUtils.hh"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"

#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/TrkFitFlag.hh"
#include "Offline/RecoDataProducts/inc/TriggerInfo.hh"

// #include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
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


namespace mu2e {


  class CaloLikelihood : public art::EDFilter {

  public:

    enum {
      kN1DVar    = 10,
      kN2DVar    = 10,
      kNCorHist  = 10
    };

    virtual ~CaloLikelihood() { }

    virtual void beginJob();
    virtual void endJob  ();
    virtual bool filter  (art::Event& event) override;
    virtual bool endRun( art::Run& run ) override;

    explicit CaloLikelihood(const fhicl::ParameterSet& PSet);

  private:

    typedef art::Ptr< CaloHit> CaloHitPtr;

    int                     _diagLevel;
    int                     _nProcess;
    int                     _nPass;
    ClusterUtils::Cogtype   _cogType;
    art::InputTag           _clTag;
    std::string             _signalTemplateFile;
    std::string             _bkgTemplateFile;
    bool                    _dropSecondDisk;
    std::string             _signalTemplates;
    std::string             _bkgTemplates;
    double                  _minClEnergy, _clEStep;
    double                  _minRDist   , _rDistStep;
    std::vector<double>     _minLH;

  //Histograms need to load the templates for the signal and background hypothesis
    TH1F*       _signalHist1D[2][kN1DVar];
    TH2F*       _signalHist2D[2][kN2DVar];

    TH1F*       _bkgHist1D   [2][kN1DVar];
    TH2F*       _bkgHist2D   [2][kN2DVar];

    TH1F*       _signalCorrHist1D[2][kN2DVar][kNCorHist];
    TH1F*       _bkgCorrHist1D   [2][kN2DVar][kNCorHist];

    void        buildTemplateHist(TH2F*Input, TH1F**FinalHist,  TString Label, double MinX, double StepSize);

    //initialize calclate prob in private
    double      calculateProb    (double&variable, TH1* templates);
    double      calculate2DProb  (double& Ref, double&Variable, TH1F** Template, double MinX, double Step);
  };


  CaloLikelihood::CaloLikelihood(const fhicl::ParameterSet & pset) :
    art::EDFilter{pset},
    _diagLevel                   (pset.get<int>("diagLevel",0)),
    _nProcess                    (0),
    _nPass                       (0),
    _cogType                     (ClusterUtils::Linear),
    _clTag                       (pset.get<art::InputTag> ("CaloClusterModuleLabel")),
    _signalTemplateFile          (pset.get<std::string>   ("SignalTemplates")),
    _bkgTemplateFile             (pset.get<std::string>   ("BackgroundTemplates")),
    _dropSecondDisk              (pset.get<bool>          ("DropSecondDisk"       , false)),
    _minClEnergy                 (pset.get<double>        ("MinClusterEnergy"     ,   50.)),   // MeV
    _clEStep                     (pset.get<double>        ("ClusterEnergyStep"    ,   10.)),   // MeV
    _minRDist                    (pset.get<double>        ("MinClusterRadialDist" ,  350.)),   // mm
    _rDistStep                   (pset.get<double>        ("ClusterRadialDistStep",   50.)),   // mm
    _minLH                       (pset.get<std::vector<double>>("MinLikelihoodCut"     , std::vector<double>{1.,1.})){   // likelihood threshold

    produces<TriggerInfo>();

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
    _signalHist2D[0][3]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk0_e2e1RDist2");
    _signalHist2D[0][4]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk0_e3RDist2");
    _signalHist2D[0][5]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk0_e4RDist2");
    _signalHist2D[0][6]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk0_nCrRDist2");
    //disk 1
    _signalHist2D[1][0]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk1_eRDist2");
    _signalHist2D[1][1]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk1_e1RDist2");
    _signalHist2D[1][2]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk1_e2RDist2");
    _signalHist2D[1][3]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk1_e2e1RDist2");
    _signalHist2D[1][4]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk1_e3RDist2");
    _signalHist2D[1][5]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk1_e4RDist2");
    _signalHist2D[1][6]   = (TH2F*)signalFile->Get("ceAna/cluster_electron/clDisk1_nCrRDist2");

    //background
    //disk 0
    _bkgHist2D   [0][0]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk0_eRDist0");
    _bkgHist2D   [0][1]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk0_e1RDist0");
    _bkgHist2D   [0][2]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk0_e2RDist0");
    _bkgHist2D   [0][3]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk0_e2e1RDist0");
    _bkgHist2D   [0][4]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk0_e3RDist0");
    _bkgHist2D   [0][5]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk0_e4RDist0");
    _bkgHist2D   [0][6]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk0_nCrRDist0");
    //disk 1
    _bkgHist2D   [1][0]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk1_eRDist0");
    _bkgHist2D   [1][1]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk1_e1RDist0");
    _bkgHist2D   [1][2]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk1_e2RDist0");
    _bkgHist2D   [1][3]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk1_e2e1RDist0");
    _bkgHist2D   [1][4]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk1_e3RDist0");
    _bkgHist2D   [1][5]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk1_e4RDist0");
    _bkgHist2D   [1][6]   = (TH2F*)bkgFile->Get("ceAna/cluster_all/clDisk1_nCrRDist0");

    //make correlation templates
    //make correlation templates
    int     nCaloDisks(2);
    for (int i=0; i<nCaloDisks; ++i){
      //cluster energy vs cluster radial distance
      buildTemplateHist(_signalHist2D[i][0], _signalCorrHist1D[i][0], Form("SignalDisk%i_eRDist", i),  _minRDist, _rDistStep);
      buildTemplateHist(_bkgHist2D   [i][0], _bkgCorrHist1D   [i][0], Form("BkgDisk%i_eRDist"   , i),  _minRDist, _rDistStep);

      //E1 (seedHitEnergy/cluster_energy) vs cluster radial distance
      buildTemplateHist(_signalHist2D[i][1], _signalCorrHist1D[i][1], Form("SignalDisk%i_e1RDist", i), _minRDist, _rDistStep);
      buildTemplateHist(_bkgHist2D   [i][1], _bkgCorrHist1D   [i][1], Form("BkgDisk%i_e1RDist"   , i), _minRDist, _rDistStep);

      //E2 (2nd_HitEnergy/cluster_energy) vs cluster radial distance
      buildTemplateHist(_signalHist2D[i][2], _signalCorrHist1D[i][2], Form("SignalDisk%i_e2RDist", i), _minRDist, _rDistStep);
      buildTemplateHist(_bkgHist2D   [i][2], _bkgCorrHist1D   [i][2], Form("BkgDisk%i_e2RDist"   , i), _minRDist, _rDistStep);

      //E2/E1 vs cluster radial distance
      buildTemplateHist(_signalHist2D[i][3], _signalCorrHist1D[i][3], Form("SignalDisk%i_e2e1RDist", i), _minRDist, _rDistStep);
      buildTemplateHist(_bkgHist2D   [i][3], _bkgCorrHist1D   [i][3], Form("BkgDisk%i_e2e1RDist"   , i), _minRDist, _rDistStep);

      //E2 (2nd_HitEnergy/cluster_energy) vs cluster radial distance
      buildTemplateHist(_signalHist2D[i][4], _signalCorrHist1D[i][4], Form("SignalDisk%i_e3RDist", i), _minRDist, _rDistStep);
      buildTemplateHist(_bkgHist2D   [i][4], _bkgCorrHist1D   [i][4], Form("BkgDisk%i_e3RDist"   , i), _minRDist, _rDistStep);

      //E2 (2nd_HitEnergy/cluster_energy) vs cluster radial distance
      buildTemplateHist(_signalHist2D[i][5], _signalCorrHist1D[i][5], Form("SignalDisk%i_e4RDist", i), _minRDist, _rDistStep);
      buildTemplateHist(_bkgHist2D   [i][5], _bkgCorrHist1D   [i][5], Form("BkgDisk%i_e4RDist"   , i), _minRDist, _rDistStep);

      //nCrystals vs cluster radial distance
      buildTemplateHist(_signalHist2D[i][6], _signalCorrHist1D[i][6], Form("SignalDisk%i_nCrRDist", i), _minRDist, _rDistStep);
      buildTemplateHist(_bkgHist2D   [i][6], _bkgCorrHist1D   [i][6], Form("BkgDisk%i_nCrRDist"   , i), _minRDist, _rDistStep);

    }
  }

  //--------------------------------------------------------------------------------
  // routine function used to produce the template histograms
  // from the 2D distributions
  //--------------------------------------------------------------------------------
  void   CaloLikelihood::buildTemplateHist(TH2F*Input, TH1F**FinalHist,  TString Label, double MinX, double StepSize){
    double    binsize = Input->GetXaxis()->GetBinWidth(1);
    double    binlow  = Input->GetXaxis()->GetBinLowEdge(1);

    for (int i=0; i<kNCorHist; ++i){
      double  start      = (MinX + i*StepSize);
      int     binstart   = (start - binlow)/binsize;
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
      std::cout << moduleDescription().moduleLabel() << " passed " <<  _nPass << " events out of " << _nProcess << " for a ratio of " << float(_nPass)/float(_nProcess) << std::endl;
    }
    return true;
  }

  //--------------------------------------------------------------------------------
  // Follow the body of the Filter logic
  //--------------------------------------------------------------------------------
  bool CaloLikelihood::filter(art::Event& event) {

    ++_nProcess;
    if (_nProcess%10==0 && _diagLevel > 0) std::cout<<"Processing event from CaloLikelihood =  "<<_nProcess  <<std::endl;

    std::unique_ptr<TriggerInfo> triginfo(new TriggerInfo);
    bool   retval(false);

    //Get handle to the calorimeter
    art::ServiceHandle<GeometryService> geom;
    if( ! geom->hasElement<Calorimeter>() )       return 0;
    Calorimeter const & cal = *(GeomHandle<Calorimeter>());

    //Get calo cluster collection
    auto  clH = event.getValidHandle<CaloClusterCollection>(_clTag);
    const CaloClusterCollection*  caloClusters = clH.product();

    const CaloHit* crystalHit(0);
    const mu2e::CaloHitPtrVector* caloClusterHits(0);

    //for loop over the clusters in the calorimeter
    for(auto icl = caloClusters->begin();icl != caloClusters->end(); ++icl){
      auto const& cluster = *icl;

      int clSection = cluster.diskID();
      if (_dropSecondDisk && (clSection == 1))        continue;

      double clEnergy = cluster.energyDep();
      if ( clEnergy < _minClEnergy)                   continue;

      ClusterUtils cluUtil(cal,cluster,_cogType);
      const auto cog = cluUtil.cog3Vector();

      double                   xpos     = cog(0);
      double                   ypos     = cog(1);
      double                   rDist    = sqrt(xpos*xpos+ypos*ypos);

      caloClusterHits = &cluster.caloHitsPtrVector();

      int       nCrystalHits = caloClusterHits->size();
      double    maxECrystal(0);
      double    secondmaxECrystal(-1), thirdmaxECrystal(-1)   , fourthmaxECrystal(-1);
      double    indexMaxECrystal(0)  , index2ndMaxECrystal(-1), index3rdMaxECrystal(-1);

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
        if (j == indexMaxECrystal)              continue;
        crystalHit = &(*caloClusterHits->at(j));

        double crystalEnergy = crystalHit->energyDep();

        if (crystalEnergy > secondmaxECrystal){
          secondmaxECrystal   = crystalEnergy;
          index2ndMaxECrystal = j;
        }
      }

      for (int j=0; j<nCrystalHits; ++j){
        if (j == indexMaxECrystal ||
            j == index2ndMaxECrystal  )         continue;
        crystalHit = &(*caloClusterHits->at(j));

        double crystalEnergy = crystalHit->energyDep();

        if (crystalEnergy > thirdmaxECrystal){
          thirdmaxECrystal    = crystalEnergy;
          index3rdMaxECrystal = j;
        }
      }

      for (int j=0; j<nCrystalHits; ++j){
        if (j == indexMaxECrystal  ||
            j == index2ndMaxECrystal ||
            j == index3rdMaxECrystal )          continue;
        crystalHit = &(*caloClusterHits->at(j));

        double crystalEnergy = crystalHit->energyDep();

        if (crystalEnergy > fourthmaxECrystal){
          fourthmaxECrystal    = crystalEnergy;
        }
      }

      double   e1Ratio   = maxECrystal/clEnergy;
      double   e2Ratio   = secondmaxECrystal/clEnergy;
      double   e3Ratio   = thirdmaxECrystal/clEnergy;
      double   e4Ratio   = fourthmaxECrystal/clEnergy;
      double   e2e1Ratio = secondmaxECrystal/maxECrystal;

      //calculate probability for likelihood
      double  clusterSize             = (double) nCrystalHits;

      double  signalRDistProb         = calculateProb(rDist, _signalHist1D[clSection][3]);
      double  bkgRDistProb            = calculateProb(rDist, _bkgHist1D   [clSection][3]);
      double  logRDistRatio           = log(signalRDistProb/bkgRDistProb);

      double  signalEnergy2DProb      = calculate2DProb(rDist, rDist, _signalCorrHist1D[clSection][0], _minClEnergy, _clEStep);
      double  bkgEnergy2DProb         = calculate2DProb(rDist, rDist, _bkgCorrHist1D   [clSection][0], _minClEnergy, _clEStep);
      double  logEnergy2DRatio        = log(signalEnergy2DProb/bkgEnergy2DProb);

      double  signalE12DProb          = calculate2DProb(rDist, e1Ratio, _signalCorrHist1D[clSection][1], _minRDist, _rDistStep);
      double  bkgE12DProb             = calculate2DProb(rDist, e1Ratio, _bkgCorrHist1D   [clSection][1], _minRDist, _rDistStep);
      double  logE12DRatio            = log(signalE12DProb/bkgE12DProb);

      double  signalE22DProb          = calculate2DProb(rDist, e2Ratio, _signalCorrHist1D[clSection][2], _minRDist, _rDistStep);
      double  bkgE22DProb             = calculate2DProb(rDist, e2Ratio, _bkgCorrHist1D   [clSection][2], _minRDist, _rDistStep);
      double  logE22DRatio            = log(signalE22DProb/bkgE22DProb);

      double  signalE2E12DProb        = calculate2DProb(rDist, e2e1Ratio, _signalCorrHist1D[clSection][3], _minRDist, _rDistStep);
      double  bkgE2E12DProb           = calculate2DProb(rDist, e2e1Ratio, _bkgCorrHist1D   [clSection][3], _minRDist, _rDistStep);
      double  logE2E12DRatio          = log(signalE2E12DProb/bkgE2E12DProb);

      double  signalE32DProb          = calculate2DProb(rDist, e3Ratio, _signalCorrHist1D[clSection][4], _minRDist, _rDistStep);
      double  bkgE32DProb             = calculate2DProb(rDist, e3Ratio, _bkgCorrHist1D   [clSection][4], _minRDist, _rDistStep);
      double  logE32DRatio            = log(signalE32DProb/bkgE32DProb);

      double  signalE42DProb          = calculate2DProb(rDist, e4Ratio, _signalCorrHist1D[clSection][5], _minRDist, _rDistStep);
      double  bkgE42DProb             = calculate2DProb(rDist, e4Ratio, _bkgCorrHist1D   [clSection][5], _minRDist, _rDistStep);
      double  logE42DRatio            = log(signalE42DProb/bkgE42DProb);

      double  signalClSize2DProb      = calculate2DProb(rDist, clusterSize, _signalCorrHist1D[clSection][6], _minRDist, _rDistStep);
      double  bkgClSize2DProb         = calculate2DProb(rDist, clusterSize , _bkgCorrHist1D  [clSection][6], _minRDist, _rDistStep);
      double  logClSize2DRatio        = log(signalClSize2DProb/bkgClSize2DProb);

      double  lhValue  =  logRDistRatio + logEnergy2DRatio + logE12DRatio + logRDistRatio + logClSize2DRatio;
      if (nCrystalHits>=2) lhValue += logE22DRatio  + logE2E12DRatio;
      if (nCrystalHits>=3) lhValue += logE32DRatio;
      if (nCrystalHits>=4) lhValue += logE42DRatio;

      if ( lhValue > _minLH[clSection]) {
        retval = true;
        ++_nPass;
        // Fill the trigger info object
        // associate to the caloCluster which triggers.  Note there may be other caloClusters which also pass the filter
        // but filtering is by event!
        size_t index = std::distance(caloClusters->begin(), icl);
        triginfo->_caloClusters.push_back(art::Ptr<CaloCluster>(clH, index));

        if(_diagLevel > 1){
          std::cout << moduleDescription().moduleLabel() << " passed event " << event.id() << std::endl;
        }
      }


    }

    event.put(std::move(triginfo));
    return retval;
  }

  //--------------------------------------------------------------------------------//
 //--------------------------------------------------------------------------------
  double   CaloLikelihood::calculate2DProb(double& Ref, double&Variable, TH1F** Template, double MinX, double Step){
    double     thisprob = 1e-17;

    for (int j=0; j<kNCorHist; ++j){
      double   ref_min = (MinX + (double)j*Step);
      double   ref_max = (MinX + (double)(j+1)*Step);
      if ( ( Ref >= ref_min) && (Ref < ref_max) ){
        thisprob = calculateProb(Variable,  Template[j]);
        break;
      }
    }

    return thisprob;
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
      //      std::cout << "whoops out of range" << std::endl;
    } else{
      thisprob = Template->GetBinContent(binIndex);
    }

    if (thisprob < 1e-10) thisprob = 1e-17;

    return thisprob;

  }


}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::CaloLikelihood)


