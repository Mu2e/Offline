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
    std::string             _templateFile;
    bool                    _dropSecondDisk;
    std::string             _templates;
    double                  _minClEnergy, _clEStep;
    double                  _minRDist   , _rDistStep;
    double                  _minLH;

  //Histograms need to load the templates for the signal and background hypothesis
    TH1F*       _signalHist1D[kN1DVar];
    TH2F*       _signalHist2D[kN2DVar];

    TH1F*       _bkgHist1D   [kN1DVar];
    TH2F*       _bkgHist2D   [kN2DVar];

    TH1F*       _signalCorrHist1D[kN2DVar][kNCorHist];
    TH1F*       _bkgCorrHist1D   [kN2DVar][kNCorHist];
    
    //initialize calclate prob in private
    double calculateProb(double &variable, TH1* templates);
  };


  CaloLikelihood::CaloLikelihood(const fhicl::ParameterSet & pset) :
    _diagLevel                   (pset.get<int>("diagLevel",0)),
    _nProcess                    (0),      _signalCorrHist1D[0][i]->Scale(1/Area);
    }


    for (int j=0; j<kNCorHist; ++j){
      double Estart = _minClEnergy + j*_clEStep;
      int binstart = (Estart - binlowbkg)/binsizebkg + 1;
      int binend   = (binstart + (_clEStep/binsizebkg));
      // _bkgEnergyRad[j]= _bkgHist2D[0]->ProjectionX(Form("bkgEnRad_%i",j), binstart, binend);
      _bkgCorrHist1D[0][j] = (TH1F*)_bkgHist2D[0]->ProjectionX(Form("bkgEnRad_%i",j), binstart, binend);
      double Area = _bkgCorrHist1D[0][j]->Integral();
      _bkgCorrHist1D[0][j]->Scale(1/Area);
    }
    
    double binsizephoton2 = _signalHist2D[1]->GetYaxis()->GetBinWidth(1);
    double binsizebkg2    = _bkgHist2D   [1]->GetYaxis()->GetBinWidth(1);
    
    double binlowphoton2  = _signalHist2D[1]->GetYaxis()->GetBinLowEdge(1);
    double binlowbkg2     = _bkgHist2D   [1]->GetYaxis()->GetBinLowEdge(1);

    for (int i=0; i<kNCorHist; ++i){
      double  Estart   = (_minRDist + i*_rDistStep);
      int     binstart = (Estart - binlowphoton2)/binsizephoton2 + 1;
      int     binend   = binstart + (_rDistStep/binsizephoton2);
      _signalCorrHist1D[1][i]= (TH1F*)_signalHist2D[1]->ProjectionX(Form("photonRadRatio_%i",i), binstart, binend);
      double  Area     = _signalCorrHist1D[1][i]->Integral();
      _signalCorrHist1D[1][i]->Scale(1/Area);
    }

    for (int j=0; j<kNCorHist; ++j){
      double  Estart   = (_minRDist + j*_rDistStep);
      int     binstart = (Estart - binlowbkg2)/binsizebkg2 + 1;
      int     binend   = binstart + (_rDistStep/binsizebkg2);
      _bkgCorrHist1D[1][j]= (TH1F*)_bkgHist2D[1]->ProjectionX(Form("bkgRadRatio_%i",j), binstart, binend);
      double  Area     =_bkgCorrHist1D[1][j]->Integral();
      _bkgCorrHist1D[1][j]->Scale(1/Area);
    }
    
    double   binsizephoton3 = _signalHist2D[2]->GetYaxis()->GetBinWidth(1);
    double   binsizebkg3    = _bkgHist2D[2]->GetYaxis()->GetBinWidth(1);
    
    double   binlowphoton3  = _signalHist2D[2]->GetYaxis()->GetBinLowEdge(1);
    double   binlowbkg3     = _bkgHist2D[2]->GetYaxis()->GetBinLowEdge(1);

  
    for (int i=0; i<kNCorHist; ++i){
      double Estart   = (_minRDist + i*_rDistStep);
      int    binstart = (Estart - binlowphoton3)/binsizephoton3 + 1;
      int    binend   = binstart + (_rDistStep/binsizephoton3);
      _signalCorrHist1D[2][i] = (TH1F*)_signalHist2D[2]->ProjectionX(Form("photonRad2ndRatio_%i",i), binstart, binend);
      double Area     = _signalCorrHist1D[2][i]->Integral();
      _signalCorrHist1D[2][i]->Scale(1/Area);
    }


    for (int j=0; j<kNCorHist; ++j){
      double  Estart   = (_minRDist + j*_rDistStep);
      int     binstart = (Estart - binlowbkg3)/binsizebkg3 + 1;
      int     binend   = binstart + (_rDistStep/binsizebkg3);
      _bkgCorrHist1D[2][j] = (TH1F*)_bkgHist2D[2]->ProjectionX(Form("bkgRad2ndRatio_%i",j), binstart, binend);
      double  Area     = _bkgCorrHist1D[2][j]->Integral();
      _bkgCorrHist1D[2][j]->Scale(1/Area);
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
      
      double  signalClusterSizeProb = calculateProb(clusterSize, _signalHist1D[2]);
      double  bkgClusterSizeProb    = calculateProb(clusterSize, _bkgHist1D[2]);
      double  logClusterSizeRatio   = log(signalClusterSizeProb/bkgClusterSizeProb);
      
      double  signalRDistProb       = calculateProb(rDist, _signalHist1D[3]);
      double  bkgRDistProb          = calculateProb(rDist, _bkgHist1D[3]);
      double  logRDistRatio         = log(signalRDistProb/bkgRDistProb);

      double  signalE2E1Prob        = calculateProb(e2e1Ratio, _signalHist1D[4]);
      double  bkgE2E1Prob           = calculateProb(e2e1Ratio, _bkgHist1D[4]);
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

    double thisprob = 1;
    //not sure how to loop over variable value   i should be event number
    double     binSize   = Template->GetBinWidth(1);
    double     binIndex = (Variable - Template->GetBinLowEdge(1))/binSize + 1;
    //what is the probability for this variable
    int           templateNBins = Template->GetNbinsX();
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


