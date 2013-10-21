//
// identification and track parameter extraction modules
//
// $Id: TrackReco_module.cc,v 1.7 2013/10/21 21:01:22 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/10/21 21:01:22 $
//
// Original author G. Tassielli
//

// C++ includes.
#include <iostream>
#include <string>
#include <memory>
#include <map>
#include <utility>
#include <limits>
#include <cmath>
//#include <unordered_set>
#include <set>

#include <boost/shared_ptr.hpp>


// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


//#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/GenericFunctions/Erf.hh"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
//#include "ITrackerGeom/inc/Cell.hh"
//#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
//#include "MCDataProducts/inc/GenParticleCollection.hh"
//#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
//#include "MCDataProducts/inc/SimParticleCollection.hh"
//#include "MCDataProducts/inc/StepPointMCCollection.hh"
//#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
//#include "ITrackerGeom/inc/ITracker.hh"
#include "TTrackerGeom/inc/TTracker.hh"
//#include "FastPatternReco/inc/TTHitPerTrackData.hh"
//#include "FastPatternReco/inc/GenTrackData.hh"
//#include "MCDataProducts/inc/GenId.hh"
//#include "MCDataProducts/inc/VisibleGenElTrack.hh"
//#include "MCDataProducts/inc/VisibleGenElTrackCollection.hh"
#include "RecoDataProducts/inc/TrackerHitTimeCluster.hh"
#include "RecoDataProducts/inc/TrackerHitTimeClusterCollection.hh"
#include "RecoDataProducts/inc/SectorStationCluster.hh"
#include "RecoDataProducts/inc/SctrSttnClusterGroup.hh"
#include "RecoDataProducts/inc/ZRotStrawHitMap.hh"
#include "RecoDataProducts/inc/ZRotStrawHitMapCollection.hh"
#include "RecoDataProducts/inc/SctrSttnClusterGroupCollection.hh"

// Root includes.
#include "TApplication.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
//#include "TClonesArray.h"
//#include "TObjArray.h"
//#include "TLine.h"
//#include "TArrow.h"
//#include "TEllipse.h"
//#include "TMarker.h"
//#include "TBox.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TF1.h"
#include "TSpectrum2.h"
#include "TLatex.h"
#include "TTree.h"
#include "TPolyMarker.h"
#include "TVirtualFitter.h"

using namespace std;

//Double_t InvG(Double_t *x, Double_t *p){
//
//  Double_t xP = x[0]+p[0];
//  Double_t xM = x[0]-p[0];
//  //cout<<"x0 "<<x[0]<<" x1 "<<x[1]<<endl;
//  //return 0.0;
//
//  Double_t outVal;
//  //outVal = InvGInt(xP)-InvGInt(xM);
//  //outVal = TMath::Sqrt(-TMath::Log(xM))-TMath::Sqrt(-TMath::Log(xP));
//
//  ///if (x[0]>p[3]) {
//    Double_t x1 = (x[0]-p[3])/(1.0-p[3]);
//    if (x1<1.0){
//      outVal = p[0]/(2.0*x1*TMath::Sqrt(-TMath::Log(x1)));
//      outVal *= 2.0*p[1]*1.414213562/p[2];
//    }
//    else outVal=0.0;//x[1];
//  //}
//  //else outVal=x[1];
//
//  //Double_t x1 = (1.0+p[6]/p[1])*x[0];
//  //if (x1<1.0){
//  //  outVal += 2.0*p[6]*1.414213562/p[2]*p[0]/(2.0*x1*TMath::Sqrt(-TMath::Log(x1)));
//  //}
//  //outVal *= p[3];
//
//  //outVal += p[4]/p[5]*TMath::Exp(-x[0]/p[5]);
//  //outVal += p[6]/p[7]*TMath::Exp(-x[0]/p[7]);
//
//  //xM = 1.0-xM;
//  //xP = 1.0-xP;
//  //outVal = TMath::Erf( TMath::Sqrt(-TMath::Log(xM)) ) - TMath::Erf( TMath::Sqrt(-TMath::Log(xP)) );
//  //outVal *=p[1];
//
//  return outVal;
//}
//
//TF1 invG("invG",InvG,0,1,4);

namespace mu2e {

  double zo=-1500.0*CLHEP::mm;
  double targetRadMax=100.0*CLHEP::mm;
  const int maxNpeaks=100;

  typedef art::Ptr<StrawHit> StrawHitPtr;
  //typedef std::map<unsigned int, pair<double, double> > absSectAngleTable;  //key = absSectId, val=( cos(theta), sin(theta) )
  //absSectAngleTable _absSectAngleTablep;
  double _wireLengthStep=2.0*10.0*CLHEP::mm;

  double removeAnglePeriod(double &angle){
          // remove the 2pi period and return the angle in the range [0,2pi[
          double intpart, fractpart;
          double tmp=angle/CLHEP::twopi;
          fractpart=std::modf (tmp , &intpart);
          fractpart*=CLHEP::twopi;
          if (fractpart<0.0) fractpart+=CLHEP::twopi;
          //if (fractpart<-CLHEP::pi) fractpart+=CLHEP::twopi;
          //if (fractpart> CLHEP::pi) fractpart-=CLHEP::twopi;
          //fractpart+=CLHEP::pi;

          //cout<<"Rescaling angle in: "<<angle/CLHEP::degree<<" out: "<<fractpart/CLHEP::degree<<endl;
          return fractpart;
  }

  struct HitData{
          const StrawHitPtr _shit;
          const CLHEP::Hep3Vector _direct;
          const CLHEP::Hep3Vector _midp;
          double _radius;
          double _halfLength;
          double _theta;   // [0,2pi[
          unsigned int _absZId;
          unsigned int _absSectId;
          vector< pair<double, double> > _wireXYTable;
          /*HitData():
                  _direct(0),
                  _midp(0),
                  _halfLength(0){
          }*/
          HitData(const StrawHitPtr &shit, const CLHEP::Hep3Vector & direct, const CLHEP::Hep3Vector & midp, double halfLength, unsigned int absZId, unsigned int absSectId):
                  _shit(shit),
                  _direct(direct),
                  _midp(midp),
                  _halfLength(halfLength),
                  _absZId(absZId),
                  _absSectId(absSectId){

                  _radius = sqrt( pow(_midp.getX(),2) + pow(_midp.getY(),2) );
                  _theta=_direct.getPhi()-CLHEP::halfpi;
                  //if ( _theta<=-CLHEP::pi ) _theta+=CLHEP::twopi;
                  //if ( _theta>  CLHEP::pi ) _theta-=CLHEP::twopi;
                  //_theta+=CLHEP::pi;
                  _theta=removeAnglePeriod(_theta);

                  //absSectAngleTable::iterator _absSectAngleTablep_it = _absSectAngleTablep.find(_absSectId);
                  //if ( _absSectAngleTablep_it==_absSectAngleTablep.end() ) {
                  //        //double theta=_direct.getPhi();
                  //        //_absSectAngleTablep_it = (_absSectAngleTablep.insert( absSectAngleTable::value_type( _absSectId, pair<double, double>(cos(theta),sin(theta) ) ) ) ).first;
                  //        _absSectAngleTablep_it = (_absSectAngleTablep.insert( absSectAngleTable::value_type( _absSectId, pair<double, double>(_direct.getX(),_direct.getY() ) ) ) ).first;
                  //}
                  _wireXYTable.clear();
                  //_wireXYTable.push_back( pair<double, double>(_midp.getX(), _midp.getY()) );
                  //double tmpL, tmpX, tmpY;
                  //for (int j=1; j<=(int)floor(_halfLength/_wireLengthStep); ++j ) {
                  //        tmpL = ((double) j)*_wireLengthStep;
                  //        tmpX = tmpL*_direct.getX()/*_absSectAngleTablep_it->second.first*/;
                  //        tmpY = tmpL*_direct.getY()/*_absSectAngleTablep_it->second.second*/;
                  //        _wireXYTable.push_back( pair<double, double>(_midp.getX()+tmpX, _midp.getY()+tmpY) );
                  //        _wireXYTable.push_back( pair<double, double>(_midp.getX()-tmpX, _midp.getY()-tmpY) );
                  //}
          }

          void RecalcWXYTab (double stepSize) {
                  //absSectAngleTable::iterator _absSectAngleTablep_it = _absSectAngleTablep.find(_absSectId);
                  _wireLengthStep=stepSize;
                  _wireXYTable.clear();
                  _wireXYTable.push_back( pair<double, double>(_midp.getX(), _midp.getY()) );
                  double tmpL, tmpX, tmpY;
                  for (int j=1; j<=(int)floor(_halfLength/_wireLengthStep); ++j ) {
                          tmpL = ((double) j)*_wireLengthStep;
                          tmpX = tmpL*_direct.getX()/*_absSectAngleTablep_it->second.first*/;
                          tmpY = tmpL*_direct.getY()/*_absSectAngleTablep_it->second.second*/;
                          _wireXYTable.push_back( pair<double, double>(_midp.getX()+tmpX, _midp.getY()+tmpY) );
                          _wireXYTable.push_back( pair<double, double>(_midp.getX()-tmpX, _midp.getY()-tmpY) );
                  }
          }
  };

  typedef std::map<int, HitData *> strawData;      //key = StrawIdx.asInt()
  //typedef std::unordered_set<int> strawList;       //list of StawIdx
  typedef std::set<int> strawList;                 //list of StawIdx
  typedef std::map<int, strawList> voteArrHitMap;  //map Straws that vote for a bin, key = binId

  typedef std::multimap<int, std::pair< int, std::pair<int, int> > > voteArrComMap; //map Straws Triple Combination that votes for a bin, key = binId

  void computeHistoProf(TH2F *hin, float &meanHeight, float &sigmaHeight ) {

          //cout<<"I'm computing mean and sigma of the heights of the histog "<<hin->GetName()<<endl;
          //cout<<"Entries in it "<<hin->GetEntries()<<" nBinX "<<hin->GetNbinsX()<<" nBinY "<<hin->GetNbinsY()<<endl;
          size_t tmpBin, ibin;
          int nBinXeff = hin->GetNbinsX()+2;
          Float_t *hin_arr = hin->GetArray();

          float nHits;
          nHits=meanHeight=sigmaHeight=0.000000;

          for (int nBinY=1; nBinY<=hin->GetNbinsY(); nBinY++ ){
                  tmpBin=nBinY*nBinXeff;
                  for (int nBinX=1; nBinX<=hin->GetNbinsX(); nBinX++ ){
                          ibin=tmpBin+nBinX;
                          if (hin_arr[ibin]>0.00000) {
                                  nHits+=1.000000;
                                  meanHeight+=hin_arr[ibin];
                                  sigmaHeight+=pow(hin_arr[ibin],2);
                                  //cout<<"x "<<nBinX<<" y "<<nBinY<<" val "<<hin_arr[ibin]<<endl;
                          }
                  }
          }
          //cout<<"In Histogr, non zero hit "<<nHits<<" sum "<<meanHeight<<" sum2 "<<sigmaHeight<<endl;

          if (nHits>1.0) {
                  sigmaHeight=sqrt( ( nHits*sigmaHeight - pow(meanHeight,2) )/( nHits*( nHits - 1.0 ) ) );
                  meanHeight/=nHits;
          }
          else { sigmaHeight=0.000000; }
          //cout<<"Mean "<<meanHeight<<" Sigma "<<sigmaHeight<<endl;
  }


  class TrackReco : public art::EDAnalyzer {
  public:

    explicit TrackReco(fhicl::ParameterSet const& pset);
    virtual ~TrackReco() {
            if (_plotCanvas)       delete _plotCanvas;
            if (_plotCanvas_1)     delete _plotCanvas_1;
            if (_plotCanvas_Cal)   delete _plotCanvas_Cal;
            if (_fakeCanvas)       delete _fakeCanvas;
            if (_iXHelStep!=0)     delete [] _iXHelStep;
            if (_iXHelStep_DR!=0)  delete [] _iXHelStep_DR;
            if (_iXPhi0_DR!=0)     delete [] _iXPhi0_DR;
            if (blur_coeffX!=0)    delete [] blur_coeffX;
            if (blur_coeffY!=0)    delete [] blur_coeffY;
    }

    virtual void beginJob();
    void endJob();

    // This is called for each event.
    void analyze(art::Event const& e);

  private:

    // Start: run time parameters

//    // The module label of this module.
//    std::string _moduleLabel;
//
//    // Label of the G4 module
//    std::string _g4ModuleLabel;
//
//    // Name of the tracker StepPoint collection
//    std::string _trackerStepPoints;

    // Label of the module that made the hits.
//    std::string _makerModuleLabel;

//    // Label of the generator.
//    std::string _generatorModuleLabel;

    // Label of the module that made the hits.
//    std::string _extractElectronsData;

    // Label of the module that made the hits.
//    std::string _timeRejecterModuleLabel;

    // Label of the module that made the hits.
//    std::string _geomRejecterModuleLabel;

    // Label of the module that made the hits.
    std::string _remappingModuleLabel;

    bool _doDisplay;
    bool _doCalib;

    // End: run time parameters
    //double removeAnglePeriod(double &angle);

    std::vector< std::pair<int, int> > _hitsCouplings; //hits couples by StrawIdx.asInt()
    void computeCombination( double minHStep=0.0, double maxHStep=0.0 );

    std::map< int, std::multimap<int, int> > _hitsTripleCouplings; //hits couples by StrawIdx.asInt()
    void computeTripleCombination( double minHStep=0.0, double maxHStep=0.0 );
    void computeTripleCombination_DeltaRay();
    void measureRbyTripleCombination( double HStep, double HPhi0 );

    int butterflyFilter( float *in, float *out, int &nBinX, int &nBinY, float minCountCut );
    int butterflyFilterRot( float *in, float *out, int &nBinX, int &nBinY, float minCountCut );
    int butterflyFilterRot45( float *in, float *out, int &nBinX, int &nBinY, float minCountCut );

    int blur( float *in, float *out, int &nBinX, int &nBinY, float sigmaX, float sigmaY, float binSizeX, float binSizeY );
    bool  blur_firstTime;
    int   blur_nBinXHalfWdt;
    int   blur_nBinYHalfWdt;
    int   blur_nBinXHalfWdt1;
    int   blur_nBinYHalfWdt1;
    float *blur_coeffX;
    float *blur_coeffY;

    int smooth( float *in, float *out, int &nBinX, int &nBinY, float sigmaX, float sigmaY, float binSizeX, float binSizeY, float mincut=0.0 );

    int peakFinder(TH2F *inHist, float peakExpSigma, float peakThreshold, float* peakPositionX,float* peakPositionY, int nIter=10 );

    double _minR;
    double _maxR;
    double _stepR;
    double _minLambda;
    double _maxLambda;
    double _stepLambda;
    double _minHelStep;
    double _maxHelStep;
    double _stepHelpStep;
    double _minMom;
    double _maxMom;
    double _rStraw;
    int    _nBinR;
    int    _nBinTanLambda;
    int    _nBinHelpStep;
    double _Bfield;

    double *_iXHelStep;
    strawData *strdat;
    voteArrHitMap voteArr_Bin_Hit_rel;

    int    _nBinHelpStep_DR;
    double *_iXHelStep_DR;
    int    _nBinPhi0_DR;
    double *_iXPhi0_DR;
    voteArrHitMap voteArr_Bin_Hit_rel_DR;
    voteArrComMap voteArr_Bin_Comb_rel_DR;
    voteArrHitMap voteArr_Phi0Bin_Hit_rel_DR;


    TCanvas*      _plotCanvas;
    TCanvas*      _plotCanvas_1;
    TCanvas*      _plotCanvas_Cal;
    TCanvas*      _fakeCanvas;

    TH2F*         _hRHelStep;
    TH1F*         _hR_1;
    TH1F*         _hR;
    TH1F*         _hPhi0;
    TH2F*         _hPhi0HelStep_1;
    TH2F*         _hPhi0HelStep;
    double        _maxPhi0HelStep;
    TH1F*         _hPhi0HelStep_Height;
    TH1F*         _hPhi0HelStep_X;
    TH1F*         _hPhi0HelStep_Y;
    TH1F*         _hPhi0HelStep_SigmaX;
    TH1F*         _hPhi0HelStep_SigmaY;
    TH1F*         _hPhi0HelStep_Cut;
    TH2F*         _hPhi0HelStepL;
    TH2F*         _BF_hPhi0HelStepL;
    TH2F*         _BFRot_hPhi0HelStepL;

    TH2F*         _hPhi0HelStep_DeltaRay;
    TH2F*         _hRhoPhi0_DeltaRay;
    TH2F*         _hRhoPhi0_DeltaRay_Blr;
    TH2F*         _hAbsSctrAbsZ_DeltaRay;
    TH2F*         _hAbsSctrRadi_DeltaRay;
    TH2F*         _hAbsSctrAbsZ_DeltaRay_Cum;
    TH2F*         _hAbsSctrRadi_DeltaRay_Cum;

    //    // Pointers to histograms, ntuples, TGraphs.
//    TH1F*         _hNtrackableEv;
//    TH1F*         _hNhitTrackableEv;
//    TH1F*         _hNBestTimePeakEv;
//    TH1F*         _hNhitBestTimePeakEv;
//    TH2F*         _hLostHitBstTmPkEv;
//    TH1F*         _hNoiseHitBstTmPkEv;
//
//    TH1F*         _hNBestSGClInTmPkEv;
//    TH1F*         _hSGClHitInBstTmPkEv;
//    TH2F*         _hLostHitBstSGClBstTmPkEv;
//    TH2F*         _hLostHitBstSGClEv;
//    TH1F*         _hNoiseHitBstSGClEv;
//    TH2F*         _hNLoopResBstSGClEv;
//
//    TTree*        _dataEvalBkgRejec;
    unsigned int  runID, eventID;// evNHit, convElNHit, convElNLoop;
//    float         sel_ptMeV, sel_ptMeV_start, sel_ptMeV_end, sel_plMeV, sel_plMeV_start, sel_plMeV_end, convElFHitTime;
//    int           bestTPcElNHit, bestTPNHit, nPeaksFound;
//    int           bestClGcElNHit, bestClGNHit, bestClGcNCls, nPotentTracks;

    // Some ugly but necessary ROOT related bookkeeping:

    // The job needs exactly one instance of TApplication.  See note 1.
    unique_ptr<TApplication> _application;

    // Save directory from beginJob so that we can go there in endJob. See note 3.
    TDirectory* _directory;

  };

  TrackReco::TrackReco(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),

    // Run time parameters
    //_makerModuleLabel(pset.get<string>("makerModuleLabel")),
    //_extractElectronsData(pset.get<string>("elextractModuleLabel")),
    //_timeRejecterModuleLabel(pset.get<string>("tRejecterModuleLabel")),
    //_geomRejecterModuleLabel(pset.get<string>("gRejecterModuleLabel")),
    _remappingModuleLabel(pset.get<string>("reMapperSHModuleLabel")),
//    _moduleLabel(pset.get<string>("module_label")),/*@module_label*/
//    _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
//    _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
//    _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
    _doDisplay(pset.get<bool>("doDisplay",false)),
    _doCalib(pset.get<bool>("doCalib",false)),
    blur_coeffX(0),
    blur_coeffY(0),

    _minR(280.0/2.0*CLHEP::mm),
    _maxR(800.0/2.0*CLHEP::mm),
    _stepR(2.5*1.0*CLHEP::mm),
    _minLambda(30.0*CLHEP::deg),
    _maxLambda(65.0*CLHEP::deg),
    _stepLambda(5.0*CLHEP::mrad),
    _stepHelpStep(2.0*5.0/sqrt(12.0)*CLHEP::mm),
    _minMom(50.0*CLHEP::MeV),
    _maxMom(150.0*CLHEP::MeV),
    _rStraw(2.50*CLHEP::mm),
    _nBinR(0),
    _nBinTanLambda(0),
    _nBinHelpStep(0),
    _Bfield(1.0),
    _iXHelStep(0x0),
    _nBinHelpStep_DR(0),
    _iXHelStep_DR(0x0),
    _nBinPhi0_DR(0),
    _iXPhi0_DR(0x0),

    _plotCanvas(0),
    _plotCanvas_1(0),
    _plotCanvas_Cal(0),
    _fakeCanvas(0),
    _hRHelStep(0),
    _hR_1(0),
    _hR(0),
    _hPhi0(0),
    _hPhi0HelStep_1(0),
    _hPhi0HelStep(0),
    _maxPhi0HelStep(0.0),
    _hPhi0HelStep_Height(0),
    _hPhi0HelStep_X(0),
    _hPhi0HelStep_Y(0),
    _hPhi0HelStep_SigmaX(0),
    _hPhi0HelStep_SigmaY(0),
    _hPhi0HelStep_Cut(0),
    _hPhi0HelStepL(0),
    _BF_hPhi0HelStepL(0),
    _BFRot_hPhi0HelStepL(0),
    _hPhi0HelStep_DeltaRay(0),
    _hRhoPhi0_DeltaRay(0),
    _hRhoPhi0_DeltaRay_Blr(0),
    _hAbsSctrAbsZ_DeltaRay(0),
    _hAbsSctrRadi_DeltaRay(0),
    _hAbsSctrAbsZ_DeltaRay_Cum(0),
    _hAbsSctrRadi_DeltaRay_Cum(0),
//    _hNtrackableEv(0),
//    _hNhitTrackableEv(0),
//    _hNBestTimePeakEv(0),
//    _hNhitBestTimePeakEv(0),
//    _hLostHitBstTmPkEv(0),
//    _hNoiseHitBstTmPkEv(0),
//
//    _hNBestSGClInTmPkEv(0),
//    _hSGClHitInBstTmPkEv(0),
//    _hLostHitBstSGClBstTmPkEv(0),
//    _hLostHitBstSGClEv(0),
//    _hNoiseHitBstSGClEv(0),
//    _hNLoopResBstSGClEv(0),
//    _dataEvalBkgRejec(0),

    // Some ugly but necessary ROOT related bookkeeping.
    //_application(nullptr),
    _directory(0){
          runID=eventID=/*evNHit=*/0;
          _minHelStep=CLHEP::twopi*_minMom/CLHEP::GeV*sin(_minLambda)/(0.3*_Bfield);
          _minHelStep*=CLHEP::m;
          _minHelStep=floor(_minHelStep);
          _maxHelStep=CLHEP::twopi*_maxMom/CLHEP::GeV*sin(_maxLambda)/(0.3*_Bfield);
          _maxHelStep*=CLHEP::m;
          _maxHelStep=ceil(_maxHelStep);
          //sel_ptMeV=sel_ptMeV_start=sel_ptMeV_end=sel_plMeV=sel_plMeV_start=sel_plMeV_end=convElFHitTime=0.0;
          //bestTPcElNHit=bestTPNHit=nPeaksFound=0;
          //bestClGcElNHit=bestClGNHit=bestClGcNCls=nPotentTracks=0;

          blur_firstTime = true;
 }

  void TrackReco::beginJob(){

          cout<<"Starting TrackReco jos!"<<endl;

    // Get access to the TFile service and save current directory for later use.
    art::ServiceHandle<art::TFileService> tfs;

    if(_doDisplay) {
            // If needed, create the ROOT interactive environment. See note 1.
            if ( !gApplication ){
                    int    tmp_argc(0);
                    char** tmp_argv(0);
                    _application = unique_ptr<TApplication>(new TApplication( "noapplication", &tmp_argc, tmp_argv ));
            }

            gStyle->SetPalette(1);
            gROOT->SetStyle("Plain");
            if (_doCalib){
                    gStyle->SetOptFit(1);
            }

            //_peaksCanvHistos   = new TObjArray();

            _plotCanvas = new TCanvas("plots","Hough Transform plots container",1290,860);
            _plotCanvas->Divide(2,2);
            _plotCanvas_1 = new TCanvas("plots1","Hough Transform plots container",1290,860);
            _plotCanvas_1->Divide(2,2);
            if (_doCalib) {
                    _plotCanvas_Cal = new TCanvas("plots1","Sigma evaluation of the Hough Transform plots container",1290,860);
                    _plotCanvas_Cal->Divide(3,2);
            }
           _fakeCanvas = new TCanvas("canvas_Fake","double click for next event",500,100);

    }

    // Create a histogram.
    double tmpNBin = (_maxR-_minR)/_stepR;
    _nBinR = (int)floor(tmpNBin+0.5);
    tmpNBin = (_maxR-_minR)/_stepR;
    _nBinTanLambda = (int)floor(tmpNBin+0.5);
    tmpNBin = (_maxHelStep-_minHelStep)/_stepHelpStep;
    _nBinHelpStep = floor(tmpNBin+0.5);

    cout<<"------------ Helix Step: "<<_minHelStep<<" "<<_maxHelStep<<" "<<_nBinHelpStep<<endl;
    _hRHelStep              = tfs->make<TH2F>( "hRHelStep",   "R vs Helix Step", _nBinHelpStep, _minHelStep, _maxHelStep, /*1200, -100,500*/_nBinR, _minR, _maxR );
    _iXHelStep = new double[_nBinHelpStep];
    for ( int ibin=0; ibin<_nBinHelpStep; ++ibin ) {
            _iXHelStep[ibin] = _hRHelStep->GetXaxis()->GetBinCenter(ibin+1);
    }
    _hPhi0                  = tfs->make<TH1F>( "hPhi0",   "Phi0", 400/*523*//*1047*/, 0.0, CLHEP::twopi/*-CLHEP::pi, CLHEP::pi*/ );
    _hR                     = tfs->make<TH1F>( "hR",   "R", _nBinR, _minR, _maxR );
    _hR_1                   = tfs->make<TH1F>( "hR_1",   "R", _nBinR, _minR, _maxR );
    _hPhi0HelStep_1         = tfs->make<TH2F>( "hPhi0HelStep_1",   "Phi0 vs Helix Step", _nBinHelpStep, _minHelStep, _maxHelStep, 200/*160*//*523*//*1047*/, 0.0, CLHEP::twopi/*-CLHEP::pi, CLHEP::pi*/ );


    _hPhi0HelStep           = tfs->make<TH2F>( "hPhi0HelStep",   "Phi0 vs Helix Step", _nBinHelpStep, _minHelStep, _maxHelStep, 100/*80*//*523*//*1047*/, -CLHEP::halfpi, CLHEP::halfpi );
    if (_doCalib){
            _hPhi0HelStep_Height    = tfs->make<TH1F>( "hPhi0HelStep_Height",   "Entries Projection of Phi0 vs Helix Step", 99, 0.01, 1 );
            _hPhi0HelStep_X         = tfs->make<TH1F>( "hPhi0HelStep_X",   "X Projection of Phi0 vs Helix Step", _nBinHelpStep, _minHelStep, _maxHelStep );
            _hPhi0HelStep_Y         = tfs->make<TH1F>( "hPhi0HelStep_Y",   "Y Projection of Phi0 vs Helix Step", 125/*80*//*523*//*1047*/, -CLHEP::halfpi, 1.5*CLHEP::halfpi );
            _hPhi0HelStep_SigmaX    = tfs->make<TH1F>( "hPhi0HelStep_SigmaX",   "Sigma of X Projection of Phi0 vs Helix Step", 100, 0, 20 );
            _hPhi0HelStep_SigmaY    = tfs->make<TH1F>( "hPhi0HelStep_SigmaY",   "Sigma of Y Projection of Phi0 vs Helix Step", 100, 0, 1 );
            _hPhi0HelStep_Cut       = tfs->make<TH1F>( "hPhi0HelStep_Cut",   "Peak heigth % cutof Projection of Phi0 vs Helix Step", 20, 0.0, 1 );
    }
    _hPhi0HelStepL          = tfs->make<TH2F>( "hPhi0HelStepL",  "Phi0 vs Helix Step (enlarged)", _nBinHelpStep, _minHelStep, _maxHelStep, 125/*100*//*523*//*1047*/, -CLHEP::halfpi, 1.5*CLHEP::halfpi );
    _BF_hPhi0HelStepL       = tfs->make<TH2F>( "BF_hPhi0HelStepL",  "Butterfly Filter of Phi0 vs Helix Step (enlarged)", _nBinHelpStep, _minHelStep, _maxHelStep, 125/*100*//*523*//*1047*/, -CLHEP::halfpi, 1.5*CLHEP::halfpi );
    _BFRot_hPhi0HelStepL    = tfs->make<TH2F>( "BFRot_hPhi0HelStepL",  "Butterfly Filter of Phi0 vs Helix Step (enlarged)", _nBinHelpStep, _minHelStep, _maxHelStep, 125/*100*//*523*//*1047*/, -CLHEP::halfpi, 1.5*CLHEP::halfpi );

    _nBinHelpStep_DR = floor((15000.0)/100.0+0.5);
    _hPhi0HelStep_DeltaRay  = tfs->make<TH2F>( "hPhi0HelStep_DeltaRay",   "Phi0 vs Helix Step for DeltaRay", _nBinHelpStep_DR, 5000.0, 20000.0, 125, -CLHEP::halfpi, 1.5*CLHEP::halfpi );
    _iXHelStep_DR = new double[_nBinHelpStep_DR];
    for ( int ibin=0; ibin<_nBinHelpStep_DR; ++ibin ) {
            _iXHelStep_DR[ibin] = _hPhi0HelStep_DeltaRay->GetXaxis()->GetBinCenter(ibin+1);
    }

    //const Tracker& tracker = getTrackerOrThrow();
    //const TTracker &ttr = static_cast<const TTracker&>( tracker );

    //int nbinRho=floor( (ttr.getTrackerEnvelopeParams().outerRadius()-ttr.getTrackerEnvelopeParams().innerRadius())/ttr.strawRadius()*sqrt(12.0)+0.5 );
    int nbinRho=floor( 410.0/40.0*sqrt(12.0)+0.5 );
    //_hRhoPhi0_DeltaRay      = tfs->make<TH2F>( "hRhoPhi0_DeltaRay",   "Rho vs Phi0 for low radius DeltaRay", 200, 0, CLHEP::twopi, nbinRho, ttr.getTrackerEnvelopeParams().innerRadius(), ttr.getTrackerEnvelopeParams().outerRadius() );
    _nBinPhi0_DR=floor( CLHEP::twopi/0.1*sqrt(12.0)+0.5 );
    _hRhoPhi0_DeltaRay      = tfs->make<TH2F>( "hRhoPhi0_DeltaRay",   "Rho vs Phi0 for low radius DeltaRay", _nBinPhi0_DR, 0, CLHEP::twopi, nbinRho, 390.0, 800.0 );
    _iXPhi0_DR = new double[_nBinPhi0_DR];
    for ( int ibin=0; ibin<_nBinPhi0_DR; ++ibin ) {
            _iXPhi0_DR[ibin] = _hRhoPhi0_DeltaRay->GetXaxis()->GetBinCenter(ibin+1);
    }
    _hRhoPhi0_DeltaRay_Blr      = tfs->make<TH2F>( "hRhoPhi0_DeltaRay_Blr",   "Rho vs Phi0 for low radius DeltaRay", _nBinPhi0_DR, 0, CLHEP::twopi, nbinRho, 390.0, 800.0 );
    _hAbsSctrAbsZ_DeltaRay      = tfs->make<TH2F>( "hAbsSctrAbsZ_DeltaRay",   "DeltaRay hit distribution in AbsSector and AbsZ", 15, 0, 15, 72, 0, 72 );
    _hAbsSctrRadi_DeltaRay      = tfs->make<TH2F>( "hAbsSctrRadi_DeltaRay",   "DeltaRay hit distribution in AbsSector and Radii", 12, 0, 12, nbinRho, 390.0, 800.0 );
    _hAbsSctrAbsZ_DeltaRay_Cum  = tfs->make<TH2F>( "hAbsSctrAbsZ_DeltaRay_Cum",   "DeltaRay hit distribution in AbsSector and AbsZ", 15, 0, 15, 72, 0, 72 );
    _hAbsSctrRadi_DeltaRay_Cum  = tfs->make<TH2F>( "hAbsSctrRadi_DeltaRay_Cum",   "DeltaRay hit distribution in AbsSector and Radii", 12, 0, 12, nbinRho, 390.0, 800.0 );

    _hRHelStep->SetXTitle("Step [mm]");
    _hRHelStep->SetYTitle("R [mm]");
    _hPhi0->SetXTitle("#phi_{0} [rad]");
    _hPhi0HelStepL->SetXTitle("Step [mm]");
    _hPhi0HelStepL->SetYTitle("#phi_{0} [rad]");
    _BF_hPhi0HelStepL->SetXTitle("Step [mm]");
    _BF_hPhi0HelStepL->SetYTitle("#phi_{0} [rad]");
    _hR_1->SetXTitle("R [mm]");
    _hR_1->SetYTitle("Entries");
    _hPhi0HelStep_DeltaRay->SetXTitle("Step [mm]");
    _hPhi0HelStep_DeltaRay->SetYTitle("#phi_{0} [rad]");
    _hRhoPhi0_DeltaRay->SetXTitle("#phi_{0} [rad]");
    _hRhoPhi0_DeltaRay->SetYTitle("Rho [mm]");
    _hRhoPhi0_DeltaRay_Blr->SetXTitle("#phi_{0} [rad]");
    _hRhoPhi0_DeltaRay_Blr->SetYTitle("Rho [mm]");


//    _hNtrackableEv       = tfs->make<TH1F>( "hNtrackableEv",   "N of trackable signal electrons", 111, 49.75, 105.25  );
//    _hNhitTrackableEv    = tfs->make<TH1F>( "hNhitTrackableEv",   "N of hits for each signal electrons track", 111, 49.75, 105.25  );
//    _hNBestTimePeakEv    = tfs->make<TH1F>( "hNBestTimePeakEv",   "N of peak time that has the best agreement with trackable signal electrons", 111, 49.75, 105.25  );
//    _hNhitBestTimePeakEv = tfs->make<TH1F>( "hNhitBestTimePeakEv",   "N of hits of the signal electrons track that are found in the best time peak", 111, 49.75, 105.25  );
//    _hLostHitBstTmPkEv   = tfs->make<TH2F>( "hLostHitBstTmPkEv",   "N of lost hits of the signal electrons track that in the best time peak", 111, 49.75, 105.25, 20, 0, 20  );
//    _hNoiseHitBstTmPkEv  = tfs->make<TH1F>( "hNoiseHitBstTmPkEv",   "N of hits of noise present in the the best time peak over signal electrons track hits", 111, 49.75, 105.25  );
//
//    _hNtrackableEv      ->SetXTitle("pt [MeV]");
//    _hNhitTrackableEv   ->SetXTitle("pt [MeV]");
//    _hNBestTimePeakEv   ->SetXTitle("pt [MeV]");
//    _hNhitBestTimePeakEv->SetXTitle("pt [MeV]");
//    _hLostHitBstTmPkEv  ->SetXTitle("pt [MeV]");
//    _hNoiseHitBstTmPkEv ->SetXTitle("pt [MeV]");
//
//    _hNBestSGClInTmPkEv  = tfs->make<TH1F>( "hNBestSGClInTmPkEv",   "N of Best Selected Clusters Group in the best peak time", 111, 49.75, 105.25  );
//    _hSGClHitInBstTmPkEv = tfs->make<TH1F>( "hSGClHitInBstTmPkEv",   "N of signal electron hit selected by best geom clustering algo from the best time peak", 111, 49.75, 105.25  );
//    _hLostHitBstSGClBstTmPkEv = tfs->make<TH2F>( "hLostHitBstSGClBstTmPkEv",   "N of lost hits of the signal electrons track that in the best geom cluster from the best time peak", 111, 49.75, 105.25, 20, 0, 20  );
//    _hLostHitBstSGClEv   = tfs->make<TH2F>( "hLostHitBstSGClEv",   "N of lost hits of the signal electrons track that in the best geom cluster", 111, 49.75, 105.25, 20, 0, 20  );
//    _hNoiseHitBstSGClEv  = tfs->make<TH1F>( "hNoiseHitBstSGClEv",   "N of hits of noise present in the the best geom cluster over signal electrons track hits", 111, 49.75, 105.25  );
//    _hNLoopResBstSGClEv  = tfs->make<TH2F>( "hNLoopResBstSGClEv",   "Res of N of track loop measured by the best geom cluster", 111, 49.75, 105.25, 21, -10.5, 10.5  );
//
//    _hNBestSGClInTmPkEv  ->SetXTitle("pt [MeV]");
//    _hSGClHitInBstTmPkEv ->SetXTitle("pt [MeV]");
//    _hLostHitBstSGClBstTmPkEv ->SetXTitle("pt [MeV]");
//    _hLostHitBstSGClEv ->SetXTitle("pt [MeV]");
//    _hNoiseHitBstSGClEv ->SetXTitle("pt [MeV]");
//    _hNLoopResBstSGClEv ->SetXTitle("pt [MeV]");
//
//    _dataEvalBkgRejec = tfs->make<TTree>( "dataEvalBkgRejec", "data for Bkg Rejection performance evaluation" );
//    _dataEvalBkgRejec->Branch("Run",&runID,"runID/i");
//    _dataEvalBkgRejec->Branch("Event",&eventID,"eventID/i");
//    _dataEvalBkgRejec->Branch("EvTotNHit",&evNHit,"evNHit/i");
//    _dataEvalBkgRejec->Branch("ConvElNHit",&convElNHit,"convElNHit/i");
//    _dataEvalBkgRejec->Branch("ConvElNLoop",&convElNLoop,"convElNLoop/i");
//    _dataEvalBkgRejec->Branch("ConvEl_ptMeV_Tracker",&sel_ptMeV,"sel_ptMeV/F");
//    _dataEvalBkgRejec->Branch("ConvEl_ptMeV_start",&sel_ptMeV_start,"sel_ptMeV_start/F");
//    _dataEvalBkgRejec->Branch("ConvEl_ptMeV_end",&sel_ptMeV_end,"sel_ptMeV_end/F");
//    _dataEvalBkgRejec->Branch("ConvEl_plMeV",&sel_plMeV,"sel_plMeV/F");
//    _dataEvalBkgRejec->Branch("ConvEl_plMeV_start",&sel_plMeV_start,"sel_plMeV_start/F");
//    _dataEvalBkgRejec->Branch("ConvEl_plMeV_end",&sel_plMeV_end,"sel_plMeV_end/F");
//    _dataEvalBkgRejec->Branch("ConvEl_FrstHit_Time",&convElFHitTime,"convElFHitTime/F");
//    _dataEvalBkgRejec->Branch("BestTPcElNHit",&bestTPcElNHit,"bestTPcElNHit/I");
//    _dataEvalBkgRejec->Branch("BestTPNHit",&bestTPNHit,"bestTPNHit/I");
//    _dataEvalBkgRejec->Branch("TotNPeaksFound",&nPeaksFound,"nPeaksFound/I");
//    _dataEvalBkgRejec->Branch("BestClGcElNHit",&bestClGcElNHit,"bestClGcElNHit/I");
//    _dataEvalBkgRejec->Branch("BestClGNHit",&bestClGNHit,"bestClGNHit/I");
//    _dataEvalBkgRejec->Branch("BestClGcNCls",&bestClGcNCls,"bestClGcNCls/I");
//    _dataEvalBkgRejec->Branch("NPotentTracksFound",&nPotentTracks,"nPotentTracks/I");


    //_fakeCanvas = new TCanvas("canvas_Fake","double click for next event",300,100);

    // See note 3.
    _directory = gDirectory;

  }

  void TrackReco::analyze(art::Event const& event ) {


    const Tracker& tracker = getTrackerOrThrow();
    const TTracker &ttr = static_cast<const TTracker&>( tracker );
    const std::vector<Device> ttrdev = ttr.getDevices();

//    art::Handle<StrawHitCollection> pdataHandle;
//    event.getByLabel(_makerModuleLabel,pdataHandle);
//    StrawHitCollection const* hits = pdataHandle.product();
//
//    art::Handle<VisibleGenElTrackCollection> genEltrksHandle;
//    event.getByLabel(_extractElectronsData,genEltrksHandle);
//    VisibleGenElTrackCollection const* genEltrks = genEltrksHandle.product();
//
//    art::Handle<TrackerHitTimeClusterCollection> tclustHandle;
//    event.getByLabel(_timeRejecterModuleLabel,tclustHandle);
//    TrackerHitTimeClusterCollection const* tclusts = tclustHandle.product();
//
//    art::Handle<SctrSttnClusterGroupCollection> gclusgtHandle;
//    event.getByLabel(_geomRejecterModuleLabel,gclusgtHandle);
//    SctrSttnClusterGroupCollection const* gclustgs = gclusgtHandle.product();
//
//    std::vector<mu2e::VisibleGenElTrack>::const_iterator genEltrk_it;

    art::Handle<ZRotStrawHitMapCollection> rmapdataHandle;
    event.getByLabel(_remappingModuleLabel,rmapdataHandle);
    ZRotStrawHitMapCollection const* mhits = rmapdataHandle.product();

    cout<<"--------------------------- Reco Track  ------------------------"<<endl;
    cout<<"Event: "<<event.id().event()<<endl;
    cout<<"----------------------------------------------------------------"<<endl;

    //StrawId sid;

    runID=event.run();
    eventID=event.event();
    //evNHit=hits->size();

    cout<<"N of group of Ttracker cluster of Hit found that could be tracks: "<<mhits->size()<<endl;
    //nPotentTracks=mhits->size();

    ZRotStrawHitMapCollection::const_iterator mhits_it;
    ZSectStrawHitMap::const_iterator zhitmap_tmp_it;
    AbsSectStrawHitMap::const_iterator scthitmap_tmp_it;

    ZStrawHitMap::const_iterator zInscthitmap_tmp1_it, zInscthitmap_tmp2_it;

    //int maxStationDist=5, maxContSect=5;
    //unsigned int negativeSectOver=2;
    //int snegativeSectOver=-((int)negativeSectOver);
    //int tmpSecDist;
    //bool goodCoupling;

    unsigned int nCoupling, nCouplingInGroup/*, checkNCoupling*/;
    std::vector<unsigned int> nCouplingForGroups;

    //unsigned int firstStation, tmpStation, lastStation, lastGroupStation;

    /*strawData **/strdat = new strawData;
    //StrawIndex firstStrawIdx, secondStrawIdx;

    //float binSizeX = _hPhi0HelStepL->GetXaxis()->GetBinWidth(1);
    //float binSizeY = _hPhi0HelStepL->GetYaxis()->GetBinWidth(1);
    float hPhi0HelS_HMean, hPhi0HelS_HSigma;
    //float peakSigmaX, peakSigmaY, peakSigma, peakSigmaBins;
    //peakSigmaX = 7.4/*2.0*binSizeX*/;
    //peakSigmaY = 0.04/*1.0*binSizeY*/;
    //peakSigma  = sqrt(pow(peakSigmaX,2)+pow(peakSigmaY,2));
    //peakSigmaBins  = sqrt(pow(peakSigmaX/binSizeX,2)+pow(peakSigmaY/binSizeY,2));
    //float peakPosX[maxNpeaks], peakPosY[maxNpeaks];
    //int npeaks;
    float peakPosX1[maxNpeaks], peakPosY1[maxNpeaks];
    int npeaks1;
    bool potentialDRfound, hitToSkip;
    float tmpDRRho, tmpiBinHit;
    size_t binx, biny, tmpBin, nBinsXEff, tmpbinx, tmpbiny;
    int nBinsX, nBinsY;
    float maxVote, tmpVote;
    float drPeakVote, tmpDrPeakVote, tmpBinDrPeak, tmpBinDRPeak, equalVoteCut;
    float sumAbsZ, sum2AbsZ, sumAbsSec, sum2AbsSec, nHitToRejec, nExtraHitInAbsSec;
    float meanAbsZ, sdAbsZ, meanAbsSec, sdAbsSec;
    float minAbsZ, maxAbsZ, minAbsSec, maxAbsSec;
    std::map<int, float* > sumsRvsAbsSec;  //key= AbsSectID, float[3] = {nHitR, sumR, sim2R}
    float meanR_AbsSec[12], sdR_AbsSec[12];
    float minR_AbsSec[12], maxR_AbsSec[12];
    bool noDrHitremoved, noDrHitremovedForZ, noDrHitremovedForS, noDrHitremovedForR;
    for (int iSec=0; iSec<12; iSec++) {
            sumsRvsAbsSec.insert( std::pair<int, float*>(iSec, new float[3]) );
    }

    std::set<int> equalMaxVote;
    strawList allHitInPeak, skipHit, peakSkipHit, potDrAlreadyDisp;
    int nHitDRcut=11;
    int nBinsX_DR, nBinsY_DR, nBinsXEff_DR, iterDR;
    nBinsX_DR = _hRhoPhi0_DeltaRay->GetNbinsX();
    nBinsY_DR = _hRhoPhi0_DeltaRay->GetNbinsY();
    float hRhoPhi0_HMean, hRhoPhi0_HSigma;

    for ( mhits_it=mhits->begin(); mhits_it!=mhits->end(); ++mhits_it ) {

            cout<<*mhits_it;
            //for ( ZSectStrawHitMap::const_iterator zhitmap_it = mhits_it->_zsctTrackerHits.begin(); zhitmap_it != mhits_it->_zsctTrackerHits.end(); ++zhitmap_it ) {
            //        for ( AbsSectStrawHitMap::const_iterator scthitmap_it =  zhitmap_it->second.begin(); scthitmap_it !=  zhitmap_it->second.end(); ++scthitmap_it ) {
            //                cout<<"Hit (z,sec): "<<zhitmap_it->first<<" "<<scthitmap_it->first<<endl;
            //        }
            //}
            //cout<<"---------------------------------"<<endl;

            for ( strawData::iterator strdat_it=strdat->begin(); strdat_it!=strdat->end(); ++strdat_it ) {
                    delete strdat_it->second;
            }
            strdat->clear();

            //----- straight line Delta Ray extraction ---------

            skipHit.clear();
            equalMaxVote.clear();
            potDrAlreadyDisp.clear();

            potentialDRfound=true;
            _hAbsSctrAbsZ_DeltaRay_Cum->Reset();
            _hAbsSctrRadi_DeltaRay_Cum->Reset();

            iterDR=0;
            while (potentialDRfound){
                    iterDR++;
                    voteArr_Phi0Bin_Hit_rel_DR.clear();
                    _hRhoPhi0_DeltaRay->Reset();
                    _hRhoPhi0_DeltaRay_Blr->Reset();
                    //_hAbsSctrAbsZ_DeltaRay->Reset();
                    maxVote=-1;

                    cout<<"---------------- N tot of hit to skip "<<skipHit.size()<<endl;
                    for ( SectZStrawHitMap::const_iterator secthitmap_it = mhits_it->_sctZTrackerHits.begin(); secthitmap_it != mhits_it->_sctZTrackerHits.end(); ++secthitmap_it ) {
                            for ( ZStrawHitMap::const_iterator zInscthitmap_it =  secthitmap_it->second.begin(); zInscthitmap_it !=  secthitmap_it->second.end(); ++zInscthitmap_it ) {
                                    StrawIndex const&firstStrawIdx=((StrawHitPtr) zInscthitmap_it->second)->strawIndex();
                                    if ( skipHit.find(firstStrawIdx.asInt())!=skipHit.end() ) continue;
                                    const Straw & fstr = ttr.getStraw(firstStrawIdx);
                                    strawData::iterator strdat_it;
                                    if (strdat->count(firstStrawIdx.asInt())==0) {
                                            strdat_it=strdat->insert(
                                                            strawData::value_type(
                                                                            firstStrawIdx.asInt(),
                                                                            new HitData(zInscthitmap_it->second, fstr.getDirection(), fstr.getMidPoint(), fstr.getHalfLength(), zInscthitmap_it->first, secthitmap_it->first)
                                                            )
                                            ).first;
                                    }
                                    else {
                                            strdat_it=strdat->find(firstStrawIdx.asInt());
                                    }
                                    for (int ibin=0; ibin<_nBinPhi0_DR; ibin++){
                                            tmpDRRho = strdat_it->second->_radius/cos( _iXPhi0_DR[ibin]-strdat_it->second->_theta );
                                            tmpiBinHit = _hRhoPhi0_DeltaRay->Fill(_iXPhi0_DR[ibin],tmpDRRho);
                                            if (tmpiBinHit>0) {
                                                    voteArrHitMap::iterator voteArr_Phi0Bin_Hit_rel_DR_it = voteArr_Phi0Bin_Hit_rel_DR.find(tmpiBinHit);
                                                    if ( voteArr_Phi0Bin_Hit_rel_DR_it==voteArr_Phi0Bin_Hit_rel_DR.end() ) {
                                                            voteArr_Phi0Bin_Hit_rel_DR_it = voteArr_Phi0Bin_Hit_rel_DR.insert( voteArrHitMap::value_type( tmpiBinHit, strawList()/* inBinHitList*/ ) ).first;
                                                            voteArr_Phi0Bin_Hit_rel_DR_it->second.clear();
                                                    }
                                                    voteArr_Phi0Bin_Hit_rel_DR_it->second.insert(strdat_it->first);
                                                    tmpVote = _hRhoPhi0_DeltaRay->GetBinCenter(tmpiBinHit);
                                                    if (tmpVote>maxVote) maxVote=tmpVote;
                                            }

                                    }
                            }

                    }
                    potentialDRfound=false;
                    if (maxVote>=nHitDRcut) {

                            //Float_t *hRhoPhi0_DeltaRay_arr = _hRhoPhi0_DeltaRay->GetArray();
                            //Float_t *hRhoPhi0_DeltaRay_Blr_arr = _hRhoPhi0_DeltaRay_Blr->GetArray();
                            //int nBinY_dr =  _hRhoPhi0_DeltaRay->GetNbinsY();
                            //_hRhoPhi0_DeltaRay_Blr->SetEntries( blur( hRhoPhi0_DeltaRay_arr, hRhoPhi0_DeltaRay_Blr_arr, _nBinPhi0_DR, nBinY_dr, 1.0, 1.0, 1.0, 1.0 ) );
                            //npeaks1 = peakFinder(_hRhoPhi0_DeltaRay_Blr, 4.2, 60, peakPosX1, peakPosY1, 2 );

                            computeHistoProf( _hRhoPhi0_DeltaRay, hRhoPhi0_HMean, hRhoPhi0_HSigma );
                            cout<<"hRhoPhi0_DeltaRay Prof val : mean "<<hRhoPhi0_HMean<<" sigma "<<hRhoPhi0_HSigma<<endl;


                            npeaks1 = peakFinder(_hRhoPhi0_DeltaRay, 4.2, 70, peakPosX1, peakPosY1, 10 );
                            nBinsXEff_DR = nBinsX_DR+2;
                            for(int ip=0; ip<npeaks1; ip++) {
                                    binx   = (int)floor(peakPosX1[ip]+0.5)+1;
                                    biny   = (int)floor(peakPosY1[ip]+0.5)+1;
                                    if ( (binx<1 || binx>((size_t)nBinsX_DR) ) || (biny<1 || biny>((size_t)nBinsY_DR)) ) continue;
                                    tmpBin = biny*nBinsXEff_DR + binx;
                                    drPeakVote = _hRhoPhi0_DeltaRay->GetBinContent(tmpBin);
                                    equalVoteCut = 0.9*drPeakVote;

                                    cout<<"Peak at pos "<<binx<<" = "<<_hRhoPhi0_DeltaRay->GetXaxis()->GetBinCenter(binx)<<" - "<<biny<<" = "<<_hRhoPhi0_DeltaRay->GetYaxis()->GetBinCenter(biny)<<" straws in peak:"<<endl;
                                    voteArrHitMap::iterator voteArr_Phi0Bin_Hit_rel_DR_it = voteArr_Phi0Bin_Hit_rel_DR.find(tmpBin);
                                    cout<<"Hits in peak "<<voteArr_Phi0Bin_Hit_rel_DR_it->second.size()<<endl;
                                    tmpBinDrPeak=tmpBin;
                                    for (int ixbin=-1; ixbin<=1; ixbin++){
                                            tmpbinx = binx+ixbin;
                                            if (tmpbinx<1) continue;
                                            else if (tmpbinx>((size_t)_nBinPhi0_DR)) continue;
                                            for (int iybin=-1; iybin<=1; iybin++){
                                                    tmpbiny = biny+iybin;
                                                    if (tmpbiny<1) continue;
                                                    else if (tmpbiny>((size_t)nBinsY_DR)) continue;
                                                    tmpBinDRPeak = tmpbiny*nBinsXEff_DR + tmpbinx;
                                                    cout<<"--- tmpbinx "<<tmpbinx<<" tmpbiny "<<tmpbiny<<" tmpBinDRPeak "<<tmpBinDRPeak<<endl;
                                                    tmpDrPeakVote = _hRhoPhi0_DeltaRay->GetBinContent(tmpBinDRPeak);
                                                    if (tmpDrPeakVote>drPeakVote) {
                                                            drPeakVote = tmpDrPeakVote;
                                                            equalVoteCut = 0.9*drPeakVote;
                                                            tmpBinDrPeak = tmpBinDRPeak;
                                                            equalMaxVote.clear();
                                                    }
                                                    else if (tmpDrPeakVote>=equalVoteCut/*tmpDrPeakVote==drPeakVote*/) {
                                                            equalMaxVote.insert(tmpBinDRPeak);
                                                    }
                                            }
                                    }
                                    if (tmpBinDrPeak!=tmpBin) {
                                            tmpBin=tmpBinDrPeak;
                                            cout<<"Higher peak in "<<tmpBin<<endl;
                                            voteArr_Phi0Bin_Hit_rel_DR_it = voteArr_Phi0Bin_Hit_rel_DR.find(tmpBin);
                                            cout<<"Hits in peak "<<voteArr_Phi0Bin_Hit_rel_DR_it->second.size()<<endl;
                                    }

                                    allHitInPeak.clear();
                                    peakSkipHit.clear();
                                    sumAbsZ=sum2AbsZ=sumAbsSec=sum2AbsSec=nHitToRejec=nExtraHitInAbsSec=0.00000000;
                                    std::map<int, float* >::iterator sumsRvsAbsSec_it;
                                    for ( sumsRvsAbsSec_it=sumsRvsAbsSec.begin(); sumsRvsAbsSec_it!=sumsRvsAbsSec.end(); ++sumsRvsAbsSec_it ){
                                            sumsRvsAbsSec_it->second[0]=0.00000000;
                                            sumsRvsAbsSec_it->second[1]=0.00000000;
                                            sumsRvsAbsSec_it->second[2]=0.00000000;
                                    }
                                    /*
                            for (int ixbin=-1; ixbin<=1; ixbin++){
                                    for (int iybin=-1; iybin<=1; iybin++){
                                            tmpBin = (biny+iybin)*nBinsXEff_DR + (binx+ixbin);
                                            voteArrHitMap::iterator voteArr_Phi0Bin_Hit_rel_DR_it2 = voteArr_Phi0Bin_Hit_rel_DR.find(tmpBin);
                                            for ( strawList::iterator inpeakStraws_it=voteArr_Phi0Bin_Hit_rel_DR_it2->second.begin(); inpeakStraws_it!=voteArr_Phi0Bin_Hit_rel_DR_it2->second.end(); ++inpeakStraws_it) {
                                                    allHitInPeak.insert(*inpeakStraws_it);
                                            }
                                    }
                            }
                                     */

                                    _hAbsSctrAbsZ_DeltaRay->Reset();
                                    _hAbsSctrRadi_DeltaRay->Reset();

                                    for ( strawList::iterator inpeakStraws_it=voteArr_Phi0Bin_Hit_rel_DR_it->second.begin(); inpeakStraws_it!=voteArr_Phi0Bin_Hit_rel_DR_it->second.end(); ++inpeakStraws_it) {
                                            allHitInPeak.insert(*inpeakStraws_it);
                                            peakSkipHit.insert(*inpeakStraws_it);
                                            HitData *drhit  = strdat->find(*inpeakStraws_it)->second;
                                            sumAbsZ+=drhit->_absZId;
                                            sum2AbsZ+=pow(drhit->_absZId,2);
                                            sumAbsSec+=drhit->_absSectId;
                                            sum2AbsSec+=pow(drhit->_absSectId,2);
                                            if (drhit->_absSectId<3) {
                                                    sumAbsSec+=(drhit->_absSectId+12);
                                                    sum2AbsSec+=pow(drhit->_absSectId+12,2);
                                                    nExtraHitInAbsSec+=1.000000;
                                            }
                                            sumsRvsAbsSec_it=sumsRvsAbsSec.find(drhit->_absSectId);
                                            sumsRvsAbsSec_it->second[0]+=1.000000;
                                            sumsRvsAbsSec_it->second[1]+=drhit->_radius;
                                            sumsRvsAbsSec_it->second[2]+=pow(drhit->_radius,2);
                                            if ( potDrAlreadyDisp.insert(*inpeakStraws_it).second ) {
                                                    _hAbsSctrAbsZ_DeltaRay->Fill(drhit->_absSectId,drhit->_absZId);
                                                    _hAbsSctrRadi_DeltaRay->Fill(drhit->_absSectId,drhit->_radius);
                                                    _hAbsSctrAbsZ_DeltaRay_Cum->Fill(drhit->_absSectId,drhit->_absZId);
                                                    _hAbsSctrRadi_DeltaRay_Cum->Fill(drhit->_absSectId,drhit->_radius);
                                                    if (drhit->_absSectId<3) {
                                                            _hAbsSctrAbsZ_DeltaRay->Fill(drhit->_absSectId+12,drhit->_absZId);
                                                            _hAbsSctrAbsZ_DeltaRay_Cum->Fill(drhit->_absSectId+12,drhit->_absZId);
                                                    }
                                            }
                                    }
                                    for ( std::set<int>::iterator equalMaxVote_it = equalMaxVote.begin(); equalMaxVote_it != equalMaxVote.end(); ++equalMaxVote_it ) {
                                            voteArrHitMap::iterator voteArr_Phi0Bin_Hit_rel_DR_it2 = voteArr_Phi0Bin_Hit_rel_DR.find(*equalMaxVote_it);
                                            for ( strawList::iterator inpeakStraws_it=voteArr_Phi0Bin_Hit_rel_DR_it2->second.begin(); inpeakStraws_it!=voteArr_Phi0Bin_Hit_rel_DR_it2->second.end(); ++inpeakStraws_it) {
                                                    allHitInPeak.insert(*inpeakStraws_it);
                                                    hitToSkip = peakSkipHit.insert(*inpeakStraws_it).second;
                                                    if (hitToSkip) {
                                                            HitData *drhit  = strdat->find(*inpeakStraws_it)->second;
                                                            sumAbsZ+=drhit->_absZId;
                                                            sum2AbsZ+=pow(drhit->_absZId,2);
                                                            sumAbsSec+=drhit->_absSectId;
                                                            sum2AbsSec+=pow(drhit->_absSectId,2);
                                                            if (drhit->_absSectId<3) {
                                                                    sumAbsSec+=(drhit->_absSectId+12);
                                                                    sum2AbsSec+=pow(drhit->_absSectId+12,2);
                                                                    nExtraHitInAbsSec+=1.000000;
                                                            }
                                                            sumsRvsAbsSec_it=sumsRvsAbsSec.find(drhit->_absSectId);
                                                            sumsRvsAbsSec_it->second[0]+=1.000000;
                                                            sumsRvsAbsSec_it->second[1]+=drhit->_radius;
                                                            sumsRvsAbsSec_it->second[2]+=pow(drhit->_radius,2);
                                                            if ( potDrAlreadyDisp.insert(*inpeakStraws_it).second ) {
                                                                    _hAbsSctrAbsZ_DeltaRay->Fill(drhit->_absSectId,drhit->_absZId);
                                                                    _hAbsSctrRadi_DeltaRay->Fill(drhit->_absSectId,drhit->_radius);
                                                                    _hAbsSctrAbsZ_DeltaRay_Cum->Fill(drhit->_absSectId,drhit->_absZId);
                                                                    _hAbsSctrRadi_DeltaRay_Cum->Fill(drhit->_absSectId,drhit->_radius);
                                                                    if (drhit->_absSectId<3) {
                                                                            _hAbsSctrAbsZ_DeltaRay->Fill(drhit->_absSectId+12,drhit->_absZId);
                                                                            _hAbsSctrAbsZ_DeltaRay_Cum->Fill(drhit->_absSectId+12,drhit->_absZId);
                                                                    }
                                                            }
                                                    }
                                            }
                                    }
                                    cout<<"Hits in peak in all equal value neighbors "<<allHitInPeak.size()<<endl;

                                    noDrHitremoved=true;
                                    while (noDrHitremoved) {
                                            noDrHitremoved=false;
                                            nHitToRejec=peakSkipHit.size();
                                            if (nHitToRejec>=nHitDRcut) {
                                                    potentialDRfound=true;
                                                    noDrHitremovedForZ=true;
                                                    while (noDrHitremovedForZ) {
                                                            noDrHitremovedForZ=false;
                                                            nHitToRejec = peakSkipHit.size();
                                                            meanAbsZ    = sumAbsZ/nHitToRejec;
                                                            sdAbsZ      = sqrt( ( nHitToRejec*sum2AbsZ - pow(sumAbsZ,2) )/( nHitToRejec*( nHitToRejec - 1.0 ) ) );
                                                            minAbsZ     = 1.732050808*sdAbsZ;
                                                            maxAbsZ     = floor(meanAbsZ+minAbsZ+0.5);
                                                            minAbsZ     = floor(meanAbsZ-minAbsZ+0.5);
                                                            cout<<"Z- Mean "<<meanAbsZ<<" sigma "<<sdAbsZ<<" min "<<minAbsZ<<" max "<<maxAbsZ<<endl;
                                                            for ( strawList::iterator inpeakStraws_it=peakSkipHit.begin(); inpeakStraws_it!=peakSkipHit.end(); ++inpeakStraws_it) {
                                                                    HitData *drhit  = strdat->find(*inpeakStraws_it)->second;
                                                                    if ( drhit->_absZId<minAbsZ || drhit->_absZId>maxAbsZ ) {
                                                                            noDrHitremoved=true;
                                                                            noDrHitremovedForZ=true;
                                                                            peakSkipHit.erase(inpeakStraws_it);
                                                                            --inpeakStraws_it;
                                                                            sumAbsZ-=drhit->_absZId;
                                                                            sum2AbsZ-=pow(drhit->_absZId,2);
                                                                            sumAbsSec-=drhit->_absSectId;
                                                                            sum2AbsSec-=pow(drhit->_absSectId,2);
                                                                            if (drhit->_absSectId<3) {
                                                                                    sumAbsSec-=(drhit->_absSectId+12);
                                                                                    sum2AbsSec-=pow(drhit->_absSectId+12,2);
                                                                                    nExtraHitInAbsSec-=1.000000;
                                                                            }
                                                                            sumsRvsAbsSec_it=sumsRvsAbsSec.find(drhit->_absSectId);
                                                                            sumsRvsAbsSec_it->second[0]-=1.000000;
                                                                            sumsRvsAbsSec_it->second[1]-=drhit->_radius;
                                                                            sumsRvsAbsSec_it->second[2]-=pow(drhit->_radius,2);
                                                                    }
                                                            }
                                                    }
                                                    noDrHitremovedForS=true;
                                                    while (noDrHitremovedForS) {
                                                            noDrHitremovedForS=false;
                                                            if (nExtraHitInAbsSec<0) cerr<<"!!!!!!!!!!!!!!! Wrong nExtraHitInAbsSec !!!!!!!!!!!!!!! "<<nExtraHitInAbsSec<<endl;
                                                            nHitToRejec = peakSkipHit.size()+nExtraHitInAbsSec;
                                                            meanAbsSec  = sumAbsSec/nHitToRejec;
                                                            sdAbsSec    = sqrt( ( nHitToRejec*sum2AbsSec - pow(sumAbsSec,2) )/( nHitToRejec*( nHitToRejec - 1.0 ) ) );
                                                            minAbsSec   = 1.732050808*sdAbsSec;
                                                            maxAbsSec   = floor(meanAbsSec+minAbsSec+0.5);
                                                            minAbsSec   = floor(meanAbsSec-minAbsSec+0.5);
                                                            cout<<"S- Mean "<<meanAbsSec<<" sigma "<<sdAbsSec<<" min "<<minAbsSec<<" max "<<maxAbsSec<<endl;
                                                            for ( strawList::iterator inpeakStraws_it=peakSkipHit.begin(); inpeakStraws_it!=peakSkipHit.end(); ++inpeakStraws_it) {
                                                                    HitData *drhit  = strdat->find(*inpeakStraws_it)->second;
                                                                    if ( drhit->_absSectId<minAbsSec || drhit->_absSectId>maxAbsSec ) {
                                                                            noDrHitremoved=true;
                                                                            noDrHitremovedForS=true;
                                                                            peakSkipHit.erase(inpeakStraws_it);
                                                                            --inpeakStraws_it;
                                                                            sumAbsZ-=drhit->_absZId;
                                                                            sum2AbsZ-=pow(drhit->_absZId,2);
                                                                            sumAbsSec-=drhit->_absSectId;
                                                                            sum2AbsSec-=pow(drhit->_absSectId,2);
                                                                            if (drhit->_absSectId<3) {
                                                                                    sumAbsSec-=(drhit->_absSectId+12);
                                                                                    sum2AbsSec-=pow(drhit->_absSectId+12,2);
                                                                                    nExtraHitInAbsSec-=1.000000;
                                                                            }
                                                                            sumsRvsAbsSec_it=sumsRvsAbsSec.find(drhit->_absSectId);
                                                                            sumsRvsAbsSec_it->second[0]-=1.000000;
                                                                            sumsRvsAbsSec_it->second[1]-=drhit->_radius;
                                                                            sumsRvsAbsSec_it->second[2]-=pow(drhit->_radius,2);
                                                                    }
                                                            }
                                                    }
                                                    noDrHitremovedForR=true;
                                                    while (noDrHitremovedForR) {
                                                            noDrHitremovedForR=false;
                                                            for ( sumsRvsAbsSec_it=sumsRvsAbsSec.begin(); sumsRvsAbsSec_it!=sumsRvsAbsSec.end(); ++sumsRvsAbsSec_it ){
                                                                    nHitToRejec = sumsRvsAbsSec_it->second[0];
                                                                    if (nHitToRejec>1.0) {
                                                                            meanR_AbsSec[sumsRvsAbsSec_it->first] = sumsRvsAbsSec_it->second[1]/nHitToRejec;
                                                                            sdR_AbsSec[sumsRvsAbsSec_it->first]   = sqrt( ( nHitToRejec*sumsRvsAbsSec_it->second[2] - pow(sumsRvsAbsSec_it->second[1],2) )/( nHitToRejec*( nHitToRejec - 1.0 ) ) );
                                                                            minR_AbsSec[sumsRvsAbsSec_it->first]  = 1.732050808*sdR_AbsSec[sumsRvsAbsSec_it->first];
                                                                            maxR_AbsSec[sumsRvsAbsSec_it->first]  = meanR_AbsSec[sumsRvsAbsSec_it->first]+minR_AbsSec[sumsRvsAbsSec_it->first];
                                                                            minR_AbsSec[sumsRvsAbsSec_it->first]  = meanR_AbsSec[sumsRvsAbsSec_it->first]-minR_AbsSec[sumsRvsAbsSec_it->first];
                                                                            cout<<"R_"<<sumsRvsAbsSec_it->first<<" - Mean "<<meanR_AbsSec[sumsRvsAbsSec_it->first]<<" sigma "<<sdR_AbsSec[sumsRvsAbsSec_it->first]<<" min "<<minR_AbsSec[sumsRvsAbsSec_it->first]<<" max "<<maxR_AbsSec[sumsRvsAbsSec_it->first]<<endl;
                                                                    }
                                                            }
                                                            for ( strawList::iterator inpeakStraws_it=peakSkipHit.begin(); inpeakStraws_it!=peakSkipHit.end(); ++inpeakStraws_it) {
                                                                    HitData *drhit  = strdat->find(*inpeakStraws_it)->second;
                                                                    if ( drhit->_radius<minR_AbsSec[drhit->_absSectId] || drhit->_radius>maxR_AbsSec[drhit->_absSectId] ) {
                                                                            noDrHitremoved=true;
                                                                            noDrHitremovedForR=true;
                                                                            peakSkipHit.erase(inpeakStraws_it);
                                                                            --inpeakStraws_it;
                                                                            sumAbsZ-=drhit->_absZId;
                                                                            sum2AbsZ-=pow(drhit->_absZId,2);
                                                                            sumAbsSec-=drhit->_absSectId;
                                                                            sum2AbsSec-=pow(drhit->_absSectId,2);
                                                                            if (drhit->_absSectId<3) {
                                                                                    sumAbsSec-=(drhit->_absSectId+12);
                                                                                    sum2AbsSec-=pow(drhit->_absSectId+12,2);
                                                                                    nExtraHitInAbsSec-=1.000000;
                                                                            }
                                                                            sumsRvsAbsSec_it=sumsRvsAbsSec.find(drhit->_absSectId);
                                                                            sumsRvsAbsSec_it->second[0]-=1.000000;
                                                                            sumsRvsAbsSec_it->second[1]-=drhit->_radius;
                                                                            sumsRvsAbsSec_it->second[2]-=pow(drhit->_radius,2);
                                                                    }
                                                            }
                                                    }
                                            }
                                            else {
                                                    potentialDRfound=false;
                                                    peakSkipHit.clear();
                                            }
                                    }

                                    cout<<"Hits selected to be a DR "<<peakSkipHit.size()<<endl;
                                    for ( strawList::iterator inpeakStraws_it=peakSkipHit.begin(); inpeakStraws_it!=peakSkipHit.end(); ++inpeakStraws_it) {
                                            ///hitToSkip = skipHit.insert(*inpeakStraws_it).second;
                                            skipHit.insert(*inpeakStraws_it);
                                            HitData *drhit  = strdat->find(*inpeakStraws_it)->second;
                                            cout<<"\t"<<*inpeakStraws_it<<"\t"<<drhit->_midp<<endl;
                                    }

                                    /*
                            if (voteArr_Phi0Bin_Hit_rel_DR_it->second.size()>=15) {
                                    for ( strawList::iterator inpeakStraws_it=voteArr_Phi0Bin_Hit_rel_DR_it->second.begin(); inpeakStraws_it!=voteArr_Phi0Bin_Hit_rel_DR_it->second.end(); ++inpeakStraws_it) {
                                            HitData *drhit  = strdat->find(*inpeakStraws_it)->second;
                                            _hAbsSctrAbsZ_DeltaRay->Fill(drhit->_absSectId,drhit->_absZId);
                                            cout<<"\t"<<*inpeakStraws_it<<"\t"<<drhit->_midp<<endl;
                                    }
                            }
                                     */

                                    if (_doDisplay) {

                                            _plotCanvas_1->cd(1);
                                            _hRhoPhi0_DeltaRay->Draw("col z");
                                            _plotCanvas_1->cd(2);
                                            _hRhoPhi0_DeltaRay_Blr->Draw("col z");
                                            //_hPhi0HelStep_DeltaRay->Draw("col z");
                                            _plotCanvas_1->cd(3);
                                            _hAbsSctrAbsZ_DeltaRay->Draw("col z");
                                            _plotCanvas_1->cd(4);
                                            _hAbsSctrRadi_DeltaRay->Draw("col z");
                                            _plotCanvas_1->Update();

                                            cerr << "Double click in the canvas_Fake to continue:"<<endl ;
                                            _fakeCanvas->Clear();
                                            _fakeCanvas->cd();
                                            TLatex *printEvN = new TLatex(0.15,0.4,Form("Current Event: %d _DR iter %d _peak %d",event.id().event(),iterDR,ip));
                                            printEvN->SetTextFont(62);
                                            printEvN->SetTextSizePixels(180);
                                            printEvN->Draw();
                                            _fakeCanvas->Update();
                                            _fakeCanvas->WaitPrimitive();
                                            cerr << endl;
                                            delete printEvN;
                                    }

                            }

                    }

                    //Float_t *hRhoPhi0_DeltaRay_arr = _hRhoPhi0_DeltaRay->GetArray();
                    //Float_t *hRhoPhi0_DeltaRay_Blr_arr = _hRhoPhi0_DeltaRay_Blr->GetArray();
                    //int nBinY_dr =  _hRhoPhi0_DeltaRay->GetNbinsY();
                    //_hRhoPhi0_DeltaRay_Blr->SetEntries( blur( hRhoPhi0_DeltaRay_arr, hRhoPhi0_DeltaRay_Blr_arr, _nBinPhi0_DR, nBinY_dr, 1.0, 1.0, 1.0, 1.0 ) );
                    //npeaks1 = peakFinder(_hRhoPhi0_DeltaRay_Blr, 4.2, 60, peakPosX1, peakPosY1, 2 );

            }

            //----- end straight line Delta Ray extraction ---------
            cout<<"---------------- N tot of hits assigned into DRs "<<skipHit.size()<<endl;


            _hitsCouplings.clear();
            _hitsTripleCouplings.clear();

            _hRHelStep->Reset();
            _hPhi0->Reset();
            _hR->Reset();
            _hR_1->Reset();
            _hPhi0HelStep_1->Reset();

            _hPhi0HelStep->Reset();
            _hPhi0HelStepL->Reset();
            _BF_hPhi0HelStepL->Reset();
            _BFRot_hPhi0HelStepL->Reset();

            _hPhi0HelStep_DeltaRay->Reset();

            voteArr_Bin_Hit_rel.clear();
            voteArr_Bin_Hit_rel_DR.clear();
            voteArr_Bin_Comb_rel_DR.clear();
    //if(false){
            nCoupling=0;
            nCouplingForGroups.clear();
            for ( SectZStrawHitMap::const_iterator secthitmap_it = mhits_it->_sctZTrackerHits.begin(); secthitmap_it != mhits_it->_sctZTrackerHits.end(); ++secthitmap_it ) {
                    nCouplingInGroup=0;
                    for ( ZStrawHitMap::const_iterator zInscthitmap_it =  secthitmap_it->second.begin(); zInscthitmap_it !=  secthitmap_it->second.end(); ++zInscthitmap_it ) {
                            StrawIndex const&firstStrawIdx=((StrawHitPtr) zInscthitmap_it->second)->strawIndex();
                            if ( skipHit.find(firstStrawIdx.asInt())!=skipHit.end() ) continue;
                            const Straw & fstr = ttr.getStraw(firstStrawIdx);
                            if (strdat->count(firstStrawIdx.asInt())==0) strdat->insert(
                                            strawData::value_type(
                                                            firstStrawIdx.asInt(),
                                                            new HitData(zInscthitmap_it->second, fstr.getDirection(), fstr.getMidPoint(), fstr.getHalfLength(), zInscthitmap_it->first, secthitmap_it->first)
                                            )
                            );

                            std::multimap<int, int> tmpCoupling;
                            zInscthitmap_tmp1_it = zInscthitmap_it;
                            ++zInscthitmap_tmp1_it;
                            for ( ; zInscthitmap_tmp1_it !=  secthitmap_it->second.end(); ++zInscthitmap_tmp1_it ) {
                                    StrawIndex const&secondStrawIdx=((StrawHitPtr) zInscthitmap_tmp1_it->second)->strawIndex();
                                    if ( skipHit.find(secondStrawIdx.asInt())!=skipHit.end() ) continue;
                                    const Straw & sstr = ttr.getStraw(secondStrawIdx);
                                    if ( sstr.getMidPoint().getZ() < fstr.getMidPoint().getZ() ) continue;
                                    if (strdat->count(secondStrawIdx.asInt())==0) strdat->insert(
                                                    strawData::value_type(
                                                                    secondStrawIdx.asInt(),
                                                                    new HitData(zInscthitmap_tmp1_it->second, sstr.getDirection(), sstr.getMidPoint(), sstr.getHalfLength(), zInscthitmap_tmp1_it->first, secthitmap_it->first)
                                                    )
                                    );
                                    _hitsCouplings.push_back(std::pair<int,int>(firstStrawIdx.asInt(),secondStrawIdx.asInt()));

                                    zInscthitmap_tmp2_it = zInscthitmap_tmp1_it;
                                    ++zInscthitmap_tmp2_it;

                                    for ( ; zInscthitmap_tmp2_it !=  secthitmap_it->second.end(); ++zInscthitmap_tmp2_it ) {
                                            StrawIndex const&thirdStrawIdx=((StrawHitPtr) zInscthitmap_tmp2_it->second)->strawIndex();
                                            if ( skipHit.find(thirdStrawIdx.asInt())!=skipHit.end() ) continue;
                                            const Straw & tstr = ttr.getStraw(thirdStrawIdx);
                                            if ( tstr.getMidPoint().getZ() < sstr.getMidPoint().getZ() ) continue;
                                            if (strdat->count(thirdStrawIdx.asInt())==0) strdat->insert(
                                                            strawData::value_type(
                                                                            thirdStrawIdx.asInt(),
                                                                            new HitData(zInscthitmap_tmp2_it->second, tstr.getDirection(), tstr.getMidPoint(), tstr.getHalfLength(), zInscthitmap_tmp2_it->first, secthitmap_it->first)
                                                            )
                                            );

                                            tmpCoupling.insert( std::pair<int,int>( secondStrawIdx.asInt(),thirdStrawIdx.asInt() ) );
                                            ++nCouplingInGroup;
                                            ++nCoupling;
                                    }

                            }
                            _hitsTripleCouplings.insert( std::pair< int,std::multimap<int, int> >( firstStrawIdx.asInt(), tmpCoupling) );
                    }
                    if ( nCouplingInGroup>0 ) {
                            cout<<"N Coupling in previos group "<<nCouplingInGroup<<endl;
                            nCouplingForGroups.push_back(nCouplingInGroup);
                            nCouplingForGroups.push_back(nCouplingInGroup);
                    }
            }

            cout<<"Tot N Coupling "<<nCoupling<<endl;
	    
	    if (nCoupling<20) continue;

	    //computeTripleCombination_DeltaRay();
	    //Float_t *hPhi0HelStep_DR_arr = _hPhi0HelStep_DeltaRay->GetArray();
            //nBinsX  = _hPhi0HelStep_DeltaRay->GetNbinsX();
            //nBinsY  = _hPhi0HelStep->GetNbinsY();//_hPhi0HelStep_DeltaRay->GetNbinsY();
            //int nBinsYCopy = 0.5*CLHEP::halfpi/_hPhi0HelStep_DeltaRay->GetYaxis()->GetBinWidth(1);
            //cout<<"------ nBinsYCopy "<<nBinsYCopy<<endl;
            //int newNEntries = _hPhi0HelStep_DeltaRay->GetEntries()+(int)_hPhi0HelStep_DeltaRay->Integral(1,nBinsX,1,nBinsYCopy+1);
            //memcpy(hPhi0HelStep_DR_arr+((nBinsY+1)*(nBinsX+2)),hPhi0HelStep_DR_arr+1*(nBinsX+2),nBinsYCopy*(nBinsX+2)*sizeof(Float_t));
            //_hPhi0HelStep_DeltaRay->SetEntries(newNEntries);

            //computeCombination( mhits_it->_min_HStep, mhits_it->_max_HStep);

            computeTripleCombination( mhits_it->_min_HStep, mhits_it->_max_HStep);
            cout<<"Maximum in hPhi0HelStep = "<<_maxPhi0HelStep<<endl;
            computeHistoProf( _hPhi0HelStep, hPhi0HelS_HMean, hPhi0HelS_HSigma );
            cout<<"hPhi0HelStep Prof val : mean "<<hPhi0HelS_HMean<<" sigma "<<hPhi0HelS_HSigma<<endl;

            Float_t *hPhi0HelStep_arr = _hPhi0HelStep->GetArray();
            Float_t *hPhi0HelStepL_arr = _hPhi0HelStepL->GetArray();
            /*int */nBinsX  = _hPhi0HelStep->GetNbinsX();
            /*int */nBinsY  = _hPhi0HelStep->GetNbinsY();
            int nBinsYL = _hPhi0HelStepL->GetNbinsY();
            memcpy(hPhi0HelStepL_arr+(1*(nBinsX+2)),hPhi0HelStep_arr+(1*(nBinsX+2)),nBinsY*(nBinsX+2)*sizeof(Float_t));
            memcpy(hPhi0HelStepL_arr+((nBinsY+1)*(nBinsX+2)),hPhi0HelStep_arr+1*(nBinsX+2),(nBinsYL-nBinsY)*(nBinsX+2)*sizeof(Float_t));
            _hPhi0HelStepL->SetEntries(_hPhi0HelStep->GetEntries()+(int)_hPhi0HelStep->Integral(1,nBinsX,1,(nBinsYL-nBinsY)));
    //}
            if (_doCalib) {
                    //TCanvas *tmpCanv = new TCanvas();
                    //tmpCanv->cd();
                    //_fakeCanvas->cd();
                    //_hPhi0HelStep->Draw("col z");
                    //tmpCanv->Update();
                    //_fakeCanvas->Update();
                    /*
                    float maxHeight = _hPhi0HelStep->GetMaximum();
                    float tmpBinVal;
                    if (maxHeight<100){
                            for (int jx=1; jx<=_hPhi0HelStep->GetNbinsX(); ++jx) {
                                    for (int jy=1; jy<=_hPhi0HelStep->GetNbinsY(); ++jy) {
                                            tmpBinVal=_hPhi0HelStep->GetBinContent(jx,jy);
                                            if (tmpBinVal>maxHeight) maxHeight=tmpBinVal;
                                    }
                            }
                    }
                    //tmpCanv->WaitPrimitive();
                    //tmpCanv->Delete();
                    cout<<"maximum "<<maxHeight<<endl;
                    */
                    double integralX, integralY;

                    Double_t amin0/*,edm0,errdef0*/, amin0Y;
                    //Int_t nvpar0,nparx0;
                    double maxLK_X=1.0e+10, sigma_X, bestSigma_X=0.0, bestCut=0.0;
                    double maxLK_Y=1.0e+10, sigma_Y, bestSigma_Y=0.0;
                    float minCut;
                    float tmpPhi0, phioPerCut=-0.5*CLHEP::halfpi;
                    float /*tmpRange,*/ tmpBinVal;
                    TF1 fitFX("fitFX","gaus",0,3000);
                    //TF1 fitFX("fitFX","[0]+[1]*((1.0-[4])*TMath::Gaus(x,[2],[3],1)+[4]*TMath::Gaus(x,[5],[6],1))",0,3000);
                    TF1 fitFY("fitFY","[0]*((1.0-[3])*TMath::Gaus(x,[1],[2],1)+[3]*TMath::Gaus(x,[4],[5],1))",-1.6,2.5);
                    //TF1 fitF("fitF","[0]+[1]*(TMath::Gaus(x,[2],[3],1)+TMath::Gaus(x,[4],[5],1))",0,3000);
                    TH1F *bestProf_X=0x0, *bestProf_Y=0x0;

                    fitFX.SetParLimits(1,_hPhi0HelStep_X->GetXaxis()->GetXmin(),_hPhi0HelStep_X->GetXaxis()->GetXmax());
                    fitFX.SetParLimits(2,0,100);

                    fitFY.SetParLimits(1,_hPhi0HelStep_Y->GetXaxis()->GetXmin(),_hPhi0HelStep_Y->GetXaxis()->GetXmax());
                    fitFY.SetParLimits(2,0,100);
                    fitFY.SetParLimits(3,0,1);
                    fitFY.SetParLimits(4,_hPhi0HelStep_Y->GetXaxis()->GetXmin(),_hPhi0HelStep_Y->GetXaxis()->GetXmax());
                    fitFY.SetParLimits(5,0,1000);

                    for (double cut=0.8; cut>=0.2; cut-=0.05){
                            minCut = floor ( cut*_maxPhi0HelStep );

                            cout<<"Cutting at "<<cut<<endl;

                            _hPhi0HelStep_X->Reset();
                            _hPhi0HelStep_Y->Reset();
                            integralX=0.0;
                            integralY=0.0;
                            for (int jx=1; jx<=_hPhi0HelStep->GetNbinsX(); ++jx) {
                                    for (int jy=1; jy<=_hPhi0HelStep->GetNbinsY(); ++jy) {
                                            tmpBinVal=_hPhi0HelStep->GetBinContent(jx,jy);
                                            //_hPhi0HelStep_Height->Fill(tmpBinVal/_maxPhi0HelStep);
                                            tmpBinVal-=minCut;
                                            if ( tmpBinVal>0.0 ) {
                                                    _hPhi0HelStep_X->Fill(_hPhi0HelStep_X->GetBinCenter(jx),tmpBinVal);
                                                    integralX+=tmpBinVal;
                                                    tmpPhi0 = _hPhi0HelStep_Y->GetBinCenter(jy);
                                                    _hPhi0HelStep_Y->Fill(tmpPhi0,tmpBinVal);
                                                    integralY+=tmpBinVal;
                                                    if (tmpPhi0<=phioPerCut) {
                                                            _hPhi0HelStep_Y->Fill(tmpPhi0+CLHEP::pi,tmpBinVal);
                                                            integralY+=tmpBinVal;
                                                    }

                                            }
                                    }
                            }

                            //_hPhi0HelStep_X->Fit("fitFtemp");
                            //fitFX.SetParameter(0,1);
                            //tmpRange=_hPhi0HelStep_X->GetMaximum();
                            //fitFX.SetParLimits(0,0.0,tmpRange);
                            fitFX.SetParameter(1,_hPhi0HelStep_X->GetMean(1));
                            fitFX.SetParameter(2,_hPhi0HelStep_X->GetRMS(1)*0.1);
                            //fitFX.SetParameter(4,0.05);
                            //fitFX.SetParameter(5,_hPhi0HelStep_X->GetMean(1));
                            //fitFX.SetParameter(6,_hPhi0HelStep_X->GetRMS(1));

                            _hPhi0HelStep_X->Fit("fitFX");
                            _hPhi0HelStep_X->Fit("fitFX");
                            //_hPhi0HelStep_X->Fit("fitFX");
                            //_hPhi0HelStep_X->Fit("fitFX","L");
                            //TVirtualFitter *fitter_X = TVirtualFitter::Fitter(_hPhi0HelStep_X);
                            //fitter_X->GetStats(amin0,edm0,errdef0,nvpar0,nparx0);
                            //sigma_X=fitter_X->GetParameter(3);
                            amin0=abs(1.0-fitFX.GetChisquare()/((double)fitFX.GetNDF()));
                            sigma_X=fitFX.GetParameter(2);

                            //fitFY.SetParameter(0,1);
                            //tmpRange=_hPhi0HelStep_Y->GetMaximum();
                            //fitFY.SetParLimits(0,0.0,tmpRange);
                            integralY*=_hPhi0HelStep_Y->GetBinWidth(1);
                            fitFY.SetParLimits(0,0.0,integralY+10.0);
                            fitFY.SetParameter(0,integralY);
                            fitFY.SetParameter(1,_hPhi0HelStep_Y->GetMean(1));
                            fitFY.SetParameter(2,_hPhi0HelStep_Y->GetRMS(1)*0.1);
                            fitFY.SetParameter(3,0.1);
                            fitFY.SetParameter(4,_hPhi0HelStep_Y->GetMean(1));
                            fitFY.SetParameter(5,_hPhi0HelStep_Y->GetRMS(1));
                            //fitFY.FixParameter(4,0);
                            //fitFY.FixParameter(5,0);
                            //fitFY.FixParameter(6,1);
                            _hPhi0HelStep_Y->Fit("fitFY");
                            _hPhi0HelStep_Y->Fit("fitFY");
                            //_hPhi0HelStep_Y->Fit("fitFY");
                            //_hPhi0HelStep_Y->Fit("fitFY","L");
                            //TVirtualFitter *fitter_Y = TVirtualFitter::Fitter(_hPhi0HelStep_Y);
                            //fitter_Y->GetStats(amin0Y,edm0,errdef0,nvpar0,nparx0);
                            //sigma_Y=fitter_Y->GetParameter(3);
                            amin0Y=abs(1.0-fitFY.GetChisquare()/((double)fitFY.GetNDF()));
                            if (fitFY.GetParameter(3)>0.6) {
                                    sigma_Y=fitFY.GetParameter(5);
                            }
                            else if (fitFY.GetParameter(3)<0.4) {
                                    sigma_Y=fitFY.GetParameter(2);
                            }
                            else {
                                    sigma_Y=sqrt( pow((1.0-fitFY.GetParameter(3))*fitFY.GetParameter(2),2)+pow(fitFY.GetParameter(3)*fitFY.GetParameter(5),2) );
                            }

                            cout<<"---------------"<<endl;
                            cout << " |1-Chi2/NDF| value Fitting " << amin0 <<" "<<amin0Y<< endl;
                            cout<<"---------------"<<endl;
                            //amin0 *=-1.0;
                            //amin0Y*=-1.0;
                            if ( amin0<maxLK_X && amin0Y<maxLK_Y ){
                                    maxLK_X=amin0;
                                    maxLK_Y=amin0Y;
                                    bestCut=cut;
                                    bestSigma_X=sigma_X;
                                    bestSigma_Y=sigma_Y;
                                    if (bestProf_X!=0x0) { bestProf_X->Delete(); }
                                    if (bestProf_Y!=0x0) { bestProf_Y->Delete(); }
                                    bestProf_X = (TH1F*) _hPhi0HelStep_X->Clone("bestProf_X");
                                    bestProf_Y = (TH1F*) _hPhi0HelStep_Y->Clone("bestProf_Y");
                            }
                            //delete fitter_X;
                            //delete fitter_Y;
                            cout<<"End cutting at "<<cut<<endl;
                    }
                    cout<<"Best Fitted: Likelihood value "<<maxLK_X<<" "<<maxLK_Y<<" @ cut "<<bestCut<<" with sigma "<<bestSigma_X<<" "<<bestSigma_Y<<endl;
                    _hPhi0HelStep_SigmaX->Fill(bestSigma_X);
                    _hPhi0HelStep_SigmaY->Fill(bestSigma_Y);
                    _hPhi0HelStep_Cut->Fill(bestCut);

                    if (_doDisplay) {
                            _hPhi0HelStep_X->Reset();
                            _hPhi0HelStep_Y->Reset();
                            for (int iPbin=1; iPbin<=bestProf_X->GetNbinsX(); ++iPbin) {
                                    _hPhi0HelStep_X->SetBinContent(iPbin, bestProf_X->GetBinContent(iPbin));
                            }
                            for (int iPbin=1; iPbin<=bestProf_Y->GetNbinsX(); ++iPbin) {
                                    _hPhi0HelStep_Y->SetBinContent(iPbin, bestProf_Y->GetBinContent(iPbin));
                            }
                            //delete bestProf_X;
                            //delete bestProf_Y;
                            bestProf_X->Delete();
                            bestProf_Y->Delete();
                            fitFX.SetParameter(3,bestSigma_X);
                            _hPhi0HelStep_X->Fit("fitFX");//,"L");
                            _hPhi0HelStep_X->Fit("fitFX");//,"L");
                            //_hPhi0HelStep_X->Fit("fitFX","L");
                            fitFY.SetParameter(3,bestSigma_Y);
                            _hPhi0HelStep_Y->Fit("fitFY");//,"L");
                            _hPhi0HelStep_Y->Fit("fitFY");//,"L");
                            //_hPhi0HelStep_Y->Fit("fitFY","L");
                    }


//                    TH2F *hPhi0HelStep_temp = (TH2F *) _hPhi0HelStep->Clone("hPhi0HelStep_temp");
//                    hPhi0HelStep_temp->Rebin2D(2);
//                    _hPhi0HelStep_Height->Reset();
//                    for (int jx=1; jx<=hPhi0HelStep_temp->GetNbinsX(); ++jx) {
//                            for (int jy=1; jy<=hPhi0HelStep_temp->GetNbinsY(); ++jy) {
//                                    tmpBinVal=hPhi0HelStep_temp->GetBinContent(jx,jy);
//                                    _hPhi0HelStep_Height->Fill(tmpBinVal/_maxPhi0HelStep);
//                            }
//                    }
//
//                    invG.FixParameter(0,_hPhi0HelStep_Height->GetBinWidth(1));
//                    invG.FixParameter(2,hPhi0HelStep_temp->GetXaxis()->GetBinWidth(1)*hPhi0HelStep_temp->GetYaxis()->GetBinWidth(1) );
//                    invG.SetParameter(1,1);
//                    invG.SetParLimits(1,1.0e-3,100);
//                    //invG.FixParameter(3,1);
//                    //invG.SetParameter(3,0.5);
//                    //invG.SetParLimits(3,0,0.9);
//
//
////                    invG.SetParameter(4,0.5);
////                    invG.SetParLimits(4,1.0e-3,1000);
////                    invG.SetParameter(5,0.1);
////                    invG.SetParLimits(5,1.0e-4,1000);
//                    //invG.FixParameter(4,0);
//                    //invG.FixParameter(5,1);
//
////                    //invG.SetParameter(6,0.1);
////                    //invG.SetParLimits(6,1.0e-3,1000);
////                    //invG.SetParameter(7,0.0125);
////                    //invG.SetParLimits(7,1.0e-4,1000);
////                    invG.FixParameter(6,0);
////                    invG.FixParameter(7,1);
//
//                    //invG.SetParameter(6,20);
//                    //invG.SetParLimits(6,1.0e-3,1000);
//
//                    Double_t amin0,edm0,errdef0;
//                    Int_t nvpar0,nparx0;
//                    //std::vector< std::pair<double, double> > firRes;
//                    double maxLK=0.0, bestCut, bestSigma;
//                    int maxPos=-1;
//                    for (double cut=0.9; cut>=0.0; cut-=0.05){
//                            invG.FixParameter(3,cut);
//                            _hPhi0HelStep_Height->Fit("invG","L","",cut,1.0);
//                            TVirtualFitter *fitter0 = TVirtualFitter::Fitter(_hPhi0HelStep_Height);
//                            fitter0->GetStats(amin0,edm0,errdef0,nvpar0,nparx0);
//                            //firRes.push_back( std::pair<amin0, fitter0->GetParameter(1)> );
//                            if (amin0>maxLK) {
//                                    maxLK=amin0;
//                                    bestCut=cut;
//                                    bestSigma=fitter0->GetParameter(1);
//                            }
//                            cout<<"---------------"<<endl;
//                            cout << " Likelihood value Fitting " << amin0 << endl;
//                            cout<<"---------------"<<endl;
//                            delete fitter0;
//                    }
//                    cout<<"Best Fitted: Likelihood value "<<maxLK<<" @ cut "<<bestCut<<" with sigma "<<bestSigma<<endl;
//                    invG.FixParameter(3,bestCut);
//                    _hPhi0HelStep_Height->Fit("invG","L","",bestCut,1.0);
//
//
//                    delete hPhi0HelStep_temp;
//                    _hPhi0HelStep_X->Fit("gaus");
//                    _hPhi0HelStep_SigmaX->Fill( _hPhi0HelStep_X->GetFunction("gaus")->GetParameter(2) );
//                    _hPhi0HelStep_Y->Fit("gaus");
//                    _hPhi0HelStep_SigmaY->Fill( _hPhi0HelStep_Y->GetFunction("gaus")->GetParameter(2) );
            }
            else {

                    /*int binMaxX, binMaxY, maxVal;
	            _hPhi0HelStep->GetMaximumBin(binMaxX,binMaxY,maxVal);
	            maxVal=(int)_hPhi0HelStep->GetBinContent(binMaxX,binMaxY);
	            cout<<"Maximun in the Histo:"<<maxVal<<" at bin "<<binMaxX-1<<" "<<binMaxY-1<<endl;
	            measureRbyTripleCombination( _hPhi0HelStep->GetXaxis()->GetBinCenter(binMaxX), _hPhi0HelStep->GetYaxis()->GetBinCenter(binMaxY) );*/

                    //int binx, biny, tmpBin, nBinsXEff;

//                    //------------------- Delta Ray extraction
//                    cout<<"Form Delta Ray extraction"<<endl;
//                    npeaks1 = peakFinder(_hPhi0HelStep_DeltaRay, 6/*peakSigma*/, 50, peakPosX1, peakPosY1 );
//                    nBinsXEff = _hPhi0HelStep_DeltaRay->GetNbinsX()+2;
//                    for(int ip=0; ip<npeaks1; ip++) {
//                            binx   = (int)floor(peakPosX1[ip]+0.5)+1;
//                            biny   = (int)floor(peakPosY1[ip]+0.5)+1;
//                            tmpBin = biny*nBinsXEff + binx;
//
//                            cout<<"Peak at pos "<<binx<<" = "<<_hPhi0HelStep_DeltaRay->GetXaxis()->GetBinCenter(binx)<<" - "<<biny<<" = "<<_hPhi0HelStep_DeltaRay->GetYaxis()->GetBinCenter(biny)<<" straws in peak:"<<endl;
//                            if (biny<=nBinsY) {
//                                    voteArrHitMap::iterator voteArr_Bin_Hit_rel_DR_it = voteArr_Bin_Hit_rel_DR.find(tmpBin);
//                                    cout<<"Hits in peak "<<voteArr_Bin_Hit_rel_DR_it->second.size()<<endl;
//                                    for ( strawList::iterator inpeakStraws_it=voteArr_Bin_Hit_rel_DR_it->second.begin(); inpeakStraws_it!=voteArr_Bin_Hit_rel_DR_it->second.end(); ++inpeakStraws_it) {
//                                            cout<<"\t"<<*inpeakStraws_it<<endl;
//                                    }
//                                    //std::pair<voteArrComMap::iterator,voteArrComMap::iterator> voteArr_Bin_Comb_rel_DR_range = voteArr_Bin_Comb_rel_DR.equal_range(tmpBin);
//                                    //cout<<"Triple combination in peak: "<<endl;
//                                    //for (voteArrComMap::iterator voteArr_Bin_Comb_rel_DR_it=voteArr_Bin_Comb_rel_DR_range.first; voteArr_Bin_Comb_rel_DR_it!=voteArr_Bin_Comb_rel_DR_range.second; ++voteArr_Bin_Comb_rel_DR_it) {
//                                    //        cout<<"\t"<<voteArr_Bin_Comb_rel_DR_it->second.first<<" "<<voteArr_Bin_Comb_rel_DR_it->second.second.first<<" "<<voteArr_Bin_Comb_rel_DR_it->second.second.second<<endl;
//                                    //}
//                           }
//                    }
//
//                    //------------------- End Delta Ray extraction


                    int maxInBin=-1, tmpNHitInBin;
                    std::vector<int> binWithMaxs;
                    //int binx, biny, tmpBin, nBinsXEff;
                    nBinsXEff = nBinsX+2;
                    for (voteArrHitMap::iterator voteArr_Bin_Hit_rel_it = voteArr_Bin_Hit_rel.begin(); voteArr_Bin_Hit_rel_it != voteArr_Bin_Hit_rel.end(); ++voteArr_Bin_Hit_rel_it) {
                            tmpNHitInBin=voteArr_Bin_Hit_rel_it->second.size();
                            //cout<<"bin "<<voteArr_Bin_Hit_rel_it->first<<" nHit "<<tmpNHitInBin<<endl;

                            _BFRot_hPhi0HelStepL->SetBinContent(voteArr_Bin_Hit_rel_it->first,tmpNHitInBin);

                            if (tmpNHitInBin>maxInBin) {
                                    binWithMaxs.clear();
                                    maxInBin=tmpNHitInBin;
                                    binWithMaxs.push_back(voteArr_Bin_Hit_rel_it->first);
                            }
                            else if (tmpNHitInBin==maxInBin) {
                                    binWithMaxs.push_back(voteArr_Bin_Hit_rel_it->first);
                            }
                    }

                    cout<<"Max Hit in a bin "<<maxInBin<<endl;
                    for (std::vector<int>::iterator binWithMaxs_it=binWithMaxs.begin(); binWithMaxs_it!=binWithMaxs.end(); ++binWithMaxs_it) {
                            tmpBin=*binWithMaxs_it;
                            biny = tmpBin/nBinsXEff;
                            binx = tmpBin-nBinsXEff*biny;
                            cout<<"\t Bin with max "<<tmpBin<<" step= "<< _hPhi0HelStepL->GetXaxis()->GetBinCenter(binx) <<" phi0= "<< _hPhi0HelStepL->GetYaxis()->GetBinCenter(biny) <<endl;
                    }

                    //npeaks = peakFinder(_hPhi0HelStepL, peakSigmaBins/*peakSigma*/, 50, peakPosX, peakPosY, 5 );

                    //int binxMax, binyMax, tmpBinMax;
                    //for(int ip=0; ip<npeaks; ip++) {
                    //        binx   = (int)floor(peakPosX[ip]+0.5)+1;
                    //        biny   = (int)floor(peakPosY[ip]+0.5)+1;
                    //        tmpBin = biny*nBinsXEff + binx;

                    //        cout<<"Peak at pos "<<binx<<" = "<<_hPhi0HelStepL->GetXaxis()->GetBinCenter(binx)<<" - "<<biny<<" = "<<_hPhi0HelStepL->GetYaxis()->GetBinCenter(biny)<<" straws in peak:"<<endl;
                    //        if (biny<=nBinsY) {
                    //                voteArrHitMap::iterator voteArr_Bin_Hit_rel_it = voteArr_Bin_Hit_rel.find(tmpBin);
                    //                cout<<"Hits in peak "<<voteArr_Bin_Hit_rel_it->second.size()<<endl;
                    //                for ( strawList::iterator inpeakStraws_it=voteArr_Bin_Hit_rel_it->second.begin(); inpeakStraws_it!=voteArr_Bin_Hit_rel_it->second.end(); ++inpeakStraws_it) {
                    //                        cout<<"\t"<<*inpeakStraws_it<<endl;
                    //                }
                    //        }
                    //}

                    Float_t *BF_hPhi0HelStepL_arr = _BF_hPhi0HelStepL->GetArray();
                    //Float_t *BFRot_hPhi0HelStepL_arr = _BFRot_hPhi0HelStepL->GetArray();
                    //computeHistoProf( _BFRot_hPhi0HelStepL, hPhi0HelS_HMean, hPhi0HelS_HSigma );
                    //cout<<"hPhi0HelStep Prof val : mean "<<hPhi0HelS_HMean<<" sigma "<<hPhi0HelS_HSigma<<endl;
                    //_BF_hPhi0HelStepL->SetEntries( smooth( BFRot_hPhi0HelStepL_arr, BF_hPhi0HelStepL_arr, nBinsX, nBinsYL, 0.0, 0.0, 1.0, 1.0, hPhi0HelS_HMean+4.0*hPhi0HelS_HSigma ) );
                    _BF_hPhi0HelStepL->SetEntries( smooth( hPhi0HelStepL_arr, BF_hPhi0HelStepL_arr, nBinsX, nBinsYL, 0.0, 0.0, 1.0, 1.0, hPhi0HelS_HMean+5.0*hPhi0HelS_HSigma ) );

                    /*
                    float maxHeight = _hPhi0HelStep->GetMaximum();
                    float tmpBinVal;
                    if (maxHeight<200){
                            for (int jx=1; jx<=_hPhi0HelStep->GetNbinsX(); ++jx) {
                                    for (int jy=1; jy<=_hPhi0HelStep->GetNbinsY(); ++jy) {
                                            tmpBinVal=_hPhi0HelStep->GetBinContent(jx,jy);
                                            if (tmpBinVal>maxHeight) maxHeight=tmpBinVal;
                                    }
                            }
                    }
                    cout<<"Maximum in hPhi0HelStep = "<<maxHeight<<endl;
                    */

                    //_BF_hPhi0HelStepL->SetEntries( smooth( hPhi0HelStepL_arr, BF_hPhi0HelStepL_arr, nBinsX, nBinsYL, peakSigmaX, peakSigmaY, binSizeX, binSizeY, 0.5*_maxPhi0HelStep ) );
                    //_hPhi0HelStepL->Scale(1.0/_hPhi0HelStepL->GetEntries());
                    //_BF_hPhi0HelStepL->Scale(1.0/_hPhi0HelStepL->GetEntries());
                    //_BF_hPhi0HelStepL->SetEntries( blur( hPhi0HelStepL_arr, BF_hPhi0HelStepL_arr, nBinsX, nBinsYL, peakSigmaX, peakSigmaY, binSizeX, binSizeY ) );
//                    _BF_hPhi0HelStepL->SetEntries( butterflyFilterRot45( hPhi0HelStepL_arr, BF_hPhi0HelStepL_arr, nBinsX, nBinsYL, 0.0/*30.0*/ ) );
                    //_BFRot_hPhi0HelStepL->SetEntries( butterflyFilterRot( hPhi0HelStepL_arr, BFRot_hPhi0HelStepL_arr, nBinsX, nBinsYL, 30.0 ) );


                    //cout<<"Form Hit second wiev"<<endl;
                    //npeaks1 = peakFinder(_BF_hPhi0HelStepL, peakSigmaBins/*peakSigma*/, 5, peakPosX1, peakPosY1, 5 );
                    //for(int ip=0; ip<npeaks1; ip++) {
                    //        binx   = (int)floor(peakPosX1[ip]+0.5)+1;
                    //        biny   = (int)floor(peakPosY1[ip]+0.5)+1;
                    //        tmpBin = biny*nBinsXEff + binx;

                    //        cout<<"Peak at pos "<<binx<<" = "<<_BF_hPhi0HelStepL->GetXaxis()->GetBinCenter(binx)<<" - "<<biny<<" = "<<_BF_hPhi0HelStepL->GetYaxis()->GetBinCenter(biny)<<" straws in peak:"<<endl;
                    //        if (biny<=nBinsY) {
                    //                voteArrHitMap::iterator voteArr_Bin_Hit_rel_it = voteArr_Bin_Hit_rel.find(tmpBin);
                    //                cout<<"Hits in peak "<<voteArr_Bin_Hit_rel_it->second.size()<<endl;
                    //                for ( strawList::iterator inpeakStraws_it=voteArr_Bin_Hit_rel_it->second.begin(); inpeakStraws_it!=voteArr_Bin_Hit_rel_it->second.end(); ++inpeakStraws_it) {
                    //                        cout<<"\t"<<*inpeakStraws_it<<endl;
                    //                }
                    //        }
                    //}



                    //            _BFRot_hPhi0HelStepL->SetEntries( blur( BF_hPhi0HelStepL_arr, BFRot_hPhi0HelStepL_arr, nBinsX, nBinsYL, peakSigmaX, peakSigmaY, binSizeX, binSizeY ) );
                    //            cout<<"Form Hit multiplicity wiev"<<endl;
                    //            npeaks1 = peakFinder(_BFRot_hPhi0HelStepL, peakSigma, 70, peakPosX1, peakPosY1 );
                    //            for(int ip=0; ip<npeaks1; ip++) {
                    //                    binx   = (int)floor(peakPosX1[ip]+0.5)+1;
                    //                    biny   = (int)floor(peakPosY1[ip]+0.5)+1;
                    //                    tmpBin = biny*nBinsXEff + binx;
                    //
                    //                    cout<<"Peak at pos "<<binx<<" = "<<_BFRot_hPhi0HelStepL->GetXaxis()->GetBinCenter(binx)<<" - "<<biny<<" = "<<_BFRot_hPhi0HelStepL->GetYaxis()->GetBinCenter(biny)<<" straws in peak:"<<endl;
                    //                    if (biny<=nBinsY) {
                    //                            voteArrHitMap::iterator voteArr_Bin_Hit_rel_it = voteArr_Bin_Hit_rel.find(tmpBin);
                    //                            cout<<"Hits in peak "<<voteArr_Bin_Hit_rel_it->second.size()<<endl;
                    //                            for ( strawList::iterator inpeakStraws_it=voteArr_Bin_Hit_rel_it->second.begin(); inpeakStraws_it!=voteArr_Bin_Hit_rel_it->second.end(); ++inpeakStraws_it) {
                    //                                    cout<<"\t"<<*inpeakStraws_it<<endl;
                    //                            }
                    //                    }
                    //            }




                    int binMaxX, binMaxY, maxVal;
                    _BF_hPhi0HelStepL->GetMaximumBin(binMaxX,binMaxY,maxVal);
                    maxVal=(int)_BF_hPhi0HelStepL->GetBinContent(binMaxX,binMaxY);
                    cout<<"Maximun in the Histo:"<<maxVal<<" at bin "<<binMaxX-1<<" "<<binMaxY-1<<endl;
                    measureRbyTripleCombination( _BF_hPhi0HelStepL->GetXaxis()->GetBinCenter(binMaxX), _BF_hPhi0HelStepL->GetYaxis()->GetBinCenter(binMaxY) );

                    //TSpectrum2 peaksSearcher(10);
                    //int nfound =  peaksSearcher.Search( _hPhi0HelStepL, 2.0, "", 0.1 );
                    //TSpectrum2 peaksSearcherBF(10);
                    //int nfoundBF =  peaksSearcherBF.Search( _BF_hPhi0HelStepL, 2.0, "", 0.1 );

                    //cout<<"Peaks found in Voting Array "<<nfound<<" in BFiltered "<<nfoundBF<<endl;

                    //TH1D *_BF_hPhi0HelStepL_Y = _BF_hPhi0HelStepL->ProjectionY();
            }


//            _hRHelStep->Reset();

//            nCoupling=0;
//            nCouplingInGroup=0;
//            lastGroupStation=0;
//            nCouplingForGroups.clear();
//            for ( ZSectStrawHitMap::const_iterator zhitmap_it = mhits_it->_zsctTrackerHits.begin(); zhitmap_it != mhits_it->_zsctTrackerHits.end(); ++zhitmap_it ) {
//                    firstStation = (unsigned int)zhitmap_it->first/4;
//                    //cout<<"******** firstStation "<<firstStation<<" lastGroupStation "<<lastGroupStation<<endl;
//                    if ( nCoupling>0 && ((int)(firstStation-lastGroupStation))>0 ) {
//                            cout<<"N Coupling in previos group "<<nCouplingInGroup<<endl;
//                            nCouplingForGroups.push_back(nCouplingInGroup);
//                            nCouplingInGroup=0;
//                    }
//                    //goodCoupling=false;
//                    for ( AbsSectStrawHitMap::const_iterator scthitmap_it =  zhitmap_it->second.begin(); scthitmap_it !=  zhitmap_it->second.end(); ++scthitmap_it ) {
//                            //cout<<"Starting hit at z "<<zhitmap_it->first<<" sec "<<scthitmap_it->first<<endl;
//                            lastStation=firstStation;
//                            //if (lastStation>lastGroupStation) lastGroupStation=lastStation;
//                            zhitmap_tmp_it  = zhitmap_it;
//                            scthitmap_tmp_it = scthitmap_it;
//                            ++scthitmap_tmp_it;
//                            goodCoupling=false;
//                            //cout<<"Coupling Z: "<<zhitmap_it->first<<" - "<<zhitmap_tmp_it->first<<":"<<endl;
//                            for ( ; scthitmap_tmp_it != zhitmap_tmp_it->second.end(); ++scthitmap_tmp_it) {
//                                    if (scthitmap_tmp_it->first != scthitmap_it->first) break;
//                                    //tmpSecDist = scthitmap_tmp_it->first - scthitmap_it->first;
//                                    //if ( tmpSecDist<snegativeSectOver || tmpSecDist>=maxContSect ) break;
//                                    //cout<<"\t Sect: "<<scthitmap_it->first<<" - "<<scthitmap_tmp_it->first<<endl;
//                                    goodCoupling=true;
//                                    nCoupling++;
//                                    nCouplingInGroup++;
//
//                                    StrawIndex const&firstStrawIdx=((StrawHitPtr) scthitmap_it->second)->strawIndex();
//                                    StrawIndex const&secondStrawIdx=((StrawHitPtr) scthitmap_tmp_it->second)->strawIndex();
//                                    _hitsCouplings.push_back(std::pair<int,int>(firstStrawIdx.asInt(),secondStrawIdx.asInt()));
//                                    const Straw & fstr = ttr.getStraw(firstStrawIdx);
//                                    if (strdat->count(firstStrawIdx.asInt())==0) strdat->insert(
//                                                    strawData::value_type(
//                                                                    firstStrawIdx.asInt(),
//                                                                    new HitData(scthitmap_it->second, fstr.getDirection(), fstr.getMidPoint(), fstr.getHalfLength(), zhitmap_it->first, scthitmap_it->first)
//                                                    )
//                                    );
//                                    const Straw & sstr = ttr.getStraw(secondStrawIdx);
//                                    if (strdat->count(secondStrawIdx.asInt())==0) strdat->insert(
//                                                    strawData::value_type(
//                                                                    secondStrawIdx.asInt(),
//                                                                    new HitData(scthitmap_it->second, sstr.getDirection(), sstr.getMidPoint(), sstr.getHalfLength(), zhitmap_tmp_it->first, scthitmap_tmp_it->first)
//                                                    )
//                                    );
//
//                            }
//                            ++zhitmap_tmp_it;
//
//                            for ( ; zhitmap_tmp_it != mhits_it->_zsctTrackerHits.end(); ++zhitmap_tmp_it ) {
//                                    tmpStation = (unsigned int)zhitmap_tmp_it->first/4;
//                                    if ( ((int)(tmpStation - lastStation))>1 || ((int)(tmpStation - firstStation))>maxStationDist ) break;
//                                    goodCoupling=false;
//                                    //cout<<"Coupling Z: "<<zhitmap_it->first<<" - "<<zhitmap_tmp_it->first<<":"<<endl;
//                                    //cout<<"("<<"tmpStation "<<tmpStation<<" lastStation "<<lastStation<<" firstStation "<<firstStation<<")"<<endl;
//
//                                    if (scthitmap_it->first<negativeSectOver) {
//                                            scthitmap_tmp_it = zhitmap_tmp_it->second.lower_bound( 0 );
//                                            for ( ; scthitmap_tmp_it != zhitmap_tmp_it->second.end(); ++scthitmap_tmp_it) {
//                                                    tmpSecDist = scthitmap_tmp_it->first - scthitmap_it->first;
//                                                    if ( tmpSecDist<snegativeSectOver || tmpSecDist>=maxContSect ) break;
//                                                    //cout<<"\t Sect: "<<scthitmap_it->first<<" - "<<scthitmap_tmp_it->first<<endl;
//                                                    goodCoupling=true;
//                                                    nCoupling++;
//                                                    nCouplingInGroup++;
//
//                                                    StrawIndex const&firstStrawIdx=((StrawHitPtr) scthitmap_it->second)->strawIndex();
//                                                    StrawIndex const&secondStrawIdx=((StrawHitPtr) scthitmap_tmp_it->second)->strawIndex();
//                                                    _hitsCouplings.push_back(std::pair<int,int>(firstStrawIdx.asInt(),secondStrawIdx.asInt()));
//                                                    const Straw & fstr = ttr.getStraw(firstStrawIdx);
//                                                    if (strdat->count(firstStrawIdx.asInt())==0) strdat->insert(
//                                                                    strawData::value_type(
//                                                                                    firstStrawIdx.asInt(),
//                                                                                    new HitData(scthitmap_it->second, fstr.getDirection(), fstr.getMidPoint(), fstr.getHalfLength(), zhitmap_it->first, scthitmap_it->first)
//                                                                    )
//                                                    );
//                                                    const Straw & sstr = ttr.getStraw(secondStrawIdx);
//                                                    if (strdat->count(secondStrawIdx.asInt())==0) strdat->insert(
//                                                                    strawData::value_type(
//                                                                                    secondStrawIdx.asInt(),
//                                                                                    new HitData(scthitmap_it->second, sstr.getDirection(), sstr.getMidPoint(), sstr.getHalfLength(), zhitmap_tmp_it->first, scthitmap_tmp_it->first)
//                                                                    )
//                                                    );
//
//                                           }
//                                            scthitmap_tmp_it = zhitmap_tmp_it->second.lower_bound( 12-negativeSectOver );
//                                            for ( ; scthitmap_tmp_it != zhitmap_tmp_it->second.end(); ++scthitmap_tmp_it) {
//                                                    tmpSecDist = scthitmap_it->first /*+1*/-(scthitmap_tmp_it->first-12);
//                                                    if ( tmpSecDist<snegativeSectOver || tmpSecDist>=maxContSect ) break;
//                                                    //cout<<"\t Sect: "<<scthitmap_it->first<<" - "<<scthitmap_tmp_it->first<<endl;
//                                                    goodCoupling=true;
//                                                    nCoupling++;
//                                                    nCouplingInGroup++;
//
//                                                    StrawIndex const&firstStrawIdx=((StrawHitPtr) scthitmap_it->second)->strawIndex();
//                                                    StrawIndex const&secondStrawIdx=((StrawHitPtr) scthitmap_tmp_it->second)->strawIndex();
//                                                    _hitsCouplings.push_back(std::pair<int,int>(firstStrawIdx.asInt(),secondStrawIdx.asInt()));
//                                                    const Straw & fstr = ttr.getStraw(firstStrawIdx);
//                                                    if (strdat->count(firstStrawIdx.asInt())==0) strdat->insert(
//                                                                    strawData::value_type(
//                                                                                    firstStrawIdx.asInt(),
//                                                                                    new HitData(scthitmap_it->second, fstr.getDirection(), fstr.getMidPoint(), fstr.getHalfLength(), zhitmap_it->first, scthitmap_it->first)
//                                                                    )
//                                                    );
//                                                    const Straw & sstr = ttr.getStraw(secondStrawIdx);
//                                                    if (strdat->count(secondStrawIdx.asInt())==0) strdat->insert(
//                                                                    strawData::value_type(
//                                                                                    secondStrawIdx.asInt(),
//                                                                                    new HitData(scthitmap_it->second, sstr.getDirection(), sstr.getMidPoint(), sstr.getHalfLength(), zhitmap_tmp_it->first, scthitmap_tmp_it->first)
//                                                                    )
//                                                    );
//
//                                           }
//                                    }
//                                    else {
//                                            scthitmap_tmp_it = zhitmap_tmp_it->second.lower_bound( scthitmap_it->first-negativeSectOver );
//                                            for ( ; scthitmap_tmp_it != zhitmap_tmp_it->second.end(); ++scthitmap_tmp_it) {
//                                                    tmpSecDist = scthitmap_tmp_it->first - scthitmap_it->first;
//                                                    if ( tmpSecDist<snegativeSectOver || tmpSecDist>=maxContSect ) break;
//                                                    //cout<<"\t Sect: "<<scthitmap_it->first<<" - "<<scthitmap_tmp_it->first<<endl;
//                                                    goodCoupling=true;
//                                                    nCoupling++;
//                                                    nCouplingInGroup++;
//
//                                                    StrawIndex const&firstStrawIdx=((StrawHitPtr) scthitmap_it->second)->strawIndex();
//                                                    StrawIndex const&secondStrawIdx=((StrawHitPtr) scthitmap_tmp_it->second)->strawIndex();
//                                                    _hitsCouplings.push_back(std::pair<int,int>(firstStrawIdx.asInt(),secondStrawIdx.asInt()));
//                                                    const Straw & fstr = ttr.getStraw(firstStrawIdx);
//                                                    if (strdat->count(firstStrawIdx.asInt())==0) strdat->insert(
//                                                                    strawData::value_type(
//                                                                                    firstStrawIdx.asInt(),
//                                                                                    new HitData(scthitmap_it->second, fstr.getDirection(), fstr.getMidPoint(), fstr.getHalfLength(), zhitmap_it->first, scthitmap_it->first)
//                                                                    )
//                                                    );
//                                                    const Straw & sstr = ttr.getStraw(secondStrawIdx);
//                                                    if (strdat->count(secondStrawIdx.asInt())==0) strdat->insert(
//                                                                    strawData::value_type(
//                                                                                    secondStrawIdx.asInt(),
//                                                                                    new HitData(scthitmap_it->second, sstr.getDirection(), sstr.getMidPoint(), sstr.getHalfLength(), zhitmap_tmp_it->first, scthitmap_tmp_it->first)
//                                                                    )
//                                                    );
//
//                                            }
//                                    }
//                                    if ( scthitmap_it->first>(12-maxStationDist) ) {
//                                            scthitmap_tmp_it = zhitmap_tmp_it->second.lower_bound( 0 );
//                                            for ( ; scthitmap_tmp_it != zhitmap_tmp_it->second.end() && (scthitmap_tmp_it->first - (scthitmap_it->first-12))<maxContSect ; ++scthitmap_tmp_it) {
//                                                    //cout<<"\t Sect: "<<scthitmap_it->first<<" - "<<scthitmap_tmp_it->first<<endl;
//                                                    goodCoupling=true;
//                                                    nCoupling++;
//                                                    nCouplingInGroup++;
//
//                                                    StrawIndex const&firstStrawIdx=((StrawHitPtr) scthitmap_it->second)->strawIndex();
//                                                    StrawIndex const&secondStrawIdx=((StrawHitPtr) scthitmap_tmp_it->second)->strawIndex();
//                                                    _hitsCouplings.push_back(std::pair<int,int>(firstStrawIdx.asInt(),secondStrawIdx.asInt()));
//                                                    const Straw & fstr = ttr.getStraw(firstStrawIdx);
//                                                    if (strdat->count(firstStrawIdx.asInt())==0) strdat->insert(
//                                                                    strawData::value_type(
//                                                                                    firstStrawIdx.asInt(),
//                                                                                    new HitData(scthitmap_it->second, fstr.getDirection(), fstr.getMidPoint(), fstr.getHalfLength(), zhitmap_it->first, scthitmap_it->first)
//                                                                    )
//                                                    );
//                                                    const Straw & sstr = ttr.getStraw(secondStrawIdx);
//                                                    if (strdat->count(secondStrawIdx.asInt())==0) strdat->insert(
//                                                                    strawData::value_type(
//                                                                                    secondStrawIdx.asInt(),
//                                                                                    new HitData(scthitmap_it->second, sstr.getDirection(), sstr.getMidPoint(), sstr.getHalfLength(), zhitmap_tmp_it->first, scthitmap_tmp_it->first)
//                                                                    )
//                                                    );
//
//                                            }
//                                   }
//                                   if (goodCoupling) {
//                                           lastStation=tmpStation;
//                                           if (lastStation>lastGroupStation) lastGroupStation=lastStation;
//                                           //cout<<"--------- lastStation "<<lastStation<<" lastGroupStation "<<lastGroupStation<<endl;
//                                   }
//                            }
//                    }
//                    //if (goodCoupling) lastStation=firstStation;
//            }
//            if ( nCouplingInGroup>0 ) {
//                    cout<<"N Coupling in previos group "<<nCouplingInGroup<<endl;
//                    nCouplingForGroups.push_back(nCouplingInGroup);
//            }
//
//            cout<<endl<<"N of total hit coupling for hits gorup: "<<nCoupling<<" = to n poit "<<(1+sqrt(1+8*nCoupling))/2<<endl<<endl;
//            checkNCoupling=0;
//            for (std::vector<unsigned int>::iterator nCouplingForGroups_it=nCouplingForGroups.begin(); nCouplingForGroups_it!=nCouplingForGroups.end(); ++nCouplingForGroups_it){
//                    checkNCoupling+=* nCouplingForGroups_it;
//            }
//            if (nCouplingForGroups.size()>1 && checkNCoupling==nCoupling) cout<<"Good separtion in groups!!!!"<<endl;
//
//            //for ( strawData::iterator  strdat_it=strdat->begin(); strdat_it!=strdat->end(); ++strdat_it ) {
//            //        cout<<strdat_it->first <<" absZ "<<strdat_it->second->_absZId<<" absSect "<<strdat_it->second->_absSectId<<" dir: "<<strdat_it->second->_direct<<" mid "<<strdat_it->second->_midp<<endl;
//            //        cout<<"wire points (x,y):"<<endl;
//            //        for ( vector< pair<double, double> >::iterator wxy_it = strdat_it->second->_wireXYTable.begin(); wxy_it != strdat_it->second->_wireXYTable.end(); ++wxy_it ) {
//            //                cout<<"\t"<<wxy_it->first<<" , "<<wxy_it->second<<endl;
//            //        }
//            //}
//
//            computeCombination( mhits_it->_min_HStep, mhits_it->_max_HStep);
//            int binMaxX, binMaxY, maxVal;
//            _hRHelStep->GetMaximumBin(binMaxX,binMaxY,maxVal);
//            maxVal=(int)_hRHelStep->GetBinContent(binMaxX,binMaxY);
//            cout<<"Maximun in the Histo:"<<maxVal<<" at bin "<<binMaxX-1<<" "<<binMaxY-1<<endl;

            if (_doDisplay) {

                    if (mhits_it->_min_HStep>0.0 && mhits_it->_max_HStep>0.0) {
                            _hRHelStep->GetXaxis()->SetRangeUser( mhits_it->_min_HStep, mhits_it->_max_HStep);
                            _hPhi0HelStep_1->GetXaxis()->SetRangeUser( mhits_it->_min_HStep, mhits_it->_max_HStep);

                            _hPhi0HelStep->GetXaxis()->SetRangeUser( mhits_it->_min_HStep, mhits_it->_max_HStep);
                            _hPhi0HelStepL->GetXaxis()->SetRangeUser( mhits_it->_min_HStep, mhits_it->_max_HStep);
                            _BF_hPhi0HelStepL->GetXaxis()->SetRangeUser( mhits_it->_min_HStep, mhits_it->_max_HStep);
                            _BFRot_hPhi0HelStepL->GetXaxis()->SetRangeUser( mhits_it->_min_HStep, mhits_it->_max_HStep);
                    }

                    //_plotCanvas->cd();
                    //_hRHelStep->Draw("col z");
                    _plotCanvas->cd(1);
                    //_hRHelStep->Draw("col z");
                    //_hPhi0HelStep->Draw("col z");
                    _hPhi0HelStepL->Draw("col z");
                    _plotCanvas->cd(2);
                    //_hPhi0HelStepL->Draw("col z");
                    _BF_hPhi0HelStepL->Draw("col z");
                    _plotCanvas->cd(3);
                    _hR_1->Draw();
                    //_hPhi0->Draw();
                    //_BF_hPhi0HelStepL_Y->Draw();
                    //_BF_hPhi0HelStepL->Draw("col z");
                    _plotCanvas->cd(4);
                    //_BF_hPhi0HelStepL->Draw("col z");
                    _BFRot_hPhi0HelStepL->Draw("col z");
                    _plotCanvas->Update();

                    _plotCanvas_1->cd(1);
                    _hRhoPhi0_DeltaRay->Draw("col z");
                    _plotCanvas_1->cd(2);
                    _hRhoPhi0_DeltaRay_Blr->Draw("col z");
                    //_hPhi0HelStep_DeltaRay->Draw("col z");
                    _plotCanvas_1->cd(3);
                    _hAbsSctrAbsZ_DeltaRay_Cum->Draw("col z");
                    _plotCanvas_1->cd(4);
                    _hAbsSctrRadi_DeltaRay_Cum->Draw("col z");
                    _plotCanvas_1->Update();

                    if (_doCalib) {
                            _plotCanvas_Cal->cd(1);
                            _hPhi0HelStep_Height->Draw();
                            _plotCanvas_Cal->cd(2);
                            _hPhi0HelStep_X->Draw();
                            //_hPhi0HelStepL_X->Draw();
                            //_hRHelStep->Draw("col z");
                            _plotCanvas_Cal->cd(3);
                            _hPhi0HelStep_Y->Draw();
                            //_hPhi0HelStepL_Y->Draw();
                            //_hPhi0HelStep_1->Draw("col z");
                            _plotCanvas_Cal->cd(4);
                            _hPhi0HelStep_Cut->Draw();
                            _plotCanvas_Cal->cd(5);
                            _hPhi0HelStep_SigmaX->Draw();
                            //_hR->Draw();
                            _plotCanvas_Cal->cd(6);
                            _hPhi0HelStep_SigmaY->Draw();
                            //_hPhi0->Draw();
                            _plotCanvas_Cal->Update();
                    }

                    cerr << "Double click in the canvas_Fake to continue:" ;
                    _fakeCanvas->Clear();
                    _fakeCanvas->cd();
                    TLatex *printEvN = new TLatex(0.15,0.4,Form("Current Event: %d",event.id().event()));
                    printEvN->SetTextFont(62);
                    printEvN->SetTextSizePixels(180);
                    printEvN->Draw();
                    _fakeCanvas->Update();
                    _fakeCanvas->WaitPrimitive();
                    cerr << endl;
                    delete printEvN;
            }


    }

    for ( strawData::iterator strdat_it=strdat->begin(); strdat_it!=strdat->end(); ++strdat_it ) {
            delete strdat_it->second;
    }
    strdat->clear();
    delete strdat;


    cout<<"--------------------------- End of Reco Track ------------------"<<endl;
    cout<<"Event: "<<event.id().event()<<endl;
    cout<<"----------------------------------------------------------------"<<endl;


  } // end analyze

  void TrackReco::endJob(){

    if (_doCalib){
            _hPhi0HelStep_SigmaX->Fit("gaus");
            _hPhi0HelStep_SigmaY->Fit("gaus");
            _hPhi0HelStep_Cut->Fit("gaus");
    }

    // cd() to correct root directory. See note 3.
    TDirectory* save = gDirectory;
    _directory->cd();

    // Write canvas.  See note 4.
//    _canvas->Write();

    // cd() back to where we were.  See note 3.
    save->cd();

  }

  void TrackReco::computeCombination( double minHStep, double maxHStep ){
          cout<<"In TrackReco::computeCombination"<<endl;
          cout<<"-------- minHStep "<<minHStep<<" maxHStep "<<maxHStep<<endl;
          int firstXbin=0, lastXBin=_nBinHelpStep;
          if ( minHStep>_minHelStep) {
                  for (int ibin=0; ibin<_nBinHelpStep; ++ibin){
                          if (minHStep<=_iXHelStep[ibin]) {
                                  firstXbin=ibin;
                                  break;
                          }
                  }
          }
          if ( maxHStep>0.0 && maxHStep<_maxHelStep ) {
                  lastXBin=firstXbin;
                  for (int ibin=firstXbin; ibin<_nBinHelpStep; ++ibin){
                          if (maxHStep<=_iXHelStep[ibin]) {
                                  break;
                          }
                          ++lastXBin;
                  }
          }
          cout<<"-------- firstXbin "<<firstXbin<<" lastXBin "<<lastXBin<<endl;
          //firstXbin=0;


          //double tmpR1, tmpDen1, tmpNum1, tmpR2, tmpDen2, tmpNum2;
          //double tmpXo, tmpYo, tmpXC, tmpYC, tmpCdist, tmpCspread;

          double DeltaZ1, DeltaZ2/*, DelataX12*/;
          double Phi1, Phi2, Phi0_1, Phi0_2/*, Phi0_1_1, Phi0_2_1*/;
          //double cosPhi1, sinPhi1, cosPhi2, sinPhi2;
          double cosTheta, sinTheta/*, cosPhi1lTheta, sinPhi1lTheta, cosPhi2lTheta, sinPhi2lTheta*/;
          //double A, B, C, parA[5], parB[5], parC[5];
          //double sign1, sign2;
          //bool   goodPhi0;

          double threeOverTwoPi = 3.0*CLHEP::halfpi;

          double tmpVal/*, tmpDelta*/, DiffAngle1, DiffAngle2, tmpPhi0;
          //double tmpNBin, stepXY, maxXY = _maxR-targetRadMax;
          //tmpNBin = 2.0*maxXY/_stepR;
          //int nXYStep = (int)floor(tmpNBin+0.5);
          //stepXY=2.0*maxXY/((double)nXYStep);
          //double tmpTwopiTanLambda1, tmpTwopiTanLambda2, twopiTanLambdaMin = CLHEP::twopi*tan(_minLambda), twoPitanLambdaMax = CLHEP::twopi*tan(_maxLambda);

          for ( vector< std::pair<int, int> >::iterator hitsCouplings_it = _hitsCouplings.begin(); hitsCouplings_it != _hitsCouplings.end(); ++hitsCouplings_it ) {
                  cout<<"I'm analyzing the couple: "<<hitsCouplings_it->first<<" "<<hitsCouplings_it->second<<endl;

                  HitData *first  = strdat->find(hitsCouplings_it->first)->second;
                  HitData *second = strdat->find(hitsCouplings_it->second)->second;
                  DeltaZ1   =  first->_midp.getZ() - zo;
                  DeltaZ2   = second->_midp.getZ() - zo;
                  cosTheta  =  first->_direct.getY();
                  sinTheta  = -first->_direct.getX();
                  //DelataX12 = second->_radius-first->_radius;
                  //cout<<"AbsSect "<<first->_absSectId<<" mid "<<first->_midp<<" theta "<<first->_theta<<endl;
                  for ( int ibin=firstXbin; ibin<lastXBin; ++ibin ) {
                          tmpVal        = CLHEP::twopi/_iXHelStep[ibin];

                          Phi1          = tmpVal*DeltaZ1;
                          Phi1          = removeAnglePeriod(Phi1);
                          //cosPhi1       = cos(Phi1);
                          //sinPhi1       = sin(Phi1);
                          //cosPhi1lTheta = cosPhi1*cosTheta+sinPhi1*sinTheta;
                          //sinPhi1lTheta = sinPhi1*cosTheta-cosPhi1*sinTheta;
                          Phi2          = tmpVal*DeltaZ2;
                          Phi2          = removeAnglePeriod(Phi2);
                          //cosPhi2       = cos(Phi2);
                          //sinPhi2       = sin(Phi2);
                          //cosPhi2lTheta = cosPhi2*cosTheta+sinPhi2*sinTheta;
                          //sinPhi2lTheta = sinPhi2*cosTheta-cosPhi2*sinTheta;

                          tmpPhi0 = atan( (sinTheta-cosTheta*cosTheta)/(sinTheta*cosTheta));
                          tmpPhi0 = removeAnglePeriod(tmpPhi0);
                          Phi0_1  = tmpPhi0-CLHEP::halfpi;
                          DiffAngle1 = Phi0_1-first->_theta;
                          DiffAngle1 = removeAnglePeriod( DiffAngle1 );
                          if (DiffAngle1>threeOverTwoPi) {
                                  Phi0_1 -= Phi1;
                                  Phi0_1 = removeAnglePeriod(Phi0_1);
                                  _hPhi0HelStep_1->Fill(_iXHelStep[ibin],Phi0_1);
                          }
                          tmpPhi0+= CLHEP::pi;
                          tmpPhi0 = removeAnglePeriod(tmpPhi0);
                          Phi0_1  = tmpPhi0-CLHEP::halfpi;
                          DiffAngle1 = Phi0_1-first->_theta;
                          DiffAngle1 = removeAnglePeriod( DiffAngle1 );
                          if (DiffAngle1<CLHEP::halfpi) {
                                  Phi0_1 -= Phi1;
                                  Phi0_1 = removeAnglePeriod(Phi0_1);
                                  _hPhi0HelStep_1->Fill(_iXHelStep[ibin],Phi0_1);
                          }

                          tmpPhi0 = atan(-(sinTheta+cosTheta*cosTheta)/(sinTheta*cosTheta));
                          tmpPhi0 = removeAnglePeriod(tmpPhi0);
                          Phi0_2  = tmpPhi0-CLHEP::halfpi;
                          DiffAngle2 = Phi0_2-second->_theta;
                          DiffAngle2 = removeAnglePeriod( DiffAngle2 );
                          if (DiffAngle2<CLHEP::halfpi) {
                                  Phi0_2 -= Phi2;
                                  Phi0_2 = removeAnglePeriod(Phi0_2);
                                  _hPhi0HelStep_1->Fill(_iXHelStep[ibin],Phi0_2);
                          }
                          tmpPhi0+= CLHEP::pi;
                          tmpPhi0 = removeAnglePeriod(tmpPhi0);
                          Phi0_2  = tmpPhi0-CLHEP::halfpi;
                          DiffAngle2 = Phi0_2-second->_theta;
                          DiffAngle2 = removeAnglePeriod( DiffAngle2 );
                          if (DiffAngle2>threeOverTwoPi) {
                                  Phi0_2 -= Phi2;
                                  Phi0_2 = removeAnglePeriod(Phi0_2);
                                  _hPhi0HelStep_1->Fill(_iXHelStep[ibin],Phi0_2);
                          }


                          /*
                          parA[0] =  cosTheta*cosPhi1lTheta*cosPhi2lTheta;
                          parA[1] = -cosPhi1*cosPhi2lTheta;
                          parA[2] =  cosPhi2*cosPhi1lTheta;
                          parA[3] = -sinTheta*sinPhi1*cosPhi2lTheta;
                          parA[4] =  sinTheta*sinPhi2*cosPhi1lTheta;

                          parB[0] =  cosTheta*(sinPhi1*cosPhi2+cosPhi1*sinPhi2);
                          parB[1] = -(sinPhi1*cosPhi2lTheta+cosPhi1*sinPhi2lTheta);
                          parB[2] = -parB[1];
                          parB[3] =  sinTheta*(cosPhi1*cosPhi2lTheta-sinPhi1*sinPhi2lTheta);
                          parB[4] = -parB[3];

                          parC[0] =  cosTheta*sinPhi1lTheta*sinPhi2lTheta;
                          parC[1] = -sinPhi1*sinPhi2lTheta;
                          parC[2] =  sinPhi2*sinPhi1lTheta;
                          parC[3] =  sinTheta*cosPhi1*sinPhi2lTheta;
                          parC[4] = -sinTheta*cosPhi2*sinPhi1lTheta;

                          sign1=1.0;
                          for (int is1=0; is1<2; is1++){
                                  sign1+=-2*is1;
                                  sign2=1.0;
                                  for (int is2=0; is2<2; is2++) {
                                          sign2+=-2*is2;
                                          cout<<"analyzing s1 "<<sign1<<" s2 "<<sign2<<endl;
                                          tmpVal = DelataX12+(sign1-sign2);//sinTheta
                                          A = tmpVal*parA[0] + sign1*(parA[1] + parA[3]) + sign2*(parA[2] + parA[4]);
                                          B = tmpVal*parB[0] + sign1*(parB[1] + parB[3]) + sign2*(parB[2] + parB[4]);
                                          C = tmpVal*parC[0] + sign1*(parC[1] + parC[3]) + sign2*(parC[2] + parC[4]);
                                          tmpDelta = B*B-4.0*A*C;
                                          if ( tmpDelta<0.0 ) {
                                                  cout<<"Negative discriminant!!!!!!!!"<<endl;
                                                  continue;
                                          }
                                          else if ( tmpDelta==0.00000 ) {
                                                  Phi0_1 = atan( -B/(2.0*A) );
                                                  Phi0_1 = removeAnglePeriod( Phi0_1 );
                                                  Phi0_2 = Phi0_1;
                                          }
                                          else {
                                                  tmpDelta=sqrt(tmpDelta);
                                                  Phi0_1 = atan( (-B+tmpDelta)/(2.0*A) );
                                                  Phi0_1 = removeAnglePeriod( Phi0_1 );
                                                  Phi0_2 = atan( (-B-tmpDelta)/(2.0*A) );
                                                  Phi0_2 = removeAnglePeriod( Phi0_2 );
                                          }
                                          goodPhi0   = true;
                                          DiffAngle1 = Phi1+Phi0_1;
                                          DiffAngle1 = removeAnglePeriod( DiffAngle1 );
                                          DiffAngle1-= first->_theta;
                                          DiffAngle1 = removeAnglePeriod( DiffAngle1 );
                                          DiffAngle2 = Phi2+Phi0_1;
                                          DiffAngle2 = removeAnglePeriod( DiffAngle2 );
                                          DiffAngle2-= second->_theta;
                                          DiffAngle2 = removeAnglePeriod( DiffAngle2 );
                                          if ( sign1>0.0 ) { if ( DiffAngle1<threeOverTwoPi ) { goodPhi0 = false; } }
                                          else if ( DiffAngle1>CLHEP::halfpi ) { goodPhi0 = false; }
                                          if (goodPhi0) {
                                                  if ( sign2>0.0 ) { if ( DiffAngle2<threeOverTwoPi ) { goodPhi0 = false; } }
                                                  else if ( DiffAngle2>CLHEP::halfpi ) { goodPhi0 = false; }
                                          }
                                          if (goodPhi0) {
                                                //...
                                                _hPhi0HelStep_1->Fill(_iXHelStep[ibin],Phi0_1);
                                                cout<<"For Step "<<_iXHelStep[ibin]<<" : s1 "<<sign1<<" s2 "<<sign2<<" phi01 "<<Phi0_1<<endl;
                                          }

                                          goodPhi0   = true;
                                          Phi0_1+=CLHEP::pi;
                                          Phi0_1=removeAnglePeriod( Phi0_1 );
                                          DiffAngle1 = Phi1+Phi0_1;
                                          DiffAngle1 = removeAnglePeriod( DiffAngle1 );
                                          DiffAngle1-= first->_theta;
                                          DiffAngle1 = removeAnglePeriod( DiffAngle1 );
                                          DiffAngle2 = Phi2+Phi0_1;
                                          DiffAngle2 = removeAnglePeriod( DiffAngle2 );
                                          DiffAngle2-= second->_theta;
                                          DiffAngle2 = removeAnglePeriod( DiffAngle2 );
                                          if ( sign1>0.0 ) { if ( DiffAngle1<threeOverTwoPi ) { goodPhi0 = false; } }
                                          else if ( DiffAngle1>CLHEP::halfpi ) { goodPhi0 = false; }
                                          if (goodPhi0) {
                                                  if ( sign2>0.0 ) { if ( DiffAngle2<threeOverTwoPi ) { goodPhi0 = false; } }
                                                  else if ( DiffAngle2>CLHEP::halfpi ) { goodPhi0 = false; }
                                          }
                                          if (goodPhi0) {
                                                //...
                                                _hPhi0HelStep_1->Fill(_iXHelStep[ibin],Phi0_1);
                                                cout<<"For Step "<<_iXHelStep[ibin]<<" : s1 "<<sign1<<" s2 "<<sign2<<" phi01+pi "<<Phi0_1<<endl;
                                          }

                                          goodPhi0   = true;
                                          DiffAngle1 = Phi1+Phi0_2;
                                          DiffAngle1 = removeAnglePeriod( DiffAngle1 );
                                          DiffAngle1-= first->_theta;
                                          DiffAngle1 = removeAnglePeriod( DiffAngle1 );
                                          DiffAngle2 = Phi2+Phi0_2;
                                          DiffAngle2 = removeAnglePeriod( DiffAngle2 );
                                          DiffAngle2-= second->_theta;
                                          DiffAngle2 = removeAnglePeriod( DiffAngle2 );
                                          if ( sign1>0.0 ) { if ( DiffAngle1<threeOverTwoPi ) { goodPhi0 = false; } }
                                          else if ( DiffAngle1>CLHEP::halfpi ) { goodPhi0 = false; }
                                          if (goodPhi0) {
                                                  if ( goodPhi0 && sign2>0.0 ) { if ( DiffAngle2<threeOverTwoPi ) { goodPhi0 = false; } }
                                                  else if ( DiffAngle2>CLHEP::halfpi ) { goodPhi0 = false; }
                                          }
                                          if (goodPhi0) {
                                                //...
                                                _hPhi0HelStep_1->Fill(_iXHelStep[ibin],Phi0_2);
                                                cout<<"For Step "<<_iXHelStep[ibin]<<" : s1 "<<sign1<<" s2 "<<sign2<<" phi02 "<<Phi0_2<<endl;
                                          }

                                          goodPhi0   = true;
                                          Phi0_2+=CLHEP::pi;
                                          Phi0_2=removeAnglePeriod( Phi0_2 );
                                          DiffAngle1 = Phi1+Phi0_2;
                                          DiffAngle1 = removeAnglePeriod( DiffAngle1 );
                                          DiffAngle1-= first->_theta;
                                          DiffAngle1 = removeAnglePeriod( DiffAngle1 );
                                          DiffAngle2 = Phi2+Phi0_2;
                                          DiffAngle2 = removeAnglePeriod( DiffAngle2 );
                                          DiffAngle2-= second->_theta;
                                          DiffAngle2 = removeAnglePeriod( DiffAngle2 );
                                          if ( sign1>0.0 ) { if ( DiffAngle1<threeOverTwoPi ) { goodPhi0 = false; } }
                                          else if ( DiffAngle1>CLHEP::halfpi ) { goodPhi0 = false; }
                                          if (goodPhi0) {
                                                  if ( goodPhi0 && sign2>0.0 ) { if ( DiffAngle2<threeOverTwoPi ) { goodPhi0 = false; } }
                                                  else if ( DiffAngle2>CLHEP::halfpi ) { goodPhi0 = false; }
                                          }
                                          if (goodPhi0) {
                                                //...
                                                _hPhi0HelStep_1->Fill(_iXHelStep[ibin],Phi0_2);
                                                cout<<"For Step "<<_iXHelStep[ibin]<<" : s1 "<<sign1<<" s2 "<<sign2<<" phi02+pi "<<Phi0_2<<endl;
                                          }
                                  }
                          }*/

//                          ////if (Phi1==0.0 || /*Phi1==-CLHEP::pi ||*/Phi1==CLHEP::twopi || Phi1==CLHEP::pi) {Phi0lTheta1 = Phi1+CLHEP::halfpi; cout<<"No der"<<endl;}
//                          //if (Phi1==0.0 || /*Phi1==-CLHEP::pi ||*/Phi1==CLHEP::twopi ) Phi0lTheta1 = CLHEP::halfpi;
//                          //else if (Phi1==CLHEP::pi ) Phi0lTheta1 = -CLHEP::halfpi;
//                          //else {
//                          //        //Phi0lTheta1  = std::atan(2.0*std::tan(Phi1+CLHEP::halfpi));
//                          //        //if (Phi1>=0.0 && Phi1<CLHEP::halfpi) Phi0lTheta1+=CLHEP::pi;
//                          //        //if (Phi1>CLHEP::halfpi && Phi1<CLHEP::pi) Phi0lTheta1-=CLHEP::pi;
//                          //        Phi0lTheta1  = std::atan(2.0*std::tan(Phi1+CLHEP::halfpi));
//                          //        if (Phi1<=CLHEP::pi) Phi0lTheta1+=CLHEP::pi;
//                          //}
//                          //Phi0lTheta1=removeAnglePeriod(Phi0lTheta1);
//                          tmpVal = Phi1+Phi0lTheta1;
//                          tmpVal = removeAnglePeriod(tmpVal);
//                          if (tmpVal>=3.0*CLHEP::halfpi) {
//                                  if (Phi1>CLHEP::pi/*Phi1>=0.0&&Phi1<CLHEP::pi*/) continue;
//                                  Phi0lTheta1=0.0;
//                          }
//                          else if (tmpVal>=CLHEP::halfpi) continue;
//                          Phi01 = Phi0lTheta1+first->_theta;
//                          Phi01 = removeAnglePeriod(Phi01);
//                          _hPhi0HelStep_1->Fill(_iXHelStep[ibin],Phi01);
//                          cosPhi1 = cos(Phi01);
//                          sinPhi1 = sin(Phi01);
//
//                          Phi2        = tmpVal*DeltaZ2;
//                          Phi2        = removeAnglePeriod(Phi2);
//                          //if (Phi2==0.0 || /*Phi2==-CLHEP::pi ||*/Phi2==CLHEP::twopi || Phi2==CLHEP::pi) {Phi0lTheta2 = Phi2+CLHEP::halfpi; cout<<"No der"<<endl;}
//                          if (Phi2==0.0 || /*Phi1==-CLHEP::pi ||*/Phi2==CLHEP::twopi ) Phi0lTheta2 = CLHEP::halfpi;
//                          else if (Phi2==CLHEP::pi ) Phi0lTheta2 = -CLHEP::halfpi;
//                          else {
//                                  //Phi0lTheta2  = std::atan(2.0*std::tan(Phi2+CLHEP::halfpi));
//                                  //if (Phi2>=0.0 && Phi2<CLHEP::halfpi) Phi0lTheta2+=CLHEP::pi;
//                                  //if (Phi2>CLHEP::halfpi && Phi2<CLHEP::pi) Phi0lTheta2-=CLHEP::pi;
//                                  Phi0lTheta2  = std::atan(2.0*std::tan(Phi2+CLHEP::halfpi));
//                                  if (Phi2<=CLHEP::pi) Phi0lTheta2+=CLHEP::pi;
//                          }
//                          Phi0lTheta2=removeAnglePeriod(Phi0lTheta2);
//                          tmpVal = Phi2+Phi0lTheta2;
//                          tmpVal = removeAnglePeriod(tmpVal);
//                          if (tmpVal>=3.0*CLHEP::halfpi) {
//                                  if (Phi2>CLHEP::pi/*Phi2>=0.0*/) continue;
//                                  Phi0lTheta2=0.0;
//                          }
//                          else if (tmpVal>=CLHEP::halfpi) continue;
//                          Phi02 = Phi0lTheta2+second->_theta;
//                          Phi02 = removeAnglePeriod(Phi02);
//                          _hPhi0HelStep_1->Fill(_iXHelStep[ibin],Phi02);
//                          cosPhi2 = cos(Phi02);
//                          sinPhi2 = sin(Phi02);
//
//                          tmpDen1 = 1.0/( cos(Phi0lTheta1)*(cos(Phi1)-1.0) - sin(Phi0lTheta1)*sin(Phi1) );
//                          tmpDen2 = 1.0/( cos(Phi0lTheta2)*(cos(Phi2)-1.0) - sin(Phi0lTheta2)*sin(Phi2) );
//                          tmpXo=-maxXY+0.5*stepXY;



/*                          for ( int ix=0; ix<nXYStep; ++ix ) {
                                  tmpYo=-maxXY+0.5*stepXY;
                                  for ( int iy=0; iy<nXYStep; ++iy ) {
                                          tmpVal = tmpXo*first->_direct.getX() + tmpYo*first->_direct.getX();
                                          tmpNum1 = first->_radius - tmpVal;
                                          tmpR1   = tmpNum1*tmpDen1;
                                          if (tmpR1<=0.0) continue;
                                          tmpTwopiTanLambda1 =  _iXHelStep[ibin]/tmpR1;
                                          if ( tmpTwopiTanLambda1>twoPitanLambdaMax || (tmpTwopiTanLambda1>-twopiTanLambdaMin && tmpTwopiTanLambda1<twopiTanLambdaMin) || tmpTwopiTanLambda1<-twoPitanLambdaMax ) continue;
                                          tmpXC = tmpXo - tmpR1*cosPhi1;
                                          tmpYC = tmpYo - tmpR1*sinPhi1;
                                          tmpCdist = sqrt( pow(tmpXC,2) + pow(tmpYC,2) );
                                          tmpCspread = tmpCdist - tmpR1;
                                          if ( tmpCspread>-targetRadMax && tmpCspread<targetRadMax ) {
                                                  _hRHelStep->Fill(_iXHelStep[ibin],tmpR1);
                                                  _hPhi0->Fill(Phi01);
                                                  _hR->Fill(tmpR1);
                                                  _hPhi0HelStep_1->Fill(_iXHelStep[ibin],Phi01);

                                                  tmpNum2 = second->_radius - tmpVal;
                                                  tmpR2   = tmpNum2*tmpDen2;
                                                  if (tmpR2<=0.0) continue;
                                                  tmpTwopiTanLambda2 =  _iXHelStep[ibin]/tmpR2;
                                                  if ( tmpTwopiTanLambda2>twoPitanLambdaMax || (tmpTwopiTanLambda2>-twopiTanLambdaMin && tmpTwopiTanLambda2<twopiTanLambdaMin) || tmpTwopiTanLambda2<-twoPitanLambdaMax ) continue;
                                                  tmpXC = tmpXo - tmpR2*cosPhi2;
                                                  tmpYC = tmpYo - tmpR2*sinPhi2;
                                                  tmpCdist = sqrt( pow(tmpXC,2) + pow(tmpYC,2) );
                                                  tmpCspread = tmpCdist - tmpR2;
                                                  if ( tmpCspread>-targetRadMax && tmpCspread<targetRadMax ) {
                                                          _hRHelStep->Fill(_iXHelStep[ibin],tmpR2);
                                                          _hPhi0->Fill(Phi02);
                                                          _hR->Fill(tmpR2);
                                                          _hPhi0HelStep_1->Fill(_iXHelStep[ibin],Phi02);

                                                  }
                                          }
                                          tmpYo+=stepXY;
                                  }
                                  tmpXo+=stepXY;
                          }*/
                  }
          }
  }

//  void TrackReco::computeCombination( double minHStep, double maxHStep ){
//          cout<<"-------- minHStep "<<minHStep<<" maxHStep "<<maxHStep<<endl;
//          int firstXbin=0, lastXBin=_nBinHelpStep;
//          if ( minHStep>_minHelStep) {
//                  for (int ibin=0; ibin<_nBinHelpStep; ++ibin){
//                          if (minHStep<=_iXHelStep[ibin]) {
//                                  firstXbin=ibin;
//                                  break;
//                          }
//                  }
//          }
//          if ( maxHStep>0.0 && maxHStep<_maxHelStep ) {
//                  lastXBin=firstXbin;
//                  for (int ibin=firstXbin; ibin<_nBinHelpStep; ++ibin){
//                          if (maxHStep<=_iXHelStep[ibin]) {
//                                  break;
//                          }
//                          ++lastXBin;
//                  }
//          }
//          cout<<"-------- firstXbin "<<firstXbin<<" lastXBin "<<lastXBin<<endl;
//          //firstXbin=0;
//
//
//          double DeltaZ, chord, tmpR;
//          for ( vector< std::pair<int, int> >::iterator hitsCouplings_it = _hitsCouplings.begin(); hitsCouplings_it != _hitsCouplings.end(); ++hitsCouplings_it ) {
//                  //cout<<"I'm analyzing the couple: "<<hitsCouplings_it->first<<" "<<hitsCouplings_it->second<<endl;
//
//                  HitData *first  = strdat->find(hitsCouplings_it->first)->second;
//                  HitData *second = strdat->find(hitsCouplings_it->second)->second;
//                  DeltaZ = second->_midp.getZ() - first->_midp.getZ();
//                  for ( vector< pair<double, double> >::iterator fWXY_it = first->_wireXYTable.begin(); fWXY_it != first->_wireXYTable.end(); ++fWXY_it ) {
//                          for ( vector< pair<double, double> >::iterator sWXY_it = second->_wireXYTable.begin(); sWXY_it != second->_wireXYTable.end(); ++sWXY_it ) {
//                                  chord = sqrt( pow( sWXY_it->first - fWXY_it->first ,2) + pow( sWXY_it->second - fWXY_it->second ,2) );
//                                  for ( int ibin=firstXbin; ibin<lastXBin; ++ibin ) {
//                                          tmpR=2.0*sin(CLHEP::pi*DeltaZ/_iXHelStep[ibin]);
//                                          tmpR=chord/tmpR;
//                                          _hRHelStep->Fill(_iXHelStep[ibin],tmpR);                                  }
//                          }
//                  }
//          }
//  }

  void TrackReco::computeTripleCombination( double minHStep, double maxHStep ){
          cout<<"In TrackReco::computeTripleCombination"<<endl;
          cout<<"-------- minHStep "<<minHStep<<" maxHStep "<<maxHStep<<endl;
          int firstXbin=0, lastXBin=_nBinHelpStep;
          if ( minHStep>_minHelStep) {
                  for (int ibin=0; ibin<_nBinHelpStep; ++ibin){
                          if (minHStep<=_iXHelStep[ibin]) {
                                  firstXbin=ibin;
                                  break;
                          }
                  }
          }
          if ( maxHStep>0.0 && maxHStep<_maxHelStep ) {
                  lastXBin=firstXbin;
                  for (int ibin=firstXbin; ibin<_nBinHelpStep; ++ibin){
                          if (maxHStep<=_iXHelStep[ibin]) {
                                  break;
                          }
                          ++lastXBin;
                  }
          }
          cout<<"-------- firstXbin "<<firstXbin<<" lastXBin "<<lastXBin<<endl;
          //firstXbin=0;

          double DeltaZ21, DeltaZ31, DeltaZ1;
          double DeltaPhi21, cosDPhi21_1, sinDPhi21, DeltaPhi31, cosDPhi31_1, sinDPhi31;
          double Phi1, tanPhi1;
          double X2X1, X3X1, XsRatio;
          double tmpVal;
          double turnsCut = 5.0*CLHEP::twopi;
          int tmpiBinHit;
          double tmpBinVal;
          _maxPhi0HelStep=0.0;

          for ( std::map< int, std::multimap<int, int> >::iterator hitsTripleCouplings_it = _hitsTripleCouplings.begin();
                  hitsTripleCouplings_it != _hitsTripleCouplings.end(); ++hitsTripleCouplings_it ) 
          {
                  HitData *first  = strdat->find(hitsTripleCouplings_it->first)->second;
                  DeltaZ1 = first->_midp.getZ() - zo;
                  for ( std::multimap<int, int>::iterator secondHitsPairs_it = hitsTripleCouplings_it->second.begin();
                          secondHitsPairs_it != hitsTripleCouplings_it->second.end(); ++secondHitsPairs_it )
                  {
                          HitData *second = strdat->find(secondHitsPairs_it->first)->second;
                          HitData *third  = strdat->find(secondHitsPairs_it->second)->second;
                          DeltaZ21 = second->_midp.getZ() - first->_midp.getZ();
                          DeltaZ31 = third->_midp.getZ() - first->_midp.getZ();
                          X2X1 = second->_radius - first->_radius;
                          X3X1 = third->_radius - first->_radius;
                          XsRatio = X3X1/X2X1;
                          for ( int ibin=firstXbin; ibin<lastXBin; ++ibin ) {
                                  tmpVal      = CLHEP::twopi/_iXHelStep[ibin];
                                  DeltaPhi21  = tmpVal*DeltaZ21;
                                  DeltaPhi31  = tmpVal*DeltaZ31;
                                  Phi1        = tmpVal*DeltaZ1;
                                  if (Phi1>turnsCut || DeltaPhi21>turnsCut|| DeltaPhi31>turnsCut) continue;
                                  cosDPhi21_1 = cos(DeltaPhi21) - 1.0;
                                  cosDPhi31_1 = cos(DeltaPhi31) - 1.0;
                                  sinDPhi21   = sin(DeltaPhi21);
                                  sinDPhi31   = sin(DeltaPhi31);
                                  tanPhi1     = tan(Phi1);

                                  tmpVal = XsRatio*( first->_direct.getY()*cosDPhi21_1 -first->_direct.getX()*sinDPhi21 );
                                  tmpVal -= ( first->_direct.getY()*cosDPhi31_1 -first->_direct.getX()*sinDPhi31 );
                                  tmpVal /= ( XsRatio*( first->_direct.getX()*cosDPhi21_1 +first->_direct.getY()*sinDPhi21 )
                                                  -( first->_direct.getX()*cosDPhi31_1 +first->_direct.getY()*sinDPhi31 ) );

                                  tmpiBinHit=_hPhi0HelStep->Fill( _iXHelStep[ibin], /*atan2( tmpVal-tanPhi1, 1.0+tmpVal*tanPhi1 )*/ atan( (tmpVal-tanPhi1)/(1.0+tmpVal*tanPhi1) ) );
                                  if (tmpiBinHit>0) {
                                          tmpBinVal=_hPhi0HelStep->GetBinContent(tmpiBinHit);
                                          if (tmpBinVal>_maxPhi0HelStep) _maxPhi0HelStep=tmpBinVal;
                                          voteArrHitMap::iterator voteArr_Bin_Hit_rel_it = voteArr_Bin_Hit_rel.find(tmpiBinHit);
                                          if ( voteArr_Bin_Hit_rel_it==voteArr_Bin_Hit_rel.end() ) {
                                                  voteArr_Bin_Hit_rel_it = voteArr_Bin_Hit_rel.insert( voteArrHitMap::value_type( tmpiBinHit, strawList()/* inBinHitList*/ ) ).first;
                                                  voteArr_Bin_Hit_rel_it->second.clear();
                                          }
                                          voteArr_Bin_Hit_rel_it->second.insert(hitsTripleCouplings_it->first);
                                          voteArr_Bin_Hit_rel_it->second.insert(secondHitsPairs_it->first);
                                          voteArr_Bin_Hit_rel_it->second.insert(secondHitsPairs_it->second);
                                  }
                          }
                  }
          }
  }

  void TrackReco::computeTripleCombination_DeltaRay(){
          cout<<"In TrackReco::computeTripleCombination_DeltaRay"<<endl;
          int firstXbin=0, lastXBin=_nBinHelpStep_DR;

          double DeltaZ21, DeltaZ31, DeltaZ1;
          double DeltaPhi21, cosDPhi21_1, sinDPhi21, DeltaPhi31, cosDPhi31_1, sinDPhi31;
          double Phi1, tanPhi1;
          double X2X1, X3X1, XsRatio;
          double tmpVal;
          double turnsCut = 5.0*CLHEP::twopi;
          int tmpiBinHit;

          for ( std::map< int, std::multimap<int, int> >::iterator hitsTripleCouplings_it = _hitsTripleCouplings.begin();
                  hitsTripleCouplings_it != _hitsTripleCouplings.end(); ++hitsTripleCouplings_it )
          {
                  HitData *first  = strdat->find(hitsTripleCouplings_it->first)->second;
                  DeltaZ1 = first->_midp.getZ() - zo;
                  for ( std::multimap<int, int>::iterator secondHitsPairs_it = hitsTripleCouplings_it->second.begin();
                          secondHitsPairs_it != hitsTripleCouplings_it->second.end(); ++secondHitsPairs_it )
                  {
                          HitData *second = strdat->find(secondHitsPairs_it->first)->second;
                          HitData *third  = strdat->find(secondHitsPairs_it->second)->second;
                          DeltaZ21 = second->_midp.getZ() - first->_midp.getZ();
                          DeltaZ31 = third->_midp.getZ() - first->_midp.getZ();
                          X2X1 = second->_radius - first->_radius;
                          X3X1 = third->_radius - first->_radius;
                          XsRatio = X3X1/X2X1;
                          for ( int ibin=firstXbin; ibin<lastXBin; ++ibin ) {
                                  tmpVal      = CLHEP::twopi/_iXHelStep[ibin];
                                  DeltaPhi21  = tmpVal*DeltaZ21;
                                  DeltaPhi31  = tmpVal*DeltaZ31;
                                  Phi1        = tmpVal*DeltaZ1;
                                  if (Phi1>turnsCut || DeltaPhi21>turnsCut|| DeltaPhi31>turnsCut) continue;
                                  cosDPhi21_1 = cos(DeltaPhi21) - 1.0;
                                  cosDPhi31_1 = cos(DeltaPhi31) - 1.0;
                                  sinDPhi21   = sin(DeltaPhi21);
                                  sinDPhi31   = sin(DeltaPhi31);
                                  tanPhi1     = tan(Phi1);

                                  tmpVal = XsRatio*( first->_direct.getY()*cosDPhi21_1 -first->_direct.getX()*sinDPhi21 );
                                  tmpVal -= ( first->_direct.getY()*cosDPhi31_1 -first->_direct.getX()*sinDPhi31 );
                                  tmpVal /= ( XsRatio*( first->_direct.getX()*cosDPhi21_1 +first->_direct.getY()*sinDPhi21 )
                                                  -( first->_direct.getX()*cosDPhi31_1 +first->_direct.getY()*sinDPhi31 ) );

                                  tmpiBinHit=_hPhi0HelStep_DeltaRay->Fill( _iXHelStep_DR[ibin], /*atan2( tmpVal-tanPhi1, 1.0+tmpVal*tanPhi1 )*/ atan( (tmpVal-tanPhi1)/(1.0+tmpVal*tanPhi1) ) );
                                  if (tmpiBinHit>0) {
                                          voteArr_Bin_Comb_rel_DR.insert( voteArrComMap::value_type( tmpiBinHit, std::pair<int, std::pair<int,int> >(hitsTripleCouplings_it->first, std::pair<int,int>(secondHitsPairs_it->first,secondHitsPairs_it->second) ) ) );

                                          voteArrHitMap::iterator voteArr_Bin_Hit_rel_DR_it = voteArr_Bin_Hit_rel_DR.find(tmpiBinHit);
                                          if ( voteArr_Bin_Hit_rel_DR_it==voteArr_Bin_Hit_rel_DR.end() ) {
                                                  voteArr_Bin_Hit_rel_DR_it = voteArr_Bin_Hit_rel_DR.insert( voteArrHitMap::value_type( tmpiBinHit, strawList()/* inBinHitList*/ ) ).first;
                                                  voteArr_Bin_Hit_rel_DR_it->second.clear();
                                          }
                                          voteArr_Bin_Hit_rel_DR_it->second.insert(hitsTripleCouplings_it->first);
                                          voteArr_Bin_Hit_rel_DR_it->second.insert(secondHitsPairs_it->first);
                                          voteArr_Bin_Hit_rel_DR_it->second.insert(secondHitsPairs_it->second);
                                  }
                          }
                  }
          }
  }

  void TrackReco::measureRbyTripleCombination( double HStep, double HPhi0 ){
          double DeltaZ21, DeltaZ31, DeltaZ1;
          double DeltaPhi21, cosDPhi21_1, sinDPhi21, DeltaPhi31, cosDPhi31_1, sinDPhi31;
          double Phi1/*, tanPhi1*/;
          double X2X1, X3X1/*, XsRatio*/;
          double tmpR21, tmpR31;
          double tmpVal = CLHEP::twopi/HStep;
          double cosPhi01, sinPhi01;
          double tmpAngle, HPhi0compl;
          HPhi0compl=HPhi0+CLHEP::pi;

          for ( std::map< int, std::multimap<int, int> >::iterator hitsTripleCouplings_it = _hitsTripleCouplings.begin();
                  hitsTripleCouplings_it != _hitsTripleCouplings.end(); ++hitsTripleCouplings_it )
          {
                  HitData *first  = strdat->find(hitsTripleCouplings_it->first)->second;
                  DeltaZ1 = first->_midp.getZ() - zo;
                  for ( std::multimap<int, int>::iterator secondHitsPairs_it = hitsTripleCouplings_it->second.begin();
                          secondHitsPairs_it != hitsTripleCouplings_it->second.end(); ++secondHitsPairs_it )
                  {
                          HitData *second = strdat->find(secondHitsPairs_it->first)->second;
                          HitData *third  = strdat->find(secondHitsPairs_it->second)->second;
                          DeltaZ21 = second->_midp.getZ() - first->_midp.getZ();
                          DeltaZ31 = third->_midp.getZ() - first->_midp.getZ();
                          X2X1 = second->_radius - first->_radius;
                          X3X1 = third->_radius - first->_radius;
                          DeltaPhi21  = tmpVal*DeltaZ21;
                          DeltaPhi31  = tmpVal*DeltaZ31;
                          Phi1        = tmpVal*DeltaZ1;
                          cosDPhi21_1 = cos(DeltaPhi21) - 1.0;
                          cosDPhi31_1 = cos(DeltaPhi31) - 1.0;
                          sinDPhi21   = sin(DeltaPhi21);
                          sinDPhi31   = sin(DeltaPhi31);

                          tmpAngle=HPhi0+Phi1;
                          cosPhi01=cos(tmpAngle);
                          sinPhi01=sin(tmpAngle);
                          tmpR21 = X2X1/( first->_direct.getY()*( cosDPhi21_1*cosPhi01 - sinDPhi21*sinPhi01 ) -first->_direct.getX()*( cosDPhi21_1*sinPhi01 + sinDPhi21*cosPhi01 ) );
                          tmpR31 = X3X1/( first->_direct.getY()*( cosDPhi31_1*cosPhi01 - sinDPhi31*sinPhi01 ) -first->_direct.getX()*( cosDPhi31_1*sinPhi01 + sinDPhi31*cosPhi01 ) );
                          _hR_1->Fill(tmpR21);
                          _hR_1->Fill(tmpR31);
                          tmpAngle=HPhi0compl+Phi1;
                          cosPhi01=cos(tmpAngle);
                          sinPhi01=sin(tmpAngle);
                          tmpR21 = X2X1/( first->_direct.getY()*( cosDPhi21_1*cosPhi01 - sinDPhi21*sinPhi01 ) -first->_direct.getX()*( cosDPhi21_1*sinPhi01 + sinDPhi21*cosPhi01 ) );
                          tmpR31 = X3X1/( first->_direct.getY()*( cosDPhi31_1*cosPhi01 - sinDPhi31*sinPhi01 ) -first->_direct.getX()*( cosDPhi31_1*sinPhi01 + sinDPhi31*cosPhi01 ) );
                          _hR_1->Fill(tmpR21);
                          _hR_1->Fill(tmpR31);
                  }
          }
  }

  int TrackReco::butterflyFilter( float *in, float *out, int &nBinX, int &nBinY, float minCountCut ){
          // the 3x3 mask for this application is:
          //    0 -3  0
          //    2  2  2
          //    0 -3  0

          float maskY = -3.0;
          // the mask on X is {2.0, 2.0, 2.0} but for optimization is better to directly sum the 3 bins and after multiply by 2

          int nOutEntries=0;
          int tmpRow, tmpRowUp, tmpRowDown, tmpCol;
          int nBinXeff = nBinX+2;
          int nBinYeff = nBinY+2;

          float *buff = new float [nBinXeff*nBinYeff];
          for (int ib=0; ib<nBinXeff*nBinYeff; ib++) buff[ib]=0.00000;

          for (int j=1; j<(nBinY-1); j++){
                  tmpRowDown = j*nBinXeff;
                  tmpRow     = tmpRowDown + nBinXeff; //(j+1)*nBinXeff
                  tmpRowUp   = tmpRow + nBinXeff;     //(j+2)*nBinXeff
                  for (int i=0; i<nBinX; i++){
                          tmpCol=i+1;
                          //buff[tmpRow+tmpCol] = maskY[0]*in[tmpRowDown+tmpCol] /*+ maskY[1]*in[tmpRow+tmpCol]*/ + maskY[2]*in[tmpRowUp+tmpCol];  //the center element will be computed during the sum on rows
                          buff[tmpRow+tmpCol] = maskY*(in[tmpRowDown+tmpCol] + in[tmpRowUp+tmpCol]);  //the center element will be computed during the sum on rows
                  }
          }

          int tmpColR, tmpColL;
          for (int i=1; i<(nBinX-1); i++){
                  tmpColL = i;
                  tmpCol  = tmpColL+1;
                  tmpColR = tmpCol+1;
                  for (int j=1; j<(nBinY-1); j++){
                          tmpRow     = (j+1)*nBinXeff;
                          out[tmpRow+tmpCol] = 2.0*( in[tmpRow+tmpColL] + in[tmpRow+tmpCol] + in[tmpRow+tmpColR] ) + buff[tmpRow+tmpCol];
                          if (out[tmpRow+tmpCol]<minCountCut) out[tmpRow+tmpCol]=0.0;
                          else nOutEntries+=(int)out[tmpRow+tmpCol];
                  }
          }

          delete [] buff;
          return nOutEntries;
  }

  int TrackReco::butterflyFilterRot( float *in, float *out, int &nBinX, int &nBinY, float minCountCut ){
          // the 3x3 mask for this application is:
          //    0  2  0
          //   -3  2 -3
          //    0  2  0

          float maskX = -3.0;
          // the mask on Y is {2.0, 2.0, 2.0} but for optimization is better to directly sum the 3 bins and after multiply by 2

          int nOutEntries=0;
          int tmpRow, tmpRowUp, tmpRowDown, tmpCol;
          int nBinXeff = nBinX+2;
          int nBinYeff = nBinY+2;

          float *buff = new float [nBinXeff*nBinYeff];
          for (int ib=0; ib<nBinXeff*nBinYeff; ib++) buff[ib]=0.00000;

          for (int j=1; j<(nBinY-1); j++){
                  tmpRowDown = j*nBinXeff;
                  tmpRow     = tmpRowDown + nBinXeff; //(j+1)*nBinXeff
                  tmpRowUp   = tmpRow + nBinXeff;     //(j+2)*nBinXeff
                  for (int i=0; i<nBinX; i++){
                          tmpCol=i+1;
                          buff[tmpRow+tmpCol] = 2.0*(in[tmpRowDown+tmpCol] +in[tmpRow+tmpCol] + in[tmpRowUp+tmpCol]);  //the center element will be computed during the sum on rows
                  }
          }

          int tmpColR, tmpColL;
          for (int i=1; i<(nBinX-1); i++){
                  tmpColL = i;
                  tmpCol  = tmpColL+1;
                  tmpColR = tmpCol+1;
                  for (int j=1; j<(nBinY-1); j++){
                          tmpRow     = (j+1)*nBinXeff;
                          out[tmpRow+tmpCol] = maskX*(in[tmpRow+tmpColL] + in[tmpRow+tmpColR] ) + buff[tmpRow+tmpCol];
                          if (out[tmpRow+tmpCol]<minCountCut) out[tmpRow+tmpCol]=0.0;
                          else nOutEntries+=(int)out[tmpRow+tmpCol];
                  }
          }

          delete [] buff;
          return nOutEntries;
  }

  int TrackReco::butterflyFilterRot45( float *in, float *out, int &nBinX, int &nBinY, float minCountCut ){
          // the 3x3 mask for this application is:
          //  (  0  1  0  )
          //  (  1  1  1  )/5
          //  (  0  1  0  )
          ////    0  0 -1
          ////    0  2  0
          ////   -1  0  0
          ////   -3  0  2
          ////    0  2  0
          ////    2  0 -3
          ////    2  0 -3
          ////    0  2  0
          ////   -3  0  2

          //float maskD = -3.0;//0;
          // the mask on Y is {2.0, 2.0, 2.0} but for optimization is better to directly sum the 3 bins and after multiply by 2

          int nOutEntries=0;
          int tmpRow, tmpRowUp, tmpRowDown, tmpCol;
          int nBinXeff = nBinX+2;
          int nBinYeff = nBinY+2;
          int tmpColR, tmpColL;

          float *buff = new float [nBinXeff*nBinYeff];
          for (int ib=0; ib<nBinXeff*nBinYeff; ib++) buff[ib]=0.00000;

//          for (int j=1; j<(nBinY-1); j++){
//                  tmpRowDown = j*nBinXeff;
//                  tmpRow     = tmpRowDown + nBinXeff; //(j+1)*nBinXeff
//                  tmpRowUp   = tmpRow + nBinXeff;     //(j+2)*nBinXeff
//                  for (int i=1; i<nBinX; i++){
//                          tmpColL = i;
//                          tmpCol  = tmpColL+1;
//                          tmpColR = tmpCol+1;
//                          buff[tmpRow+tmpCol] = -(in[tmpRowDown+tmpColL] + in[tmpRowUp+tmpColR]);  //the center element will be computed during the sum on rows
//                          //buff[tmpRow+tmpCol] = 2.0*(in[tmpRowDown+tmpColL] +in[tmpRow+tmpCol] + in[tmpRowUp+tmpColR]);  //the center element will be computed during the sum on rows
//                          //buff[tmpRow+tmpCol] = 2.0*(in[tmpRowDown+tmpColR] +in[tmpRow+tmpCol] + in[tmpRowUp+tmpColL]);  //the center element will be computed during the sum on rows
//                  }
//          }

          for (int i=1; i<nBinX; i++){
                  tmpColL = i;
                  tmpCol  = tmpColL+1;
                  tmpColR = tmpCol+1;
                  for (int j=1; j<(nBinY-1); j++){
                          tmpRowDown = j*nBinXeff;
                          tmpRow     = tmpRowDown + nBinXeff; //(j+1)*nBinXeff
                          tmpRowUp   = tmpRow + nBinXeff;     //(j+2)*nBinXeff

                          buff[tmpRow+tmpCol] = in[tmpRow+tmpColL]+in[tmpRow+tmpCol]+in[tmpRow+tmpColR];
                          out[tmpRow+tmpCol]  = in[tmpRowDown+tmpCol] + in[tmpRowUp+tmpCol] + buff[tmpRow+tmpCol];
                          out[tmpRow+tmpCol] *= 0.2;

                          //out[tmpRow+tmpCol] = 2.0*in[tmpRow+tmpCol] + buff[tmpRow+tmpCol];
                          //out[tmpRow+tmpCol] = maskD*(in[tmpRowUp+tmpColL] + in[tmpRowDown+tmpColR] ) + buff[tmpRow+tmpCol];
                          //cout<<"Buff "<<buff[tmpRow+tmpCol]<<" out "<<out[tmpRow+tmpCol]<<endl;
                          //out[tmpRow+tmpCol] = maskD*(in[tmpRowUp+tmpColR] + in[tmpRowDown+tmpColL] ) + buff[tmpRow+tmpCol];
                          if (out[tmpRow+tmpCol]<minCountCut) out[tmpRow+tmpCol]=0.0;
                          else nOutEntries+=(int)out[tmpRow+tmpCol];
                  }
          }

          delete [] buff;
          return nOutEntries;
  }

  int TrackReco::blur( float *in, float *out, int &nBinX, int &nBinY, float sigmaX, float sigmaY, float binSizeX, float binSizeY ){
          // the blur is a circular (2D) gaussian smoothing but it is equivalent to apply a 1D gaussian smoothing in X and after in Y

          if (blur_firstTime) {
                  blur_nBinXHalfWdt  = (int)floor(3.0*sigmaX/binSizeX);  //the window size on witch apply the smoothing is 2*nBin+1
                  blur_nBinYHalfWdt  = (int)floor(3.0*sigmaY/binSizeY);
                  blur_nBinXHalfWdt1 = blur_nBinXHalfWdt+1;
                  blur_nBinYHalfWdt1 = blur_nBinYHalfWdt+1;
                  blur_coeffX     = new float [blur_nBinXHalfWdt1];
                  blur_coeffY     = new float [blur_nBinYHalfWdt1];

                  Genfun::Erf erff;
                  float lastPVal, fact;

                  fact      = 3.0/((float)blur_nBinXHalfWdt + 0.5);
                  lastPVal  = erff(0.5*fact*0.707106781);                //1/sqrt(2)=0.707106781
                  blur_coeffX[0] = lastPVal;
                  for (int icf=1; icf<blur_nBinXHalfWdt1; icf++) {
                          blur_coeffX[icf] = lastPVal;
                          lastPVal    = erff((0.5+(float)icf)*fact*0.707106781);
                          blur_coeffX[icf] = 0.5*(lastPVal-blur_coeffX[icf]);
                  }
                  fact      = 3.0/((float)blur_nBinYHalfWdt + 0.5);
                  lastPVal  = erff(0.5*fact*0.707106781);                //1/sqrt(2)=0.707106781
                  blur_coeffY[0] = lastPVal;
                  for (int icf=1; icf<blur_nBinYHalfWdt1; icf++) {
                          blur_coeffY[icf] = lastPVal;
                          lastPVal    = erff((0.5+(float)icf)*fact*0.707106781);
                          blur_coeffY[icf] = 0.5*(lastPVal-blur_coeffY[icf]);
                  }
                  blur_firstTime=false;
          }


          int nOutEntries=0;
          int tmpRow, tmpBin, tmpBin1, tmpBin2;
          int nBinXeff = nBinX+2;
          int nBinYeff = nBinY+2;

          float *buff = new float [nBinXeff*nBinYeff];
          for (int ib=0; ib<nBinXeff*nBinYeff; ib++) buff[ib]=0.00000;

          for (int j=1; j<=nBinY; j++){
                  tmpRow     = j*nBinXeff;
                  for (int i=1; i<=nBinX; i++){
                          tmpBin       = tmpRow+i;
                          buff[tmpBin] = blur_coeffX[0]*in[tmpBin];
                          for (int kb=1; kb<blur_nBinXHalfWdt1; kb++){
                                  tmpBin1 = tmpBin - ( ((i-kb)<1) ? -kb: kb );
                                  tmpBin2 = tmpBin + ( ((i+kb)>=nBinX) ? -kb: kb );
                                  buff[tmpBin] += blur_coeffX[kb]*(in[tmpBin1] + in[tmpBin2]);
                          }
                  }
          }

          for (int i=1; i<=nBinX; i++){
                  for (int j=1; j<=nBinY; j++){
                          tmpRow       = j*nBinXeff;
                          tmpBin       = tmpRow+i;
                          out[tmpBin]  = blur_coeffY[0]*buff[tmpBin];
                          for (int kb=1; kb<blur_nBinYHalfWdt1; kb++) {
                                  tmpBin1 = tmpBin - ( ((j-kb)<1) ? -kb: kb )*nBinXeff;
                                  tmpBin2 = tmpBin + ( ((j+kb)>=nBinY) ? -kb: kb )*nBinXeff;
                                  out[tmpBin] += blur_coeffY[kb]*(buff[tmpBin1] + buff[tmpBin2]);
                          }
                          if (out[tmpBin]>0.0) {
                                  out[tmpBin]=floor(out[tmpBin]+0.5);
                                  nOutEntries+=(int)out[tmpBin];
                          }
                  }
          }

          delete [] buff;

          return nOutEntries;
  }

  int TrackReco::smooth( float *in, float *out, int &nBinX, int &nBinY, float sigmaX, float sigmaY, float binSizeX, float binSizeY, float mincut ){
          // (2D) gaussian smoothing but it is equivalent to apply a 1D gaussian smoothing in X and after in Y

          int smooth_nBinXHalfWdt  = (int)floor(3.0*sigmaX/binSizeX);  //the window size on witch apply the smoothing is 2*nBin+1
          int smooth_nBinYHalfWdt  = (int)floor(3.0*sigmaY/binSizeY);
          int smooth_nBinXHalfWdt1 = smooth_nBinXHalfWdt+1;
          int smooth_nBinYHalfWdt1 = smooth_nBinYHalfWdt+1;

          int nOutEntries=0;
          int tmpRow, tmpBin, tmpBin1, tmpBin2;
          int nBinXeff = nBinX+2;
          int nBinYeff = nBinY+2;

          float *buff = new float [nBinXeff*nBinYeff];
          for (int ib=0; ib<nBinXeff*nBinYeff; ib++) buff[ib]=0.00000;

          for (int j=1; j<=nBinY; j++){
                  tmpRow     = j*nBinXeff;
                  for (int i=1; i<=nBinX; i++){
                          tmpBin       = tmpRow+i;
                          buff[tmpBin] = (in[tmpBin]>=mincut) ? in[tmpBin] : 0.0;
                          for (int kb=1; kb<smooth_nBinXHalfWdt1; kb++){
                                  if ((i-kb)>=1) {
                                          tmpBin1 = tmpBin - kb;
                                          buff[tmpBin] += ((in[tmpBin1]>=mincut) ? in[tmpBin1] : 0.0);
                                  }

                                  if ((i+kb)<=nBinX) {
                                          tmpBin2 = tmpBin + kb;
                                          buff[tmpBin] += ((in[tmpBin2]>=mincut) ? in[tmpBin2] : 0.0);
                                  }
                          }
                  }
          }

          for (int i=1; i<=nBinX; i++){
                  for (int j=1; j<=nBinY; j++){
                          tmpRow       = j*nBinXeff;
                          tmpBin       = tmpRow+i;
                          out[tmpBin]  = buff[tmpBin];
                          for (int kb=1; kb<smooth_nBinYHalfWdt1; kb++) {
                                  if ((j-kb)>=1) {
                                          tmpBin1 = tmpBin - kb*nBinXeff;
                                          out[tmpBin] += buff[tmpBin1];
                                  }
                                  if ((j+kb)<=nBinY) {
                                          tmpBin2 = tmpBin + kb*nBinXeff;
                                          out[tmpBin] += buff[tmpBin2];
                                  }
                          }
                          if (out[tmpBin]>0.0) {
                                  out[tmpBin]=floor(out[tmpBin]+0.5);
                                  nOutEntries+=(int)out[tmpBin];
                          }
                  }
          }

          delete [] buff;

          return nOutEntries;
  }

  int TrackReco::peakFinder(TH2F *inHist, float peakExpSigma, float peakThreshold, float* peakPositionX, float* peakPositionY, int nIter ) {

          int i, j, nfound;
          int nbinsx = inHist->GetNbinsX();
          int nbinsy = inHist->GetNbinsY();

          //double xmin = inHist->GetXaxis()->GetXmin();
          //double xmax = inHist->GetXaxis()->GetXmax();
          //double ymin = inHist->GetYaxis()->GetXmin();
          //double ymax = inHist->GetYaxis()->GetXmax();

          Float_t ** source = new Float_t *[nbinsx];
          for (i=0;i<nbinsx;i++) source[i]=new float[nbinsy];

          Float_t ** dest = new Float_t *[nbinsx];
          for (i=0;i<nbinsx;i++) dest[i]=new float[nbinsy];

          //TSpectrum2 *s = new TSpectrum2();
          TSpectrum2 finder;

          //Float_t *histArr = inHist->GetArray();
          for (i = 0; i < nbinsx; i++){
            for (j = 0; j < nbinsy; j++){
                       source[i][j] = inHist->GetBinContent(i + 1,j + 1);
                    }
          }

          nfound = finder.SearchHighRes(source, dest, nbinsx, nbinsy, peakExpSigma, peakThreshold, kTRUE, nIter, kFALSE, 0);

          //TH2F *dec = new TH2F("dec","deconvoluted hist",nbinsx,xmin,xmax,nbinsy,ymin,ymax);
          //unsigned int nEntries=0;
          //int tmpval;
          //for (i = 0; i < nbinsx; i++){
          //  for (j = 0; j < nbinsy; j++){
          //             tmpval=floor(dest[i][j]+0.5);
          //             dec->SetBinContent(i+1,j+1,tmpval);
          //             nEntries+=tmpval;
          //          }
          //}
          //dec->SetEntries(nEntries);
          //Searching->cd(2);
          //dec->Draw("col z");

          if (nfound>maxNpeaks) nfound=maxNpeaks;

          memcpy( peakPositionX, finder.GetPositionX(), nfound*sizeof(Float_t) );
          memcpy( peakPositionY, finder.GetPositionY(), nfound*sizeof(Float_t) );

          cout<<"Found candidate peaks "<<nfound<<endl;
          if (!nfound) return 0;
          if (_doDisplay) {
                  TPolyMarker * pm = (TPolyMarker*)inHist->GetListOfFunctions()->FindObject("TPolyMarker");
                  if (pm) {
                          inHist->GetListOfFunctions()->Remove(pm);
                          delete pm;
                  }
                  float *posX = new float[nfound];
                  float *posY = new float[nfound];
                  int binx, biny;
                  for(i=0;i<nfound;i++) {
                          binx = (int)floor(peakPositionX[i]+0.5)+1;
                          biny = (int)floor(peakPositionY[i]+0.5)+1;
                          posX[i]=inHist->GetXaxis()->GetBinCenter(binx);
                          posY[i]=inHist->GetYaxis()->GetBinCenter(biny);
                          cout<<"posx= "<<peakPositionX[i]<<" = "<<binx<<" posy= "<<peakPositionY[i]<<" = "<<biny<<" value= "<<inHist->GetBinContent(binx,biny)<<endl;//source[(int)(peakPositionX[i]+0.5)][(int)(peakPositionY[i]+0.5)]<<endl;
                  }
                  pm = new TPolyMarker(nfound, posX, posY);
                  pm->SetMarkerStyle(23);
                  pm->SetMarkerColor(kBlue-1);
                  pm->SetMarkerSize(1.3);
                  inHist->GetListOfFunctions()->Add(pm);
          }

          for (i=0;i<nbinsx;i++) {
                  delete [] dest[i];
                  delete [] source[i];
          }
          delete [] dest;
          delete [] source;

          return nfound;

  }



}  // end namespace mu2e

using mu2e::TrackReco;
DEFINE_ART_MODULE(TrackReco);
