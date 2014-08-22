//
// TTracker time peak finder
//
// $Id: TimePeakFinder_module.cc,v 1.1 2014/08/22 16:10:41 tassiell Exp $
// $Author: tassiell $
// $Date: 2014/08/22 16:10:41 $
//
// Original author D. Brown and G. Tassielli
//

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
//#include "art/Framework/Services/Optional/TFileService.h"
// conditions
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
// Mu2e
#include "TrkPatRec/inc/TrkPatRecUtils.hh"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/TrackerHitTimeCluster.hh"
// root
#include "TFile.h"
#include "TMVA/Reader.h"
// C++
#include <memory>
using namespace std; 


namespace mu2e {

  class TimePeakFinder : public art::EDProducer {
  public:

    explicit TimePeakFinder(fhicl::ParameterSet const& pset);
    virtual ~TimePeakFinder();

    virtual void beginJob();
    void endJob();

    // This is called for each event.
    void produce(art::Event & e);

  private:

    // Start: run time parameters

    int _debug;
    int _printfreq;

    unsigned _iev;

    // event object labels
    std::string _shLabel;
    std::string _shpLabel;
    std::string _shfLabel;
    StrawHitFlag _tsel, _hsel;
    StrawHitFlag _tbkg, _hbkg;
    double _maxdt;
    // time spectrum parameters
    bool _findtpeak;
    unsigned _maxnpeak;
    unsigned _minnhits;
    bool _cleanpeaks;
    double _minpeakmva, _maxpeakdt, _maxpeakdphi;
    std::string _PMVAType; // type of MVA
    std::string _PMVAWeights; // file of MVA weights
    double _tmin;
    double _tmax;
    double _tbin;
    unsigned _nbins;
    double _ymin;
    double _1dthresh,_tssigma;

    // cache of event objects
    const StrawHitCollection* _shcol;
    const StrawHitFlagCollection* _shfcol;
    StrawHitFlagCollection* _flags;
    const StrawHitPositionCollection* _shpcol;

    // cache of time peaks
    std::vector<TrkTimePeak> _tpeaks;
    //string _iname; // data instance name

    void initializeReaders();
    TMVA::Reader *_peakMVA; // MVA for peak cleaning
    TimePeakMVA _pmva; // input variables to TMVA for peak cleaning

    art::Handle<mu2e::StrawHitCollection> _strawhitsH;

    bool findData(const art::Event& e);

  };

  TimePeakFinder::~TimePeakFinder() {
    if (_peakMVA!=NULL) { delete _peakMVA; }
  }
  
  TimePeakFinder::TimePeakFinder(fhicl::ParameterSet const& pset) :
                    _debug(pset.get<int>("debugLevel",0)),
                    _printfreq(pset.get<int>("printFrequency",101)),
                    _shLabel(pset.get<std::string>("StrawHitCollectionLabel","makeSH")),
                    _shpLabel(pset.get<std::string>("StrawHitPositionCollectionLabel","MakeStereoHits")),
                    _shfLabel(pset.get<std::string>("StrawHitFlagCollectionLabel","FlagBkgHits")),
                    _tsel(pset.get<std::vector<std::string> >("TimeSelectionBits")),
                    _hsel(pset.get<std::vector<std::string> >("HelixFitSelectionBits")),
                    _tbkg(pset.get<vector<string> >("TimeBackgroundBits",vector<string>{"DeltaRay","Isolated"})),
                    _hbkg(pset.get<vector<string> >("HelixFitBackgroundBits",vector<string>{"DeltaRay","Isolated"})),
                    _maxdt(pset.get<double>("DtMax",30.0)),
                    _findtpeak(pset.get<bool>("FindTimePeaks",true)),
                    _maxnpeak(pset.get<unsigned>("MaxNPeaks",50)),
                    _minnhits(pset.get<unsigned>("MinNHits",0)),
                    _cleanpeaks(pset.get<bool>("CleanTimePeaks",true)),
                    _minpeakmva(pset.get<double>("MinTimePeakMVA",0.5)),
                    _maxpeakdt(pset.get<double>("MaxTimePeakDeltat",25.0)),
                    _maxpeakdphi(pset.get<double>("MaxTimePeakDeltaPhi",1.0)),
                    _PMVAType(pset.get<std::string>("TimePeakMVAType","MLP method")),
                    _tmin(pset.get<double>("tmin",500.0)),
                    _tmax(pset.get<double>("tmax",1700.0)),
                    _tbin(pset.get<double>("tbin",20.0)),
                    _ymin(pset.get<double>("ymin",8.0)),
                    _1dthresh(pset.get<double>("OneDPeakThreshold",5.0)),
                    _peakMVA(NULL)
  {
          // tag the data product instance by the direction and particle type found by this fitter
          produces<StrawHitFlagCollection>();
          // set # bins for time spectrum plot
          _nbins = (unsigned)rint((_tmax-_tmin)/_tbin);
          // location-independent files
          ConfigFileLookupPolicy configFile;
          std::string weights = pset.get<std::string>("PeakMVAWeights","TrkPatRec/test/TimePeak.weights.xml");
          _PMVAWeights = configFile(weights);

          // Tell the framework what we make.
          produces<TrackerHitTimeClusterCollection>();

  }

  void TimePeakFinder::beginJob(){

    std::cout<<"Time peak finder jos!"<<std::endl;
//    initializeReaders(_peakMVA, _pmva, _PMVAType, _PMVAWeights);
    initializeReaders();

    // Get access to the TFile service and save current directory for later use.
//    art::ServiceHandle<art::TFileService> tfs;
    _iev = 0;
  }

  void TimePeakFinder::produce(art::Event & event ) {

          _iev=event.id().event();
          if((_iev%_printfreq)==0)cout<<"TimePeakFinder: event="<<_iev<<endl;
          // find the data
          if(!findData(event)){
            cout << "TimePeakFinder: No straw hits found, event="<<_iev << endl;
            return;
          }
          // copy in the existing flags
          _flags = new StrawHitFlagCollection(*_shfcol);
          std::unique_ptr<StrawHitFlagCollection> flags(_flags );

          // find the time peaks in the time spectrum of selected hits.  Otherwise, take all
          // selected hits as a peak
          _tpeaks.clear();
          if(_findtpeak){
            findTimePeaks( _shcol, _shpcol, _flags,
                           _tsel, _hsel, _tbkg, _hbkg,
                           _tpeaks, _maxnpeak,
                           _nbins, _tmin, _tmax,
                           _1dthresh, _ymin, _maxdt, _minnhits,
                           _cleanpeaks, _pmva, _peakMVA, _PMVAType, 
                           _minpeakmva, _maxpeakdt, _maxpeakdphi );

          } else {
            createTimePeak( _shcol, _flags,
                            _tsel, _hsel, _tbkg, _hbkg,
                            _tpeaks, _minnhits );

          }

          std::unique_ptr<TrackerHitTimeClusterCollection> thcc(new TrackerHitTimeClusterCollection);

          for (std::vector<TrkTimePeak>::iterator itp=_tpeaks.begin(); itp!=_tpeaks.end(); ++itp) {
                  thcc->push_back(TrackerHitTimeCluster());
                  thcc->back()._meanTime=itp->_tpeak;
                  thcc->back()._peakmax=itp->_peakmax;

                  thcc->back()._minHitTime = thcc->back()._maxHitTime = _shcol->at( itp->_trkptrs.begin()->_index ).time();
                  double sumOfSqr(0.0), sum(0.0);
                  for (std::vector<hitIndex>::iterator hittpit=itp->_trkptrs.begin(); hittpit!=itp->_trkptrs.end(); ++hittpit) {
                          sum+=_shcol->at(hittpit->_index).time();
                          sumOfSqr+=_shcol->at(hittpit->_index).time()*_shcol->at(hittpit->_index).time();
                          thcc->back()._selectedTrackerHits.push_back( art::Ptr<StrawHit> ( _strawhitsH, hittpit->_index ) );
                          if ( _shcol->at(hittpit->_index).time() < thcc->back()._minHitTime ) { thcc->back()._minHitTime=_shcol->at(hittpit->_index).time(); }
                          if ( _shcol->at(hittpit->_index).time() > thcc->back()._maxHitTime ) { thcc->back()._maxHitTime=_shcol->at(hittpit->_index).time(); }
                          flags->at(hittpit->_index).merge(StrawHitFlag::timesel);
                  }

                  if (itp->_trkptrs.size()>1) {
                          thcc->back()._sigma=sumOfSqr-sum*sum/((double)itp->_trkptrs.size());
                          if (thcc->back()._sigma>0.0) {
                                  thcc->back()._sigma=sqrt(thcc->back()._sigma/((double)(itp->_trkptrs.size()-1.0)));
                                  thcc->back()._nominalWidth=3.0*thcc->back()._sigma;
                          } else {
                                  thcc->back()._sigma=0.0;
                                  thcc->back()._nominalWidth = thcc->back()._maxHitTime-thcc->back()._minHitTime;
                          }
                  } else {
                          thcc->back()._nominalWidth = thcc->back()._maxHitTime-thcc->back()._minHitTime;
                  }

          }

          int iPeak(0);
          if (_debug>0) {
                  std::cout<<"n Peaks "<<thcc->size()<<std::endl;
                  for (TrackerHitTimeClusterCollection::iterator tpeak_it=thcc->begin(); tpeak_it!=thcc->end(); ++tpeak_it ) {
                          std::cout<<"iPeak "<<std::endl;
                          std::cout<<*tpeak_it<<std::endl;
                          ++iPeak;
                  }
          }

          event.put(std::move(thcc));

          event.put(std::move(flags));

  } // end produce

  void TimePeakFinder::endJob(){
  }

  // find the input data objects
   bool TimePeakFinder::findData(const art::Event& evt){
     _shcol = 0;
     _shfcol = 0;
     _shpcol = 0;

     if(evt.getByLabel(_shLabel,_strawhitsH))
       _shcol = _strawhitsH.product();
     art::Handle<mu2e::StrawHitPositionCollection> shposH;
     if(evt.getByLabel(_shpLabel,shposH))
       _shpcol = shposH.product();
     art::Handle<mu2e::StrawHitFlagCollection> shflagH;
     if(evt.getByLabel(_shfLabel,shflagH))
       _shfcol = shflagH.product();
 // don't require stereo hits: they are only used for diagnostics
     return _shcol != 0 && _shfcol != 0 && _shpcol != 0;
   }

   void TimePeakFinder::initializeReaders() {
     _peakMVA = new TMVA::Reader();
     _peakMVA->AddVariable("_dt",&_pmva._dt);
     _peakMVA->AddVariable("_dphi",&_pmva._dphi);
     _peakMVA->AddVariable("_rho",&_pmva._rho);
     _peakMVA->BookMVA(_PMVAType,_PMVAWeights);
   }

}  // end namespace mu2e

using mu2e::TimePeakFinder;
DEFINE_ART_MODULE(TimePeakFinder);
