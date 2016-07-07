//
// TTracker time peak finder
//
// $Id: TimePeakDiag_module.cc,v 1.3 2014/08/25 12:08:29 tassiell Exp $
// $Author: tassiell $
// $Date: 2014/08/25 12:08:29 $
//
// Original author D. Brown and G. Tassielli
//

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
// Mu2e
#include "TrkPatRec/inc/TrkPatRec.hh"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/TrackSeed.hh"
#include "RecoDataProducts/inc/TrackSeedCollection.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
// root
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TGMsgBox.h"
#include "TTree.h"
#include "TMarker.h"
// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/mean.hpp>
// C++
#include <memory>
using namespace std; 
using namespace boost::accumulators;

namespace mu2e {

  class TimePeakDiag : public art::EDAnalyzer {
  public:

    explicit TimePeakDiag(fhicl::ParameterSet const& pset);
    virtual ~TimePeakDiag();

    virtual void beginJob();
    void endJob();

    // This is called for each event.
    void analyze(art::Event & e);

  private:

    // Start: run time parameters

    int           _diag;
    bool	  _mcdiag;

    unsigned      _iev;

    // event object labels
    std::string   _shLabel;
    std::string   _shpLabel;
    std::string   _shfLabel;
    std::string   _mcdigislabel;

    // time spectrum parameters
    double        _tmin;
    double        _tmax;
    double        _tbin;
    unsigned      _nbins;
    // cache of event objects
    const StrawHitCollection*             _shcol;
    const StrawHitFlagCollection*         _shfcol;
    const StrawHitPositionCollection*     _shpcol;
    const StrawDigiMCCollection*          _mcdigis;

    void fillTimeDiag     ();
    void fillPeakDiag     (size_t ip, TrkTimePeak const& tp);

// time peak diag variables
    TTree*                                _tpdiag;
    Int_t                                 _tpeventid, _peakid, _pmax, _nphits, _ncphits, _nchits;
    Float_t                               _ptime, _pdtimemax, _ctime, _cdtimemax;
    Float_t                               _pphi, _cphi, _cphirange, _pdphimax, _cdphimax;
    vector<TimePeakHitInfo>               _tphinfo;
  };

  TimePeakDiag::~TimePeakDiag() {
  }
  
  TimePeakDiag::TimePeakDiag(fhicl::ParameterSet const& pset) :
    _diag              (pset.get<int>("diagLevel",0)),
    _mcdiag            (pset.get<int>("MCdiag",true)),
    _shLabel           (pset.get<std::string>("StrawHitCollectionLabel","makeSH")),
    _shpLabel          (pset.get<std::string>("StrawHitPositionCollectionLabel","MakeStereoHits")),
    _shfLabel          (pset.get<std::string>("StrawHitFlagCollectionLabel","FlagBkgHits")),
    _mcdigislabel      (pset.get<string>("StrawDigiMCLabel")),
    _tmin              (pset.get<double>("tmin",500.0)),
    _tmax              (pset.get<double>("tmax",1700.0)),
    _tbin              (pset.get<double>("tbin",20.0)),
  {
    // set # bins for time spectrum plot
    _nbins = (unsigned)rint((_tmax-_tmin)/_tbin);
  }

  void TimePeakDiag::beginJob(){
    if(_diag > 0)createDiagnostics();
    _iev = 0;
  }

  void TimePeakDiag::analyze(art::Event & event ) {

    _iev=event.id().event();
    // find the data
    if(!findData(event)){
      throw cet::exception("RECO")<<"mu2e::TimePeakDiag: data missing or incomplete"<< endl;
      return;
    }
    fillTimeDiag();

    for (std::vector<TrkTimePeak>::iterator itp=_tpeaks.begin(); itp!=_tpeaks.end(); ++itp) {
      double    meanTime(0), sigma(0), minHitTime(0), maxHitTime(0), nominalWidth(0);
      double    sumOfSqr(0), sum(0);
      double    t0(0), errt0(0);

      TrackSeed tmpSeed;
      
      meanTime   = itp->_tpeak;
      minHitTime = _shcol->at( itp->_trkptrs.begin()->_index ).time();
      maxHitTime = _shcol->at( itp->_trkptrs.begin()->_index ).time();

      for (std::vector<hitIndex>::iterator hittpit=itp->_trkptrs.begin(); hittpit!=itp->_trkptrs.end(); ++hittpit) {
	double htime = _shcol->at(hittpit->_index).time();

	sum         += htime;
	sumOfSqr    += htime*htime;
	
	tmpSeed._timeCluster._strawHitIdxs.push_back( mu2e::hitIndex( hittpit->_index, hittpit->_ambig) );

	if ( htime < minHitTime ) { minHitTime = htime; }
	if ( htime > maxHitTime ) { maxHitTime = htime; }
	
	flags->at(hittpit->_index).merge(StrawHitFlag::timesel);
      }

    }

  } // end analyze

  void TimePeakDiag::endJob(){
  }

  // find the input data objects
   bool TimePeakDiag::findData(const art::Event& evt){
     _shcol  = 0;
     _shfcol = 0;
     _shpcol = 0;
    art::Handle<mu2e::StrawHitCollection> strawhitsH;
     if(evt.getByLabel(_shLabel,strawhitsH))
       _shcol = strawhitsH.product();
     art::Handle<mu2e::StrawHitPositionCollection> shposH;
     if(evt.getByLabel(_shpLabel,shposH))
       _shpcol = shposH.product();
     art::Handle<mu2e::StrawHitFlagCollection> shflagH;
     if(evt.getByLabel(_shfLabel,shflagH))
       _shfcol = shflagH.product();

    if(_diag > 0){
      art::Handle<StrawDigiMCCollection> mcdigisHandle;
      if(evt.getByLabel(_mcdigislabel,"StrawHitMC",mcdigisHandle))
	_mcdigis = mcdigisHandle.product();
    }

     return _shcol != 0 && _shfcol != 0 && _shpcol != 0;
   }

  void TimePeakDiag::fillPeakDiag(size_t ip,TrkTimePeak const& tp) {
    _tpeventid = _iev;
    _peakid = ip;
    _pmax = tp._peakmax;
    _nphits = tp._trkptrs.size();
    _ncphits = 0;
    _pdtimemax = 0.0;
    _pdphimax = 0.0;
    _cdtimemax = 0.0;
    _cdphimax = 0.0;
    _tphinfo.clear();
    accumulator_set<double, stats<tag::mean > > facc;
    accumulator_set<double, stats<tag::mean > > tacc;
    for(size_t iph=0;iph<tp._trkptrs.size();++iph){
      unsigned ish = tp._trkptrs[iph]._index;
      StrawHit const& sh = _shcol->at(ish);
      StrawHitPosition const& shp = _shpcol->at(ish);
      double time = sh.time();
      tacc(time);
      CLHEP::Hep3Vector const& pos = shp.pos();
      double phi = pos.phi();
      if(extract_result<tag::count>(facc) > 0){
	double dphi = phi - extract_result<tag::mean>(facc);
	if(dphi > M_PI){
	  phi -= 2*M_PI;
	} else if(dphi < -M_PI){
	  phi += 2*M_PI;
	}
      }
      facc(phi);
    }
    _pphi =extract_result<tag::mean>(facc);
    _ptime =extract_result<tag::mean>(tacc);
    
    for(size_t iph=0;iph<tp._trkptrs.size();++iph){
      unsigned ish = tp._trkptrs[iph]._index;
      StrawHitPosition const& shp = _shpcol->at(ish);
      CLHEP::Hep3Vector const& pos = shp.pos();
      double phi = pos.phi();
      double dphi = phi - _pphi;
      if(dphi > M_PI)
	dphi -= 2*M_PI;
      else if(dphi < -M_PI)
	dphi += 2*M_PI;
      if(fabs(dphi) > _pdphimax)_pdphimax = fabs(dphi);

      double dt = _shcol->at(ish).time() - _ptime;
      if(fabs(dt) > _pdtimemax)_pdtimemax=fabs(dt);

      StrawDigiMC const& mcdigi = _mcdigis->at(ish);
      // use TDC channel 0 to define the MC match
      StrawDigi::TDCChannel itdc = StrawDigi::zero;
      if(!mcdigi.hasTDC(StrawDigi::one)) itdc = StrawDigi::one;
      art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
      art::Ptr<SimParticle> const& spp = spmcp->simParticle();
      int gid(-1);
      if(spp->genParticle().isNonnull())
	gid = spp->genParticle()->generatorId().id();

      bool conversion = (spp->pdgId() == 11 && gid == 2 && spmcp->momentum().mag()>90.0);

      if(conversion){
	_ncphits++;
	if(fabs(dt) > _cdtimemax)_cdtimemax = fabs(dt);
      }
      if(conversion && dphi > _cdphimax)_cdphimax = dphi;

      _pmva._dt = dt;
      _pmva._dphi = dphi;
      _pmva._rho = pos.perp();
      double mvaout = _peakMVA.evalMVA(_pmva._pars);

      TimePeakHitInfo tph;
      tph._dt = dt;
      tph._dphi = dphi;
      tph._rho = pos.perp();
      tph._mva = mvaout;
      tph._mcpdg = spp->pdgId();
      tph._mcgen = gid;
      tph._mcproc = spp->creationCode();
      _tphinfo.push_back(tph);
    }
    _tpdiag->Fill();
//    if(_ncphits > 0.5*_nchits)_icepeak = ip;
  }

  void TimePeakDiag::fillTimeDiag() {
    art::ServiceHandle<art::TFileService> tfs;
    TH1F *ctsp, *rtsp, *ttsp, *ltsp, *tdtsp, *ptsp;

    char rsname[100];
    char csname[100];
    char tsname[100];
    char lsname[100];
    char tdsname[100];
    char tpsname[100];
    snprintf(rsname,100,"rawtspectrum%i",_iev);
    snprintf(csname,100,"convtspectrum%i",_iev);
    snprintf(tsname,100,"tighttspectrum%i",_iev);
    snprintf(lsname,100,"loosetspectrum%i",_iev);
    snprintf(tdsname,100,"tightnodeltatspectrum%i",_iev);
    snprintf(tpsname,100,"protontspectrum%i",_iev);
    ttsp = tfs->make<TH1F>(tsname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    ttsp->SetLineColor(kCyan);
    ltsp = tfs->make<TH1F>(lsname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    ltsp->SetLineColor(kGreen);
    rtsp = tfs->make<TH1F>(rsname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    rtsp->SetLineColor(kBlue);
    ctsp = tfs->make<TH1F>(csname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    ctsp->SetLineColor(kRed);
    ctsp->SetFillColor(kRed);
    ptsp = tfs->make<TH1F>(tpsname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    ptsp->SetLineColor(kBlack);
    ptsp->SetFillColor(kBlack);
    tdtsp = tfs->make<TH1F>(tdsname,"time spectrum;nsec",_nbins,_tmin,_tmax);
    tdtsp->SetLineColor(kOrange);

    unsigned nstrs = _shcol->size();
    for(unsigned istr=0; istr<nstrs;++istr){
      double time = _shcol->at(istr).time();
      bool conversion(false);
      bool proton(false);
      // summarize the MC truth for this strawhit
      if(_mcdigis != 0) {
	StrawDigiMC const& mcdigi = _mcdigis->at(istr);
	// use TDC channel 0 to define the MC match
	StrawDigi::TDCChannel itdc = StrawDigi::zero;
	if(!mcdigi.hasTDC(StrawDigi::one)) itdc = StrawDigi::one;
	art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
	art::Ptr<SimParticle> const& spp = spmcp->simParticle();
	int gid(-1);
	if(spp->genParticle().isNonnull())
	  gid = spp->genParticle()->generatorId().id();

	conversion = (spp->pdgId() == 11 && gid == 2 && spmcp->momentum().mag()>90.0);
	proton = spp->pdgId()==2212;
      }
      // fill plots
      rtsp->Fill(time);
      if(_flags->at(istr).hasAllProperties(_hsel)){
	ttsp->Fill(time);
      }
      if(_flags->at(istr).hasAllProperties(_hsel) && !_flags->at(istr).hasAnyProperty(_hbkg)){
	tdtsp->Fill(time);
      }
      if(_flags->at(istr).hasAllProperties(_hsel) && !_flags->at(istr).hasAnyProperty(_hbkg)){
	ltsp->Fill(time);
	if(proton)
	  ptsp->Fill(time);
      }
      if(conversion)
	ctsp->Fill(time);

    }
    // plot time peaks
    TList* flist = tdtsp->GetListOfFunctions();
    for(auto ipeak=_tpeaks.begin();ipeak!=_tpeaks.end();++ipeak){
      TMarker* smark = new TMarker(ipeak->_tpeak,ipeak->_peakmax,23);
      smark->SetMarkerColor(kRed);
      smark->SetMarkerSize(1.5);
      flist->Add(smark);
    }
  }

  void TimePeakDiag::createDiagnostics() {
    art::ServiceHandle<art::TFileService> tfs;
   // time peak diagnostics
    _tpdiag=tfs->make<TTree>("tpdiag","time peak diagnostics");
    _tpdiag->Branch("eventid",&_tpeventid,"eventid/I");
    _tpdiag->Branch("peakid",&_peakid,"peakid/I");
    _tpdiag->Branch("pmax",&_pmax,"pmax/I");
    _tpdiag->Branch("nphits",&_nphits,"nphits/I");
    _tpdiag->Branch("ncphits",&_ncphits,"ncphits/I");
    _tpdiag->Branch("nchits",&_nchits,"nchits/I");
    _tpdiag->Branch("ptime",&_ptime,"ptime/F");
    _tpdiag->Branch("ctime",&_ctime,"ctime/F");
    _tpdiag->Branch("pdtimemax",&_pdtimemax,"pdtimemax/F");
    _tpdiag->Branch("cdtimemax",&_cdtimemax,"cdtimemax/F");
    _tpdiag->Branch("pphi",&_pphi,"pphi/F");
    _tpdiag->Branch("cphi",&_cphi,"cphi/F");
    _tpdiag->Branch("cphirange",&_cphirange,"cphirange/F");
    _tpdiag->Branch("pdphimax",&_pdphimax,"pdphimax/F");
    _tpdiag->Branch("cdphimax",&_cdphimax,"cdphimax/F");
    _tpdiag->Branch("tphinfo",&_tphinfo);
 }



}  // end namespace mu2e

using mu2e::TimePeakDiag;
DEFINE_ART_MODULE(TimePeakDiag);
