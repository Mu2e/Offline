//
// TTracker commons for Kalman Fit
//
// $Id: TrkRecBase.hh,v 1.1 2014/08/22 16:10:41 tassiell Exp $
// $Author: tassiell $
// $Date: 2014/08/22 16:10:41 $
//
//
// Original author D. Brown and G. Tassielli
//
#ifndef TrkRecBase_HH
#define TrkRecBase_HH

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "GeometryService/inc/DetectorSystem.hh"
#include "art/Framework/Services/Optional/TFileService.h"
// conditions
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StereoHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/KalRepPayload.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
// Mu2e
#include "KalmanTests/inc/KalFit.hh"
#include "KalmanTests/inc/KalFitMC.hh"
#include "KalmanTests/inc/KalRepCollection.hh"
#include "KalmanTests/inc/KalFitResult.hh"
#include "TrkPatRec/inc/TrkHitFilter.hh"
#include "TrkPatRec/inc/StrawHitInfo.hh"
#include "TrkPatRec/inc/HelixFit.hh"
#include "TrkPatRec/inc/TrkPatRecUtils.hh"
#include "TrkPatRec/inc/PayloadSaver.hh"
//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// root 
#include "TFile.h"
#include "TTree.h"
// C++
#include <string>
#include <memory>
#include <functional>
#include <vector>
//using namespace std;

namespace mu2e 
{
  class TrkRecBase
  {
    public:
      enum fitType {helixFit=0,seedFit,kalFit};
      TrkRecBase(fhicl::ParameterSet const&);
      ~TrkRecBase();
      void bgnJob();
    protected:
      unsigned _iev;
      // configuration parameters
      int _diag,_debug;
      int _printfreq;
      bool _addhits; 
      // event object labels
      string _shLabel;
      string _shpLabel;
      string _stLabel;
      string _shfLabel;
      StrawHitFlag _tsel, _hsel, _addsel;
      StrawHitFlag _tbkg, _hbkg, _addbkg;
      double _maxdt, _maxdtmiss;
      // time spectrum parameters
      bool _findtpeak;
      unsigned _maxnpeak;
      unsigned _minnhits;
      bool _cleanpeaks;
      double _minpeakmva, _maxpeakdt, _maxpeakdphi;
      std::string _PMVAType; // type of MVA
      std::string _PMVAWeights; // file of MVA weights
      double _maxphirange;
      double _tmin;
      double _tmax;
      double _tbin;
      unsigned _nbins;
      double _ymin;
      double _1dthresh,_tssigma;
      // outlier cuts
      double _maxseeddoca,_maxhelixdoca,_maxadddoca, _maxaddchi;
      TrkParticle _tpart; // particle type being searched for
      TrkFitDirection _fdir;  // fit direction in search
      // cache of event objects
      const StrawHitCollection* _shcol;
      const StrawHitFlagCollection* _shfcol;
      StrawHitFlagCollection* _flags;
      const StrawHitPositionCollection* _shpcol;
      const StereoHitCollection* _stcol;
      // Kalman fitters.  Seed fit has a special configuration
      KalFit _seedfit, _kfit;
      // robust helix fitter
      HelixFit _hfit;
      // cache of time peaks
      vector<TrkTimePeak> _tpeaks;
      string _iname; // data instance name
      //
      PayloadSaver _payloadSaver;
      // helper functions
      bool findData(const art::Event& e);
//      void filterOutliers(TrkDef& mytrk,Trajectory const& traj,double maxdoca,vector<TrkHitFilter>& thfvec);
      void findMissingHits(KalFitResult& kalfit, vector<hitIndex>& indices);
      void createDiagnostics();
      void fillStrawDiag();
      void fillTimeDiag();
      void fillPeakDiag(size_t ip, TrkTimePeak const& tp);
      void fillFitDiag(int ipeak, HelixFitResult const& helixfit,
	  KalFitResult const& seedfit,KalFitResult const& kalfit);
      void fillStrawHitInfo(size_t ish, StrawHitInfo& shinfo) const;
      void initializeReaders();

      TMVA::Reader *_peakMVA; // MVA for peak cleaning
      TimePeakMVA _pmva; // input variables to TMVA for peak cleaning

      // MC tools
      KalFitMC _kfitmc;
      // strawhit tuple variables
      TTree *_shdiag;
      Int_t _eventid;
      threevec _shp;
      Float_t _edep;
      Float_t _time, _deltat, _rho;
      Int_t _nmcsteps;
      Int_t _mcnunique,_mcnmax;
      Int_t _mcpdg,_mcgen,_mcproc;
      Int_t _mcppdg,_mcpproc;
      Int_t _mcgid, _mcgpdg;
      Float_t _mcge, _mcgt;
      threevec _mcshp, _mcop, _mcpop, _mcgpos;
      Float_t _mcoe, _mcpoe, _mcom, _mcpom;
      Float_t _mcshlen,_mcshd;
      Float_t _mcedep,_mcemax;
      Float_t _pdist,_pperp,_pmom;
      Float_t _mctime, _mcptime;
      Int_t _esel,_rsel, _timesel,  _delta, _stereo, _isolated;
      Int_t _device, _sector, _layer, _straw;
      Float_t _shpres, _shrres, _shchisq, _shdt, _shdist;
      Float_t _shmct0, _shmcmom, _shmctd;
      Bool_t _xtalk;
      // time peak diag variables
      TTree* _tpdiag;
      Int_t _tpeventid, _peakid, _pmax, _nphits, _ncphits, _nchits;
      Float_t _ptime, _pdtimemax, _ctime, _cdtimemax;;
      Float_t _pphi, _cphi, _cphirange, _pdphimax, _cdphimax;
      vector<TimePeakHitInfo> _tphinfo;

      // fit tuple variables
      Int_t _nadd,_ipeak;
      Float_t _hcx, _hcy, _hr, _hdfdz, _hfz0;
      Float_t _mccx, _mccy, _mcr, _mcdfdz, _mcfz0;
      Int_t _helixfail,_seedfail,_kalfail;
      helixpar _hpar,_spar;
      helixpar _hparerr,_sparerr;
      Int_t _snhits, _snactive, _sniter, _sndof, _snweediter;
      Float_t _schisq, _st0;
      Int_t _nchit;
      Int_t _npeak, _nmc;
      Float_t _peakmax, _tpeak;
      // hit filtering tuple variables
      vector<TrkHitFilter> _sfilt, _hfilt;
      // flow diagnostic
      TH1F* _cutflow, *_ccutflow;
      int _icepeak;
  };
}

#endif /* TrkRecBase_HH */

