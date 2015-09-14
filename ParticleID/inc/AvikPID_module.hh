///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __ParticleID_AvikPID__
#define __ParticleID_AvikPID__
// C++ includes.
#include <iostream>
#include <string>
#include <sstream>

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

//ROOTs
#include "TH1F.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"

#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "BTrk/TrkBase/TrkHotList.hh"
#include "BTrk/TrkBase/TrkHitOnTrk.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/PIDProduct.hh"
#include "RecoDataProducts/inc/PIDProductCollection.hh"
#include "KalmanTests/inc/Doublet.hh"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "KalmanTests/inc/TrkFitDirection.hh"

#include "ParticleID/inc/PIDUtilities.hh"
#include "RecoDataProducts/inc/AvikPIDProductCollection.hh"

namespace mu2e {
  class AvikPID : public art::EDProducer {

  private:

    enum { nbounds = 11 };
    static float _pathbounds[nbounds];

    AvikPIDProduct _pid;

    int    _debugLevel;
    int    _verbosity;
    int    _diagLevel;

    int    _processed_events;
    int    _evtid;

    string _trkPatRecDemModuleLabel;
    string _trkPatRecDmmModuleLabel;

    string _eleTemplates;
    string _muoTemplates;

    string _eleDedxTemplateFile;
    string _muoDedxTemplateFile;

    TrkParticle     _tpart;
    TrkFitDirection _fdir;
    std::string     _iname; // data instance name

    TH1D* _heletemp[nbounds];
    TH1D* _hmuotemp[nbounds];

    int   _templatesnbins ;
    float _templateslastbin ;
    float _templatesbinsize ;

    const KalRepPtrCollection* _listOfEleTracks;
    const KalRepPtrCollection* _listOfMuoTracks;

    int    _trkid;
    double _trkmom;

    double _logDedxProbEle;
    double _logDedxProbMuo;
//-----------------------------------------------------------------------------
// Vadim's global slopes
//-----------------------------------------------------------------------------
    double _drdsVadimEle;
    double _drdsVadimEleErr;

    double _drdsVadimMuo;
    double _drdsVadimMuoErr;
//-----------------------------------------------------------------------------
// Avik's sums
//-----------------------------------------------------------------------------
    double _sumAvikEle;
    double _sumAvikMuo;

    int    _nMatched;
    int    _nMatchedAll;
    double _sq2AvikEle;
    double _sq2AvikMuo;
//-----------------------------------------------------------------------------
// same-sign slopes
//-----------------------------------------------------------------------------
    int    _ele_nusedSs;
    int    _muo_nusedSs;
    double _drdsSsEle;
    double _drdsSsEleErr;
    double _drdsSsMuo;
    double _drdsSsMuoErr;
    double _logRatioSs;
//-----------------------------------------------------------------------------
// Avik doesn't calculate OS slopes! ...
//-----------------------------------------------------------------------------
    double _drdsOsEle;
    double _drdsOsEleErr;
    double _drdsOsMuo;
    double _drdsOsMuoErr;
    double _logRatioOs;

    double _ele_resSumOs;
    double _muo_resSumOs;
    int    _ele_nusedOs;
    int    _muo_nusedOs;

    TTree *       _pidtree;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
    explicit AvikPID(fhicl::ParameterSet const& pset);
    virtual ~AvikPID() { }
    void beginJob();
    void beginRun(art::Run &run);
    void beginSubRun(art::SubRun & lblock );
    virtual void produce(art::Event& event);
    void endJob();


    static  void myfcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t);
    static  int  findlowhist(float d);

    bool calculateVadimSlope(std::vector<double> vresd, 
			     std::vector<double> vflt,
			     std::vector<double> evresd,
			     std::vector<double> evflt,  
			     double              *slope, 
			     double              *eslope);

    double calculateDedxProb(std::vector<double>gaspaths, std::vector<double>edeps, TH1D** templates);

    void   doubletMaker(const KalRep* ele_Trk, const KalRep* muo_Trk);

    //    double calculateAvikSums();
    
    int    CalculateSlope(vector<double>& Fltlen, vector<double>& Resid, 
			  double&         Slope , double&         SlopeErr);

    int    AddHits(const Doublet* Multiplet, vector<double>& Fltlen, vector<double>& Resid);

    int    AddSsMultiplets(const vector<Doublet>* ListOfDoublets,
			   vector<double>&        Fltlen        , 
			   vector<double>&        Resid         );

    void   calculateSsSums(const vector<Doublet>* ele_LOD, const vector<Doublet>* muo_LOD);

    void   calculateOsSums(const vector<Doublet>* ele_LOD, const vector<Doublet>* muo_LOD);

    double weightedResidual(double R);

    // Save directory from beginJob so that we can go there in endJob. 
    //    TDirectory* _directory;


  };
}
#endif
