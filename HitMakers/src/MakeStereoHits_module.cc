//
// A module to create simple stereo hits out of StrawHits.  This can work
// with either tracker.  StrawHit selection is done by flagging in an upstream module
//
// $Id: MakeStereoHits_module.cc,v 1.19 2014/05/31 03:23:52 mwang Exp $
// $Author: mwang $
// $Date: 2014/05/31 03:23:52 $
// 
//  Original Author: David Brown, LBNL
//  

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StereoHitCollection.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "KalmanTests/inc/KalFitMC.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "HitMakers/inc/MVATools.hh"

// art includes.
#include "art/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"

// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// root
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TMarker.h"
#include "TList.h"
#include "TLegend.h"
#include "TTree.h"
// C++ includes.
#include <iostream>
#include <float.h>

using namespace std;

namespace mu2e {
  // used to sort sectors by dphi
  struct SecPhi {
    int isec;
    double dphi;
  };

  // comparison functor for sorting by dphi
  struct sortSec : public std::binary_function<SecPhi,SecPhi,bool> {
    bool operator()(SecPhi const& v1, SecPhi const& v2) { return (v1.dphi < v2.dphi); }
  };

  class MakeStereoHits : public art::EDProducer {

  public:
    explicit MakeStereoHits(fhicl::ParameterSet const& pset);
    // Accept compiler written d'tor.

    void produce( art::Event& e);
    void beginJob();

  private:

    // Diagnostics level.
    int _diagLevel;
    // Debug level.
    int _debugLevel;
    // Name of the StrawHit collection
    string _shLabel, _shfLabel;
    string _mcdigislabel;
  // Parameters
    StrawHitFlag _shsel; // flag selection
    double _maxDt; // maximum time separation between hits
    double _maxDE; // maximum deposited energy deference: this excludes inconsistent hits
    double _maxDZ; // maximum longitudinal separation
    double _maxDPerp; // maximum transverse separation
    double _minDdot; // minimum dot product of straw directions
    double _minDL; // minimum distance from end of active straw;
    double _maxChi; // maximum # of TimeDivision sigmas past the active edge of a wire to allow making stereo hits
    double _maxChisq; // maximum # of TimeDivision consistency chisquared to allow making stereo hits
    double _minMVA; // minimum MVA output
    bool _writepairs; // write out the stereo pairs
    StrawHitFlag _stmask; // selection 
    MVATools _mvatool;
    std::string _MVAType; // type of MVA
    std::string _MVAWeights; // file of MVA weights
    vector<double> _vmva; // input variables to TMVA for stereo selection
    // for optimized Stereo Hit finding
    size_t _nsta;
    size_t _nsec;
    vector <vector<int> >_imap;             // remap sector indices to new scheme
    vector <vector<int> >_dosec;            // list of overlapping sectors in "road" to search in
    void genMap(const Tracker& tracker);    // function to generate imap from tracker geometry
    double deltaPhi(double phi1,double phi2) const;
    // diagnostics
    TH1F* _nhits;
    TH1F* _deltat;
    TH1F* _deltaE;
    TH1F* _deltaz;
    TH1F* _fsep;
    TH1F* _dL;
    TH1F* _mva;
    vector<TH2F*> _stations;
    const StrawDigiMCCollection *_mcdigis;
    TTree *_spdiag, *_sdiag;
    Float_t _shphi, _stphi, _mcshphi;
    Float_t _shrho, _strho, _mcshrho;
    Float_t _de, _dt, _dist, _dperp, _dz, _rho, _dl1, _dl2, _dc1, _dc2, _chi2, _mvaout, _ddot;
    Float_t _schi2, _smvaout, _sddot, _sdist, _sdz;
    Float_t _mcdist;
    Int_t _stereo, _fs, _sfs, _mcr, _mcrel, _mcpdg, _mcgen, _mcproc;
 };

  MakeStereoHits::MakeStereoHits(fhicl::ParameterSet const& pset) :
    // Parameters
    _diagLevel(pset.get<int>("diagLevel",0)),
    _debugLevel(pset.get<int>("debugLevel",0)),
    _shLabel(pset.get<string>("StrawHitCollectionLabel","makeSH")),
    _shfLabel(pset.get<std::string>("StrawHitFlagCollectionLabel","FSHPreStereo")),
    _mcdigislabel(pset.get<string>("StrawHitMCLabel","makeSH")),
    _shsel(pset.get<vector<string> >("StrawHitSelectionBits",vector<string>{"EnergySelection","TimeSelection"} )),
    _maxDt(pset.get<double>("maxDt",40.0)), // nsec
    _maxDE(pset.get<double>("maxDE",0.99)), // dimensionless, from 0 to 1
    _maxDZ(pset.get<double>("maxDZ",40.)), // mm, maximum longitudinal distance between straws
    _maxDPerp(pset.get<double>("maxDPerp",500.)), // mm, maximum perpendicular distance between time-division points
    _minDdot(pset.get<double>("minDdot",0.6)), // minimum angle between straws
    _minDL(pset.get<double>("minDL",-20.0)), // extent along straw
    _maxChisq(pset.get<double>("maxChisquared",80.0)), // position matching
    _minMVA(pset.get<double>("minMVA",0.7)), // MVA cut
    _writepairs(pset.get<bool>("WriteStereoPairs",false)),
    _mvatool(pset.get<fhicl::ParameterSet>("MVATool",fhicl::ParameterSet())),
    _MVAType(pset.get<std::string>("MVAType","MLP method")),
    _nhits(0),_deltat(0),_deltaE(0),_deltaz(0),_fsep(0),_dL(0),_mva(0),
    _mcdigis(0),_sdiag(0)
  {
    _maxChi = sqrt(_maxChisq);
    // location-independent files
    ConfigFileLookupPolicy configFile;
    std::string weights = pset.get<std::string>("MVAWeights","HitMakers/test/StereoMVA.weights.xml");
    _MVAWeights = configFile(weights);
    // Tell the framework what we make.
    if(_writepairs)produces<StereoHitCollection>();
    produces<StrawHitPositionCollection>();
  }

  void MakeStereoHits::beginJob(){
  // initialize MVA
    _mvatool.initMVA();
    _mvatool.showMVA();
    _vmva.resize(4);
    // create diagnostics if requested
    if(_diagLevel > 0){
      art::ServiceHandle<art::TFileService> tfs;
      _nhits = tfs->make<TH1F>("nhits","NHits",500,0,5000);
      _deltat = tfs->make<TH1F>("deltat","#Delta t;ns",100,-200.0,200.0);
      _deltaE = tfs->make<TH1F>("deltaE","#Delta E/#Sigma E;Ratio",100,-1.0,1.0);
      _deltaz = tfs->make<TH1F>("deltaz","#Delta d;mm",120,0.0,120.0);
      _fsep = tfs->make<TH1F>("sep","Face separation",6,-0.5,5.5);
      _dL = tfs->make<TH1F>("dL","Length Difference;mm",100,-200.0,100.0);
      _mva = tfs->make<TH1F>("mva","MVA output",100,-0.05,1.05);
      if( _diagLevel > 1){
	// detailed diagnostics
	_spdiag=tfs->make<TTree>("spdiag","stereo position diagnostics");
	_spdiag->Branch("shphi",&_shphi,"shphi/F");
	_spdiag->Branch("shrho",&_shrho,"shrho/F");
	_spdiag->Branch("stphi",&_stphi,"stphi/F");
	_spdiag->Branch("strho",&_strho,"strho/F");
	_spdiag->Branch("stereo",&_stereo,"stereo/I");
	_spdiag->Branch("chisq",&_schi2,"chisq/F");
	_spdiag->Branch("mvaout",&_smvaout,"mvaout/F");
	_spdiag->Branch("dist",&_sdist,"dist/F");
	_spdiag->Branch("dz",&_sdz,"dz/F");
	_spdiag->Branch("ddot",&_sddot,"ddot/F");
	_spdiag->Branch("mcrel",&_mcr,"mcr/I");
	_spdiag->Branch("mcpdg",&_mcpdg,"mcpdg/I");
	_spdiag->Branch("mcgen",&_mcgen,"mcgen/I");
	_spdiag->Branch("mcproc",&_mcproc,"mcproc/I");
	_spdiag->Branch("mcshphi",&_mcshphi,"mcshphi/F");
	_spdiag->Branch("mcshrho",&_mcshrho,"mcshrho/F");
	_spdiag->Branch("fs",&_sfs,"fs/I");
	if(_diagLevel > 2){
	  _sdiag=tfs->make<TTree>("sdiag","stereo diagnostics");
	  _sdiag->Branch("fs",&_fs,"fs/I");
	  _sdiag->Branch("de",&_de,"de/F");
	  _sdiag->Branch("dt",&_dt,"dt/F");
	  _sdiag->Branch("dz",&_dz,"dz/F");
	  _sdiag->Branch("rho",&_rho,"rho/F");
	  _sdiag->Branch("dist",&_dist,"dist/F");
	  _sdiag->Branch("dperp",&_dperp,"dperp/F");
	  _sdiag->Branch("dl1",&_dl1,"dl1/F");
	  _sdiag->Branch("dl2",&_dl2,"dl2/F");
	  _sdiag->Branch("dc1",&_dc1,"dc1/F");
	  _sdiag->Branch("dc2",&_dc2,"dc2/F");
	  _sdiag->Branch("chi2",&_chi2,"chi2/F");
	  _sdiag->Branch("mvaout",&_mvaout,"mvaout/F");
	  _sdiag->Branch("ddot",&_ddot,"ddot/F");
	  _sdiag->Branch("mcrel",&_mcrel,"mcrel/I");
	  _sdiag->Branch("mcdist",&_mcdist,"mcdist/F");
	}
      }
    }
  }

  void MakeStereoHits::produce(art::Event& event) {
    // Get a reference to one of the L or T trackers.
    // Throw exception if not successful.
    const Tracker& tracker = getTrackerOrThrow();

    // setup imap, dosec, and hdx
    static bool genmap(true);
    if(genmap){
      genmap = false;
      genMap(tracker);
    }

    static bool first(true);
    if(_diagLevel >0 && first){
      first = false;
      const TTracker& tt = dynamic_cast<const TTracker&>(tracker);
      art::ServiceHandle<art::TFileService> tfs;
      unsigned nsta = tt.nDevices()/2;
      for(unsigned ista=0;ista<nsta;++ista){
	char name[100];
	snprintf(name,100,"station%i",ista);
	_stations.push_back( tfs->make<TH2F>(name,name,100,-700,700,100,-700,700));
	_stations[ista]->SetStats(false);
	TList* flist = _stations[ista]->GetListOfFunctions();
	TLegend* sleg = new TLegend(0.1,0.6,0.3,0.9);
	flist->Add(sleg);
	for(int idev=0;idev<2;++idev){
	  const Device& dev = tt.getDevice(2*ista+idev);
	  const vector<Sector>& sectors = dev.getSectors();
	  for(size_t isec=0;isec<sectors.size();++isec){
	    int iface = isec%2;
	    const Sector& sec = sectors[isec];
	    CLHEP::Hep3Vector spos = sec.straw0MidPoint();
	    CLHEP::Hep3Vector sdir = sec.straw0Direction();
	    CLHEP::Hep3Vector end0 = spos - 100.0*sdir;
	    CLHEP::Hep3Vector end1 = spos + 100.0*sdir;
	    TLine* sline = new TLine(end0.x(),end0.y(),end1.x(),end1.y());
	    sline->SetLineColor(isec+1);
	    sline->SetLineStyle(2*idev+iface+1);
	    flist->Add(sline);
	    TMarker* smark = new TMarker(end0.x(),end0.y(),8);
	    smark->SetMarkerColor(isec+1);
	    smark->SetMarkerSize(2);
	    flist->Add(smark);
	    char label[80];
	    snprintf(label,80,"dev %i sec %i",idev,(int)isec);
	    sleg->AddEntry(sline,label,"l");
	  }
	}
      }
    }

    // Handle to the conditions service
    ConditionsHandle<TrackerCalibrations> tcal("ignored");

    art::Handle<mu2e::StrawHitCollection> strawhitsH; 
    const StrawHitCollection* strawhits(0);
    if(event.getByLabel(_shLabel,strawhitsH))
      strawhits = strawhitsH.product();
    if(strawhits == 0){
      throw cet::exception("RECO")<<"mu2e::MakeStereoHits: No StrawHit collection found for label " <<  _shLabel << endl;
    }
    const StrawHitFlagCollection* shfcol(0);
    art::Handle<mu2e::StrawHitFlagCollection> shflagH;
    if(event.getByLabel(_shfLabel,shflagH))
      shfcol = shflagH.product();
    if(shfcol == 0){
      throw cet::exception("RECO")<<"mu2e::MakeStereoHits: No StrawHitFlag collection found for label " <<  _shfLabel << endl;
    }
    if(shfcol->size() != strawhits->size()){
      throw cet::exception("RECO")<<"mu2e::MakeStereoHits: StrawHitFlag collection size " <<  shfcol->size() <<
       " doesn't match StrawHit collection size " << strawhits->size() << endl;
    }
    // create a collection of StrawHitPosition, and intialize them using the time division
    size_t nsh = strawhits->size();
    if(_diagLevel > 0)_nhits->Fill(nsh);
    if(_diagLevel > 1){
      art::Handle<StrawDigiMCCollection> mcHandle;
      if(event.getByLabel(_mcdigislabel,"StrawHitMC",mcHandle))
	_mcdigis = mcHandle.product();
    }
    unique_ptr<StrawHitPositionCollection> shpos(new StrawHitPositionCollection);
    shpos->reserve(2*nsh);

    size_t nres = max(size_t(100),nsh/20);
    vector <vector<vector<int> > >hdx(_nsta,vector<vector<int> >(_nsec*2));
    for(vector<vector<int>> rhdx : hdx){
      for(vector<int>  chdx : rhdx){
        chdx.reserve(nres);
      }
    }

    for(size_t ish=0;ish<nsh;++ish){
      StrawHit const& hit = strawhits->at(ish);
      Straw const& straw = tracker.getStraw(hit.strawIndex());

      int idev = straw.id().getDevice();
      int isec = straw.id().getSector();
      int jdev = idev/2;
      int jsec = _imap[idev][isec];
      hdx[jdev][jsec].push_back(ish);

      SHInfo shinfo;
      tcal->StrawHitInfo(straw,hit,shinfo);
      shpos->push_back(StrawHitPosition(hit,straw,shinfo));
    }
    // create the stereo hits
    StereoHitCollection stereohits;
    stereohits.reserve(3*nsh);
    vector<double> maxMVA(nsh,-FLT_MAX);
    vector<int> minsep(nsh,SectorId::apart);
    vector<size_t> ibest(nsh);
    // double loop over selected straw hits

    for(vector<vector<int>> rhdx : hdx){                                // loop over stations
      for(vector<int>::size_type isec=0; isec!=rhdx.size(); isec++){    // loop over sectors
        if(!_dosec[isec].empty()){
          for(int ish : rhdx[isec]){                                    // loop over hit1
            if(shfcol->at(ish).hasAllProperties(_shsel)){
              StrawHit const& sh1 = strawhits->at(ish);
              Straw const& straw1 = tracker.getStraw(sh1.strawIndex());
              StrawHitPosition const& shp1 = (*shpos)[ish];
              for( int jsec : _dosec[isec]){                            // loop over overlapping sectors
                for(int jsh : rhdx[jsec]){                              // loop over hit2
	          if(shfcol->at(jsh).hasAllProperties(_shsel)){
	            StrawHit const& sh2 = strawhits->at(jsh);
	            Straw const& straw2 = tracker.getStraw(sh2.strawIndex());
	            StrawHitPosition const& shp2 = (*shpos)[jsh];
	            double ddot = straw1.direction().dot(straw2.direction());
	            SectorId::isep sep = straw1.id().getSectorId().separation(straw2.id().getSectorId());
	            double de = min((float)1.0,fabs((sh1.energyDep() - sh2.energyDep())/(sh1.energyDep()+sh2.energyDep())));
	            CLHEP::Hep3Vector dp = shp1.pos()-shp2.pos();
	            double dist = dp.mag();
	            double dperp = dp.perp();
	            double dz = fabs(dp.z());
	            double dt = fabs(sh1.time()-sh2.time()); 
	            if( sep != SectorId::same && sep < SectorId::apart // hits are in the same station but not the same sector
	                && (sep <= minsep[ish] || sep <= minsep[jsh]) // this separation is at least as good as the current best for one of the hits
	                && ddot > _minDdot // negative crosings are in opposite quadrants
	                && dt < _maxDt // hits are close in time
	                && de < _maxDE   // deposited energy is roughly consistent (should compare dE/dx but have no path yet!)
	                && dz < _maxDZ // longitudinal separation isn't too big
	                && dperp < _maxDPerp) { // transverse separation isn't too big
	              // tentative stereo hit: this solves for the POCA
	              StereoHit sth(*strawhits,tracker,ish,jsh);
	              double dl1 = straw1.getDetail().activeHalfLength()-fabs(sth.wdist1());
	              double dl2 = straw2.getDetail().activeHalfLength()-fabs(sth.wdist2());
	              if(_diagLevel > 1 ) {
	                _dL->Fill(dl1);
	                _dL->Fill(dl2);
	              }
	              if( dl1 > _minDL && dl2 > _minDL) {
	                // stereo point is inside the active length
	                // compute difference between stereo points and TD prediction
	                double chi1 = (shp1.wireDist()-sth.wdist1())/shp1.posRes(StrawHitPosition::phi);	      
	                double chi2 = (shp2.wireDist()-sth.wdist2())/shp2.posRes(StrawHitPosition::phi);
	                if(fabs(chi1) <_maxChi && fabs(chi2) < _maxChi)
	                {
	                  // compute chisquared
	                  double chisq = chi1*chi1+chi2*chi2; 
	                  if(chisq < _maxChisq){
		            sth.setChisquared(chisq);
                            // compute MVA
		            _vmva[0] = de;
		            _vmva[1] = dt;
		            _vmva[2] = chisq;
		            _vmva[3] = sth.pos().perp();
		            double mvaout = _mvatool.evalMVA(_vmva);
                            //double mvaout=0.;
		            if(mvaout > _minMVA){
		              stereohits.push_back(sth);
		              // choose the best pair as:
		              // 1) take the pair with the minimum plane separation
		              // 2) otherwise, take the pair with the maximum MVA output
		              if(sep < minsep[ish] || (sep == minsep[ish] && mvaout > maxMVA[ish])){
		                minsep[ish] = sep;
		                maxMVA[ish] = mvaout;
		                ibest[ish] = stereohits.size()-1;
		              }
		              if(sep < minsep[jsh] || (sep == minsep[jsh] && mvaout > maxMVA[jsh])){
		                minsep[jsh] = sep;
		                maxMVA[jsh] = mvaout;
		                ibest[jsh] = stereohits.size()-1;
		              }
		              if(_diagLevel > 2){
		                _fs = sep;
		                _de = de;
		                _dt = dt;
		                _dz = dz;
		                _rho = sth.pos().perp();
		                _dist = dist;
		                _dperp = dperp;
		                _dl1 = dl1;
		                _dl2 = dl2;
		                _dc1 = chi1;
		                _dc2 = chi2;
		                _chi2 = chisq;
		                _mvaout = mvaout;
		                _ddot = ddot;
		                _mcdist = -1.0;
		                if(_mcdigis != 0){
		                  StrawDigiMC const& mcd1 = _mcdigis->at(ish);
		                  StrawDigiMC const& mcd2 = _mcdigis->at(jsh);
		                  _mcrel = KalFitMC::relationship(mcd1,mcd2);
		                  if(mcd1.stepPointMC(StrawDigi::zero).isNonnull() &&
			             mcd2.stepPointMC(StrawDigi::zero).isNonnull() )
			            _mcdist = (mcd1.stepPointMC(StrawDigi::zero)->position() -
			                mcd2.stepPointMC(StrawDigi::zero)->position()).mag();
		                }
		                _sdiag->Fill();
		              } // _diagLevel > 2
		            }
	                  }
	                }
	              }
	            }
	          } // shfcol->at(jsh).hasAllProperties(_shsel)
                } // loop over hit2: jsh
              } // loop over overlapping sectors: jsec
            } // shfcol->at(ish).hasAllProperties(_shsel)
          } // loop over hit1: ish
        } // dosec not empty
      } // loop over 2nd dim: isec, sectors
    } // loop over 1st dim: rhdx, stations
// now, overwrite the positions for those hits which have stereosresolve the stereo hits to find the best position for each hit that particpates.  The algorithm is:
    for(size_t ish=0; ish<nsh;++ish){
      bool stereo(false);
      if(minsep[ish] < SectorId::apart){
	shpos->at(ish) = StrawHitPosition(stereohits,ibest[ish],ish);
	stereo = true;
      }
      if(_diagLevel > 0){
	_fsep->Fill(minsep[ish]);
	_mva->Fill(maxMVA[ish]);
      }
      if(_diagLevel > 1){
	StrawHit const& hit = strawhits->at(ish);
	Straw const& straw = tracker.getStraw(hit.strawIndex());
	SHInfo shinfo;
	tcal->StrawHitInfo(straw,hit,shinfo);
	StrawHitPosition shp(hit,straw,shinfo);
	_shphi = shp.pos().phi();
	_shrho = shp.pos().perp();
	_stphi = shpos->at(ish).pos().phi();
	_strho = shpos->at(ish).pos().perp();
	_stereo = stereo;
	_mcr = _sfs = -1;
	_schi2 = _sdist = _sdz = _sddot = -1.0;
	if(stereo){
	  StereoHit const& sh = stereohits.at(ibest[ish]);
	  _schi2 = sh.chisq();
	  _smvaout = maxMVA[ish];
	  _sfs = sh.sectorSeparation();
	  StrawHit const h1 = strawhits->at(sh.hitIndex1());
	  StrawHit const h2 = strawhits->at(sh.hitIndex2());
	  Straw const& straw1 = tracker.getStraw(h1.strawIndex());
	  Straw const& straw2 = tracker.getStraw(h2.strawIndex());
	  SHInfo shi1, shi2;
	  tcal->StrawHitInfo(straw1,h1,shi1);
	  tcal->StrawHitInfo(straw2,h2,shi2);
	  StrawHitPosition shp1(h1,straw1,shi1);
	  StrawHitPosition shp2(h1,straw1,shi2);
	  _sdist = (shp1.pos()-shp2.pos()).mag();
	  _sddot = straw1.direction().dot(straw2.direction());
	  _sdz = fabs(straw1.getMidPoint().z()-straw2.getMidPoint().z());
	  StrawDigiMC const& mcd1 = _mcdigis->at(sh.hitIndex1());
	  StrawDigiMC const& mcd2 = _mcdigis->at(sh.hitIndex2());
	  _mcr = KalFitMC::relationship(mcd1,mcd2);
	}
	StrawDigiMC const& mcd = _mcdigis->at(ish);
	StrawDigi::TDCChannel itdc = StrawDigi::zero;
	_mcshphi = _mcshrho = -1000.0;
	_mcpdg = _mcproc = _mcgen = 0;
	if(mcd.hasTDC(itdc) && mcd.stepPointMC(itdc).isNonnull() ){
	  _mcshphi = mcd.stepPointMC(itdc)->position().phi();
	  _mcshrho = mcd.stepPointMC(itdc)->position().perp();
	  if(mcd.stepPointMC(itdc)->simParticle().isNonnull()){
	    _mcpdg = mcd.stepPointMC(itdc)->simParticle()->pdgId();
	    _mcproc = mcd.stepPointMC(itdc)->simParticle()->creationCode();
	    if(mcd.stepPointMC(itdc)->simParticle()->genParticle().isNonnull()){
	      _mcgen = mcd.stepPointMC(itdc)->simParticle()->genParticle()->generatorId().id();
	    }
	  }
	}
	_spdiag->Fill();
      }
    }
    if(_writepairs){
      unique_ptr<StereoHitCollection> sthits(new StereoHitCollection(stereohits));
      event.put(move(sthits));
    }
    event.put(move(shpos));
  } // end MakeStereoHits::produce.

  void MakeStereoHits::genMap(const Tracker& tracker) {

    vector <vector<int> >jdosec
                    { { 3, 5, 6, 8, 9,11},{ 3, 4, 6, 7, 9,10},{ 4, 5, 7, 8,10,11},
                      { 6, 8, 9,11},{ 6, 7, 9,10},{ 7, 8,10,11},{ 9,11},{ 9,10},{10,11},{},{},{} };
    const TTracker& tt = dynamic_cast<const TTracker&>(tracker);
    size_t ndev = tt.nDevices();
    _nsta = ndev/2;
    _nsec = tt.getDevice(0).getSectors().size();
    vector <vector<int> >jmap(ndev,vector<int>(_nsec));
    vector <vector<int> >newisec(2,vector<int>(_nsec));
    // clockwise sequence of sector indices per layer in remapped scheme
    //vector <vector<int> >newisec{ { 0, 3, 1, 4, 2, 5},{ 6, 9, 7,10, 8,11} };
    for(size_t i=0; i<2; i++){
      for(size_t j=0; j<_nsec-1; j+=2){
        int k=j/2+i*_nsec;
        newisec[i][j]=k;
        newisec[i][j+1]=k+3;
      }
    }

    for(size_t ista=0;ista<_nsta;++ista){
      CLHEP::Hep3Vector spos0;
      for(int idev=0;idev<2;++idev){
        vector<SecPhi>vsecphi;
        int jdev=2*ista+idev;
        const Device& dev = tt.getDevice(jdev);
        const vector<Sector>& sectors = dev.getSectors();
        // for 1st device in station (even numbered dev)
        if(idev==0){
          for(size_t isec=0;isec<_nsec;++isec){
            const Sector& sec = sectors[isec];
            CLHEP::Hep3Vector spos = sec.straw0MidPoint();
            // calculate deltaphi of all sectors relative to sector0
            if(isec==0)spos0 = spos;
            double dphi = spos0.deltaPhi(spos);
            if(dphi<0)dphi += 2*M_PI;
            SecPhi secphi;
            secphi.isec=isec;
            secphi.dphi=dphi;
            vsecphi.push_back(secphi);
          }
          // sort all sectors of device by deltaphi
          sort(vsecphi.begin(), vsecphi.end(),sortSec());
        // for 2nd device in station (odd numbered devices)
        }else{
          vector<SecPhi>vtmp;
          for(int ilyr=0;ilyr<2;ilyr++){
            for(size_t isec=ilyr;isec<_nsec;isec+=2){
              const Sector& sec = sectors[isec];
              CLHEP::Hep3Vector spos = sec.straw0MidPoint();
              double dphi = spos0.deltaPhi(spos);
              if(dphi<0)dphi += 2*M_PI;
              SecPhi secphi;
              secphi.isec=isec;
              secphi.dphi=dphi;
              vtmp.push_back(secphi);
            }
            // sort sectors of each layer separately by deltaphi
            sort(vtmp.begin(), vtmp.end(),sortSec());
            if(ilyr==0){
              // append front layer sector info to existing vector from 1st device
              vsecphi.insert(vsecphi.end(),vtmp.begin(),vtmp.end());
            }else{
              // insert rear layer sector info between elements of front layer
              int msec=_nsec/2;
              for(int i=0;i<msec;i++){
                vsecphi.insert(vsecphi.end()-msec+1+i,vtmp[i]);
              }
            }  
            vtmp.clear();
          }
        }
        for(vector<SecPhi>::size_type isec=0; isec!=vsecphi.size(); isec++){    // loop over sectors
          int sp_isec = vsecphi[isec].isec;
          jmap[jdev][sp_isec]=newisec[idev][isec];
        }
      } // loop over front and back
    } // loop over stations
    _imap=jmap;
    _dosec=jdosec;
    if(_debugLevel >0){
      // print out the generated map
      const string stars(16,'*');
      ostringstream sndev,snsec;
      sndev << ndev;
      snsec << _nsec;
      const string label = " Generated imap[" + sndev.str() + "][" + snsec.str() + "] ";
      cout << stars << label << stars << endl;
      for(size_t i=0;i<ndev;i++){
        if(i%2==0)cout << "Station " << i/2 << ": ";
        for(size_t j=0;j<_nsec;j++){
          if(j==0)cout << "{ ";
          cout << _imap[i][j];
          if(j<5){
            cout << ", ";
          }else{
            if(i%2==0){
              cout << " },";
            }else{
              cout << " }\n";
            }
          }
        }
      }
      const string morestars(32+label.size(),'*');
      cout << morestars << endl;
    }
  }

  double MakeStereoHits::deltaPhi(double phi1,double phi2) const {
    double dphi = phi2-phi1;
    if(dphi > M_PI){
      dphi -= 2*M_PI;
    } else if(dphi < -M_PI){
      dphi += 2*M_PI;
    }
    return dphi;
  }

} // end namespace mu2e

using mu2e::MakeStereoHits;
DEFINE_ART_MODULE(MakeStereoHits)

