//
// module to convert G4 steps into straw digis.
// It also builds the truth match
//
// Original author David Brown, LBNL
//
// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Core/EDProducer.h"
#include "GeometryService/inc/DetectorSystem.hh"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "SeedService/inc/SeedService.hh"
#include "cetlib_except/exception.h"
// conditions
#include "ProditionsService/inc/ProditionsHandle.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "TrackerConditions/inc/DeadStraw.hh"
#include "TrackerConditions/inc/StrawElectronics.hh"
#include "TrackerConditions/inc/StrawPhysics.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "BTrk/BField/BField.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
// utiliities
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "DataProducts/inc/TrkTypes.hh"
// data
#include "DataProducts/inc/EventWindowMarker.hh"
#include "DataProducts/inc/StrawId.hh"
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
// MC structures
#include "TrackerMC/inc/StrawClusterSequencePair.hh"
#include "TrackerMC/inc/StrawWaveform.hh"
#include "TrackerMC/inc/IonCluster.hh"
//CLHEP
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Vector/LorentzVector.h"
// root
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TMarker.h"
#include "TTree.h"
// C++
#include <map>
#include <algorithm>
#include <array>
#include <iostream>
#include <limits>
using namespace std;
using CLHEP::Hep3Vector;
namespace mu2e {
  namespace TrackerMC {
    using namespace TrkTypes;

    struct WireCharge { // charge at the wire after drift
      double _charge; // charge at the wire, in units of pC
      double _time; // relative time at the wire, relative to ionzation time (ns)
      double _dd;  // transverse distance drifted to the wrie
      double _phi; //JB: angle between E and B at ionization event
      double _wpos; // position long the wire, WRT the wire center, signed by the wire direction

    };

    struct WireEndCharge { // charge at one end of the wire after propagation
      double _charge; // charge at the wire, in units of pC
      double _time; // time at the wire end, relative to the time the charge reached the wire (ns)
      double _wdist; // propagation distance from the point of collection to the end
    };

    class StrawDigisFromStepPointMCs : public art::EDProducer {

    public:
      typedef map<StrawId,StrawClusterSequencePair> StrawClusterMap;  // clusts by straw
      typedef vector<art::Ptr<StepPointMC> > StrawSPMCPV; // vector of associated StepPointMCs for a single straw/particle
      // work with pairs of waveforms, one for each straw end
      typedef std::array<StrawWaveform,2> SWFP;
      typedef std::array<WFX,2> WFXP;
      typedef list<WFXP> WFXPList;
      typedef WFXPList::const_iterator WFXPI;

      explicit StrawDigisFromStepPointMCs(fhicl::ParameterSet const& pset);
      // Accept compiler written d'tor.

    private:

      void beginJob() override;
      void beginRun(art::Run& run) override;
      void produce(art::Event& e) override;

      // Diagnostics level.
      int _debug, _diag, _printLevel;
      unsigned _maxhist;
      bool  _xtalkhist;
      unsigned _minnxinghist;
      double _tstep, _nfall;
      // Limit on number of events for which there will be full printout.
      int _maxFullPrint;

      // Parameters
      bool   _addXtalk; // should we add cross talk hits?
      double _ctMinCharge; // minimum charge to add cross talk (for performance issues)
      bool   _addNoise; // should we add noise hits?
      double _preampxtalk, _postampxtalk; // x-talk parameters; these should come from conditions, FIXME!!
      double _bgcut; // cut dividing 'min-ion' particles from highly-ionizing
      double _minstepE; // minimum step energy to simulate
      art::InputTag _ewMarkerTag; // name of the module that makes eventwindowmarkers
      double _mbtime; // period of 1 microbunch
      double _mbbuffer; // buffer on that for ghost clusts (for waveform)
      double _adcbuffer; // time buffer for ADC
      double _steptimebuf; // buffer for MC step point times
      double _tdcbuf; // buffer for TDC jitter
      uint16_t _allStraw; // minimum straw # to read all hits
      std::vector<uint16_t> _allPlanes; // planes in which to read all hits
      // models of straw response to stimuli
      ProditionsHandle<StrawPhysics> _strawphys_h;
      ProditionsHandle<StrawElectronics> _strawele_h;
      SimParticleTimeOffset _toff;
      StrawElectronics::Path _diagpath; // electronics path for waveform diagnostics
      // Random number distributions
      art::RandomNumberGenerator::base_engine_t& _engine;
      CLHEP::RandGaussQ _randgauss;
      CLHEP::RandFlat _randflat;
      CLHEP::RandExponential _randexp;
      CLHEP::RandPoisson _randP;
      // A category for the error logger.
      const string _messageCategory;
      // Give some informationation messages only on the first event.
      bool _firstEvent;
      // record the BField direction at the tracker center
      Hep3Vector _bdir;
      // minimum pt (perp to bfield) to assume straight trajectory in a starw
      double _ptmin, _ptfac;
      // max # clusters for modeling non-minion steps
      unsigned _maxnclu;
      bool _sort; // sort cluster sizes before filling energy
      // List of dead straws as a parameter set; needed at beginRun time.
      ProditionsHandle<DeadStraw> _deadStraw_h;
      // StepPointMC selector
      // This selector will select only data products with the given instance name.
      // optionally exclude modules: this is a fix
      art::Selector _selector;

      // diagnostics
      TTree* _swdiag;
      Int_t _swplane, _swpanel, _swlayer, _swstraw, _ndigi;
      Float_t _hqsum[2], _vmax[2], _tvmax[2], _sesum[2];
      Int_t _wmcpdg[2], _wmcproc[2], _nxing[2], _nclu[2];
      Int_t _nsteppoint[2], _npart[2];
      Float_t _mce[2], _slen[2], _sedep[2];
      Float_t _tmin[2], _tmax[2], _txing[2], _xddist[2], _xwdist[2], _xpdist[2];
      TTree* _sddiag;
      Int_t _sdplane, _sdpanel, _sdlayer, _sdstraw;
      Int_t _ncludd[2], _iclust[2];
      Int_t _nstep;
      Float_t _ectime[2], _ecddist[2], _ecdtime[2], _ecptime[2];
      Float_t _xtime[2], _tctime[2], _charge[2], _ddist[2], _dtime[2], _ptime[2];
      Float_t _wdist[2], _vstart[2], _vcross[2];
      Float_t _phi[2]; //JB
      Float_t _mcenergy, _mctrigenergy, _mcthreshenergy;
      Double_t _mctime;
      Int_t _mcthreshpdg, _mcthreshproc, _mcnstep;
      Float_t _mcdca, _mcdcaphi, _mcdcadtime;
      Int_t _dmcpdg, _dmcproc, _dmcgen;
      Float_t _dmcmom;
      Bool_t _xtalk;
      vector<unsigned> _adc;
      Int_t _tdc[2], _tot[2];
      TTree* _sdiag;
      Float_t _steplen, _stepE, _qsum, _esum, _eesum, _qe, _partP, _steptime;
      Int_t _nclusd, _netot, _partPDG;
      vector<IonCluster> _clusters;
      Float_t _ewMarkerOffset;
      array<Float_t, StrawId::_nupanels> _ewMarkerROCdt;

      //  helper functions
      void fillClusterMap(art::Event const& event, StrawClusterMap & hmap);
      void addStep(StrawPhysics const& strawphys,
            StrawElectronics const& strawele,
            art::Ptr<StepPointMC> const& spmcptr,
            Straw const& straw, StrawClusterSequencePair& shsp);
      void divideStep(StrawPhysics const& strawphys,
                      StrawElectronics const& strawele,
                      StepPointMC const& step, vector<IonCluster>& clusters);
      void driftCluster(StrawPhysics const& strawphys, Straw const& straw,
                        IonCluster const& cluster, WireCharge& wireq);
      void propagateCharge(StrawPhysics const& strawphys, Straw const& straw,
                 WireCharge const& wireq, StrawEnd end, WireEndCharge& weq);
      double microbunchTime(StrawElectronics const& strawele, double globaltime) const;
      void addGhosts(StrawElectronics const& strawele, StrawCluster const& clust,StrawClusterSequence& shs);
      void addNoise(StrawClusterMap& hmap);
      void findThresholdCrossings(StrawElectronics const& strawele, SWFP const& swfp, WFXPList& xings);
      void createDigis(StrawPhysics const& strawphys,
                       StrawElectronics const& strawele,
                       StrawClusterSequencePair const& hsp,
                       XTalk const& xtalk,
                       StrawDigiCollection* digis, StrawDigiMCCollection* mcdigis,
                       PtrStepPointMCVectorCollection* mcptrs );
      void fillDigis(StrawPhysics const& strawphys,
                     StrawElectronics const& strawele,
                     WFXPList const& xings,SWFP const& swfp , StrawId sid,
                     StrawDigiCollection* digis, StrawDigiMCCollection* mcdigis,
                     PtrStepPointMCVectorCollection* mcptrs );
      bool createDigi(StrawElectronics const& strawele,WFXP const& xpair, SWFP const& wf, StrawId sid, StrawDigiCollection* digis);
      void findCrossTalkStraws(Straw const& straw,vector<XTalk>& xtalk);
      void fillClusterNe(StrawPhysics const& strawphys,std::vector<unsigned>& me);
      void fillClusterPositions(Straw const& straw, StepPointMC const& step, std::vector<Hep3Vector>& cpos);
      void fillClusterMinion(StrawPhysics const& strawphys, StepPointMC const& step, std::vector<unsigned>& me, std::vector<double>& cen);
      bool readAll(StrawId const& sid) const;
      // diagnostic functions
      void waveformHist(StrawElectronics const& strawele,
                        SWFP const& wf, WFXPList const& xings);
      void waveformDiag(StrawElectronics const& strawele,
                        SWFP const& wf, WFXPList const& xings);
      void digiDiag(StrawPhysics const& strawphys, SWFP const& wf, WFXP const& xpair, StrawDigi const& digi,StrawDigiMC const& mcdigi);
    };

    StrawDigisFromStepPointMCs::StrawDigisFromStepPointMCs(fhicl::ParameterSet const& pset) :
      EDProducer{pset},
    // diagnostic parameters
    _debug(pset.get<int>("debugLevel",0)),
    _diag(pset.get<int>("diagLevel",0)),
    _printLevel(pset.get<int>("printLevel",0)),
    _maxhist(pset.get<unsigned>("MaxHist",100)),
    _xtalkhist(pset.get<bool>("CrossTalkHist",false)),
    _minnxinghist(pset.get<int>("MinNXingHist",1)), // minimum # of crossings to histogram waveform
    _tstep(pset.get<double>("WaveformStep",0.1)), // ns
    _nfall(pset.get<double>("WaveformTail",10.0)),  // # of decay lambda past last signal to record waveform
    // Parameters
    _maxFullPrint(pset.get<int>("maxFullPrint",2)),
    _addXtalk(pset.get<bool>("addCrossTalk",false)),
    _ctMinCharge(pset.get<double>("xtalkMinimumCharge",0)),
    _addNoise(pset.get<bool>("addNoise",false)),
    _preampxtalk(pset.get<double>("preAmplificationCrossTalk",0.0)),
    _postampxtalk(pset.get<double>("postAmplificationCrossTalk",0.02)), // dimensionless relative coupling
    _bgcut(pset.get<double>("BetaGammaCut",0.5)), // treat particles with beta-gamma above this as minimum-ionizing
    _minstepE(pset.get<double>("minstepE",2.0e-6)), // minimum step energy depostion to turn into a straw signal (MeV)
    _ewMarkerTag(pset.get<art::InputTag>("EventWindowMarkerLabel","EWMProducer")),
    _steptimebuf(pset.get<double>("StepPointMCTimeBuffer",100.0)), // nsec
    _tdcbuf(pset.get<double>("TDCTimeBuffer",2.0)), // nsec
    _allStraw(pset.get<uint16_t>("AllHitsStraw",90)),
    _allPlanes(pset.get<std::vector<uint16_t>>("AllHitsPlanes",std::vector<uint16_t>{})), // planes to read all hits
    _toff(pset.get<fhicl::ParameterSet>("TimeOffsets", {})),
    _diagpath(static_cast<StrawElectronics::Path>(pset.get<int>("WaveformDiagPath",StrawElectronics::thresh))),
    // Random number distributions
    _engine(createEngine( art::ServiceHandle<SeedService>()->getSeed())),
    _randgauss( _engine ),
    _randflat( _engine ),
    _randexp( _engine),
    _randP( _engine),
    _messageCategory("HITS"),
    _firstEvent(true),      // Control some information messages.
    _ptfac(pset.get<double>("PtFactor", 2.0)), // factor for defining curling in a straw
    _maxnclu(pset.get<unsigned>("MaxNClusters", 10)), // max # of clusters for low-PT steps
    _sort(pset.get<bool>("SortClusterEnergy",false)),
    _selector{art::ProductInstanceNameSelector(pset.get<string>("trackerStepPoints","tracker")) &&
      !art::ModuleLabelSelector(pset.get<string>("SkipTheseStepPoints","")) } {
      // Tell the framework what we consume.
      consumesMany<StepPointMCCollection>();
      consumes<EventWindowMarker>(_ewMarkerTag);
      // Since SimParticleTimeOffset calls getValidHandle, we have to
      // declare the consumes statements here.
      auto const& toffInputs = pset.get<std::vector<std::string>>("TimeOffsets.inputs", {});
      for (auto const& tag : toffInputs) {
        consumes<SimParticleTimeMap>(tag);
      }
      // Tell the framework what we make.
      produces<StrawDigiCollection>();
      produces<PtrStepPointMCVectorCollection>();
      produces<StrawDigiMCCollection>();
    }

    void StrawDigisFromStepPointMCs::beginJob(){

      if(_diag > 0){

        art::ServiceHandle<art::TFileService> tfs;
        _sdiag =tfs->make<TTree>("sdiag","Step diagnostics");
        _sdiag->Branch("steplen",&_steplen,"steplen/F");
        _sdiag->Branch("stepE",&_stepE,"stepE/F");
        _sdiag->Branch("partP",&_partP,"partP/F");
        _sdiag->Branch("qsum",&_qsum,"qsum/F");
        _sdiag->Branch("esum",&_esum,"esum/F");
        _sdiag->Branch("eesum",&_eesum,"eesum/F");
        _sdiag->Branch("qe",&_qe,"qe/F");
        _sdiag->Branch("steptime",&_steptime,"steptime/F");
        _sdiag->Branch("nclust",&_nclusd,"nclust/I");
        _sdiag->Branch("netot",&_netot,"netot/I");
        _sdiag->Branch("partPDG",&_partPDG,"partPDG/I");
        _sdiag->Branch("clusters",&_clusters);

        _swdiag =tfs->make<TTree>("swdiag","StrawWaveform diagnostics");
        _swdiag->Branch("plane",&_swplane,"plane/I");
        _swdiag->Branch("panel",&_swpanel,"panel/I");
        _swdiag->Branch("layer",&_swlayer,"layer/I");
        _swdiag->Branch("straw",&_swstraw,"straw/I");
        _swdiag->Branch("ndigi",&_ndigi,"ndigi/I");
        _swdiag->Branch("hqsum",&_hqsum,"hqsumcal/F:hqsumhv/F");
        _swdiag->Branch("vmax",&_vmax,"vmaxcal/F:vmaxhv/F");
        _swdiag->Branch("tvmax",&_tvmax,"tvmaxcal/F:tvmaxhv/F");
        _swdiag->Branch("mcpdg",&_wmcpdg,"mcpdgcal/I:mcpdghv/I");
        _swdiag->Branch("mcproc",&_wmcproc,"mcproccal/I:mcprochv/I");
        _swdiag->Branch("mce",&_mce,"mcecal/F:mcehv/F");
        _swdiag->Branch("slen",&_slen,"slencal/F:slenhv/F");
        _swdiag->Branch("sedep",&_sedep,"sedepcal/F:sedephv/F");
        _swdiag->Branch("nxing",&_nxing,"nxingcal/I:nxinghv/I");
        _swdiag->Branch("nclust",&_nclu,"nclucal/I:ncluhv/I");
        _swdiag->Branch("nstep",&_nsteppoint,"nscal/I:nshv/I");
        _swdiag->Branch("sesum",&_sesum,"sesumcal/F:sesumhv/F");
        _swdiag->Branch("npart",&_npart,"npart/I");
        _swdiag->Branch("tmin",&_tmin,"tmincal/F:tminhv/F");
        _swdiag->Branch("tmax",&_tmax,"tmaxcal/F:tmaxhv/F");
        _swdiag->Branch("txing",&_txing,"txcal/F:txhv/F");
        _swdiag->Branch("xddist",&_xddist,"xdcal/F:xdhv/F");
        _swdiag->Branch("xwdist",&_xwdist,"xwdcal/F:xwdhv/F");
        _swdiag->Branch("xpdist",&_xpdist,"xpdcal/F:xpdhv/F");


        if(_diag > 1){
          _sddiag =tfs->make<TTree>("sddiag","StrawDigi diagnostics");
          _sddiag->Branch("plane",&_sdplane,"plane/I");
          _sddiag->Branch("panel",&_sdpanel,"panel/I");
          _sddiag->Branch("layer",&_sdlayer,"layer/I");
          _sddiag->Branch("straw",&_sdstraw,"straw/I");
          _sddiag->Branch("nstep",&_nstep,"nstep/I");
          _sddiag->Branch("xtime",&_xtime,"xtimecal/F:xtimehv/F");
          _sddiag->Branch("tctime",&_tctime,"tctimecal/F:tctimehv/F");
          _sddiag->Branch("ectime",&_ectime,"ectimecal/F:ectimehv/F");
          _sddiag->Branch("ecdtime",&_ecdtime,"ecdtimecal/F:ecdtimehv/F");
          _sddiag->Branch("ecptime",&_ecptime,"ecptimecal/F:ecptimehv/F");
          _sddiag->Branch("charge",&_charge,"chargecal/F:chargehv/F");
          _sddiag->Branch("wdist",&_wdist,"wdistcal/F:wdisthv/F");
          _sddiag->Branch("phi",&_phi,"phical/F:phihv/F");//JB
          _sddiag->Branch("ecddist",&_ecddist,"ecddistcal/F:ecddisthv/F");
          _sddiag->Branch("vstart",&_vstart,"vstartcal/F:vstarthv/F");
          _sddiag->Branch("vcross",&_vcross,"vcrosscal/F:vcrosshv/F");
          _sddiag->Branch("ddist",&_ddist,"ddistcal/F:ddisthv/F");
          _sddiag->Branch("dtime",&_dtime,"dtimecal/F:dtimehv/F");
          _sddiag->Branch("ptime",&_ptime,"ptimecal/F:ptimehv/F");
          _sddiag->Branch("nclust",&_ncludd,"nclustcal/I:nclusthv/I");
          _sddiag->Branch("iclust",&_iclust,"iclustcal/I:iclusthv/I");
          _sddiag->Branch("tdc",&_tdc,"tdccal/I:tdchv/I");
          _sddiag->Branch("tot",&_tot,"totcal/I:tothv/I");
          _sddiag->Branch("adc",&_adc);
          _sddiag->Branch("mctime",&_mctime,"mctime/D");
          _sddiag->Branch("mcenergy",&_mcenergy,"mcenergy/F");
          _sddiag->Branch("mctrigenergy",&_mctrigenergy,"mctrigenergy/F");
          _sddiag->Branch("mcthreshenergy",&_mcthreshenergy,"mcthreshenergy/F");
          _sddiag->Branch("mcthreshpdg",&_mcthreshpdg,"mcthreshpdg/I");
          _sddiag->Branch("mcthreshproc",&_mcthreshproc,"mcthreshproc/I");
          _sddiag->Branch("mcnstep",&_mcnstep,"mcnstep/I");
          _sddiag->Branch("mcdca",&_mcdca,"mcdca/F");
          _sddiag->Branch("mcdcaphi",&_mcdcaphi,"mcdcaphi/F");
          _sddiag->Branch("mcdcadtime",&_mcdcadtime,"mcdcadtime/F");
          _sddiag->Branch("mcpdg",&_dmcpdg,"mcpdg/I");
          _sddiag->Branch("mcproc",&_dmcproc,"mcproc/I");
          _sddiag->Branch("mcgen",&_dmcgen,"mcgen/I");
          _sddiag->Branch("mcmom",&_dmcmom,"mcmom/F");
          _sddiag->Branch("xtalk",&_xtalk,"xtalk/B");

        }
      }
    }

    void StrawDigisFromStepPointMCs::beginRun( art::Run& run ){
      // get field at the center of the tracker
      GeomHandle<BFieldManager> bfmgr;
      GeomHandle<DetectorSystem> det;
      Hep3Vector vpoint_mu2e = det->toMu2e(Hep3Vector(0.0,0.0,0.0));
      Hep3Vector b0 = bfmgr->getBField(vpoint_mu2e);
      // if the field is too small, don't perform any curvature-based analysis
      if ( b0.mag() < 1.0e-4 ){
        _bdir = Hep3Vector(0.0,0.0,1.0);
        _ptmin = -1.0;
      } else {
        _bdir = b0.unit();
        //compute the transverse momentum for which a particle will curl up in a straw
        const Tracker& tracker = *GeomHandle<Tracker>();
        const Straw& straw = tracker.getStraw(StrawId(0,0,0));
        double rstraw = straw.getRadius();
        _ptmin = _ptfac*BField::mmTeslaToMeVc*b0.mag()*rstraw;
      }
      if ( _printLevel > 0 ) {
        auto const& strawphys = _strawphys_h.get(run.id());
        strawphys.print(cout);
      }
    }

    void StrawDigisFromStepPointMCs::produce(art::Event& event) {
      if ( _printLevel > 1 ) cout << "StrawDigisFromStepPointMCs: produce() begin; event " << event.id().event() << endl;
      static int ncalls(0);
      ++ncalls;
      // update conditions caches.
      ConditionsHandle<AcceleratorParams> accPar("ignored");
      _mbtime = accPar->deBuncherPeriod;
      _toff.updateMap(event);
      StrawElectronics const& strawele = _strawele_h.get(event.id());
      StrawPhysics const& strawphys = _strawphys_h.get(event.id());
      art::Handle<EventWindowMarker> ewMarkerHandle;
      event.getByLabel(_ewMarkerTag, ewMarkerHandle);
      const EventWindowMarker& ewMarker(*ewMarkerHandle);
      _ewMarkerOffset = ewMarker.timeOffset();
      // calculate event window marker jitter for this microbunch for each panel
      for (size_t i=0;i<StrawId::_nupanels;i++){
        _ewMarkerROCdt.at(i) = _randgauss.fire(0,strawele.eventWindowMarkerROCJitter());
      }
      const Tracker& tracker = *GeomHandle<Tracker>();
      // make the microbunch buffer long enough to get the full waveform
      _mbbuffer = (strawele.nADCSamples() - strawele.nADCPreSamples())*strawele.adcPeriod();
      _adcbuffer = 0.01*strawele.adcPeriod();
      // Containers to hold the output information.
      unique_ptr<StrawDigiCollection> digis(new StrawDigiCollection);
      unique_ptr<StrawDigiMCCollection> mcdigis(new StrawDigiMCCollection);
      unique_ptr<PtrStepPointMCVectorCollection> mcptrs(new PtrStepPointMCVectorCollection);
      // create the StrawCluster map
      StrawClusterMap hmap;
      // fill this from the event
      fillClusterMap(event,hmap);
      // add noise clusts
      if(_addNoise)addNoise(hmap);
      // loop over the clust sequences
      for(auto ihsp=hmap.begin();ihsp!= hmap.end();++ihsp){
        StrawClusterSequencePair const& hsp = ihsp->second;
        // create primary digis from this clust sequence
        XTalk self(hsp.strawId()); // this object represents the straws coupling to itself, ie 100%
        createDigis(strawphys,strawele,hsp,self,digis.get(),mcdigis.get(),mcptrs.get());
        // if we're applying x-talk, look for nearby coupled straws
        if(_addXtalk) {
          // only apply if the charge is above a threshold
          double totalCharge = 0;
          for(auto ih=hsp.clustSequence(StrawEnd::cal).clustList().begin();ih!= hsp.clustSequence(StrawEnd::cal).clustList().end();++ih){
            totalCharge += ih->charge();
          }
          if( totalCharge > _ctMinCharge){
            vector<XTalk> xtalk;
            Straw const& straw = tracker.getStraw(hsp.strawId());
            findCrossTalkStraws(straw,xtalk);
            for(auto ixtalk=xtalk.begin();ixtalk!=xtalk.end();++ixtalk){
              createDigis(strawphys,strawele,hsp,*ixtalk,digis.get(),mcdigis.get(),mcptrs.get());
            }
          }
        }
      }
      // store the digis in the event
      event.put(move(digis));
      // store MC truth match
      event.put(move(mcdigis));
      event.put(move(mcptrs));
      if ( _printLevel > 1 ) cout << "StrawDigisFromStepPointMCs: produce() end" << endl;
      // Done with the first event; disable some messages.
      _firstEvent = false;
    } // end produce

    void StrawDigisFromStepPointMCs::createDigis(
          StrawPhysics const& strawphys,
          StrawElectronics const& strawele,
          StrawClusterSequencePair const& hsp, XTalk const& xtalk,
          StrawDigiCollection* digis, StrawDigiMCCollection* mcdigis,
          PtrStepPointMCVectorCollection* mcptrs ) {
      // instantiate waveforms for both ends of this straw
      SWFP waveforms  ={ StrawWaveform(hsp.clustSequence(StrawEnd::cal),xtalk),
        StrawWaveform(hsp.clustSequence(StrawEnd::hv),xtalk) };
      // find the threshold crossing points for these waveforms
      WFXPList xings;
      // find the threshold crossings
      findThresholdCrossings(strawele,waveforms,xings);
      // convert the crossing points into digis, and add them to the event data
      fillDigis(strawphys,strawele,xings,waveforms,xtalk._dest,digis,mcdigis,mcptrs);
      // waveform diagnostics
      if (_diag >1 && (
                       waveforms[0].clusts().clustList().size() > 0 ||
                       waveforms[1].clusts().clustList().size() > 0 ) ) {
        // waveform xing diagnostics
        _ndigi = digis->size();
        waveformDiag(strawele,waveforms,xings);
        // waveform histograms
        if(_diag > 2 )waveformHist(strawele,waveforms,xings);
      }
    }

    void StrawDigisFromStepPointMCs::fillClusterMap(art::Event const& event, StrawClusterMap & hmap){
      // get conditions
      DeadStraw const& deadStraw = _deadStraw_h.get(event.id());
      StrawPhysics const& strawphys = _strawphys_h.get(event.id());
      StrawElectronics const& strawele = _strawele_h.get(event.id());
      const Tracker& tracker = *GeomHandle<Tracker>();
      // Get all of the tracker StepPointMC collections from the event:
      typedef vector< art::Handle<StepPointMCCollection> > HandleVector;
      HandleVector stepsHandles;
      event.getMany( _selector, stepsHandles);
      // Informational message on the first event.
      if ( _firstEvent ) {
        mf::LogInfo log(_messageCategory);
        log << "StrawDigisFromStepPointMCs::fillHitMap will use StepPointMCs from: \n";
        for ( HandleVector::const_iterator i=stepsHandles.begin(), e=stepsHandles.end();
             i != e; ++i ){
          art::Provenance const& prov(*(i->provenance()));
          log  << "   " << prov.branchName() << "\n";
        }
      }
      if(stepsHandles.empty()){
        throw cet::exception("SIM")<<"mu2e::StrawDigisFromStepPointMCs: No StepPointMC collections found for tracker" << endl;
      }
      // Loop over StepPointMC collections
      for ( HandleVector::const_iterator ispmcc=stepsHandles.begin(), espmcc=stepsHandles.end();ispmcc != espmcc; ++ispmcc ){
        art::Handle<StepPointMCCollection> const& handle(*ispmcc);
        StepPointMCCollection const& steps(*handle);
        // Loop over the StepPointMCs in this collection
        for (size_t ispmc =0; ispmc<steps.size();++ispmc){
          // find straw id
          StrawId const & sid = steps[ispmc].strawId();
          // Skip dead straws, and straws that don't exist
          if(tracker.strawExists(sid)) {
            // lookup straw here, to avoid having to find the tracker for every step
            Straw const& straw = tracker.getStraw(sid);
            // Skip steps that occur in the deadened region near the end of each wire,
            // or in dead regions of the straw
            double wpos = fabs((steps[ispmc].position()-straw.getMidPoint()).dot(straw.getDirection()));

            if(wpos <  straw.activeHalfLength() &&
               deadStraw.isAlive(sid,wpos) &&
               steps[ispmc].ionizingEdep() > _minstepE){
              // create ptr to MC truth, used for references
              art::Ptr<StepPointMC> spmcptr(handle,ispmc);
              // create a clust from this step, and add it to the clust map
              addStep(strawphys,strawele,spmcptr,straw,hmap[sid]);
            }
          }
        }
      }
    }

    void StrawDigisFromStepPointMCs::addStep(StrawPhysics const& strawphys,
                StrawElectronics const& strawele,
                art::Ptr<StepPointMC> const& spmcptr,
                Straw const& straw, StrawClusterSequencePair& shsp) {
      StepPointMC const& step = *spmcptr;
      StrawId sid = straw.id();
     // get time offset for this step
      double tstep = _toff.timeWithOffsetsApplied(step);
      // test if this step point is roughly in the digitization window
      double mbtime = microbunchTime(strawele,tstep);
      if( (mbtime > strawele.flashEnd() - _steptimebuf
            && mbtime <  strawele.flashStart())
          || readAll(sid)) {
        // Subdivide the StepPointMC into ionization clusters
        _clusters.clear();
        divideStep(strawphys,strawele,step,_clusters);
        // check
        if(_debug > 1){
          double ec(0.0);
          double ee(0.0);
          double eq(0.0);
          for (auto const& cluster : _clusters) {
            ec += cluster._eion;
            ee += strawphys.ionizationEnergy(cluster._ne);
            eq += strawphys.ionizationEnergy(cluster._charge);
          }
          cout << "step with ionization edep = " << step.ionizingEdep()
            << " creates " << _clusters.size()
            << " clusters with total cluster energy = " << ec
            << " electron count energy = " << ee
            << " charge energy = " << eq << endl;
        }
        // drift these clusters to the wire, and record the charge at the wire
        for(auto iclu = _clusters.begin(); iclu != _clusters.end(); ++iclu){
          WireCharge wireq;
          driftCluster(strawphys,straw,*iclu,wireq);
          // propagate this charge to each end of the wire
          for(size_t iend=0;iend<2;++iend){
            StrawEnd end(static_cast<StrawEnd::End>(iend));
            // compute the longitudinal propagation effects
            WireEndCharge weq;
            propagateCharge(strawphys,straw,wireq,end,weq);
            // compute the total time, modulo the microbunch
            double gtime = tstep + wireq._time + weq._time;
            // convert from
            double ctime = microbunchTime(strawele,gtime);
            // create the clust
            StrawCluster clust(StrawCluster::primary,sid,end,ctime,weq._charge,wireq._dd,wireq._phi,weq._wdist,wireq._time,weq._time,
                               spmcptr,CLHEP::HepLorentzVector(iclu->_pos,mbtime)); //JB: + wireq._phi

            // add the clusts to the appropriate sequence.
            shsp.clustSequence(end).insert(clust);
            // if required, add a 'ghost' copy of this clust
            addGhosts(strawele,clust,shsp.clustSequence(end));
          }
        }
      }
    }

    void StrawDigisFromStepPointMCs::divideStep(StrawPhysics const& strawphys,
                StrawElectronics const& strawele,
                StepPointMC const& step, vector<IonCluster>& clusters) {
      // get particle charge
      double charge(0.0);
      GlobalConstantsHandle<ParticleDataTable> pdt;
      if(pdt->particle(step.simParticle()->pdgId()).isValid()){
        charge = pdt->particle(step.simParticle()->pdgId()).ref().charge();
      }

      // get tracker information
      const Tracker& tracker = *GeomHandle<Tracker>(); //JB
      const Straw& straw = tracker.getStraw(step.strawId());//JB
      // if the step length is small compared to the mean free path, or this is an
      // uncharged particle, put all the energy in a single cluster
      if (charge == 0.0 || step.stepLength() < strawphys.meanFreePath()){
        double cen = step.ionizingEdep();
        double fne = cen/strawphys.meanElectronEnergy();
        unsigned ne = std::max( static_cast<unsigned>(_randP(fne)),(unsigned)1);

        Hep3Vector cdir = (step.position()-straw.getMidPoint());//JB
        cdir -= straw.getDirection()*(cdir.dot(straw.getDirection()));//JB
        double phi = cdir.theta(); //JB
        for (size_t i=0;i<ne;i++){
          IonCluster cluster(step.position(),phi,strawphys.ionizationCharge((unsigned)1),cen/(float)ne,1); //JB + phi
          clusters.push_back(cluster);
        }
      }
      else {
        // use beta-gamma to decide if this is a min-ion particle or not
        bool minion(false);
        if(pdt->particle(step.simParticle()->pdgId()).isValid()){
          double mass = pdt->particle(step.simParticle()->pdgId()).ref().mass();
          double mom = step.momentum().mag();
          // approximate pt
          double apt = 0.;
          if ( _bdir.mag() > 0 ) apt = step.momentum().perpPart(_bdir).mag();
          double bg = mom/mass; // beta x gamma
          minion = bg > _bgcut && apt > _ptmin;
        }


        // compute the number of clusters for this step from the mean free path
        double fnc = step.stepLength()/strawphys.meanFreePath();
        // use a truncated Poisson distribution; this keeps both the mean and variance physical
        unsigned nc = std::max(static_cast<unsigned>(_randP.fire(fnc)),(unsigned)1);
        if(!minion)nc = std::min(nc,_maxnclu);
        // require clusters not exceed the energy sum required for single-electron clusters
        nc = std::min(nc,static_cast<unsigned>(floor(step.ionizingEdep()/strawphys.ionizationEnergy((unsigned)1))));
        // generate random positions for the clusters
        std::vector<Hep3Vector> cpos(nc);
        fillClusterPositions(straw,step,cpos);
        // generate electron counts and energies for these clusters: minion model is more detailed
        std::vector<unsigned> ne(nc);
        std::vector<double> cen(nc);
        if(minion){
          fillClusterMinion(strawphys,step,ne,cen);
        } else {
          // get Poisson distribution of # of electrons for the average energy
          double fne = step.ionizingEdep()/(nc*strawphys.meanElectronEnergy()); // average # of electrons/cluster for non-minion clusters
          for(unsigned ic=0;ic<nc;++ic){
            ne[ic] = static_cast<unsigned>(std::max(_randP.fire(fne),(long)1));
            cen[ic] = ne[ic]*strawphys.meanElectronEnergy(); // average energy per electron, works for large numbers of electrons
          }
        }
        // create the cluster objects
        for(unsigned ic=0;ic<nc;++ic){
          // compute phi for this cluster
          Hep3Vector cdir = (cpos[ic]-straw.getMidPoint());
          cdir -= straw.getDirection()*(cdir.dot(straw.getDirection()));
          double phi = cdir.theta(); //JB: point1 and point2 are the DCA on/to the track/wire. theta takes the angle to the z-axis
          for (size_t i=0;i<ne[ic];i++){
            IonCluster cluster(cpos[ic],phi,strawphys.ionizationCharge((unsigned)1),cen[ic]/(float)ne[ic],1);
            clusters.push_back(cluster);
          }
          //cout <<"phi1b = "<<phi<<"\n";
        }
      }
      // diagnostics
      if(_diag > 0){
        _steplen = step.stepLength();
        _stepE = step.ionizingEdep();
        _steptime = microbunchTime(strawele,_toff.timeWithOffsetsApplied(step));
        _partP = step.momentum().mag();
        _partPDG = step.simParticle()->pdgId();
        _nclusd = (int)clusters.size();
        _netot = 0;
        _qsum = _esum = _eesum = 0.0;
        for(auto iclust=clusters.begin();iclust != clusters.end();++iclust){
          _netot += iclust->_ne;
          _qsum += iclust->_charge;
          _esum += iclust->_eion;
          _eesum += strawphys.meanElectronEnergy()*iclust->_ne;
        }
        _qe = strawphys.ionizationEnergy(_qsum);
        _sdiag->Fill();
      }
    }

    void StrawDigisFromStepPointMCs::driftCluster(
                StrawPhysics const& strawphys,Straw const& straw,
                IonCluster const& cluster, WireCharge& wireq ) {
      // Compute the vector from the cluster to the wire
      Hep3Vector cpos = cluster._pos-straw.getMidPoint();
      // drift distance perp to wire, and angle WRT magnetic field (for Lorentz effect)
      double dd = min(cpos.perp(straw.getDirection()),straw.innerRadius());
      // sample the gain for this cluster
      double gain = strawphys.clusterGain(_randgauss, _randflat, cluster._ne);
      wireq._charge = cluster._charge*(gain);
      // compute drift time for this cluster
      double dt = strawphys.driftDistanceToTime(dd,cluster._phi); //JB: this is now from the lorentz corrected r-component of the drift
      wireq._phi = cluster._phi; //JB
      wireq._time = _randgauss.fire(dt,strawphys.driftTimeSpread(dd));
      wireq._dd = dd;
      // position along wire
      wireq._wpos = cpos.dot(straw.getDirection());

    }

    void StrawDigisFromStepPointMCs::propagateCharge(
          StrawPhysics const& strawphys, Straw const& straw,
          WireCharge const& wireq, StrawEnd end, WireEndCharge& weq) {
      // compute distance to the appropriate end
      double wlen = straw.halfLength(); // use the full length, not the active length
      // NB: the following assumes the straw direction points in increasing azimuth.  FIXME!
      if(end == StrawEnd::hv)
        weq._wdist = wlen - wireq._wpos;
      else
        weq._wdist = wlen + wireq._wpos;
      // split the charge
      weq._charge = 0.5*wireq._charge;
      weq._time = strawphys.propagationTime(weq._wdist);
    }

    double StrawDigisFromStepPointMCs::microbunchTime(StrawElectronics const& strawele, double globaltime) const {
      // converts time from proton beam time (StepPointMC time) to event window marker time
      // fold time relative to MB frequency
      double mbtime = fmod(globaltime - _ewMarkerOffset,_mbtime);
      // keep the microbunch time contiguous
      if(mbtime < strawele.flashStart()-_mbtime ) mbtime += _mbtime;
      return mbtime;
    }

    void StrawDigisFromStepPointMCs::addGhosts(StrawElectronics const& strawele,StrawCluster const& clust,StrawClusterSequence& shs) {
      // add enough buffer to cover both the flash blanking and the ADC waveform
      if(clust.time() < strawele.flashStart() - _mbtime + _mbbuffer)
        shs.insert(StrawCluster(clust,_mbtime));
      if(clust.time() > _mbtime - _mbbuffer) shs.insert(StrawCluster(clust,-_mbtime));
    }

    void StrawDigisFromStepPointMCs::findThresholdCrossings(StrawElectronics const& strawele, SWFP const& swfp, WFXPList& xings){
      //randomize the threshold to account for electronics noise; this includes parts that are coherent
      // for both ends (coming from the straw itself)
      // Keep track of crossings on each end to keep them in sequence
      double strawnoise = _randgauss.fire(0,strawele.strawNoise());
      // add specifics for each end
      double thresh[2] = {_randgauss.fire(strawele.threshold(swfp[0].strawId(),static_cast<StrawEnd::End>(0))+strawnoise,strawele.analogNoise(StrawElectronics::thresh)),
        _randgauss.fire(strawele.threshold(swfp[0].strawId(),static_cast<StrawEnd::End>(1))+strawnoise,strawele.analogNoise(StrawElectronics::thresh))};
      // Initialize search when the electronics becomes enabled:
      double tstart =strawele.flashEnd() - 10.0; // this buffer should be a parameter FIXME!
      // for reading all hits, make sure we start looking for clusters at the minimum possible cluster time
      // this accounts for deadtime effects from previous microbunches
      if(readAll(swfp[0].strawId()))tstart = -strawele.deadTimeAnalog();
      WFXP wfx = {WFX(swfp[0],tstart),WFX(swfp[1],tstart)};
      // search for coherent crossings on both ends
      bool crosses[2];
      for(size_t iend=0;iend<2;++iend){
        crosses[iend] = swfp[iend].crossesThreshold(strawele,thresh[iend],wfx[iend]);
      }
      // loop until we hit the end of the waveforms.  Require both in time.  Buffer to account for eventual TDC jitter
      // this is a loose pre-selection, final selection is done at digitization
      while( crosses[0] && crosses[1] && std::max(wfx[0]._time,wfx[1]._time)
          < strawele.flashStart() + strawele.electronicsTimeDelay() + _tdcbuf){
        // see if the crossings match
        if(strawele.combineEnds(wfx[0]._time,wfx[1]._time)){
          // put the pair of crossings in the crosing list
          // make sure the time is positive in case this was a hit from the 'previous' microbunch
          if(std::min(wfx[0]._time,wfx[1]._time) > 0.0 )xings.push_back(wfx);
          // search for next crossing:
          // update threshold for straw noise
          strawnoise = _randgauss.fire(0,strawele.strawNoise());
          for(unsigned iend=0;iend<2;++iend){
            // insure a minimum time buffer between crossings
            wfx[iend]._time += strawele.deadTimeAnalog();
            // skip to the next clust
            ++(wfx[iend]._iclust);
            // update threshold for incoherent noise
            thresh[iend] = _randgauss.fire(strawele.threshold(swfp[0].strawId(),static_cast<StrawEnd::End>(iend)),strawele.analogNoise(StrawElectronics::thresh));
            // find next crossing
            crosses[iend] = swfp[iend].crossesThreshold(strawele,thresh[iend],wfx[iend]);
          }
        } else {
          // skip to the next crossing on the earlier waveform
          unsigned iearly = wfx[0]._time < wfx[1]._time ? 0 : 1;
          ++(wfx[iearly]._iclust);
          wfx[iearly]._time += strawele.deadTimeAnalog();
          crosses[iearly] = swfp[iearly].crossesThreshold(strawele,thresh[iearly],wfx[iearly]);
        }
      }
    }

    void StrawDigisFromStepPointMCs::fillDigis(StrawPhysics const& strawphys,
                                               StrawElectronics const& strawele,
        WFXPList const& xings, SWFP const& wf,
        StrawId sid,
        StrawDigiCollection* digis, StrawDigiMCCollection* mcdigis,
        PtrStepPointMCVectorCollection* mcptrs ){
      // loop over crossings
      for(auto xpair : xings) {
        // create a digi from this pair.  This also performs a finial test
        // on whether the pair should make a digi
        if(createDigi(strawele,xpair,wf,sid,digis)){
          // fill associated MC truth matching. Only count the same step once
          set<art::Ptr<StepPointMC> > xmcsp;
          double wetime[2] = {-100.,-100.};
          CLHEP::HepLorentzVector cpos[2];
          art::Ptr<StepPointMC> stepMC[2];
          set<art::Ptr<StepPointMC> > spmcs;
          for (size_t iend=0;iend<2;++iend){
            StrawCluster const& sc = *(xpair[iend]._iclust);
            xmcsp.insert(sc.stepPointMC());
            wetime[iend] = sc.time();
            cpos[iend] = sc.clusterPosition();
            stepMC[iend] = sc.stepPointMC();
	    // make sure the trigter StepPoints also go in the StrawDigiMC
	    spmcs.insert(sc.stepPointMC());
          }
          // choose the minimum time from either end, as the ADC sums both
          double ptime = 1.0e10;
          for (size_t iend=0;iend<2;++iend){
            ptime = std::min(ptime,wetime[iend]);
          }
          // subtract a small buffer
          ptime -= _adcbuffer;
          // pickup all StepPointMCs associated with clusts inside the time window of the ADC digitizations (after the threshold)
          for (auto ih=wf[0].clusts().clustList().begin();ih!=wf[0].clusts().clustList().end();++ih){
            if (ih->time() >= ptime && ih->time() < ptime +
                ( strawele.nADCSamples()-strawele.nADCPreSamples())*strawele.adcPeriod())
              spmcs.insert(ih->stepPointMC());
          }
          vector<art::Ptr<StepPointMC> > stepMCs;
          stepMCs.reserve(spmcs.size());
          for(auto ispmc=spmcs.begin(); ispmc!= spmcs.end(); ++ispmc){
            stepMCs.push_back(*ispmc);
          }
          PtrStepPointMCVector mcptr;
          for(auto ixmcsp=xmcsp.begin();ixmcsp!=xmcsp.end();++ixmcsp)
            mcptr.push_back(*ixmcsp);
          mcptrs->push_back(mcptr);
          mcdigis->push_back(StrawDigiMC(sid,wetime,cpos,stepMC,stepMCs));
          // diagnostics
          if(_diag > 1)digiDiag(strawphys,wf,xpair,digis->back(),mcdigis->back());
        }
      }
    }

    bool StrawDigisFromStepPointMCs::createDigi(StrawElectronics const& strawele, WFXP const& xpair, SWFP const& waveform,
                                                StrawId sid, StrawDigiCollection* digis){
      // initialize the float variables that we later digitize
      TDCTimes xtimes = {0.0,0.0};
      TrkTypes::TOTValues tot;
      // get the ADC sample times from the electroincs.  Use the cal side time to randomize
      // the phase, this doesn't really matter
      TrkTypes::ADCTimes adctimes;
      strawele.adcTimes(xpair[0]._time,adctimes);
      //  sums voltages from both waveforms for ADC
      ADCVoltages wf[2];
      // add the jitter in the EventWindowMarker time for this Panel (constant for a whole microbunch, same for both sides)
      double dt = _ewMarkerROCdt[sid.getPanel()];
      // loop over the associated crossings
      for(size_t iend = 0;iend<2; ++iend){
        WFX const& wfx = xpair[iend];
        // record the crossing time for this end, including clock jitter  These already include noise effects
        // add noise for TDC on each side
        double tdc_jitter = _randgauss.fire(0.0,strawele.TDCResolution());
        xtimes[iend] = wfx._time+dt+tdc_jitter;
        // randomize threshold using the incoherent noise
        double threshold = _randgauss.fire(wfx._vcross,strawele.analogNoise(StrawElectronics::thresh));
        // find TOT
        tot[iend] = waveform[iend].digitizeTOT(strawele,threshold,wfx._time + dt);
        // sample ADC
        waveform[iend].sampleADCWaveform(strawele,adctimes,wf[iend]);
      }
      // uncalibrate
      strawele.uncalibrateTimes(xtimes,sid);
      // add ends and add noise
      ADCVoltages wfsum; wfsum.reserve(adctimes.size());
      for(unsigned isamp=0;isamp<adctimes.size();++isamp){
        wfsum.push_back(wf[0][isamp]+wf[1][isamp]+_randgauss.fire(0.0,strawele.analogNoise(StrawElectronics::adc)));
      }
      // digitize, and make final test.  This call includes the clock error WRT the proton pulse
      TrkTypes::TDCValues tdcs;
      bool digitize;
      if(readAll(sid))
        digitize = strawele.digitizeAllTimes(xtimes,_mbtime,tdcs);
      else
        digitize = strawele.digitizeTimes(xtimes,tdcs);

      if(digitize){
        TrkTypes::ADCWaveform adc;
        strawele.digitizeWaveform(sid,wfsum,adc);
        // create the digi from this
        digis->push_back(StrawDigi(sid,tdcs,tot,adc));
      }
      return digitize;
    }

    // find straws which couple to the given one, and record them and their couplings in XTalk objects.
    // For now, this is just a fixed number for adjacent straws,
    // the couplings and straw identities should eventually come from a database, FIXME!!!
    void StrawDigisFromStepPointMCs::findCrossTalkStraws(Straw const& straw, vector<XTalk>& xtalk) {
      StrawId selfid = straw.id();
      xtalk.clear();
      // find straws sensitive to straw-to-straw cross talk
      vector<StrawId> const& strawNeighbors = straw.nearestNeighboursById();
      // find straws sensitive to electronics cross talk
      vector<StrawId> const& preampNeighbors = straw.preampNeighboursById();
      // convert these to cross-talk
      for(auto isid=strawNeighbors.begin();isid!=strawNeighbors.end();++isid){
        xtalk.push_back(XTalk(selfid,*isid,_preampxtalk,0));
      }
      for(auto isid=preampNeighbors.begin();isid!=preampNeighbors.end();++isid){
        xtalk.push_back(XTalk(selfid,*isid,0,_postampxtalk));
      }
    }

    // functions that need implementing:: FIXME!!!!!!
    // Could also fold in beam-off random trigger hits from real data
    void StrawDigisFromStepPointMCs::addNoise(StrawClusterMap& hmap){
      // create random noise clusts and add them to the sequences of random straws.
    }
    // diagnostic functions
    void StrawDigisFromStepPointMCs::waveformHist(StrawElectronics const& strawele, SWFP const& wfs, WFXPList const& xings) {
      // histogram individual waveforms
      static unsigned nhist(0);// maximum number of histograms per job!
      for(size_t iend=0;iend<2;++iend){
        // step to the 1st cluster past the blanking time to avoid double-counting
        ClusterList const& clist = wfs[iend].clusts().clustList();
        auto icl = clist.begin();
        while(icl->time() < strawele.flashEnd())
          icl++;
        if(icl != clist.end() && nhist < _maxhist && xings.size() >= _minnxinghist &&
           ( ((!_xtalkhist) && wfs[iend].xtalk().self()) || (_xtalkhist && !wfs[iend].xtalk().self()) ) ) {
          double tstart = icl->time()-_tstep;
          double tfall = strawele.fallTime(_diagpath);
          double tend = clist.rbegin()->time() + _nfall*tfall;
          ADCTimes times;
          ADCVoltages volts;
          times.reserve(size_t(ceil(tend-tstart)/_tstep));
          volts.reserve(size_t(ceil(tend-tstart)/_tstep));
          double t = tstart;
          while(t<tend){
            times.push_back(t);
            volts.push_back(wfs[iend].sampleWaveform(strawele,_diagpath,t));
            t += _tstep;
          }
          ++nhist;
          art::ServiceHandle<art::TFileService> tfs;
          char name[60];
          char title[100];
          snprintf(name,60,"SWF%i_%i",wfs[iend].clusts().strawId().asUint16(),nhist);
          snprintf(title,100,"Electronic output for straw %i end %i path %i;time (nSec);Waveform (mVolts)",wfs[iend].clusts().strawId().asUint16(),(int)iend,_diagpath);
          TH1F* wfh = tfs->make<TH1F>(name,title,volts.size(),times.front(),times.back());
          for(size_t ibin=0;ibin<times.size();++ibin)
            wfh->SetBinContent(ibin+1,volts[ibin]);
          TList* flist = wfh->GetListOfFunctions();
          for(auto ixing=xings.begin();ixing!=xings.end();++ixing){
            TMarker* smark = new TMarker(ixing->at(iend)._time,ixing->at(iend)._vcross,8);
            smark->SetMarkerColor(kGreen);
            smark->SetMarkerSize(2);
            flist->Add(smark);
          }
        }
      }
    }

    void StrawDigisFromStepPointMCs::waveformDiag(
                     StrawElectronics const& strawele,
                     SWFP const& wfs, WFXPList const& xings) {
      const Tracker& tracker = *GeomHandle<Tracker>();
      const Straw& straw = tracker.getStraw( wfs[0].clusts().strawId() );
      _swplane = straw.id().getPlane();
      _swpanel = straw.id().getPanel();
      _swlayer = straw.id().getLayer();
      _swstraw = straw.id().getStraw();
      for(size_t iend=0;iend<2; ++iend){
        ClusterList const& clusts = wfs[iend].clusts().clustList();
        size_t nclust = clusts.size();
        set<art::Ptr<StepPointMC> > steps;
        set<art::Ptr<SimParticle> > parts;
        _nxing[iend] = 0;
        _txing[iend] = strawele.flashStart() + _mbbuffer;
        _xddist[iend] = _xwdist[iend] = _xpdist[iend] = -1.0;
        for(auto ixing=xings.begin();ixing!=xings.end();++ixing){
          ++_nxing[iend];
          _txing[iend] = min(_txing[iend],static_cast<float_t>(ixing->at(iend)._time));
          _xddist[iend] = ixing->at(iend)._iclust->driftDistance();
          _xwdist[iend] = ixing->at(iend)._iclust->wireDistance();
          // compute direction perpendicular to wire and momentum
          art::Ptr<StepPointMC> const& spp = ixing->at(iend)._iclust->stepPointMC();
          if(!spp.isNull()){
            Hep3Vector pdir = straw.getDirection().cross(spp->momentum()).unit();
            // project the differences in position to get the perp distance
            _xpdist[iend] = pdir.dot(spp->position()-straw.getMidPoint());
          }
        }
        if(_nxing[iend] == 0){
          // no xings: just take the 1st clust
          if(nclust > 0 ){
            _xddist[iend] = clusts.front().driftDistance();
            _xwdist[iend] = clusts.front().wireDistance();
            art::Ptr<StepPointMC> const& spp = clusts.front().stepPointMC();
            if(!spp.isNull()){
              Hep3Vector pdir = straw.getDirection().cross(spp->momentum()).unit();
              // project the differences in position to get the perp distance
              _xpdist[iend] = pdir.dot(spp->position()-straw.getMidPoint());
            }
          }
        }
        if(nclust > 0){
          _tmin[iend] = clusts.begin()->time();
          _tmax[iend] = clusts.rbegin()->time();
        } else {
          _tmin[iend] = _mbtime+_mbbuffer;
          _tmax[iend] = -100.0;
        }

        _hqsum[iend] = 0.0;
        _vmax[iend] = _tvmax[iend] = 0.0;
        _wmcpdg[iend] = _wmcproc[iend] = 0;
        for(auto iclu=clusts.begin();iclu!=clusts.end();++iclu){
          if(iclu->stepPointMC().isNonnull()){
            steps.insert(iclu->stepPointMC());
            parts.insert(iclu->stepPointMC()->simParticle());
            _hqsum[iend] += iclu->charge();
            double ctime = iclu->time()+strawele.maxResponseTime(_diagpath,iclu->wireDistance());
            double vout = wfs[iend].sampleWaveform(strawele,_diagpath,ctime);
            if(vout > _vmax[iend]){
              _vmax[iend] = vout;
              _tvmax[iend] = ctime;
              _wmcpdg[iend] = iclu->stepPointMC()->simParticle()->pdgId();
              _wmcproc[iend] = iclu->stepPointMC()->simParticle()->creationCode();
              _mce[iend] = iclu->stepPointMC()->simParticle()->startMomentum().e();
              _slen[iend] = iclu->stepPointMC()->stepLength();
              _sedep[iend] = iclu->stepPointMC()->ionizingEdep();
            }
          }
        }
        _nsteppoint[iend] = steps.size();
        _npart[iend] = parts.size();
        _sesum[iend] = 0.0;
        for(auto istep=steps.begin();istep!=steps.end();++istep)
          _sesum [iend]+= (*istep)->ionizingEdep();
      }
      _swdiag->Fill();
    }

    void StrawDigisFromStepPointMCs::digiDiag(StrawPhysics const& strawphys, SWFP const& wfs, WFXP const& xpair, StrawDigi const& digi,StrawDigiMC const& mcdigi) {
      const Tracker& tracker = *GeomHandle<Tracker>();
      const Straw& straw = tracker.getStraw( digi.strawId() );
      _sdplane = straw.id().getPlane();
      _sdpanel = straw.id().getPanel();
      _sdlayer = straw.id().getLayer();
      _sdstraw = straw.id().getStraw();

      for(size_t iend=0;iend<2;++iend){
        _xtime[iend] = xpair[iend]._time;
        _tctime[iend] = xpair[iend]._iclust->time();
        _charge[iend] = xpair[iend]._iclust->charge();
        _ddist[iend] = xpair[iend]._iclust->driftDistance();
        _dtime[iend] = xpair[iend]._iclust->driftTime();
        _ptime[iend] = xpair[iend]._iclust->propTime();
        _phi[iend] = xpair[iend]._iclust->phi(); //JB
        _wdist[iend] = xpair[iend]._iclust->wireDistance();
        _vstart[iend] = xpair[iend]._vstart;
        _vcross[iend] = xpair[iend]._vcross;
        _tdc[iend] = digi.TDC(xpair[iend]._iclust->strawEnd());
        _tot[iend] = digi.TOT(xpair[iend]._iclust->strawEnd());
        ClusterList const& clist = wfs[iend].clusts().clustList();
        auto ctrig = xpair[iend]._iclust;
        _ncludd[iend] = clist.size();
        // find the earliest cluster from the same particle that triggered the crossing
        auto iclu = clist.begin();
        while( iclu != clist.end() && ctrig->stepPointMC()->simParticle() != iclu->stepPointMC()->simParticle() ){
          ++iclu;
        }
        if(iclu != clist.end() ){
          _ectime[iend] = iclu->time();
          _ecddist[iend] = iclu->driftDistance();
          _ecdtime[iend] = iclu->driftTime();
          _ecptime[iend] = iclu->propTime();
          // count how many clusters till we get to the trigger cluster
          size_t iclust(0);
          while( iclu != clist.end() && iclu != ctrig){
            ++iclu;
            ++iclust;
          }
          _iclust[iend] = iclust;
        }
      }
      if(xpair[0]._iclust->stepPointMC() == xpair[1]._iclust->stepPointMC())
        _nstep = 1;
      else
        _nstep = 2;
      _adc.clear();
      for(auto iadc=digi.adcWaveform().begin();iadc!=digi.adcWaveform().end();++iadc){
        _adc.push_back(*iadc);
      }
      // mc truth information
      _dmcpdg = _dmcproc = _dmcgen = 0;
      _dmcmom = -1.0;
      _mctime = _mcenergy = _mctrigenergy = _mcthreshenergy = _mcdca = -1000.0;
      _mcthreshpdg = _mcthreshproc = _mcnstep = 0;
      art::Ptr<StepPointMC> const& spmc = xpair[0]._iclust->stepPointMC();
      if(!spmc.isNull()){
        _mctime = _toff.timeWithOffsetsApplied(*spmc);
        // compute the doca for this step
        TwoLinePCA pca( straw.getMidPoint(), straw.getDirection(),
            spmc->position(), spmc->momentum().unit() );
        _mcdca = pca.dca();

        Hep3Vector mccdir = (pca.point2()-straw.getMidPoint());
        mccdir -= straw.getDirection()*(mccdir.dot(straw.getDirection()));
        _mcdcaphi = mccdir.theta();
        _mcdcadtime = strawphys.driftDistanceToTime(_mcdca,_mcdcaphi); //JB: this is now from the lorentz corrected r-component of the drift

        if(!spmc->simParticle().isNull()){
          _dmcpdg = spmc->simParticle()->pdgId();
          _dmcproc = spmc->simParticle()->creationCode();
          if(spmc->simParticle()->genParticle().isNonnull())
            _dmcgen = spmc->simParticle()->genParticle()->generatorId().id();
          _dmcmom = spmc->momentum().mag();
        }
      }
      _mcenergy = mcdigi.energySum();
      _mctrigenergy = mcdigi.triggerEnergySum(StrawEnd::cal);
      // sum the energy from the explicit trigger particle, and find it's releationship
      _mcthreshenergy = 0.0;
      _mcnstep = mcdigi.stepPointMCs().size();
      art::Ptr<StepPointMC> threshpart = mcdigi.stepPointMC(StrawEnd::cal);
      if(threshpart.isNull()) threshpart = mcdigi.stepPointMC(StrawEnd::hv);
      for(auto imcs = mcdigi.stepPointMCs().begin(); imcs!= mcdigi.stepPointMCs().end(); ++ imcs){
        // if the SimParticle for this step is the same as the one which fired the discrim, add the energy
        if( (*imcs)->simParticle() == threshpart->simParticle() )
          _mcthreshenergy += (*imcs)->eDep();
      }
      _mcthreshpdg = threshpart->simParticle()->pdgId();
      _mcthreshproc = threshpart->simParticle()->creationCode();

      _xtalk = digi.strawId() != spmc->strawId();
      // fill the tree entry
      _sddiag->Fill();
    }//End of digiDiag




    void StrawDigisFromStepPointMCs::fillClusterPositions(Straw const& straw, StepPointMC const& step, std::vector<Hep3Vector>& cpos) {
      // basic info
      double charge(0.0);
      GlobalConstantsHandle<ParticleDataTable> pdt;
      if(pdt->particle(step.simParticle()->pdgId()).isValid()){
        charge = pdt->particle(step.simParticle()->pdgId()).ref().charge();
      }
      static const double r2 = straw.innerRadius()*straw.innerRadius();
      // decide how we step; straight or helix, depending on the Pt
      Hep3Vector const& mom = step.momentum();
      Hep3Vector mdir = mom.unit();
      // approximate pt
      double apt = step.momentum().perpPart(_bdir).mag();
      if( apt > _ptmin) { // use linear approximation
        double slen = step.stepLength();
        // make sure this linear approximation doesn't extend past the physical straw
        Hep3Vector dperp = (step.position() -straw.getMidPoint()).perpPart(straw.getDirection());
        Hep3Vector mperp = mdir.perpPart(straw.getDirection());
        double dm = dperp.dot(mperp);
        double m2 = mperp.mag2();
        double dp2 = dperp.mag2();
        double sarg = dm*dm + m2*(r2 - dp2);
        // some glancing cases fail the linear math
        if(sarg > 0.0 && m2 > 0.0){
          double smax = (-dm + sqrt(sarg))/m2;
          slen = std::min(smax,slen);
        }
        // generate random cluster positions
        for(unsigned ic=0;ic < cpos.size();++ic){
          //
          cpos[ic] = step.position() +_randflat.fire(slen) *mdir;
        }
      } else {
        // Use a helix to model particles which curl on the scale of the straws
        GeomHandle<BFieldManager> bfmgr;
        GeomHandle<DetectorSystem> det;
        // find the local field vector at this step
        Hep3Vector vpoint_mu2e = det->toMu2e(step.position());
        Hep3Vector bf = bfmgr->getBField(vpoint_mu2e);
        // compute transverse radius of particle
        double rcurl = fabs(charge*(mom.perpPart(bf).mag())/BField::mmTeslaToMeVc*bf.mag());
        // basis using local Bfield direction
        Hep3Vector bfdir = bf.unit();
        Hep3Vector qmdir = (charge*mom).unit(); // charge-signed momentum direction
        Hep3Vector rdir = qmdir.cross(bfdir).unit(); // points along magnetic force, ie towards center
        Hep3Vector pdir = bfdir.cross(rdir).unit(); // perp to this and field
        // find the center of the helix at the start of this step
        Hep3Vector hcent = step.position() + rcurl*rdir;
        // find helix parameters.  By definition, z0 = phi0 = 0
        double omega = qmdir.dot(pdir)/(rcurl*qmdir.dot(bfdir));
        // compute how far this step goes along the field direction.  This includes sign information
        double zlen = step.stepLength()*mdir.dot(bfdir);
        // loop until we've found enough valid samples, or have given up trying
        unsigned iclu(0);
        unsigned ntries(0);
        unsigned nclus = cpos.size();
        while(iclu < nclus && ntries < 10*nclus){
          double zclu = _randflat.fire(zlen);
          double phi = zclu*omega;
          // build cluster position from these
          Hep3Vector cp = hcent + rcurl*(-rdir*cos(phi) + pdir*sin(phi)) + zclu*bfdir;
          // test
          double rd2 = (cp-straw.getMidPoint()).perpPart(straw.getDirection()).mag2();
          if(rd2 - r2 < 1.0e-3){
            cpos[iclu] = cp;
            ++iclu;
          } else if (_debug > 0) {
            cout << "cluster outside straw: helix " <<  sqrt(rd2) << endl;
          }
          ++ntries;
        }
        if(iclu != nclus){
          // failed to find valid steps. put any remining energy at the step
          if(_debug > 0)cout << "mu2e::StrawDigisFromStepPointMCs: Couldn't create enough clusters : "<< iclu << " wanted " << nclus << endl;
          while(iclu < nclus){
            cpos[iclu] = step.position();
            ++iclu;
          }
        }
      }
    }

    void StrawDigisFromStepPointMCs::fillClusterMinion(StrawPhysics const& strawphys, StepPointMC const& step, std::vector<unsigned>& ne, std::vector<double>& cen) {
      // Loop until we've assigned energy + electrons to every cluster
      unsigned mc(0);
      double esum(0.0);
      double etot = step.ionizingEdep();
      unsigned nc = ne.size();
      while(mc < nc){
        std::vector<unsigned> me(nc);
        // fill an array of random# of electrons according to the measured distribution.  These are returned sorted lowest-highest.
        fillClusterNe(strawphys,me);
        for(auto ie : me) {
          // maximum energy for this cluster requires at least 1 electron for the rest of the cluster
          double emax = etot - esum - (nc -mc -1)*strawphys.ionizationEnergy((unsigned)1);
          double eele = strawphys.ionizationEnergy(ie);
          if( eele < emax){
            ne[mc] = ie;
            cen[mc] = eele;
            ++mc;
            esum += eele;
          } else {
            break;
          }
        }
      }
      // distribute any residual energy randomly to these clusters.  This models delta rays
      unsigned ns;
      do{
        unsigned me = strawphys.nePerIon(_randflat.fire());
        double emax = etot - esum;
        double eele = strawphys.ionizationEnergy(me);
        if(eele < emax){
          // choose a random cluster to assign this energy to
          unsigned mc = std::min(nc-1,static_cast<unsigned>(floor(_randflat.fire(nc))));
          ne[mc] += me;
          cen[mc] += eele;
          esum += eele;
        }
        // maximum energy for this cluster requires at least 1 electron for the rest of the cluster
        ns = static_cast<unsigned>(floor((etot-esum)/strawphys.ionizationEnergy((unsigned)1)));
      } while(ns > 0);
    }

    void StrawDigisFromStepPointMCs::fillClusterNe(StrawPhysics const& strawphys,std::vector<unsigned>& me) {
      for(size_t ie=0;ie < me.size(); ++ie){
        me[ie] = strawphys.nePerIon(_randflat.fire());
      }
      if(_sort)std::sort(me.begin(),me.end());
    }

    bool StrawDigisFromStepPointMCs::readAll(StrawId const& sid) const {

      return sid.straw() >= _allStraw &&
        (std::find(_allPlanes.begin(),_allPlanes.end(),sid.plane()) != _allPlanes.end());
    }

  } // end namespace trackermc
} // end namespace mu2e

using mu2e::TrackerMC::StrawDigisFromStepPointMCs;
DEFINE_ART_MODULE(StrawDigisFromStepPointMCs);
