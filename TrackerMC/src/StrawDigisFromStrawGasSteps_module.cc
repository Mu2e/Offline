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
#include "art/Framework/Principal/Run.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Core/EDProducer.h"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "art_root_io/TFileService.h"
#include "Offline/SeedService/inc/SeedService.hh"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
// conditions
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/TrackerConditions/inc/StrawElectronics.hh"
#include "Offline/TrackerConditions/inc/StrawPhysics.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "Offline/TrackerConditions/inc/TrackerStatus.hh"
#include "BTrk/BField/BField.hh"
// utiliities
#include "Offline/Mu2eUtilities/inc/TwoLinePCA.hh"
#include "Offline/DataProducts/inc/TrkTypes.hh"
// persistent data
#include "Offline/DataProducts/inc/EventWindowMarker.hh"
#include "Offline/MCDataProducts/inc/ProtonBunchTimeMC.hh"
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
// temporary MC structures
#include "Offline/TrackerMC/inc/StrawClusterSequencePair.hh"
#include "Offline/TrackerMC/inc/StrawWaveform.hh"
#include "Offline/TrackerMC/inc/IonCluster.hh"
#include "Offline/TrackerMC/inc/StrawPosition.hh"
#include <Offline/TrackerMC/inc/StrawDigiBundle.hh>
#include <Offline/TrackerMC/inc/StrawDigiBundleCollection.hh>
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
using EventIDCollection = std::vector<art::EventID>;

namespace mu2e {
  namespace TrackerMC {
    using namespace TrkTypes;

    struct WireCharge { // charge at the wire after drift
      float _charge; // charge at the wire, in units of pC
      float _time; // relative time at the wire, relative to ionzation time (ns)
      StrawCoordinates _pos; // cylindrical coordinates WRT straw
    };

    struct WireEndCharge { // charge at one end of the wire after propagation
      float _charge; // charge at the wire, in units of pC
      float _time; // time at the wire end, relative to the time the charge reached the wire (ns)
      float _wdist; // propagation distance from the point of collection to the end
    };

    class StrawDigisFromStrawGasSteps : public art::EDProducer {
      public:
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;

        struct Config {
          fhicl::Atom<int> debug{ Name("debugLevel"), Comment("Debug Level"), 0};
          fhicl::Atom<int> diag{ Name("diagLevel"), Comment("Diag Level"), 0};
          fhicl::Atom<int> print{ Name("printLevel"), Comment("Print Level"), 0};
          fhicl::Atom<unsigned> maxhist{ Name("MaxHist"), Comment("Maximum number of waveform histograms"), 100};
          fhicl::Atom<bool> xtalkhist{ Name("CrossTalkHist"), Comment("Histogram of cross-talk"), false};
          fhicl::Atom<int> minnxinghist{ Name("MinNXingHist"), Comment("Minimum # of crossings to histogram waveform"),1};
          fhicl::Atom<float> tstep { Name("WaveformStep"), Comment("WaveformStep (nsec)"),0.1 };
          fhicl::Atom<float> tfall{ Name("WaveformTail"), Comment("Time past last signal to record waveform (ns)"),220.0};
          fhicl::Atom<unsigned> maxnclu{ Name("MaxNClusters"), Comment("Maximum number of clusters for non-minion steps"), 20};
          fhicl::Atom<bool> addXtalk{ Name("addCrossTalk"), Comment("Should we add cross talk hits?"),false };
          fhicl::Atom<bool> drift1e{ Name("DriftSingleElectrons"), Comment("Always drift single electrons"),false };
          fhicl::Atom<bool> randrad{ Name("RandomizeRadius"), Comment("Randomize the drift by step width effects"),true };
          fhicl::Atom<float> ctMinCharge{ Name("xtalkMinimumCharge"), Comment("minimum charge to add cross talk (for performance issues)") ,0};
          fhicl::Atom<bool> addNoise{ Name("addNoise"), Comment("should we add noise hits? NOT CURRENTLY IMPLEMENTED FIXME!"),false };
          fhicl::Atom<float> preampxtalk{ Name("preAmplificationCrossTalk"), Comment("Pre-amplification (straw) X-talk coupling"), 0.0 };
          fhicl::Atom<float> postampxtalk{ Name("postAmplificationCrossTalk"), Comment("Post-amplification (board) X-talk coupling") ,0.02};
          fhicl::Atom<float> minstepE{ Name("minstepE"), Comment(" minimum step energy depostion to turn into a straw signal (MeV)"),2.0e-6 };
          fhicl::Atom<art::InputTag> ewMarkerTag{ Name("EventWindowMarker"), Comment("EventWindowMarker producer"),"EWMProducer" };
          fhicl::Atom<art::InputTag> pbtmcTag{ Name("ProtonBunchTimeMC"), Comment("ProtonBunchTimeMC producer"),"EWMProducer" };
          fhicl::Atom<float> steptimebuf{ Name("StrawGasStepTimeBuffer"), Comment("buffer for MC step point times (nsec) ") ,100.0 };
          fhicl::Atom<float> flashBuffer{ Name("FlashTimeBuffer"), Comment("buffer for flash blanking times (nsec) ") ,10.0 };
          fhicl::Atom<float> tdcbuf{ Name("TDCTimeBuffer"), Comment("buffer for TDC jitter (nsec) ") ,2.0 };
          fhicl::Atom<uint16_t> allStraw{ Name("AllHitsStraw"), Comment("minimum straw # to read all hits") ,90};
          fhicl::Sequence<uint16_t> allPlanes{ Name("AllHitsPlanes"), Comment("planes to read all hits"), std::vector<uint16_t>{} };
          fhicl::Atom<int> diagpath{ Name("DiagPath"), Comment("Digitization Path for waveform diagnostics") ,0 };
          fhicl::Atom<string> spinstance { Name("StrawGasStepInstance"), Comment("StrawGasStep Instance name"),""};
          fhicl::Sequence<art::InputTag> spmodules { Name("StrawGasStepModules"), Comment("StrawGasStep Module names")};
          fhicl::Atom<bool> usestatus { Name("UseStatus"), Comment("Use TrackerStatus when making digis"), false};
          fhicl::Atom<art::InputTag> mixedDigisTag { Name("MixedDigisTag"), Comment("Source of digis to overlay event onto"), ""};
          fhicl::Atom<bool> mixDigiMCs { Name("MixDigiMCs"), Comment("Propagate mixed StrawDigiMCs through module"), false};
          fhicl::Atom<bool> allowEmptySteps { Name("AllowEmptyStrawGasSteps"), Comment("Allow digitization to proceed even without any valid straw gas step collections"), false};
        };

        typedef art::Ptr<StrawGasStep> SGSPtr;
        typedef art::Ptr<SimParticle> SPPtr;
        typedef map<StrawId,StrawClusterSequencePair> StrawClusterMap;  // clusts by straw
        // work with pairs of waveforms, one for each straw end
        typedef std::array<StrawWaveform,2> SWFP;
        typedef std::array<WFX,2> WFXP;
        typedef list<WFXP> WFXPList;
        typedef WFXPList::const_iterator WFXPI;

        using Parameters = art::EDProducer::Table<Config>;
        explicit StrawDigisFromStrawGasSteps(const Parameters& config);
        // Accept compiler written d'tor.

      private:

        void beginJob() override;
        void beginRun(art::Run& run) override;
        void produce(art::Event& e) override;

        // Diagnostics
        int _debug, _diag, _printLevel;
        unsigned _maxhist;
        bool  _xtalkhist;
        unsigned _minnxinghist;
        double _tstep, _tfall;
        // Parameters
        bool   _addXtalk, _drift1e, _randrad;
        double _ctMinCharge;
        bool   _addNoise;
        double _preampxtalk, _postampxtalk;// these should come from conditions, FIXME!!
        double _minstepE;
        art::InputTag _ewMarkerTag;
        art::InputTag _pbtmcTag;
        double _mbtime;
        double _mbbuffer;
        double _flashbuffer;
        double _adcbuffer;
        double _steptimebuf;
        double _tdcbuf;
        uint16_t _allStraw;
        std::vector<uint16_t> _allPlanes;
        unsigned _maxnclu;
        StrawElectronics::Path _diagpath;
        bool _usestatus;
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
        // digi mixing
        const art::InputTag _mixedDigisTag;
        const bool _mixDigiMCs;
        const bool _allowEmptySteps;
        // Proditions
        ProditionsHandle<StrawPhysics> _strawphys_h;
        ProditionsHandle<StrawElectronics> _strawele_h;
        ProditionsHandle<Tracker> _alignedTrackerSim_h{"Sim"};
        ProditionsHandle<TrackerStatus> _trackerStatus_h;
        const Tracker *_tracker;
        art::Selector _selector;
        double _rstraw; // cache
        // diagnostics
        TTree* _swdiag;
        Int_t _swplane, _swpanel, _swlayer, _swstraw, _ndigi;
        Float_t _hqsum[2], _vmax[2], _tvmax[2], _sesum[2];
        Int_t _wmcpdg[2], _wmcproc[2], _nxing[2], _nclu[2];
        Int_t _ngasstep[2], _npart[2];
        Float_t _mce[2], _slen[2], _sedep[2];
        Float_t _tmin[2], _tmax[2], _txing[2], _xddist[2], _xwdist[2], _xpdist[2];
        TTree* _sddiag;
        Int_t _sdplane, _sdpanel, _sdlayer, _sdstraw;
        Int_t _ncludd[2], _iclust[2];
        Int_t _nstep;
        Float_t _ectime[2], _ecddist[2], _ecdtime[2], _ecptime[2], _ecphi[2];
        Float_t _xtime[2], _tctime[2], _charge[2], _acharge[2], _ddist[2], _dtime[2], _ptime[2];
        Float_t _wdist[2], _vstart[2], _vcross[2];
        Float_t _phi[2];
        Int_t _provenance;
        Float_t _mcenergy, _mctrigenergy, _mcthreshenergy;
        Double_t _mctime;
        Int_t _mcthreshpdg, _mcthreshproc, _mcnstep;
        Float_t _mcdca, _mcdcaphi, _mcdcadtime;
        Float_t _mcstrawdca, _mcstrawdcaphi, _mcstrawdcadtime;
        Int_t _dmcpdg, _dmcproc, _dmcgen;
        Float_t _dmcmom;
        Bool_t _xtalk;
        vector<unsigned> _adc;
        unsigned _pmp;
        Int_t _tdc[2], _tot[2];
        Int_t _sdtype;
        TTree* _sdiag;
        Float_t _sdwidth, _sdlen;
        Float_t _steplen, _stepE, _qsum, _esum, _eesum, _qe, _partP, _steptime;
        Int_t _nclust, _netot, _partPDG, _stype;
        vector<IonCluster> _clusters;
        Float_t _pbtimemc;
        array<Float_t, StrawId::_nupanels> _ewMarkerROCdt;
        double _eventWindowLength;
        TDCValue _eventWindowEndTDC;
        bool _onSpill;
        double _digitizationEndFromMarker;

        //  helper functions
        void fillClusterMap(StrawPhysics const& strawphys,
            StrawElectronics const& strawele,
            art::Event const& event, StrawClusterMap & hmap);
        void addStep(StrawPhysics const& strawphys,
            StrawElectronics const& strawele,
            Straw const& straw,
            SGSPtr const& sgsptr,
            StrawClusterSequencePair& shsp);
        void divideStep(StrawPhysics const& strawphys,
            StrawElectronics const& strawele,
            Straw const& straw,
            StrawGasStep const& step,
            vector<IonCluster>& clusters);
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
            Straw const& straw,
            StrawClusterSequencePair const& hsp,
            XTalk const& xtalk,
            StrawDigiCollection* digis, StrawDigiADCWaveformCollection* digiadcs, StrawDigiMCCollection* mcdigis);
        void fillDigis(StrawPhysics const& strawphys,
            StrawElectronics const& strawele,
            WFXPList const& xings,SWFP const& swfp , StrawId sid,
            StrawDigiCollection* digis, StrawDigiADCWaveformCollection* digiadcs, StrawDigiMCCollection* mcdigis);
        bool createDigi(StrawElectronics const& strawele,WFXP const& xpair, SWFP const& wf, StrawId sid, StrawDigiCollection* digis,
            StrawDigiADCWaveformCollection* digiadcs, double &digitization_ready_time);
        void findCrossTalkStraws(Straw const& straw,vector<XTalk>& xtalk);
        void fillClusterNe(StrawPhysics const& strawphys,std::vector<unsigned>& me);
        void fillClusterPositions(StrawGasStep const& step, Straw const& straw, std::vector<StrawCoordinates>& cpos);
        void fillClusterMinion(StrawPhysics const& strawphys, StrawGasStep const& step, std::vector<unsigned>& me, std::vector<float>& cen);
        bool readAll(StrawId const& sid) const;
        // diagnostic functions
        void waveformHist(StrawElectronics const& strawele,
            SWFP const& wf, WFXPList const& xings);
        void waveformDiag(StrawElectronics const& strawele,
            SWFP const& wf, WFXPList const& xings);
        void digiDiag(StrawPhysics const& strawphys, SWFP const& wf, WFXP const& xpair, StrawDigi const& digi, StrawDigiADCWaveform const& digiadc, StrawDigiMC const& mcdigi);
        void stepDiag(StrawPhysics const& strawphys, StrawElectronics const& strawele, StrawGasStep const& sgs);
        StrawCoordinates strawCoordinates(XYZVectorF const& cpos, Straw const& straw) const;
        XYZVectorF strawCoordinatesToXYZ(StrawCoordinates const& cpos, Straw const& straw) const;
    };

    StrawDigisFromStrawGasSteps::StrawDigisFromStrawGasSteps(const Parameters& config) :
      EDProducer(config),
      _debug(config().debug()),
      _diag(config().diag()),
      _printLevel(config().print()),
      _maxhist(config().maxhist()),
      _xtalkhist(config().xtalkhist()),
      _minnxinghist(config().minnxinghist()),
      _tstep(config().tstep()),
      _tfall(config().tfall()),
      _addXtalk(config().addXtalk()),
      _drift1e(config().drift1e()),
      _randrad(config().randrad()),
      _ctMinCharge(config().ctMinCharge()),
      _addNoise(config().addNoise()),
      _preampxtalk(config().preampxtalk()),
      _postampxtalk(config().postampxtalk()),
      _minstepE(config().minstepE()),
      _ewMarkerTag(config().ewMarkerTag()),
      _pbtmcTag(config().pbtmcTag()),
      _flashbuffer(config().flashBuffer()),
      _steptimebuf(config().steptimebuf()),
      _tdcbuf(config().tdcbuf()),
      _allStraw(config().allStraw()),
      _allPlanes(config().allPlanes()),
      _maxnclu(config().maxnclu()),
      _diagpath(static_cast<StrawElectronics::Path>(config().diagpath())),
      _usestatus(config().usestatus()),
      // Random number distributions
      _engine(createEngine( art::ServiceHandle<SeedService>()->getSeed())),
      _randgauss( _engine ),
      _randflat( _engine ),
      _randexp( _engine),
      _randP( _engine),
      _messageCategory("HITS"),
      _firstEvent(true),      // Control some information messages.
      _mixedDigisTag(config().mixedDigisTag()),
      _mixDigiMCs(config().mixDigiMCs()),
      _allowEmptySteps(config().allowEmptySteps()),
      // This selector will select only data products with the given instance name.
      _selector{ art::ProductInstanceNameSelector(config().spinstance())}
      {
        vector<art::InputTag> tags = config().spmodules();
        if (0 < tags.size()){
          _selector = art::Selector(_selector && art::InputTagListSelector(tags.begin(), tags.end(), ""));
        }
        // Tell the framework what we consume.
        for (const auto& tag: tags){
          consumes<StrawGasStepCollection>(tag);
        }
        consumes<EventWindowMarker>(_ewMarkerTag);
        consumes<ProtonBunchTimeMC>(_pbtmcTag);
        // Tell the framework what we make.
        produces<StrawDigiCollection>();
        produces<StrawDigiADCWaveformCollection>();
        produces<StrawDigiMCCollection>();
      }

    void StrawDigisFromStrawGasSteps::beginJob(){

      if(_diag > 0){

        art::ServiceHandle<art::TFileService> tfs;
        _sdiag =tfs->make<TTree>("sdiag","StrawGasStep diagnostics");
        _sdiag->Branch("steplen",&_steplen,"steplen/F");
        _sdiag->Branch("stepE",&_stepE,"stepE/F");
        _sdiag->Branch("partP",&_partP,"partP/F");
        _sdiag->Branch("qsum",&_qsum,"qsum/F");
        _sdiag->Branch("esum",&_esum,"esum/F");
        _sdiag->Branch("eesum",&_eesum,"eesum/F");
        _sdiag->Branch("qe",&_qe,"qe/F");
        _sdiag->Branch("steptime",&_steptime,"steptime/F");
        _sdiag->Branch("nclust",&_nclust,"nclust/I");
        _sdiag->Branch("netot",&_netot,"netot/I");
        _sdiag->Branch("partPDG",&_partPDG,"partPDG/I");
        _sdiag->Branch("stepType",&_stype,"stepType/I");
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
        _swdiag->Branch("nstep",&_ngasstep,"nscal/I:nshv/I");
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
          _sddiag->Branch("ewlength",&_eventWindowLength,"ewlength/D");
          _sddiag->Branch("pbtimemc",&_pbtimemc,"pbtimemc/F");
          _sddiag->Branch("xtime",&_xtime,"xtimecal/F:xtimehv/F");
          _sddiag->Branch("tctime",&_tctime,"tctimecal/F:tctimehv/F");
          _sddiag->Branch("ectime",&_ectime,"ectimecal/F:ectimehv/F");
          _sddiag->Branch("ecdtime",&_ecdtime,"ecdtimecal/F:ecdtimehv/F");
          _sddiag->Branch("ecptime",&_ecptime,"ecptimecal/F:ecptimehv/F");
          _sddiag->Branch("ecphi",&_ecphi,"ecphical/F:ecphihv/F");
          _sddiag->Branch("charge",&_charge,"chargecal/F:chargehv/F");
          _sddiag->Branch("acharge",&_acharge,"achargecal/F:achargehv/F");
          _sddiag->Branch("wdist",&_wdist,"wdistcal/F:wdisthv/F");
          _sddiag->Branch("phi",&_phi,"phical/F:phihv/F");
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
          _sddiag->Branch("sdtype",&_sdtype,"sdtype/I");
          _sddiag->Branch("sdwidth",&_sdwidth,"sdwidth/F");
          _sddiag->Branch("sdlen",&_sdlen,"sdlen/F");
          _sddiag->Branch("adc",&_adc);
          _sddiag->Branch("pmp",&_pmp);
          _sddiag->Branch("provenance",&_provenance,"provenance/I");
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
          _sddiag->Branch("mcstrawdca",&_mcstrawdca,"mcstrawdca/F");
          _sddiag->Branch("mcstrawdcaphi",&_mcstrawdcaphi,"mcstrawdcaphi/F");
          _sddiag->Branch("mcstrawdcadtime",&_mcstrawdcadtime,"mcstrawdcadtime/F");
          _sddiag->Branch("mcpdg",&_dmcpdg,"mcpdg/I");
          _sddiag->Branch("mcproc",&_dmcproc,"mcproc/I");
          _sddiag->Branch("mcgen",&_dmcgen,"mcgen/I");
          _sddiag->Branch("mcmom",&_dmcmom,"mcmom/F");
          _sddiag->Branch("xtalk",&_xtalk,"xtalk/B");

        }
      }
    }

    void StrawDigisFromStrawGasSteps::beginRun( art::Run& run ){
      const Tracker& tracker = *GeomHandle<Tracker>();
      _rstraw = tracker.strawProperties()._strawInnerRadius;
      if ( _printLevel > 0 ) {
        auto const& strawphys = _strawphys_h.get(run.id());
        strawphys.print(cout);
      }
    }

    void StrawDigisFromStrawGasSteps::produce(art::Event& event) {
      if ( _printLevel > 1 ) cout << "StrawDigisFromStrawGasSteps: produce() begin; event " << event.id().event() << endl;
      static int ncalls(0);
      ++ncalls;

      // initialize "global" collection of digis
      StrawDigiBundleCollection bundles;

      // by default, fetch conditions of simulated run
      art::EventID proditions_key = event.id();

      // (optionally) fetch preexisting digis, onto which the current event
      // will be laid; in this case, want to query conditions associated with
      // the preexisting event
      if (!_mixedDigisTag.empty()){
        auto evid_handle = event.getHandle<EventIDCollection>(_mixedDigisTag);
        if (evid_handle->size() != 1){
          std::string msg = "# of mixed EventIDs != 1: ";
          msg += std::to_string(evid_handle->size()) + " of them";
          throw cet::exception("SIM") << msg << endl;
        }
        proditions_key = evid_handle->at(0);
        auto digi_handle = event.getHandle<StrawDigiCollection>(_mixedDigisTag);
        auto adcs_handle = event.getHandle<StrawDigiADCWaveformCollection>(_mixedDigisTag);
        // bundle up preexisting digi products, optionally including DigiMCs
        if (!_mixDigiMCs){
          bundles.Append(*digi_handle, *adcs_handle);
        }
        else{
          auto dgmcs_handle = event.getHandle<StrawDigiMCCollection>(_mixedDigisTag);
          bundles.Append(*digi_handle, *adcs_handle, *dgmcs_handle);
        }
      }

      // update conditions caches, etc
      StrawPhysics const& strawphys = _strawphys_h.get(proditions_key);
      StrawElectronics const& strawele = _strawele_h.get(proditions_key);
      _tracker = _alignedTrackerSim_h.getPtr(event.id()).get();
      _mbtime = GlobalConstantsHandle<PhysicsParams>()->getNominalDRPeriod();
      art::Handle<EventWindowMarker> ewMarkerHandle;
      event.getByLabel(_ewMarkerTag, ewMarkerHandle);
      const EventWindowMarker& ewMarker(*ewMarkerHandle);
      _eventWindowLength = ewMarker.eventLength();
      // this is the maximum TDC value that makes it into this event.
      // After this time, the next marker has arrived and later hits roll over into the next event
      _eventWindowEndTDC = strawele.tdcResponse( _eventWindowLength - strawele.electronicsTimeDelay());
      _onSpill = (ewMarker.spillType() == EventWindowMarker::SpillType::onspill);
      // for offspill events, we assume we digitize for the whole event length
      _digitizationEndFromMarker = strawele.digitizationEndFromMarker();
      if (!_onSpill)
        _digitizationEndFromMarker = _eventWindowLength;
      if (strawele.digitizationEndFromMarker() > _eventWindowLength){
        throw cet::exception("SIM")<<"mu2e::StrawDigisFromStrawGasSteps: digitization window extends past next event window marker" << endl;
      }
      art::Handle<ProtonBunchTimeMC> pbtmcHandle;
      event.getByLabel(_pbtmcTag, pbtmcHandle);
      const ProtonBunchTimeMC& pbtmc(*pbtmcHandle);
      _pbtimemc = pbtmc.pbtime_;
      // calculate event window marker jitter for this microbunch for each panel
      for (size_t i=0;i<StrawId::_nupanels;i++){
        _ewMarkerROCdt.at(i) = _randgauss.fire(0,strawele.eventWindowMarkerROCJitter());
      }
      // make the microbunch buffer long enough to get the full waveform
      _mbbuffer = (strawele.nADCSamples() - strawele.nADCPreSamples())*strawele.adcPeriod();
      _adcbuffer = 0.01*strawele.adcPeriod();
      // Containers to hold the output information.
      unique_ptr<StrawDigiCollection> digis(new StrawDigiCollection);
      unique_ptr<StrawDigiADCWaveformCollection> digiadcs(new StrawDigiADCWaveformCollection);
      unique_ptr<StrawDigiMCCollection> mcdigis(new StrawDigiMCCollection);
      // create the StrawCluster map
      // this is a map from straw ids to a list of all clusters on that straw from this event
      StrawClusterMap hmap;
      // fill this from the event
      fillClusterMap(strawphys,strawele,event,hmap);
      // add noise clusts
      if(_addNoise)addNoise(hmap);
      // loop over the clust sequences (i.e. loop over straws, and for each get their list of clusters)
      for(auto ihsp=hmap.begin();ihsp!= hmap.end();++ihsp){
        StrawClusterSequencePair const& hsp = ihsp->second;
        Straw const& straw = _tracker->getStraw(hsp.strawId());
        // create primary digis from this clust sequence
        XTalk self(hsp.strawId()); // this object represents the straws coupling to itself, ie 100%
        createDigis(strawphys,strawele,straw,hsp,self,digis.get(),digiadcs.get(),mcdigis.get());
        // if we're applying x-talk, look for nearby coupled straws
        if(_addXtalk) {
          // only apply if the charge is above a threshold
          double totalCharge = 0;
          for(auto ih=hsp.clustSequence(StrawEnd::cal).clustList().begin();ih!= hsp.clustSequence(StrawEnd::cal).clustList().end();++ih){
            totalCharge += ih->charge();
          }
          if( totalCharge > _ctMinCharge){
            vector<XTalk> xtalk;
            findCrossTalkStraws(straw,xtalk);
            for(auto ixtalk=xtalk.begin();ixtalk!=xtalk.end();++ixtalk){
              createDigis(strawphys,strawele,straw,hsp,*ixtalk,digis.get(),digiadcs.get(),mcdigis.get());
            }
          }
        }
      }
      // bundle up new digis in global collection
      bundles.Append(*digis, *digiadcs, *mcdigis);

      // resolve collisions between any preexisting and new digis
      StrawDigiBundleCollection resolved;
      bundles.ResolveCollisions(strawele, resolved);
      digis = resolved.GetStrawDigiPtrs();
      digiadcs = resolved.GetStrawDigiADCWaveformPtrs();
      mcdigis = resolved.GetStrawDigiMCPtrs();

      // store the digis in the event
      event.put(move(digis));
      event.put(move(digiadcs));
      // store MC truth match
      event.put(move(mcdigis));
      if ( _printLevel > 1 ) cout << "StrawDigisFromStrawGasSteps: produce() end" << endl;
      // Done with the first event; disable some messages.
      _firstEvent = false;

    } // end produce

    void StrawDigisFromStrawGasSteps::createDigis(
        StrawPhysics const& strawphys,
        StrawElectronics const& strawele,
        Straw const& straw,
        StrawClusterSequencePair const& hsp,
        XTalk const& xtalk,
        StrawDigiCollection* digis, StrawDigiADCWaveformCollection* digiadcs,
        StrawDigiMCCollection* mcdigis) {
      // instantiate waveforms for both ends of this straw
      SWFP waveforms  ={ StrawWaveform(straw,hsp.clustSequence(StrawEnd::cal),xtalk),
        StrawWaveform(straw,hsp.clustSequence(StrawEnd::hv),xtalk) };
      // find the threshold crossing points for these waveforms
      WFXPList xings;
      // find the threshold crossings
      findThresholdCrossings(strawele,waveforms,xings);
      // convert the crossing points into digis, and add them to the event data
      fillDigis(strawphys,strawele,xings,waveforms,xtalk._dest,digis,digiadcs,mcdigis);
    }

    void StrawDigisFromStrawGasSteps::fillClusterMap(StrawPhysics const& strawphys,
        StrawElectronics const& strawele,
        art::Event const& event, StrawClusterMap & hmap){
// get status if needed
      std::shared_ptr<const TrackerStatus> trackerStatus;
      if(_usestatus) {
        trackerStatus = _trackerStatus_h.getPtr(event.id());
      }
      // Get all of the tracker StrawGasStep collections from the event:
      typedef vector< art::Handle<StrawGasStepCollection> > HandleVector;
      HandleVector stepsHandles = event.getMany<StrawGasStepCollection>( _selector);
      // Informational message on the first event.
      if ( _firstEvent ) {
        mf::LogInfo log(_messageCategory);
        log << "StrawDigisFromStrawGasSteps::fillHitMap will use StrawGasSteps from: \n";
        for ( HandleVector::const_iterator i=stepsHandles.begin(), e=stepsHandles.end();
            i != e; ++i ){
          art::Provenance const& prov(*(i->provenance()));
          log  << "   " << prov.branchName() << "\n";
        }
      }
      if(stepsHandles.empty() && (!_allowEmptySteps)){
        throw cet::exception("SIM")<<"mu2e::StrawDigisFromStrawGasSteps: No StrawGasStep collections found for tracker" << endl;
      }

      // Loop over StrawGasStep collections
      for ( auto const& sgsch : stepsHandles) {
        StrawGasStepCollection const& steps(*sgsch);
        // Loop over the StrawGasSteps in this collection
        for(size_t isgs = 0; isgs < steps.size(); isgs++){
          auto const& sgs = steps[isgs];
          // lookup straw here, to avoid having to find the tracker for every step
          StrawId const & sid = sgs.strawId();
          if ( ((!_usestatus) || (!trackerStatus->noSignal(sid))) && sgs.ionizingEdep() > _minstepE){
            Straw const& straw = _tracker->getStraw(sid);
            auto sgsptr = SGSPtr(sgsch,isgs);
            // create a clust from this step, and add it to the clust map
            addStep(strawphys,strawele,straw,sgsptr,hmap[sid]);
          } else if(_debug > 0) {
            StrawStatus stat;
            if(_usestatus) stat = trackerStatus->strawStatus(sid);
            std::cout << "Suppressing sid " << sid << " Status " << stat << " energy " << sgs.ionizingEdep() << std::endl;
          }
        }
      }
    }

    void StrawDigisFromStrawGasSteps::addStep(StrawPhysics const& strawphys,
        StrawElectronics const& strawele,
        Straw const& straw,
        SGSPtr const& sgsptr,
        StrawClusterSequencePair& shsp) {
      auto const& sgs = *sgsptr;
      StrawId sid = sgs.strawId();
      // apply time offsets, and take module with MB
      double ctime  = microbunchTime(strawele,sgs.time());
      // test if this step point is roughly in the digitization window
      if( (ctime > strawele.digitizationStartFromMarker() - strawele.electronicsTimeDelay() - _steptimebuf
            && ctime <  max(_mbtime,_digitizationEndFromMarker) - strawele.electronicsTimeDelay() + _steptimebuf) || readAll(sid)) {
        // Subdivide the StrawGasStep into ionization clusters
        _clusters.clear();
        divideStep(strawphys,strawele,straw,sgs,_clusters);
        // check
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
            // time of this cluster was produced, including offset
            // compute the time the signal arrives at the wire end
            double gtime = ctime + wireq._time + weq._time;
            // create the clust
            StrawCluster clust(StrawCluster::primary,sid,end,(float)gtime,weq._charge,weq._wdist,wireq._pos,(float)wireq._time,(float)weq._time,sgsptr,(float)ctime);
            // add the clusts to the appropriate sequence.
            shsp.clustSequence(end).insert(clust);
            // if required, add a 'ghost' copy of this clust
            if (_onSpill)
              addGhosts(strawele,clust,shsp.clustSequence(end));
          }
        }
        if(_diag > 0) stepDiag(strawphys, strawele, sgs);
      }
    }

    void StrawDigisFromStrawGasSteps::divideStep(StrawPhysics const& strawphys,
        StrawElectronics const& strawele,
        Straw const& straw,
        StrawGasStep const& sgs,
        vector<IonCluster>& clusters) {
      // single cluster
      if (sgs.stepType().shape() == StrawGasStep::StepType::point || sgs.stepLength() < strawphys.meanFreePath()){
        float cen = sgs.ionizingEdep();
        float fne = cen/strawphys.meanElectronEnergy();
        unsigned ne = std::max( static_cast<unsigned>(_randP(fne)),(unsigned)1);
        auto spos = strawCoordinates(sgs.startPosition(),straw);
        if(_drift1e){
          for (size_t i=0;i<ne;i++){
            IonCluster cluster(spos,strawphys.ionizationCharge((unsigned)1),cen/(float)ne,1);
            clusters.push_back(cluster);
          }
        } else {
          IonCluster cluster(spos,strawphys.ionizationCharge(ne),cen,ne);
          clusters.push_back(cluster);
        }
      } else {
        // compute the number of clusters for this step from the mean free path
        double fnc = sgs.stepLength()/strawphys.meanFreePath();
        // use a truncated Poisson distribution; this keeps both the mean and variance physical
        unsigned nc = std::max(static_cast<unsigned>(_randP.fire(fnc)),(unsigned)1);
        // if not minion, limit the number of steps geometrically
        bool minion = (sgs.stepType().ionization()==StrawGasStep::StepType::minion);
        if(!minion )nc = std::min(nc,_maxnclu);
        // require clusters not exceed the energy sum required for single-electron clusters
        nc = std::min(nc,static_cast<unsigned>(floor(sgs.ionizingEdep()/strawphys.ionizationEnergy((unsigned)1))));
        if(nc>0){
          // generate random positions for the clusters
          std::vector<StrawCoordinates> cposv(nc);
          fillClusterPositions(sgs,straw,cposv);
          // generate electron counts and energies for these clusters: minion model is more detailed
          std::vector<unsigned> ne(nc);
          std::vector<float> cen(nc);
          if(minion){
            fillClusterMinion(strawphys,sgs,ne,cen);
          } else {
            // get Poisson distribution of # of electrons for the average energy
            double fne = sgs.ionizingEdep()/(nc*strawphys.meanElectronEnergy()); // average # of electrons/cluster for non-minion clusters
            for(unsigned ic=0;ic<nc;++ic){
              ne[ic] = static_cast<unsigned>(std::max(_randP.fire(fne),(long)1));
              cen[ic] = ne[ic]*strawphys.meanElectronEnergy(); // average energy per electron, works for large numbers of electrons
            }
          }
          // create the cluster objects
          for(unsigned ic=0;ic<nc;++ic){
            // compute phi for this cluster
            float ee = cen[ic]/(float)ne[ic];
            // let each electron drift separately (RBonvtre)
            if(_drift1e || minion){
              for (size_t ie=0;ie<ne[ic];ie++){
                IonCluster cluster(cposv[ic],strawphys.ionizationCharge((unsigned)1),ee,1);
                clusters.push_back(cluster);
              }
            } else {
              IonCluster cluster(cposv[ic],strawphys.ionizationCharge(ne[ic]),cen[ic],ne[ic]);
              clusters.push_back(cluster);
            }
          }
        }
      }
    }

    void StrawDigisFromStrawGasSteps::driftCluster(
        StrawPhysics const& strawphys,Straw const& straw,
        IonCluster const& cluster, WireCharge& wireq ) {
      // sample the gain for this cluster
      double gain = strawphys.clusterGain(_randgauss, _randflat, cluster._ne);
      wireq._charge = cluster._charge*(gain);
      // compute drift time for this cluster
      double dt = strawphys.driftDistanceToTime(cluster._pos._wirePosition.Rho(),cluster._pos._wirePosition.Phi()); // this is now from the lorentz corrected r-component of the drift
      wireq._pos = cluster._pos;
      wireq._time = _randgauss.fire(dt,strawphys.driftTimeSpread(cluster._pos._wirePosition.Rho()));
    }

    void StrawDigisFromStrawGasSteps::propagateCharge(
        StrawPhysics const& strawphys, Straw const& straw,
        WireCharge const& wireq, StrawEnd end, WireEndCharge& weq) {
      // compute distance to the appropriate end; note that the straw always points from HV to cal (Duke convention)
      double wlen = straw.halfLength(); // use the full length, not the active length
      if(end == StrawEnd::cal)
        weq._wdist = wlen - wireq._pos._wirePosition.Z();
      else
        weq._wdist = wlen + wireq._pos._wirePosition.Z();
      // split the charge
      weq._charge = 0.5*wireq._charge;
      weq._time = strawphys.propagationTime(weq._wdist);
    }

    double StrawDigisFromStrawGasSteps::microbunchTime(StrawElectronics const& strawele, double globaltime) const {
      // converts time from proton beam time (StrawGasStep time) to event window marker time
      double mbtime = globaltime + _pbtimemc;
      // only fold if simulating onspill events
      if (_onSpill){
        // fold time relative to MB frequency
        mbtime = fmod(mbtime,_mbtime);
        // keep the microbunch time contiguous
        if(mbtime < 0 ) mbtime += _mbtime;
      }
      return mbtime;
    }

    void StrawDigisFromStrawGasSteps::addGhosts(StrawElectronics const& strawele,StrawCluster const& clust,StrawClusterSequence& shs) {
      // add enough buffer to cover both the flash blanking and the ADC waveform
      // at this point cluster times are relative to marker and wrapped at 1695 (if onspill)
      // wrap from beginning of microbunch to times > 1695 to digitize ADCs for hits near end of event window
      if(clust.time() < _mbbuffer)
        shs.insert(StrawCluster(clust,_mbtime));
      // wrap from end of microbunch to negative time to digitize ADCs for hits at tdc time=0
      if(clust.time() > _mbtime - _mbbuffer) shs.insert(StrawCluster(clust,-_mbtime));
    }

    void StrawDigisFromStrawGasSteps::findThresholdCrossings(StrawElectronics const& strawele, SWFP const& swfp, WFXPList& xings){
      //randomize the threshold to account for electronics noise; this includes parts that are coherent
      // for both ends (coming from the straw itself)
      // Keep track of crossings on each end to keep them in sequence
      double strawnoise = _randgauss.fire(0,strawele.strawNoise());
      // add specifics for each end
      double thresh[2] = {_randgauss.fire(strawele.threshold(swfp[0].straw().id(),static_cast<StrawEnd::End>(0))+strawnoise,strawele.analogNoise(StrawElectronics::thresh)),
        _randgauss.fire(strawele.threshold(swfp[0].straw().id(),static_cast<StrawEnd::End>(1))+strawnoise,strawele.analogNoise(StrawElectronics::thresh))};
      // Initialize search when the electronics becomes enabled:
      double tstart =strawele.digitizationStartFromMarker() - _flashbuffer;
      // for reading all hits, make sure we start looking for clusters at the minimum possible cluster time
      // this accounts for deadtime effects from previous microbunches
      if(readAll(swfp[0].straw().id()))tstart = -strawele.deadTimeAnalog();
      WFXP wfx = {WFX(swfp[0],tstart),WFX(swfp[1],tstart)};
      // search for coherent crossings on both ends
      bool crosses[2];
      for(size_t iend=0;iend<2;++iend){
        crosses[iend] = swfp[iend].crossesThreshold(strawele,thresh[iend],wfx[iend]);
      }
      // loop until we hit the end of the waveforms.  Require both in time.  Buffer to account for eventual TDC jitter
      // this is a loose pre-selection, final selection is done at digitization
      while( crosses[0] && crosses[1] && std::max(wfx[0]._time,wfx[1]._time)
          < max(_digitizationEndFromMarker,_mbtime) - strawele.electronicsTimeDelay() + _tdcbuf){
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
            thresh[iend] = _randgauss.fire(strawele.threshold(swfp[0].straw().id(),static_cast<StrawEnd::End>(iend)),strawele.analogNoise(StrawElectronics::thresh));
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

    void StrawDigisFromStrawGasSteps::fillDigis(StrawPhysics const& strawphys,
        StrawElectronics const& strawele,
        WFXPList const& xings, SWFP const& wf,
        StrawId sid,
        StrawDigiCollection* digis, StrawDigiADCWaveformCollection* digiadcs,
        StrawDigiMCCollection* mcdigis ) {
      //
      Straw const& straw = _tracker->getStraw(sid);
      double digitization_ready_time = -9e9; //FIXME no deadtime for first hit of a microbunch
      // loop over crossings
      for(auto xpair : xings) {
        // create a digi from this pair.  This also performs a finial test
        // on whether the pair should make a digi
        if(createDigi(strawele,xpair,wf,sid,digis,digiadcs,digitization_ready_time)){
          // fill associated MC truth matching. Only count the same step once
          StrawDigiMC::SGSPA sgspa;
          StrawDigiMC::PA cpos;
          StrawDigiMC::FA ctime, wetime;
          for (size_t iend=0;iend<2;++iend){
            StrawCluster const& sc = *(xpair[iend]._iclust);
            wetime[iend] = sc.time();
            ctime[iend] = sc.cluTime();
            cpos[iend] = strawCoordinatesToXYZ(sc.cluPos(),straw);
            sgspa[iend] = sc.strawGasStep();
          }
          // choose the minimum time from either end, as the ADC sums both
          float ptime = 1.0e10;
          for (size_t iend=0;iend<2;++iend){
            ptime = std::min(ptime,wetime[iend]);
          }
          // subtract a small buffer
          ptime -= _adcbuffer;
          mcdigis->push_back(StrawDigiMC(sid,cpos,ctime,wetime,sgspa,DigiProvenance::Simulation));
          if(_diag > 1){
            digiDiag(strawphys,wf,xpair,digis->back(),digiadcs->back(),mcdigis->back());
          }
        }
      }
      // waveform diags
      if ( _diag > 1 && (wf[0].clusts().clustList().size() > 0 ||
            wf[1].clusts().clustList().size() > 0 ) ) {
        // waveform xing diagnostics
        _ndigi = digis->size();
        waveformDiag(strawele,wf,xings);
      }
    }

    bool StrawDigisFromStrawGasSteps::createDigi(StrawElectronics const& strawele, WFXP const& xpair, SWFP const& waveform,
        StrawId sid, StrawDigiCollection* digis, StrawDigiADCWaveformCollection* digiadcs, double &digitization_ready_time){
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
      double dt = _ewMarkerROCdt[sid.uniquePanel()];
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
      double digitize_time = std::max(xtimes[0],xtimes[1]);
      if (digitize_time < digitization_ready_time)
        return false;

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
        digitize = strawele.digitizeAllTimes(xtimes,tdcs,_eventWindowEndTDC);
      else
        digitize = strawele.digitizeTimes(xtimes,tdcs,_onSpill,_eventWindowEndTDC);

      if(digitize){
        TrkTypes::ADCWaveform adc;
        TrkTypes::ADCValue pmp;
        strawele.digitizeWaveform(sid,wfsum,adc,pmp);
        // create the digi from this
        digis->push_back(StrawDigi(sid,tdcs,tot,pmp));
        digiadcs->push_back(StrawDigiADCWaveform(adc));
        // update digital deadtime for this channel
        digitization_ready_time = digitize_time + strawele.deadTimeDigital();
      }
      return digitize;
    }

    // find straws which couple to the given one, and record them and their couplings in XTalk objects.
    // For now, this is just a fixed number for adjacent straws,
    // the couplings and straw identities should eventually come from a database, FIXME!!!
    void StrawDigisFromStrawGasSteps::findCrossTalkStraws(Straw const& straw, vector<XTalk>& xtalk) {
      StrawId selfid = straw.id();
      xtalk.clear();
      // find straws sensitive to straw-to-straw cross talk
      vector<StrawId> strawNeighbors;
      vector<StrawId> preampNeighbors;
      for(uint16_t istraw=0; istraw < StrawId::_nstraws; istraw++){
        StrawId nid(selfid.plane(),selfid.panel(),istraw);
        if(nid != selfid && selfid.nearestNeighbor(nid)) strawNeighbors.push_back(nid);
        if(nid != selfid && selfid.samePreamp(nid)) preampNeighbors.push_back(nid);
      }
      // find straws sensitive to electronics cross talk
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
    void StrawDigisFromStrawGasSteps::addNoise(StrawClusterMap& hmap){
      // create random noise clusts and add them to the sequences of random straws.
    }

    void StrawDigisFromStrawGasSteps::fillClusterPositions(StrawGasStep const& sgs, Straw const& straw, std::vector<StrawCoordinates>& cposv) {
      // generate a random position between the start and end points.
      XYZVectorF path = sgs.endPosition() - sgs.startPosition();
      for(auto& cpos : cposv) {
        XYZVectorF pos = sgs.startPosition() + _randflat.fire(1.0)*path;
        // randomize the position by width.  This needs to be 2-d to avoid problems at the origin
        if(_randrad){
          XYZVectorF sdir = XYZVectorF(straw.getDirection());
          XYZVectorF p1 = path.Cross(sdir).Unit();
          XYZVectorF p2 = path.Cross(p1).Unit();
          pos += p1*_randgauss.fire()*sgs.width();
          pos += p2*_randgauss.fire()*sgs.width();
        }
        cpos = strawCoordinates(pos,straw);
      }
    }

    void StrawDigisFromStrawGasSteps::fillClusterMinion(StrawPhysics const& strawphys, StrawGasStep const& step, std::vector<unsigned>& ne, std::vector<float>& cen) {
      // Loop until we've assigned energy + electrons to every cluster
      unsigned mc(0);
      double esum(0.0);
      double etot = step.ionizingEdep();
      unsigned nc = ne.size();
      while(mc < nc){
        std::vector<unsigned> me(nc);
        // fill an array of random# of electrons according to the measured distribution.
        fillClusterNe(strawphys,me);
        // loop through these as long as there's enough energy to have at least 1 electron in each cluster.  If not, re-throw the # of electrons/cluster for the remainder
        for(auto ie : me) {
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

    void StrawDigisFromStrawGasSteps::fillClusterNe(StrawPhysics const& strawphys,std::vector<unsigned>& me) {
      for(size_t ie=0;ie < me.size(); ++ie){
        me[ie] = strawphys.nePerIon(_randflat.fire());
      }
    }

    bool StrawDigisFromStrawGasSteps::readAll(StrawId const& sid) const {
      return sid.straw() >= _allStraw &&
        (std::find(_allPlanes.begin(),_allPlanes.end(),sid.plane()) != _allPlanes.end());
    }

    // diagnostic functions
    void StrawDigisFromStrawGasSteps::waveformDiag(
        StrawElectronics const& strawele,
        SWFP const& wfs, WFXPList const& xings) {
      const Straw& straw = _tracker->getStraw( wfs[0].clusts().strawId() );
      _swplane = straw.id().getPlane();
      _swpanel = straw.id().getPanel();
      _swlayer = straw.id().getLayer();
      _swstraw = straw.id().getStraw();
      for(size_t iend=0;iend<2; ++iend){
        StrawClusterList const& clusts = wfs[iend].clusts().clustList();
        size_t nclust = clusts.size();
        set<SGSPtr > steps;
        set<SPPtr > parts;
        _nxing[iend] = 0;
        _txing[iend] = _eventWindowLength + _mbbuffer;
        _xddist[iend] = _xwdist[iend] = _xpdist[iend] = -1.0;
        for(auto ixing=xings.begin();ixing!=xings.end();++ixing){
          ++_nxing[iend];
          _txing[iend] = min(_txing[iend],static_cast<float_t>(ixing->at(iend)._time));
          _xddist[iend] = ixing->at(iend)._iclust->driftDistance();
          _xwdist[iend] = ixing->at(iend)._iclust->wireDistance();
          // compute direction perpendicular to wire and momentum
          auto const& sgs = ixing->at(iend)._iclust->strawGasStep();
          if(!sgs.isNull()){
            Hep3Vector pdir = straw.getDirection().cross(GenVector::Hep3Vec(sgs->momentum())).unit();
            // project the differences in position to get the perp distance
            _xpdist[iend] = pdir.dot(sgs->position()-straw.getMidPoint());
          }
        }
        if(_nxing[iend] == 0){
          // no xings: just take the 1st clust
          if(nclust > 0 ){
            _xddist[iend] = clusts.front().driftDistance();
            _xwdist[iend] = clusts.front().wireDistance();
            auto const& sgs = clusts.front().strawGasStep();
            if(!sgs.isNull()){
              Hep3Vector pdir = straw.getDirection().cross(GenVector::Hep3Vec(sgs->momentum())).unit();
              // project the differences in position to get the perp distance
              _xpdist[iend] = pdir.dot(sgs->position()-straw.getMidPoint());
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
          if(iclu->strawGasStep().isNonnull()){
            steps.insert(iclu->strawGasStep());
            parts.insert(iclu->strawGasStep()->simParticle());
            _hqsum[iend] += iclu->charge();
            double ctime = iclu->time()+strawele.maxResponseTime(straw.id(),_diagpath,iclu->wireDistance());
            double vout = wfs[iend].sampleWaveform(strawele,_diagpath,ctime);
            if(vout > _vmax[iend]){
              _vmax[iend] = vout;
              _tvmax[iend] = ctime;
              _wmcpdg[iend] = iclu->strawGasStep()->simParticle()->pdgId();
              _wmcproc[iend] = iclu->strawGasStep()->simParticle()->creationCode();
              _mce[iend] = iclu->strawGasStep()->simParticle()->startMomentum().e();
              _slen[iend] = iclu->strawGasStep()->stepLength();
              _sedep[iend] = iclu->strawGasStep()->ionizingEdep();
            }
          }
        }
        _ngasstep[iend] = steps.size();
        _npart[iend] = parts.size();
        _sesum[iend] = 0.0;
        for(auto istep=steps.begin();istep!=steps.end();++istep)
          _sesum [iend]+= (*istep)->ionizingEdep();
      }
      _swdiag->Fill();
      // waveform histograms
      if(_diag > 2 )waveformHist(strawele,wfs,xings);
    }

    void StrawDigisFromStrawGasSteps::waveformHist(StrawElectronics const& strawele, SWFP const& wfs, WFXPList const& xings) {
      // histogram individual waveforms
      static unsigned nhist(0);// maximum number of histograms per job!
      for(size_t iend=0;iend<2;++iend){
        // step to the 1st cluster past the blanking time to avoid double-counting
        StrawClusterList const& clist = wfs[iend].clusts().clustList();
        auto icl = clist.begin();
        while(icl->time() < strawele.digitizationStartFromMarker())
          icl++;
        if(icl != clist.end() && nhist < _maxhist && xings.size() >= _minnxinghist &&
            ( ((!_xtalkhist) && wfs[iend].xtalk().self()) || (_xtalkhist && !wfs[iend].xtalk().self()) ) ) {
          double tstart = icl->time()-_tstep;
          double tend = clist.rbegin()->time() + _tfall;
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

    void StrawDigisFromStrawGasSteps::digiDiag(StrawPhysics const& strawphys, SWFP const& wfs, WFXP const& xpair,
        StrawDigi const& digi, StrawDigiADCWaveform const& digiadc, StrawDigiMC const& mcdigi) {
      const Straw& straw = _tracker->getStraw( digi.strawId() );
      _sdplane = straw.id().getPlane();
      _sdpanel = straw.id().getPanel();
      _sdlayer = straw.id().getLayer();
      _sdstraw = straw.id().getStraw();

      for(size_t iend=0;iend<2;++iend){
        _xtime[iend] = xpair[iend]._time;
        _tctime[iend] = xpair[iend]._iclust->time();
        _charge[iend] = 0;
        _acharge[iend] = 0;
        _ddist[iend] = xpair[iend]._iclust->driftDistance();
        _dtime[iend] = xpair[iend]._iclust->driftTime();
        _ptime[iend] = xpair[iend]._iclust->propTime();
        _phi[iend] = xpair[iend]._iclust->driftPhi();
        _wdist[iend] = xpair[iend]._iclust->wireDistance();
        _vstart[iend] = xpair[iend]._vstart;
        _vcross[iend] = xpair[iend]._vcross;
        _tdc[iend] = digi.TDC(xpair[iend]._iclust->strawEnd());
        _tot[iend] = digi.TOT(xpair[iend]._iclust->strawEnd());
        StrawClusterList const& clist = wfs[iend].clusts().clustList();
        auto ctrig = xpair[iend]._iclust;
        _ncludd[iend] = clist.size();
        // find the earliest cluster from the same particle that triggered the crossing
        auto iclu = clist.begin();
        while( iclu != clist.end() && ctrig->strawGasStep()->simParticle() != iclu->strawGasStep()->simParticle() ){
          ++iclu;
        }
        if(iclu != clist.end() ){
          _ectime[iend] = iclu->time();
          _ecddist[iend] = iclu->driftDistance();
          _ecdtime[iend] = iclu->driftTime();
          _ecptime[iend] = iclu->propTime();
          _ecphi[iend] = iclu->driftPhi();
          // count how many clusters till we get to the trigger cluster
          size_t iclust(0);
          while( iclu != clist.end() && iclu != ctrig){
            _charge[iend] += iclu->charge();
            _acharge[iend] += iclu->charge();
            ++iclu;
            ++iclust;
          }
          _iclust[iend] = iclust;
          while (iclu != clist.end() && iclu->time() < _xtime[iend] + _mbbuffer){
            _acharge[iend] += iclu->charge();
            ++iclu;
          }
        }
      }
      if(xpair[0]._iclust->strawGasStep() == xpair[1]._iclust->strawGasStep())
        _nstep = 1;
      else
        _nstep = 2;
      _adc.clear();
      //for(auto iadc=digi.adcWaveform().begin();iadc!=digi.adcWaveform().end();++iadc){
      for(auto iadc=digiadc.samples().begin();iadc!=digiadc.samples().end();++iadc){
        _adc.push_back(*iadc);
      }
      _pmp = digi.PMP();
      _provenance = mcdigi.provenance();
      // mc truth information
      _dmcpdg = _dmcproc = _dmcgen = 0;
      _dmcmom = -1.0;
      _mctime = _mcenergy = _mctrigenergy = _mcthreshenergy = _mcdca = -1000.0;
      _mcthreshpdg = _mcthreshproc = _mcnstep = 0;
      auto const& sgsptr = mcdigi.earlyStrawGasStep();
      auto const& sgs = *sgsptr;
      _mctime = sgs.time() + _pbtimemc;
      // compute the doca for this step
      TwoLinePCA wirepca( straw.wirePosition(), straw.wireDirection(),
          GenVector::Hep3Vec(sgs.startPosition()), GenVector::Hep3Vec((sgs.endPosition()-sgs.startPosition()).Unit()) );
      _mcdca = wirepca.dca();
      auto spos = strawCoordinates(XYZVectorF(wirepca.point2()),straw);
      _mcdcaphi = spos._wirePosition.Phi();
      _mcdcadtime = strawphys.driftDistanceToTime(_mcdca,_mcdcaphi);

      TwoLinePCA strawpca( straw.strawPosition(), straw.strawDirection(),
           GenVector::Hep3Vec(sgs.startPosition()), GenVector::Hep3Vec(sgs.endPosition()-sgs.startPosition()) );
      _mcstrawdca = strawpca.dca();
      auto strawspos = strawCoordinates(XYZVectorF(strawpca.point2()),straw);
      _mcstrawdcaphi = strawspos._strawPosition.Phi();
      _mcstrawdcadtime = strawphys.driftDistanceToTime(_mcstrawdca,_mcstrawdcaphi);
      _dmcmom = sqrt(sgs.momentum().mag2());
      auto const& sp = *sgs.simParticle();
      _dmcpdg = sp.pdgId();
      _dmcproc = sp.creationCode();
      if(sp.genParticle().isNonnull())
        _dmcgen = sp.genParticle()->generatorId().id();
      _mcthreshpdg = sp.pdgId();
      _mcthreshproc = sp.creationCode();
      _sdtype= sgs.stepType()._stype;
      _sdwidth= sgs.width();
      _sdlen= sgs.stepLength();
      _mcenergy = mcdigi.energySum();
      _mctrigenergy = mcdigi.triggerEnergySum(StrawEnd::cal);
      // sum the energy from the explicit trigger particle, and find it's releationship
      _mcthreshenergy = 0.0; // FIXME!
      _mcnstep = 2;// FIXME!
      _xtalk = digi.strawId() != mcdigi.strawId();
      // fill the tree entry
      _sddiag->Fill();

    }//End of digiDiag

    void StrawDigisFromStrawGasSteps::stepDiag( StrawPhysics const& strawphys, StrawElectronics const& strawele,
        StrawGasStep const& sgs) {
      _steplen = sgs.stepLength();
      _stepE = sgs.ionizingEdep();
      _steptime = microbunchTime(strawele,sgs.time());
      _stype = sgs.stepType()._stype;
      _partP = sqrt(sgs.momentum().mag2());
      _partPDG = sgs.simParticle()->pdgId();
      _nclust = (int)_clusters.size();
      _netot = 0;
      _qsum = _esum = _eesum = 0.0;
      for(auto const& clust : _clusters) {
        _netot += clust._ne;
        _qsum += clust._charge;
        _esum += clust._eion;
        _eesum += strawphys.meanElectronEnergy()*clust._ne;
      }
      _qe = strawphys.ionizationEnergy(_qsum);
      _sdiag->Fill();
    }

    StrawCoordinates StrawDigisFromStrawGasSteps::strawCoordinates( XYZVectorF const& cpos, Straw const& straw) const {
      StrawCoordinates retval;
      static XYZVectorF zdir(0.0,0.0,1.0);
      XYZVectorF smid = XYZVectorF(straw.strawPosition());
      XYZVectorF delta = cpos - smid; // cluster position WRT nominal straw middle
      XYZVectorF sdir = XYZVectorF(straw.strawDirection());
      XYZVectorF pdir = sdir.Cross(zdir); // radial direction
      if(pdir.Dot(smid) < 0.0)pdir *= -1.0; // sign radially outwards
      float dw = delta.Dot(sdir);
      XYZVectorF cperp = delta - dw*sdir; // just perp part
      float phi = atan2(cperp.Dot(pdir),cperp.Dot(zdir)); // angle around wire WRT z axis in range -pi,pi
      float rho = sqrt(cperp.mag2());
      XYZVectorF truncpos = cpos; // if outside straw, move inside
      if (rho > (float)_rstraw){
        truncpos = cpos - cperp.Unit()*(rho-(float)_rstraw);
        rho = (float) _rstraw;
      }
//      float rho = min(sqrt(cperp.mag2()),(float)_rstraw); // truncate!

      retval._strawPosition = StrawPosition(rho,dw,phi);
      // now we get the wire position by calculating a position relative to the misaligned straw envelope
      //XYZVectorF simAlignedPos = strawPositionToXYZ(retval._strawPosition, simStraw);
      // calculate cylidrical coordinates relative to misaligned wire
      smid = XYZVectorF(straw.wirePosition());
      delta = truncpos - smid; // cluster position WRT wire middle
      sdir = XYZVectorF(straw.wireDirection());
      pdir = sdir.Cross(zdir); // radial direction
      if(pdir.Dot(smid) < 0.0)pdir *= -1.0; // sign radially outwards
      dw = delta.Dot(sdir);
      cperp = delta - dw*sdir; // just perp part
      phi = atan2(cperp.Dot(pdir),cperp.Dot(zdir));// angle around wire WRT Z axis in range -pi,pi
      rho = sqrt(cperp.mag2());
      retval._wirePosition = StrawPosition(rho,dw,phi);
      return retval;
    }

    XYZVectorF StrawDigisFromStrawGasSteps::strawCoordinatesToXYZ( StrawCoordinates const& cpos, Straw const& straw) const {
      static XYZVectorF zdir(0.0,0.0,1.0);
      XYZVectorF smid = XYZVectorF(straw.wirePosition());
      XYZVectorF sdir = XYZVectorF(straw.wireDirection());
      XYZVectorF pdir = sdir.Cross(zdir);
      if(pdir.Dot(smid) < 0.0)pdir *= -1.0; // sign radially outwards
      XYZVectorF cdir = cos(cpos._wirePosition.Phi())*zdir + sin(cpos._wirePosition.Phi())*pdir; // cluster direction perp to wire
      XYZVectorF retval = smid + cpos._wirePosition.Z()*sdir + cpos._wirePosition.Rho()*cdir;
      return retval;
    }

    } // end namespace trackermc
  } // end namespace mu2e

  using mu2e::TrackerMC::StrawDigisFromStrawGasSteps;
  DEFINE_ART_MODULE(StrawDigisFromStrawGasSteps)
