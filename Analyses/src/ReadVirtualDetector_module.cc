//
// Plugin to read virtual detectors data and create ntuples
//
//
// Original author Ivan Logashenko
//

#include "CLHEP/Units/SystemOfUnits.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/Mu2eUtilities/inc/fromStrings.hh"
#include "Offline/MCDataProducts/inc/G4BeamlineInfo.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "TH1F.h"
#include "TNtuple.h"
#include "TTree.h"
#include "Offline/GeometryService/inc/VirtualDetector.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <cmath>
#include <iostream>
#include <string>

using namespace std;

using CLHEP::Hep3Vector;
using CLHEP::keV;

namespace mu2e {

  const unsigned int nvdet = VirtualDetectorId::lastEnum;
  const unsigned int ntvdet = 20; // maximum number of time VD

  typedef struct {

    Int_t run;
    Int_t subrun;
    Int_t evt;
    Int_t trk;

    Int_t pdg;
    Float_t time;
    Float_t gtime;
    Float_t x;
    Float_t y;
    Float_t z;
    Float_t px;
    Float_t py;
    Float_t pz;
    Float_t p;
    Int_t code;
    Int_t creation_code;

    Bool_t isstop;
    Float_t tstop;
    Float_t gtstop;
    Float_t xstop;
    Float_t ystop;
    Float_t zstop;
    Int_t codestop;
    Float_t pxstop;
    Float_t pystop;
    Float_t pzstop;
    Float_t pstop;

    Int_t g4bl_evt;
    Int_t g4bl_trk;
    Float_t g4bl_weight;
    Float_t g4bl_time;

    Int_t parent_id;
    Int_t parent_pdg;
    Float_t parent_x;
    Float_t parent_y;
    Float_t parent_z;
    Float_t parent_px;
    Float_t parent_py;
    Float_t parent_pz;
    Float_t parent_p;
    Float_t parent_pxstop;
    Float_t parent_pystop;
    Float_t parent_pzstop;
    Float_t parent_pstop;
    Int_t parent_code;
    Float_t parent_lastke;

    Int_t nvd;
    Bool_t isvd[nvdet];
    Float_t tvd[nvdet];
    Float_t gtvd[nvdet];
    Float_t xvd[nvdet];
    Float_t yvd[nvdet];
    Float_t zvd[nvdet];
    Float_t pxvd[nvdet];
    Float_t pyvd[nvdet];
    Float_t pzvd[nvdet];
    Float_t pvd[nvdet];
    Float_t xlvd[nvdet];
    Float_t ylvd[nvdet];
    Float_t zlvd[nvdet];

    Int_t ntvd;
    Bool_t istvd[ntvdet];
    Float_t ttvd[ntvdet];
    Float_t gttvd[ntvdet];
    Float_t xtvd[ntvdet];
    Float_t ytvd[ntvdet];
    Float_t ztvd[ntvdet];
    Float_t pxtvd[ntvdet];
    Float_t pytvd[ntvdet];
    Float_t pztvd[ntvdet];
    Float_t ptvd[ntvdet];
    Int_t codetvd[ntvdet];

  } NtPartData;

  class ReadVirtualDetector : public art::EDAnalyzer {

    typedef vector<int> Vint;
    typedef vector<string> Vstr;
    typedef SimParticleCollection::key_type key_type;

    // Name of the VD and TVD StepPoint collections

    art::InputTag _vdInputTag;
    art::InputTag _tvdInputTag;
    art::InputTag _simpInputTag;
    art::InputTag _generatorInputTag;
    art::InputTag _physInputTag;

    // Control printed output.
    int _nAnalyzed;
    int _maxPrint;
    int _debugout;

    TNtuple* _ntvd;
    TNtuple* _nttvd;
    TTree* _ntpart;
    TTree* _ntpart1;
    TNtuple* _ntvdext;

    float *nt; // Need this buffer to fill TTree ntvd
    float *ntext; // Need this buffer to fill TTree ntvdext
    NtPartData ntp; // Buffer to fill particles ntuple

    bool write_ntvd;
    bool write_nttvd;
    bool write_ntpart;
    bool write_ntpart1;
    bool write_ntvdext;

    // Pointers to the physical volumes we are interested in
    // -- stopping target
    map<int,int> vid_stop;

    // List of particles of interest for the particles ntuple
    set<int> pdg_save;

    // List of particles to drop from time VD
    set<int> tvd_drop_pdg;

    // List of virtual detectors to be saved
    set<int> vd_save;

    // Module label of the g4 module that made the hits.
    std::string _g4ModuleLabel;

    // Virtual detector, which has to be crossed by particle before
    // it is saved in particles ntuple
    int _vd_required;

    // Save in the particles ntuple only those particles, which die
    // after this time (in ns)
    double _timeCut;

    // Save in the particles ntuple only particle with momentum larger than this
    double _minMomentum;

    // Save only stopped particles in the particles ntuple
    bool _stopped_only;

    // Save all particles
    bool _save_all_pdg;

    // Should we add together proper time for the whole decay chain
    bool _add_proper_time;

    // If we are analyzing output of the staged simulation, look for
    // real parent, navigating through the staged SimParticles
    bool _navigate_to_parent;

  public:

    explicit ReadVirtualDetector(fhicl::ParameterSet const& pset);
    virtual ~ReadVirtualDetector() { }

    virtual void beginJob();
    virtual void beginRun(art::Run const&);

    void analyze(const art::Event& e);

  };

  ReadVirtualDetector::ReadVirtualDetector(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _nAnalyzed(0),
    _maxPrint(pset.get<int>("maxPrint",0)),
    _debugout(pset.get<int>("debugOutput",0)),
    _ntvd(0), _nttvd(0), _ntpart(0), _ntpart1(0), _ntvdext(0),
    _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")), // obsolete, left for backward compatibility
    _vd_required(pset.get<int>("requireVD",0)),
    _timeCut(pset.get<double>("timeCut",0.0)),
    _minMomentum(pset.get<double>("minMomentum",-1.0)),
    _stopped_only(pset.get<bool>("saveStopped",false)),
    _save_all_pdg(pset.get<bool>("saveAllPDG",false)),
    _add_proper_time(pset.get<bool>("addProperTime",false)),
    _navigate_to_parent(pset.get<bool>("navigateToParent",true))
  {
    _vdInputTag = pset.get<std::string>("vdStepPoints","g4run:virtualdetector");
    _tvdInputTag = pset.get<string>("tvdStepPoints","g4run:timeVD");
    _simpInputTag = pset.get<string>("simParticleColl","g4run");
    _generatorInputTag = pset.get<std::string>("generatorModuleLabel", "generate");
    _physInputTag = pset.get<std::string>("physicsVolumeColl", "g4run");

    write_ntvd    = pset.get<bool>("writeNTVD",true);
    write_nttvd   = pset.get<bool>("writeNTTVD",true);
    write_ntpart  = pset.get<bool>("writeNTPART",true);
    write_ntpart1 = pset.get<bool>("writeNTPART1",true);
    write_ntvdext = pset.get<bool>("writeNTVDEXT",false);

    if( _debugout>0 ) cout << "ReadVirtualDetector: fill ntuples "
                           << " NTVD=" << write_ntvd
                           << " NTTVD=" << write_nttvd
                           << " NTPART=" << write_ntpart
                           << " NTPART1=" << write_ntpart1
                           << " NTVDEXT=" << write_ntvdext
                           << endl;
    if (_debugout > 1){
      std::cout << "_vd_required = " << _vd_required << std::endl;
    }
    Vstr const & pdg_names = pset.get<Vstr>("savePDG", Vstr());
    if( pdg_names.size()>0 ) {
      cout << "ReadVirtualDetector: save following particle types in the ntuple: ";
      for( size_t i=0; i<pdg_names.size(); ++i ) {
        PDGCode::enum_type id = PDGCode(pdg_names[i]);
        pdg_save.insert(int(id));
        cout << pdg_names[i] << " ("<<id<<"), ";
      }
      cout << endl;
    }

    Vstr const & tvd_drop_names = pset.get<Vstr>("tvdDropPDG", Vstr());
    if( tvd_drop_names.size()>0 ) {
      cout << "ReadVirtualDetector: drop following particle types from time VD ntuple: ";
      for( size_t i=0; i<tvd_drop_names.size(); ++i ) {
        PDGCode::enum_type id = PDGCode(tvd_drop_names[i]);
        tvd_drop_pdg.insert(int(id));
        cout << tvd_drop_names[i] << "("<<id<<"), ";
      }
      cout << endl;
    }

    Vint const & vd_ids = pset.get<Vint>("saveVD", Vint());
    if( vd_ids.size()>0 ) {
      cout << "ReadVirtualDetector: save data from the following virtual detectors: ";
      for( size_t i=0; i<vd_ids.size(); ++i ) {
        vd_save.insert(vd_ids[i]);
        cout << vd_ids[i] << ", ";
      }
      cout << endl;
    }

    nt    = new float[1000];
    ntext = new float[1000];

  }

  void ReadVirtualDetector::beginJob(){

    vid_stop.clear();

    // Get access to the TFile service.

    art::ServiceHandle<art::TFileService> tfs;

    if (write_ntvd){
      _ntvd = tfs->make<TNtuple>( "ntvd", "Virtual Detectors ntuple",
                                  "evt:trk:sid:pdg:time:x:y:z:px:py:pz:"
                                  "xl:yl:zl:pxl:pyl:pzl:gtime:"
                                  "g4bl_weight:g4bl_time:run:ke:subrun:code");
    }

    if (write_ntvdext){
      _ntvdext = tfs->make<TNtuple>( "ntvdext", "Virtual Detectors ntuple (extended)",
                                     "run:subrun:evt:trk:vdid:pdg:time:gtime:ke:"
                                     "x:y:z:px:py:pz:"
                                     "xl:yl:zl:pxl:pyl:pzl:"
                                     "code:creation_code:parent_pdg:"
                                     "originx:originy:originz");
    }

    if (write_nttvd){
      _nttvd = tfs->make<TNtuple>( "nttvd", "Time Virtual Detectors ntuple",
                                   "evt:trk:sid:pdg:time:x:y:z:px:py:pz:"
                                   "gtime:code:g4bl_weight:g4bl_time:run:ke:subrun");
    }

    // Have to use TTree here, because one cannot use more than 100 variables in TNtuple

    if (write_ntpart){
      _ntpart = tfs->make<TTree>("ntpart", "Particles ntuple");
      /*
        _ntpart->Branch("all",nt,
        "evt:trk:pdg:"
        "time:gtime:x:y:z:px:py:pz:"
        "isstop:tstop:gtstop:xstop:ystop:zstop:"
        "g4bl_evt:g4bl_trk:g4bl_weight:g4bl_time:"
        "parent_id:parent_pdg:"
        "parent_x:parent_y:parent_z:"
        "parent_px:parent_py:parent_pz:"
        "nvd:isvd[20]:"
        "tvd[20]:gtvd[20]:xvd[20]:yvd[20]:zvd[20]:"
        "pxvd[20]:pyvd[20]:pzvd[20]:"
        "xlvd[20]:ylvd[20]:zlvd[20]"
        );
      */

      _ntpart->Branch("run",        &ntp.run,        "run/I");
      _ntpart->Branch("subrun",     &ntp.subrun,     "subrun/I");
      _ntpart->Branch("evt",        &ntp.evt,        "evt/I");
      _ntpart->Branch("trk",        &ntp.trk,        "trk/I");
      _ntpart->Branch("pdg",        &ntp.pdg,        "pdg/I");
      _ntpart->Branch("time",       &ntp.time,       "time/F");
      _ntpart->Branch("gtime",      &ntp.gtime,      "gtime/F");
      _ntpart->Branch("x",          &ntp.x,          "x/F");
      _ntpart->Branch("y",          &ntp.y,          "y/F");
      _ntpart->Branch("z",          &ntp.z,          "z/F");
      _ntpart->Branch("px",         &ntp.px,         "px/F");
      _ntpart->Branch("py",         &ntp.py,         "py/F");
      _ntpart->Branch("pz",         &ntp.pz,         "pz/F");
      _ntpart->Branch("p",          &ntp.p,          "p/F");
      _ntpart->Branch("code",       &ntp.code,       "code/I");
      _ntpart->Branch("isstop",     &ntp.isstop,     "isstop/O");
      _ntpart->Branch("tstop",      &ntp.tstop,      "tstop/F");
      _ntpart->Branch("gtstop",     &ntp.gtstop,     "gtstop/F");
      _ntpart->Branch("xstop",      &ntp.xstop,      "xstop/F");
      _ntpart->Branch("ystop",      &ntp.ystop,      "ystop/F");
      _ntpart->Branch("zstop",      &ntp.zstop,      "zstop/F");
      _ntpart->Branch("codestop",   &ntp.codestop,   "codestop/I");
      _ntpart->Branch("pxstop",     &ntp.pxstop,     "pxstop/F");
      _ntpart->Branch("pystop",     &ntp.pystop,     "pystop/F");
      _ntpart->Branch("pzstop",     &ntp.pzstop,     "pzstop/F");
      _ntpart->Branch("pstop",      &ntp.pstop,      "pstop/F");
      _ntpart->Branch("g4bl_evt",   &ntp.g4bl_evt,   "g4bl_evt/I");
      _ntpart->Branch("g4bl_trk",   &ntp.g4bl_trk,   "g4bl_trk/I");
      _ntpart->Branch("g4bl_weight",&ntp.g4bl_weight,"g4bl_weight/F");
      _ntpart->Branch("g4bl_time",  &ntp.g4bl_time,  "g4bl_time/F");
      _ntpart->Branch("parent_id",  &ntp.parent_id,  "parent_id/I");
      _ntpart->Branch("parent_pdg", &ntp.parent_pdg, "parent_pdg/I");
      _ntpart->Branch("parent_x",   &ntp.parent_x,   "parent_x/F");
      _ntpart->Branch("parent_y",   &ntp.parent_y,   "parent_y/F");
      _ntpart->Branch("parent_z",   &ntp.parent_z,   "parent_z/F");
      _ntpart->Branch("parent_px",  &ntp.parent_px,  "parent_px/F");
      _ntpart->Branch("parent_py",  &ntp.parent_py,  "parent_py/F");
      _ntpart->Branch("parent_pz",  &ntp.parent_pz,  "parent_pz/F");
      _ntpart->Branch("parent_p",   &ntp.parent_p,   "parent_p/F");
      _ntpart->Branch("parent_pxstop",&ntp.parent_pxstop,"parent_pxstop/F");
      _ntpart->Branch("parent_pystop",&ntp.parent_pystop,"parent_pystop/F");
      _ntpart->Branch("parent_pzstop",&ntp.parent_pzstop,"parent_pzstop/F");
      _ntpart->Branch("parent_pstop", &ntp.parent_pstop, "parent_pstop/F");
      _ntpart->Branch("parent_code",  &ntp.parent_code,  "parent_code/I");
      _ntpart->Branch("parent_lastke",&ntp.parent_lastke,"parent_lastke/F");

      _ntpart->Branch("nvd",        &ntp.nvd,        "nvd/I");
      _ntpart->Branch("isvd",        ntp.isvd,       "isvd[nvd]/O");
      _ntpart->Branch("tvd",         ntp.tvd,        "tvd[nvd]/F");
      _ntpart->Branch("gtvd",        ntp.gtvd,       "gtvd[nvd]/F");
      _ntpart->Branch("xvd",         ntp.xvd,        "xvd[nvd]/F");
      _ntpart->Branch("yvd",         ntp.yvd,        "yvd[nvd]/F");
      _ntpart->Branch("zvd",         ntp.zvd,        "zvd[nvd]/F");
      _ntpart->Branch("pxvd",        ntp.pxvd,       "pxvd[nvd]/F");
      _ntpart->Branch("pyvd",        ntp.pyvd,       "pyvd[nvd]/F");
      _ntpart->Branch("pzvd",        ntp.pzvd,       "pzvd[nvd]/F");
      _ntpart->Branch("pvd",         ntp.pvd,        "pvd[nvd]/F");
      _ntpart->Branch("xlvd",        ntp.xlvd,       "xlvd[nvd]/F");
      _ntpart->Branch("ylvd",        ntp.ylvd,       "ylvd[nvd]/F");
      _ntpart->Branch("zlvd",        ntp.zlvd,       "zlvd[nvd]/F");

      _ntpart->Branch("ntvd",       &ntp.ntvd,       "ntvd/I");
      _ntpart->Branch("istvd",       ntp.istvd,      "istvd[ntvd]/O");
      _ntpart->Branch("ttvd",        ntp.ttvd,       "ttvd[ntvd]/F");
      _ntpart->Branch("gttvd",       ntp.gttvd,      "gttvd[ntvd]/F");
      _ntpart->Branch("xtvd",        ntp.xtvd,       "xtvd[ntvd]/F");
      _ntpart->Branch("ytvd",        ntp.ytvd,       "ytvd[ntvd]/F");
      _ntpart->Branch("ztvd",        ntp.ztvd,       "ztvd[ntvd]/F");
      _ntpart->Branch("pxtvd",       ntp.pxtvd,      "pxtvd[ntvd]/F");
      _ntpart->Branch("pytvd",       ntp.pytvd,      "pytvd[ntvd]/F");
      _ntpart->Branch("pztvd",       ntp.pztvd,      "pztvd[ntvd]/F");
      _ntpart->Branch("ptvd",        ntp.ptvd,       "ptvd[ntvd]/F");
      _ntpart->Branch("codetvd",     ntp.codetvd,    "codetvd[ntvd]/I");
    }

  }

  void ReadVirtualDetector::beginRun(art::Run const& run){
  }

  void ReadVirtualDetector::analyze(const art::Event& event) {

    ++_nAnalyzed;

    // Access virtual detectors geometry information
    // If not virtual detectors are defined, skip the rest
    if (_debugout > 1){
      std::cout << "nAnalyzed = " << _nAnalyzed << std::endl;
    }
    GeomHandle<VirtualDetector> vdg;
    if( vdg->nDet()<=0 ) return;
    GlobalConstantsHandle<ParticleDataList> pdt;

    // Ask the event to give us a "handle" to the requested hits.
    art::Handle<StepPointMCCollection> hits;
    event.getByLabel(_vdInputTag,hits);
    if( _debugout>0 ) {
      cout << "ReadVirtualDetector: hits collection is " << hits.isValid();
      if( hits.isValid() ) cout << " size=" << hits->size();
      cout << endl;
    }

    art::Handle<StepPointMCCollection> thits;
    event.getByLabel(_tvdInputTag,thits);
    if( _debugout>0 ) {
      cout << "ReadVirtualDetector: thits collection is " << thits.isValid();
      if( thits.isValid() ) cout << " size=" << thits->size();
      cout << endl;
    }

    art::Handle<SimParticleCollection> simParticles;
    event.getByLabel(_simpInputTag, simParticles);
    bool haveSimPart = simParticles.isValid();
    if ( haveSimPart ) haveSimPart = !(simParticles->empty());
    if( _debugout>0 ) {
      cout << "ReadVirtualDetector: simpart collection is " << simParticles.isValid();
      if( simParticles.isValid() ) cout << " size=" << simParticles->size();
      cout << endl;
    }

    art::Handle<G4BeamlineInfoCollection> g4beamlineData;
    event.getByLabel(_generatorInputTag, g4beamlineData);
    bool haveG4BL = g4beamlineData.isValid();
    if ( haveG4BL ) haveG4BL = (g4beamlineData->size()==1);
    if( _debugout>0 ) {
      cout << "ReadVirtualDetector: beamline collection is " << g4beamlineData.isValid();
      if( g4beamlineData.isValid() ) cout << " size=" << g4beamlineData->size();
      cout << endl;
    }

    // Fill virtual detectors ntuple

    if( haveG4BL ) {
      G4BeamlineInfo const& extra = g4beamlineData->at(0);
      nt[18] = extra.weight();
      nt[19] = extra.time();
    } else {
      nt[18] = 0;
      nt[19] = 0;
    }

    // Loop over all VD hits.
    if( hits.isValid() ) for ( size_t i=0; i<hits->size(); ++i ){

        if (_debugout > 1){
          std:: cout << "     hit number " << i << std::endl;
        }
        // Alias, used for readability.
        const StepPointMC& hit = (*hits)[i];

        // Get the hit information.

        int id = hit.volumeId();
        if (_debugout > 1){
          std::cout << "id = " << id << std::endl;
          std::cout << "vd_save.size() = " << vd_save.size() << "  vd_save.find(id) = "
                    << *(vd_save.find(id)) << "  vd_save.end() = " << *(vd_save.end()) << std::endl;
        }
        // If virtual detector id is not in the list - skip it
        if( vd_save.size()>0 && vd_save.find(id) == vd_save.end() ) continue;
        if (_debugout > 1){
          std::cout << "past save id = " << id << std::endl;
        }
        const CLHEP::Hep3Vector& pos = hit.position();
        const CLHEP::Hep3Vector& mom = hit.momentum();

        CLHEP::Hep3Vector lpos = (pos-vdg->getGlobal(id));
        CLHEP::Hep3Vector lmom = mom;
        if( vdg->getRotation(id)!=0 ) {
          lpos *= *(vdg->getRotation(id));
          lmom *= *(vdg->getRotation(id));
        }

        // Get track info
        key_type trackId = hit.trackId();
        int pdgId = 0;
        double mass(0.0);
        if ( haveSimPart ){
          if( !simParticles->has(trackId) ) {
            pdgId = 0;
          } else {
            SimParticle const& sim = simParticles->at(trackId);
            pdgId = sim.pdgId();
            // If virtual detector id is not in the list - skip it
            if (_debugout > 1){
          std::cout << "pdgId = " << pdgId << std::endl;
          std::cout << "pdg_save.size() = " << pdg_save.size() << "  pdg_save.find(id) = "
                    << *(pdg_save.find(id)) << "  pdg_save.end() = " << *(pdg_save.end()) << std::endl;
                  }
            if( _save_all_pdg || pdg_save.size() == 0 || ( pdg_save.size()>0 && pdg_save.find(pdgId) != pdg_save.end()) )
              {
                mass = pdt->particle(pdgId).mass();

                // Fill the ntuple.
                nt[0]  = event.id().event();
                nt[1]  = trackId.asInt();
                nt[2]  = hit.volumeId();
                nt[3]  = pdgId;
                nt[4]  = hit.time();
                nt[5]  = pos.x();
                nt[6]  = pos.y();
                nt[7]  = pos.z();
                nt[8]  = mom.x();
                nt[9]  = mom.y();
                nt[10] = mom.z();
                nt[11] = lpos.x();
                nt[12] = lpos.y();
                nt[13] = lpos.z();
                nt[14] = lmom.x();
                nt[15] = lmom.y();
                nt[16] = lmom.z();
                nt[17] = hit.properTime();
                nt[20] = event.id().run();
                // compute kinetic energy: this is what Geant cuts on
                nt[21] = sqrt(mom.mag2()+mass*mass)-mass;
                nt[22] = event.id().subRun();
                nt[23] =  sim.creationCode();

                if (write_ntvd){
                  if (_debugout > 1){
                    std::cout << "filling ntuple with pdg = " << pdgId << " and volume Id = " << hit.volumeId() << std::endl;
                  }
                  _ntvd->Fill(nt);
                }
              }
            // //print provenance
            // cout<<"Event:"<<event.id().event();
            // //loop back through all parents until we find one that wasn't a mu2ePrimary
            // SimParticle const* sim_parent = &sim;
            // cout<<"\tPDG="<<sim_parent->pdgId()<<" Creation code:"<<sim_parent->creationCode()<<" Origin:"<<sim_parent->startPosition()<<endl;
            // while( sim_parent && sim_parent->hasParent() ) { //&& sim_parent->creationCode()==56
            //   sim_parent = simParticles->getOrNull(sim_parent->parentId());
            //   if( sim_parent && sim_parent->endDefined() ) {
            //       //when creation code is 56, the parent is actually the same particle.
            //       //so it's real parent would be the grandparent in this hierarchy
            //       SimParticle const* sim_grandparent = simParticles->getOrNull(sim_parent->parentId());
            //       cout<<"       "<<"PDG="<<sim_parent->pdgId();
            //       if (sim_grandparent) cout<<" ParentPDG:"<< sim_grandparent->pdgId();
            //       cout<<" Creation code:"<<sim_parent->creationCode()<<" Origin:"<<sim_parent->startPosition()<<endl;
            //   } else {
            //     cout<<"       "<<"Has parent, but not defined."<<endl;
            //   }
            // }//end while


            //get the additional information for the extended VD ntuple
            int parent_pdg = 0;
            CLHEP::Hep3Vector origin(0.0,0.0,0.0);
            int this_creation_code = sim.creationCode();
            if (sim.isPrimary()){ //creation code 56=mu2ePrimary, means it was the first particle of this stage
              //loop back through all parents until we find one that wasn't a mu2ePrimary
              SimParticle const* sim_parent = &sim;
              while( sim_parent && sim_parent->hasParent() ) {
                sim_parent = simParticles->getOrNull(sim_parent->parentId());
                if( sim_parent && sim_parent->endDefined() ) {
                  if (!sim_parent->isPrimary()){
                    //when creation code is 56, the parent is actually the same particle.
                    //so it's real parent would be the grandparent in this hierarchy
                    SimParticle const* sim_grandparent = simParticles->getOrNull(sim_parent->parentId());
                    if (sim_grandparent && sim_grandparent->endDefined()) parent_pdg = sim_grandparent->pdgId();
                    //else throw cet::exception("ReadVirtualDetector")<< " (Grand)parent not defined. \n";

                    origin = sim_parent->startPosition();
                    this_creation_code = sim_parent->creationCode();
                    break;//we're just looking for the first one that's not a primary
                  }
                }
              }//end while
            } else {
              SimParticle const* sim_parent = simParticles->getOrNull(sim.parentId());
              if( sim_parent && sim_parent->endDefined() ){
                parent_pdg = sim_parent->pdgId();
              }
              origin = sim.startPosition();
              this_creation_code = sim.creationCode();
            }
            // Fill the extended ntuple.
            ntext[0]  = event.id().run();
            ntext[1]  = event.id().subRun();
            ntext[2]  = event.id().event();
            ntext[3]  = trackId.asInt();
            ntext[4]  = hit.volumeId();
            ntext[5]  = pdgId;
            ntext[6]  = hit.time();
            ntext[7]  = hit.properTime();
            ntext[8]  = sqrt(mom.mag2()+mass*mass)-mass; // compute kinetic energy: this is what Geant cuts on
            ntext[9]  = pos.x();
            ntext[10] = pos.y();
            ntext[11] = pos.z();
            ntext[12] = mom.x();
            ntext[13] = mom.y();
            ntext[14] = mom.z();
            ntext[15] = lpos.x();
            ntext[16] = lpos.y();
            ntext[17] = lpos.z();
            ntext[18] = lmom.x();
            ntext[19] = lmom.y();
            ntext[20] = lmom.z();
            ntext[21] = sim.creationCode();//geant4 creation code (might be mu2ePrimary i.e. not physics)
            ntext[22] = this_creation_code;//This is supposed to be only physics creation codes, so trace back to parent when mu2ePrimary
            ntext[23] = parent_pdg;
            ntext[24] = origin.x();
            ntext[25] = origin.y();
            ntext[26] = origin.z();
            if (write_ntvdext){
              _ntvdext->Fill(ntext);
            }


          }
        }

        if ( _nAnalyzed < _maxPrint){
          cout << "VD hit: "
               << event.id().run()   << " | "
               << event.id().subRun()<< " | "
               << event.id().event() << " | "
               << hit.volumeId()     << " "
               << pdgId              << " | "
               << hit.time()         << " "
               << lpos               << " "
               << mom.mag()
               << endl;

        }

      } // end loop over hits.


    // Fill time virtual detectors ntuple

    if( haveG4BL ) {
      G4BeamlineInfo const& extra = g4beamlineData->at(0);
      nt[13] = extra.weight();
      nt[14] = extra.time();
    } else {
      nt[13] = 0;
      nt[14] = 0;
    }

    // Loop over all time VD hits.

    if( thits.isValid() ) for ( size_t i=0; i<thits->size(); ++i ){

        // Alias, used for readability.
        const StepPointMC& hit = (*thits)[i];

        // Get the hit information.

        int id = hit.volumeId();

        const CLHEP::Hep3Vector& pos = hit.position();
        const CLHEP::Hep3Vector& mom = hit.momentum();

        // Get track info
        key_type trackId = hit.trackId();
        int pdgId = 0;
        double mass(0.0);
        if ( haveSimPart ){
          if( !simParticles->has(trackId) ) {
            pdgId = 0;
          } else {
            SimParticle const& sim = simParticles->at(trackId);
            pdgId = sim.pdgId();
            mass = pdt->particle(pdgId).mass();
          }
        }

        if( tvd_drop_pdg.size()>0 && tvd_drop_pdg.find(pdgId)!=tvd_drop_pdg.end() ) continue;

        // Fill the ntuple.
        nt[0]  = event.id().event();
        nt[1]  = trackId.asInt();
        nt[2]  = id;
        nt[3]  = pdgId;
        nt[4]  = hit.time();
        nt[5]  = pos.x();
        nt[6]  = pos.y();
        nt[7]  = pos.z();
        nt[8]  = mom.x();
        nt[9]  = mom.y();
        nt[10] = mom.z();
        nt[11] = hit.properTime();
        nt[12] = hit.endProcessCode();
        nt[15] = event.id().run();
        nt[16] = sqrt(mom.mag2()+mass*mass)-mass;
        nt[17] = event.id().subRun();
        if (write_nttvd){
          _nttvd->Fill(nt);
        }

        if ( _nAnalyzed < _maxPrint){
          cout << "TVD hit: "
               << event.id().run()   << " | "
               << event.id().subRun()<< " | "
               << event.id().event() << " | "
               << hit.volumeId()     << " "
               << pdgId              << " | "
               << hit.time()         << " "
               << pos                << " "
               << mom.mag()
               << endl;

        }

      } // end loop over hits.

    // Fill tracks ntuple
    if( haveSimPart && (pdg_save.size()>0 || _save_all_pdg) ) {

      // Go through SimParticle container and analyze one particle at a time
      for ( SimParticleCollection::const_iterator isp=simParticles->begin();
            isp!=simParticles->end(); ++isp ){
        SimParticle const& sim = isp->second;

        // It particle PDG ID is not in the list - skip it
        if( !_save_all_pdg && pdg_save.find(sim.pdgId()) == pdg_save.end() ) continue;

        // Save SimParticle header info
        ntp.run = event.id().run();      // run_id
        ntp.subrun = event.id().subRun();  // subrun_id
        ntp.evt = event.id().event();    // event_id
        ntp.trk = sim.id().asInt();      // track_id
        ntp.pdg = sim.pdgId();           // PDG id

        if (_debugout > 1){
          std::cout << "ntp.pdg = " << ntp.pdg <<std::endl;
        }
        // Calculate parent proper time
        double gtime_parent = 0.0;
        if( _add_proper_time ) {
          SimParticle const* sim_par = &sim;
          while( sim_par && sim_par->hasParent() ) {
            sim_par = simParticles->getOrNull(sim_par->parentId());
            if( sim_par && sim_par->pdgId()==ntp.pdg && sim.endDefined() ) {
              gtime_parent += sim_par->endProperTime();
            }
          }
        }

        // Parent info
        SimParticle const* sim_parent = 0; // True parent
        SimParticle const* sim_child = &sim; // First incrarnation of the current particle

        while( sim_child && sim_child->hasParent() ) {
          sim_parent = simParticles->getOrNull(sim_child->parentId());
          if(_navigate_to_parent && (sim_child->creationCode() == ProcessCode::mu2ePrimary)) {
            // Thats not really a parent, this is the same particle
            // Need to navigate one more step up the chain
            sim_child = sim_parent;
          } else {
            break;
          }
        }

        // Save SimParticle other info
        ntp.time = sim.startGlobalTime(); // start time
        ntp.gtime = gtime_parent+sim.startProperTime(); // start time
        if( sim_child ) {
          CLHEP::Hep3Vector const & pos_start = sim_child->startPosition();
          CLHEP::Hep3Vector const mom_start = sim_child->startMomentum();
          ntp.x = pos_start.x();
          ntp.y = pos_start.y();
          ntp.z = pos_start.z();
          ntp.px = mom_start.x();
          ntp.py = mom_start.y();
          ntp.pz = mom_start.z();
          ntp.p = mom_start.mag();
          ntp.code = sim_child->creationCode();
        } else {
          ntp.x = 0;
          ntp.y = 0;
          ntp.z = 0;
          ntp.px = 0;
          ntp.py = 0;
          ntp.pz = 0;
          ntp.p = 0;
          ntp.code = sim.creationCode();
        }

        // Apply momentum cut
        if( ntp.p>0 && ntp.p<_minMomentum ) continue;

        // Check id of the volume where particle dies
        if( sim.endDefined() ) {
          if( vid_stop.find(sim.endVolumeIndex()) != vid_stop.end() ) {
            ntp.isstop = true;
          } else {
            ntp.isstop = false;
          }
          ntp.tstop = sim.endGlobalTime();
          ntp.gtstop = gtime_parent+sim.endProperTime();
          CLHEP::Hep3Vector const & pos_end = sim.endPosition();
          ntp.xstop = pos_end.x();
          ntp.ystop = pos_end.y();
          ntp.zstop = pos_end.z();
          ntp.codestop = sim.stoppingCode();
          CLHEP::Hep3Vector const mom_end = sim.endMomentum();
          ntp.pxstop = mom_end.x();
          ntp.pystop = mom_end.y();
          ntp.pzstop = mom_end.z();
          ntp.pstop  = mom_end.mag();
        } else {
          ntp.isstop = false;
          ntp.tstop = 0;
          ntp.gtstop = 0;
          ntp.xstop = 0;
          ntp.ystop = 0;
          ntp.zstop = 0;
          ntp.codestop = 0;
          ntp.pxstop = 0;
          ntp.pystop = 0;
          ntp.pzstop = 0;
          ntp.pstop  = 0;
        }

        if( haveG4BL ) {
          G4BeamlineInfo const& extra = g4beamlineData->at(0);
          ntp.g4bl_evt    = extra.eventId();
          ntp.g4bl_trk    = extra.trackId();
          ntp.g4bl_weight = extra.weight();
          ntp.g4bl_time   = extra.time();
        } else {
          ntp.g4bl_evt    = 0;
          ntp.g4bl_trk    = 0;
          ntp.g4bl_weight = 0;
          ntp.g4bl_time   = 0;
        }

        if( sim_parent ) {
          ntp.parent_id = sim_parent->id().asInt();
          ntp.parent_pdg = sim_parent->pdgId();
          CLHEP::Hep3Vector const & pos_parent = sim_parent->startPosition();
          CLHEP::Hep3Vector const mom_parent = sim_parent->startMomentum();
          ntp.parent_x = pos_parent.x();
          ntp.parent_y = pos_parent.y();
          ntp.parent_z = pos_parent.z();
          ntp.parent_px = mom_parent.x();
          ntp.parent_py = mom_parent.y();
          ntp.parent_pz = mom_parent.z();
          ntp.parent_p = mom_parent.mag();
          CLHEP::Hep3Vector const endmom_parent = sim_parent->endMomentum();
          ntp.parent_pxstop = endmom_parent.x();
          ntp.parent_pystop = endmom_parent.y();
          ntp.parent_pzstop = endmom_parent.z();
          ntp.parent_pstop = endmom_parent.mag();
          ntp.parent_code = sim_parent->stoppingCode();
          ntp.parent_lastke = sim_parent->preLastStepKineticEnergy();
        } else {
          ntp.parent_id = -1;
          ntp.parent_pdg = 0;
          ntp.parent_x = 0;
          ntp.parent_y = 0;
          ntp.parent_z = 0;
          ntp.parent_px = 0;
          ntp.parent_py = 0;
          ntp.parent_pz = 0;
          ntp.parent_p = 0;
          ntp.parent_pxstop = 0;
          ntp.parent_pystop = 0;
          ntp.parent_pzstop = 0;
          ntp.parent_pstop = 0;
          ntp.parent_code = 0;
          ntp.parent_lastke = 0;
        }

        // Clear up VD data
        ntp.nvd = nvdet;
        for ( size_t i=0; i<nvdet; ++i ) {
          ntp.isvd[i]=false;
          ntp.tvd[i]=0;
          ntp.gtvd[i]=0;
          ntp.xvd[i]=0;
          ntp.yvd[i]=0;
          ntp.zvd[i]=0;
          ntp.pxvd[i]=0;
          ntp.pyvd[i]=0;
          ntp.pzvd[i]=0;
          ntp.pvd[i]=0;
          ntp.xlvd[i]=0;
          ntp.ylvd[i]=0;
          ntp.zlvd[i]=0;
        }

        // Loop over all virtual detectors and fill corresponding data
        for ( size_t i=0; i<hits->size(); ++i ){

          // Alias, used for readability.
          const StepPointMC& hit = (*hits)[i];

          // Only use hits associated with current particle
          key_type trackId = hit.trackId();
          //if( trackId != isp->first ) continue;
          bool sim_vd_found = (trackId == sim.id());
          if( (!sim_vd_found) && _navigate_to_parent ) {
            SimParticle const* sim_vd = &sim;
            while( (!sim_vd_found) && (sim_vd->creationCode() == ProcessCode::mu2ePrimary) && sim_vd->hasParent() ) {
              sim_vd = simParticles->getOrNull(sim_vd->parentId());
              if( !sim_vd ) break;
              sim_vd_found = (trackId == sim_vd->id());
            }
          }
          if( !sim_vd_found ) continue;

          // Get the hit information.

          unsigned int id = hit.volumeId();

          if( id<=0 || id>nvdet || ntp.isvd[id-1] ) continue;

          const CLHEP::Hep3Vector& pos = hit.position();
          const CLHEP::Hep3Vector& mom = hit.momentum();

          CLHEP::Hep3Vector lpos = (pos-vdg->getGlobal(id));
          if( vdg->getRotation(id)!=0 ) {
            lpos *= *(vdg->getRotation(id));
          }

          ntp.isvd[id-1] = true;
          ntp.tvd[id-1]  = hit.time();
          ntp.gtvd[id-1] = gtime_parent+hit.properTime();
          ntp.xvd[id-1]  = pos.x();
          ntp.yvd[id-1]  = pos.y();
          ntp.zvd[id-1]  = pos.z();
          ntp.pxvd[id-1] = mom.x();
          ntp.pyvd[id-1] = mom.y();
          ntp.pzvd[id-1] = mom.z();
          ntp.pvd[id-1]  = mom.mag();
          ntp.xlvd[id-1] = lpos.x();
          ntp.ylvd[id-1] = lpos.y();
          ntp.zlvd[id-1] = lpos.z();

        } // end loop over hits.

        // Clear up time VD data
        ntp.ntvd = ntvdet;
        for ( size_t i=0; i<ntvdet; ++i ) {
          ntp.istvd[i]=false;
          ntp.ttvd[i]=0;
          ntp.gttvd[i]=0;
          ntp.xtvd[i]=0;
          ntp.ytvd[i]=0;
          ntp.ztvd[i]=0;
          ntp.pxtvd[i]=0;
          ntp.pytvd[i]=0;
          ntp.pztvd[i]=0;
          ntp.ptvd[i]=0;
          ntp.codetvd[i]=0;
        }

        // Loop over all virtual detectors and fill corresponding data
        if( thits.isValid() ) for ( size_t i=0; i<thits->size(); ++i ){

            // Alias, used for readability.
            const StepPointMC& hit = (*thits)[i];

            // Only use hits associated with current particle
            key_type trackId = hit.trackId();
            //if( trackId != isp->first ) continue;
            bool sim_tvd_found = (trackId == sim.id());
            if( (!sim_tvd_found) && _navigate_to_parent ) {
              SimParticle const* sim_tvd = &sim;
              while( (!sim_tvd_found) && (sim_tvd->creationCode() == ProcessCode::mu2ePrimary) && sim_tvd->hasParent() ) {
                sim_tvd = simParticles->getOrNull(sim_tvd->parentId());
                if( !sim_tvd ) break;
                sim_tvd_found = (trackId == sim_tvd->id());
              }
            }
            if( !sim_tvd_found ) continue;

            // Get the hit information.

            unsigned int id = hit.volumeId();

            if( id<=0 || id>ntvdet || ntp.istvd[id-1] ) continue;

            const CLHEP::Hep3Vector& pos = hit.position();
            const CLHEP::Hep3Vector& mom = hit.momentum();

            ntp.istvd[id-1] = true;
            ntp.ttvd[id-1]  = hit.time();
            ntp.gttvd[id-1] = gtime_parent+hit.properTime();
            ntp.xtvd[id-1]  = pos.x();
            ntp.ytvd[id-1]  = pos.y();
            ntp.ztvd[id-1]  = pos.z();
            ntp.pxtvd[id-1] = mom.x();
            ntp.pytvd[id-1] = mom.y();
            ntp.pztvd[id-1] = mom.z();
            ntp.ptvd[id-1]  = mom.mag();
            ntp.codetvd[id-1] = hit.endProcessCode();

          } // end loop over hits.

        //--------------------

        // Keep only stopped particles
        if( _stopped_only && !ntp.isstop ) continue;

        // Keep only those particles which went through required VD

        if (_debugout > 1){
          std::cout << "vd_required = " << _vd_required << " " << ntp.isvd[_vd_required-1] << std::endl;
        }
        if( _vd_required>0 && !ntp.isvd[_vd_required-1] ) continue;

        // Keep only those particles, which die late enough
        if( _timeCut>0.1 && ntp.tstop<_timeCut ) continue;

        if (write_ntpart){
          _ntpart->Fill();
        }

      }
    }

  }

}  // end namespace mu2e

using mu2e::ReadVirtualDetector;
DEFINE_ART_MODULE(ReadVirtualDetector)
