////////////////////////////////////////////////////////////////////////
// Class:       KalSeedToTrkQualTree
// Plugin Type: analyzer (art v2_06_02)
// File:        KalSeedToTrkQualTree_module.cc
//
// Generated at Fri Apr 21 14:56:02 2017 by Andrew Edmonds using cetskelgen
// from cetlib version v2_02_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TTree.h"

#include "art_root_io/TFileService.h"

#include "BTrk/TrkBase/TrkHelixUtils.hh"

#include "DataProducts/inc/IndexMap.hh"

#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/EventWeight.hh"

#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/TrkQual.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/TrkStrawHitSeed.hh"

#include "TrkDiag/inc/helixpar.hh"
#include "TrkDiag/inc/TrkMCTools.hh"

#include "BTrk/BbrGeom/HepPoint.h"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "BFieldGeom/inc/BFieldManager.hh"

#include "CLHEP/Units/PhysicalConstants.h"

namespace mu2e {
  class KalSeedToTrkQualTree;
}


class mu2e::KalSeedToTrkQualTree : public art::EDAnalyzer {
public:
  explicit KalSeedToTrkQualTree(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  KalSeedToTrkQualTree(KalSeedToTrkQualTree const &) = delete;
  KalSeedToTrkQualTree(KalSeedToTrkQualTree &&) = delete;
  KalSeedToTrkQualTree & operator = (KalSeedToTrkQualTree const &) = delete;
  KalSeedToTrkQualTree & operator = (KalSeedToTrkQualTree &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  // art input tags for collections
  art::InputTag _kalFinalTag;
  art::InputTag _strawDigiMCTag;
  art::InputTag _stepPointMCTag;
  art::InputTag _trkQualTag;
  art::InputTag _beamWeightTag;
  art::InputTag _genWeightTag;

  // For the TTree
  TTree* _trkQualTree; // the tree we will fill

  // event weight info
  float _eventWeight; // total event weight (product of weights below)
  float _beamWeight; // proton beam intensity weight
  float _genWeight; // weight from the generator (e.g. DIO spectrum)

  // variables needed for cutting but not used in the TrkQualMVA itself
  TrkFitFlag _fitStatus; // the status of the fit
  float _fitMom; // the fit momentum at the front of the tracker
  float _mcMom; // true momentum at the front of the tracker
  float _t0; // t0 of the track
  float _tandip; // tandip of the track

  // the current TrkQual information
  std::vector<float> _trkQualVariableVals; // need a vector of floats to use TMVA::Reader so I will copy things in here
  float _trkQualValue;

  // variables that we will add will go in here
  //
  //
  //


  // For condensed collections
  bool _usingCondensedCollections; // true if we are using condensed collection and so should expect an IndexMap
  art::InputTag _indexMapTag;
  const mu2e::IndexMap* _indexMap;


  // Functions that we use to fill variables
  void fillEventInfo(const art::Event& event) {
    _genWeight = _beamWeight = _eventWeight = 1;

    art::Handle<EventWeight> beamWeightHandle;
    event.getByLabel(_beamWeightTag, beamWeightHandle);
    if(beamWeightHandle.isValid()){
      _beamWeight = beamWeightHandle->weight();
      _eventWeight *= _beamWeight;
    }

    art::Handle<EventWeight> genWeightHandle;
    event.getByLabel(_genWeightTag, genWeightHandle);
    if(genWeightHandle.isValid()){
      _genWeight = genWeightHandle->weight();
      _eventWeight *= _genWeight;
    }
  }

  void fillKalSeedInfo(const KalSeed& kal_seed) {
    _fitStatus = kal_seed.status();
    _t0 = kal_seed.t0().t0();
  }

  void fillKalSegmentInfo(const KalSegment& kal_segment) {
    _fitMom = kal_segment.mom();
  }

  void fillStepPointMCInfo(const StepPointMC& step_point) {
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    GlobalConstantsHandle<ParticleDataTable> pdt;

    _mcMom = step_point.momentum().mag();

    double charge = pdt->particle(step_point.simParticle()->pdgId()).ref().charge();
    double hflt(0.0);
    CLHEP::HepVector parvec(5,0);
    static CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(CLHEP::Hep3Vector(0.0,0.0,0.0));
    static double bz = bfmgr->getBField(vpoint_mu2e).z();
    CLHEP::Hep3Vector pos(step_point.position());
    HepPoint ppos(pos.x(),pos.y(),pos.z());
    TrkHelixUtils::helixFromMom( parvec, hflt,ppos, step_point.momentum(),charge,bz);
    helixpar hpar = helixpar(parvec);
    _tandip = hpar._td;
  }

  void fillTrkQualInfo(const TrkQual& trk_qual) {

    int n_trkqual_vars = TrkQual::n_vars;
    for (int i_trkqual_var = 0; i_trkqual_var < n_trkqual_vars; ++i_trkqual_var) {
      TrkQual::MVA_varindex i_index = TrkQual::MVA_varindex(i_trkqual_var);
      _trkQualVariableVals[i_trkqual_var] = (float) trk_qual[i_index];

    }
    _trkQualValue = trk_qual.MVAOutput();
  }
};


mu2e::KalSeedToTrkQualTree::KalSeedToTrkQualTree(fhicl::ParameterSet const & pset)
  :
  EDAnalyzer(pset),
  _kalFinalTag(pset.get<art::InputTag>("kalFinalTag")),
  _strawDigiMCTag(pset.get<art::InputTag>("strawDigiMCTag")),
  _stepPointMCTag(pset.get<art::InputTag>("stepPointMCTag")),
  _trkQualTag(pset.get<art::InputTag>("trkQualTag")),
  _beamWeightTag(pset.get<art::InputTag>("beamWeightTag", "")),
  _genWeightTag(pset.get<art::InputTag>("genWeightTag", "")),
  _usingCondensedCollections(pset.get<bool>("usingCondensedCollections", "false")),
  _indexMapTag(pset.get<art::InputTag>("indexMapTag", ""))
 // More initializers here.
{
  art::ServiceHandle<art::TFileService> tfs;
  _trkQualTree=tfs->make<TTree>("trkQualTree","Variables for Testing and Training TrkQual");

  _trkQualTree->Branch("EventWeight", &_eventWeight);
  _trkQualTree->Branch("BeamWeight", &_beamWeight);
  _trkQualTree->Branch("GenWeight", &_genWeight);
  _trkQualTree->Branch("FitStatus", &_fitStatus);
  _trkQualTree->Branch("FitMom", &_fitMom);
  _trkQualTree->Branch("MCMom", &_mcMom);
  _trkQualTree->Branch("T0", &_t0);
  _trkQualTree->Branch("TanDip", &_tandip);

  int n_trkqual_vars = TrkQual::n_vars;
  _trkQualVariableVals.reserve(n_trkqual_vars);
  for (int i_trkqual_var = 0; i_trkqual_var < n_trkqual_vars; ++i_trkqual_var) {
    TrkQual::MVA_varindex i_index =TrkQual::MVA_varindex(i_trkqual_var);
    std::string varname = TrkQual::varName(i_index);

    _trkQualTree->Branch(varname.c_str(), &_trkQualVariableVals[i_trkqual_var]);
  }
  //  _trkQualTree->Branch("trkQual", &_trkQual);
  _trkQualTree->Branch("OriginalTrkQualVal", &_trkQualValue);

  // new variables will go here
}

void mu2e::KalSeedToTrkQualTree::analyze(art::Event const & event)
{
  // The collections we'll need
  const auto& finalKals = event.getValidHandle<KalSeedCollection>(_kalFinalTag);
  const auto& strawDigiMCs = event.getValidHandle<StrawDigiMCCollection>(_strawDigiMCTag);
  const auto& vdStepPointMCs = event.getValidHandle<StepPointMCCollection>(_stepPointMCTag);

  if (_usingCondensedCollections) {
    _indexMap = event.getValidHandle<IndexMap>(_indexMapTag);
  }

  fillEventInfo(event);

  // Loop through the final KalSeed fits
  for (const auto& finalKal : *finalKals) {

    fillKalSeedInfo(finalKal);

    // Get the segment of the fit at the front of the tracker
    KalSegment const& frontSegment = finalKal.segments().front(); 
    fillKalSegmentInfo(frontSegment);

    // Find the SimParticle that this track corresponds to
    const auto& hits = finalKal.hits();
    art::Ptr<mu2e::SimParticle> primary; // will be filled

    std::vector<StrawHitIndex> hit_indices;
    for ( const auto& trk_straw_hit_seed : hits) {
      // might want the condensed hit indices
      StrawHitIndex index = trk_straw_hit_seed.index();
      if (_usingCondensedCollections) {
	index = _indexMap->getCondensedIndex(index);
      }

      hit_indices.push_back(index);
    }
    TrkMCTools::primaryParticle(primary, hit_indices, strawDigiMCs);

    // Find the StepPointMC at the front of the tracker
    const auto& productid = primary.id();
    const auto& trackid = primary->id();
    for ( const auto& i_step_point : *vdStepPointMCs) {

      if ( (i_step_point.volumeId() == VirtualDetectorId::TT_FrontHollow || i_step_point.volumeId() == VirtualDetectorId::TT_FrontPA)	   
	   && i_step_point.trackId() == trackid
	   && i_step_point.simParticle().id() == productid) {

	fillStepPointMCInfo(i_step_point);
	break;
      }
    }
  }

  // Get the TrkQual information that was used originally
  const auto& trkQuals =  event.getValidHandle<TrkQualCollection>(_trkQualTag);
  for (const auto& trkQual : *trkQuals) {
    fillTrkQualInfo(trkQual);
  }

  // Fill the tree
  if (finalKals->size()>0) {
    _trkQualTree->Fill();
  }
}

DEFINE_ART_MODULE(mu2e::KalSeedToTrkQualTree)
