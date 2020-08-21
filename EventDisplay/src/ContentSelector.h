//
// Class which manages the combo boxes and list box in the event display frame. It is able to returns the data objects associated with the selected box entries.
//
//
// Original author Ralf Ehrlich
//

#ifndef EventDisplay_src_ContentSelector_h
#define EventDisplay_src_ContentSelector_h

#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CrvRecoPulseCollection.hh"
#include "RecoDataProducts/inc/CrvDigiCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include <TGComboBox.h>
#include <TGListBox.h>
#include <iostream>
#include <vector>

using namespace CLHEP;
#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/TrkExtTrajCollection.hh"

namespace mu2e_eventdisplay
{

class ContentSelector
{
  ContentSelector();
  ContentSelector(const ContentSelector &);
  ContentSelector& operator=(const ContentSelector &);

  std::vector<art::Handle<mu2e::StepPointMCCollection> > _stepPointMCVector;
  std::vector<art::Handle<mu2e::StepPointMCCollection> > _caloStepPointMCVector;
  std::vector<art::Handle<mu2e::StrawHitCollection> > _strawHitVector;
  std::vector<art::Handle<mu2e::StrawHitFlagCollection> > _strawHitFlagVector;
  std::vector<art::Handle<mu2e::StrawHitPositionCollection> > _strawHitPositionVector;
  std::vector<art::Handle<mu2e::CaloCrystalHitCollection> > _caloCrystalHitVector;
  std::vector<art::Handle<mu2e::CaloHitCollection> > _caloHitVector;
  std::vector<art::Handle<mu2e::CrvRecoPulseCollection> > _crvRecoPulseVector;
  std::vector<art::Handle<mu2e::CrvDigiCollection> > _crvDigisVector;
  std::vector<art::Handle<mu2e::SimParticleCollection> > _simParticleVector;
  std::vector<art::Handle<mu2e::MCTrajectoryCollection> > _mcTrajectoryVector;
  std::vector<art::Handle<mu2e::KalRepCollection> > _trkRecoTrkVector;
  std::vector<art::Handle<mu2e::KalRepCollection> > _hitOnTrackVector; //Hits on Tracks are stored inside of KalRep
  std::vector<art::Handle<mu2e::KalSeedCollection> > _kalSeedTrkVector;
  std::vector<art::Handle<mu2e::KalSeedCollection> > _kalSeedHitVector;
  std::vector<art::Handle<mu2e::TrkExtTrajCollection> > _trkExtTrajVector;
  art::Handle<mu2e::PhysicalVolumeInfoMultiCollection> _physicalVolumesMulti;
  bool _hasPhysicalVolumesMulti;

  TGComboBox  *_hitBox;
  TGComboBox  *_caloHitBox;
  TGComboBox  *_crvHitBox;
  TGListBox   *_trackBox;
  std::string _g4ModuleLabel;
  std::string _physicalVolumesMultiLabel;

  public:
  struct trackInfoStruct
  {
    int classID, index;
    std::string entryText;
    std::string moduleLabel, productInstanceName;
    art::ProductID productId;
  };

  private:
  struct entryStruct
  {
    int         entryID, classID, vectorPos;
    std::string entryText;
    std::string className, moduleLabel, productInstanceName;
    bool operator==(const entryStruct& rhs) const
    {
      return ((this->entryID==rhs.entryID) && (this->entryText==rhs.entryText));
    }
  };
  std::vector<entryStruct> _hitEntries, _caloHitEntries, _crvHitEntries, _trackEntries;
  std::vector<entryStruct> _hitFlagEntries, _hitPositionEntries;

  std::string _selectedHitFlagEntry, _selectedHitPositionEntry;

  template<class CollectionType> void createNewEntries(std::vector<art::Handle<CollectionType> > &dataVector,
                                                       const art::Event &event, const std::string &className,
                                                       std::vector<entryStruct> &newEntries, int classID);

  public:
  ContentSelector(TGComboBox *hitBox, TGComboBox *caloHitBox, TGComboBox *crvHitBox, TGListBox *trackBox, 
                  std::string const &g4ModuleLabel, std::string const &physicalVolumesMultiLabel);
  void firstLoop();
  void setAvailableCollections(const art::Event& event);

  bool getSelectedHitsName(std::string &className,
                           std::string &moduleLabel,
                           std::string &productInstanceName) const;
  std::vector<trackInfoStruct> getSelectedTrackNames() const;

  template<typename CollectionType> const CollectionType* getSelectedHitCollection() const;
  template<typename CollectionType> const CollectionType* getSelectedCaloHitCollection() const;
  template<typename CollectionType> const CollectionType* getSelectedCrvHitCollection() const;
  const std::vector<art::Handle<mu2e::CrvDigiCollection> >& getSelectedCrvDigiCollection() const;
  template<typename CollectionType> std::vector<const CollectionType*> getSelectedTrackCollection(std::vector<trackInfoStruct> &v) const;
  const mu2e::PhysicalVolumeInfoMultiCollection *getPhysicalVolumeInfoMultiCollection() const;
  const mu2e::MCTrajectoryCollection *getMCTrajectoryCollection(const trackInfoStruct &t) const;

  //for filter and setup dialog
  const std::vector<entryStruct> &getStrawHitFlagEntries() const     {return _hitFlagEntries;}
  const std::vector<entryStruct> &getStrawHitPositionEntries() const {return _hitPositionEntries;}
  const std::string getSelectedStrawHitFlagEntry()      {return _selectedHitFlagEntry;}
  const std::string getSelectedStrawPositionFlagEntry() {return _selectedHitPositionEntry;}
  void setSelectedStrawHitFlagEntry(std::string selectedHitFlagEntry)
                                                         {_selectedHitFlagEntry=selectedHitFlagEntry;}
  void setSelectedStrawHitPositionEntry(std::string selectedHitPositionEntry)
                                                         {_selectedHitPositionEntry=selectedHitPositionEntry;}
  const mu2e::StrawHitFlagCollection *getStrawHitFlagCollection() const;
  const mu2e::StrawHitPositionCollection *getStrawHitPositionCollection() const;

  friend class FilterDialog;
};

}

#endif /* EventDisplay_src_ContentSelector_h */
