//
// Class which manages the combo boxes and list box in the event display frame. It is able to returns the data objects associated with the selected box entries.
//
// $Id: ContentSelector.h,v 1.7 2011/05/20 20:18:23 wb Exp $
// $Author: wb $
// $Date: 2011/05/20 20:18:23 $
//
// Original author Ralf Ehrlich
//

#ifndef EventDisplay_src_ContentSelector_h
#define EventDisplay_src_ContentSelector_h

#include "ToyDP/inc/CaloCrystalHitCollection.hh"
#include "ToyDP/inc/CaloHitCollection.hh"
#include "ToyDP/inc/PhysicalVolumeInfoCollection.hh"
#include "ToyDP/inc/PointTrajectoryCollection.hh"
#include "ToyDP/inc/SimParticleCollection.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/StrawHitCollection.hh"
#include "art/Framework/Core/Event.h"
#include "art/Framework/Core/Run.h"
#include <TGComboBox.h>
#include <TGListBox.h>
#include <iostream>
#include <vector>

namespace mu2e_eventdisplay
{

class ContentSelector
{
  ContentSelector();
  ContentSelector(const ContentSelector &);
  ContentSelector& operator=(const ContentSelector &);

  std::vector<art::Handle<mu2e::StepPointMCCollection> > _stepPointMCVector;
  std::vector<art::Handle<mu2e::StrawHitCollection> > _strawHitVector;
  std::vector<art::Handle<mu2e::CaloCrystalHitCollection> > _caloCrystalHitVector;
  std::vector<art::Handle<mu2e::CaloHitCollection> > _caloHitVector;
  std::vector<art::Handle<mu2e::SimParticleCollection> > _simParticlesVector;
  art::Handle<mu2e::PhysicalVolumeInfoCollection> _physicalVolumes;
  art::Handle<mu2e::PointTrajectoryCollection> _pointTrajectories;
  bool _hasPhysicalVolumes, _hasPointTrajectories;

  TGComboBox  *_hitBox;
  TGComboBox  *_caloHitBox;
  TGListBox   *_trackBox;
  std::string _g4ModuleLabel;

  struct entryStruct
  {
    int         entryID;
    std::string entryText;
  };

//these are information stored for the minimum hit test
  bool _hitsAreSelected;
  std::string _hitsClassName, _hitsModuleLabel, _hitsProductInstanceName;

  private:
  bool compareLists(const std::vector<entryStruct> &newEntries, const TGListBox *boxContent) const;
  template<class CollectionType> void createNewEntries(std::vector<art::Handle<CollectionType> > &dataVector,
                                                       const art::Event &event, const std::string &className,
                                                       std::vector<entryStruct> &newEntries, int entryIDStart);

  public:
  ContentSelector(TGComboBox *hitBox, TGComboBox *caloHitBox, TGListBox *trackBox, std::string const &g4ModuleLabel);
  void firstLoop();
  void setAvailableCollections(const art::Event& event);

  void setSelectedHitsName();
  bool getSelectedHitsName(std::string &className,
                           std::string &moduleLabel,
                           std::string &productInstanceName) const;

  template<typename CollectionType> const CollectionType* getSelectedHitCollection() const;
  template<typename CollectionType> const CollectionType* getSelectedCaloHitCollection() const;
  template<typename CollectionType> std::vector<const CollectionType*> getSelectedTrackCollection() const;
  const mu2e::PhysicalVolumeInfoCollection *getPhysicalVolumeInfoCollection() const;
  const mu2e::PointTrajectoryCollection *getPointTrajectoryCollection() const;
};

}

#endif /* EventDisplay_src_ContentSelector_h */
