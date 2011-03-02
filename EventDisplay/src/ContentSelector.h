//
// Class which manages the combo boxes and list box in the event display frame. It is able to returns the data objects associated with the selected box entries. 
//
// $Id: ContentSelector.h,v 1.1 2011/03/02 03:25:47 ehrlich Exp $
// $Author: ehrlich $ 
// $Date: 2011/03/02 03:25:47 $
//
// Original author Ralf Ehrlich
//

#ifndef COMPONENTSELECTOR_H
#define COMPONENTSELECTOR_H

#include <iostream>
#include <vector>
#include <TGComboBox.h>
#include <TGListBox.h>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/SimParticleCollection.hh"
#include "ToyDP/inc/PointTrajectoryCollection.hh"
#include "ToyDP/inc/CaloCrystalHitCollection.hh"
#include "ToyDP/inc/PhysicalVolumeInfoCollection.hh"

namespace mu2e_eventdisplay
{

class ContentSelector
{
  std::vector<edm::Handle<mu2e::StepPointMCCollection> > _hitsVector;
  std::vector<edm::Handle<mu2e::CaloCrystalHitCollection> > _caloHitsVector;
  std::vector<edm::Handle<mu2e::SimParticleCollection> > _simParticlesVector;
  edm::Handle<mu2e::PhysicalVolumeInfoCollection> _physicalVolumes;
  edm::Handle<mu2e::PointTrajectoryCollection> _pointTrajectories;
  bool _hasPhysicalVolumes, _hasPointTrajectories;

  TGComboBox *_hitBox;
  TGComboBox *_caloHitBox;
  TGListBox  *_trackBox;

//these are information stored for the minimum hit test
  bool _hitsAreSelected;
  std::string _hitsClassName, _hitsModuleLabel, _hitsProductInstanceName;

  private:
  bool compareLists(const std::vector<std::string> &newContent, const TGListBox *boxContent) const;

  public:
  ContentSelector(TGComboBox *hitBox, TGComboBox *caloHitBox, TGListBox *trackBox);
  void firstLoop();
  void setAvailableCollections(const edm::Event& event);

  void setSelectedHitsName();
  bool getSelectedHitsName(std::string &className, 
                           std::string &moduleLabel, 
                           std::string &productInstanceName) const;

  const mu2e::StepPointMCCollection *getSelectedHitCollection() const;
  const mu2e::CaloCrystalHitCollection *getSelectedCaloHitCollection() const;
  std::vector<const mu2e::SimParticleCollection*> getSelectedTrackCollection() const;
  const mu2e::PhysicalVolumeInfoCollection *getPhysicalVolumeInfoCollection() const;
  const mu2e::PointTrajectoryCollection *getPointTrajectoryCollection() const;
};

}

#endif
