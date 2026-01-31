#include "Offline/EventDisplay/src/ContentSelector.h"
#include "art/Framework/Principal/SubRun.h"

namespace mu2e_eventdisplay
{

ContentSelector::ContentSelector(TGComboBox *hitBox, TGComboBox *caloHitBox, TGComboBox *crvHitBox, TGListBox *trackBox, std::string const &g4ModuleLabel, std::string const &physicalVolumesMultiLabel, std::string const &protonBunchTimeLabel)
  :
  _hitBox(hitBox),
  _caloHitBox(caloHitBox),
  _crvHitBox(crvHitBox),
  _trackBox(trackBox),
  _g4ModuleLabel(g4ModuleLabel),
  _physicalVolumesMultiLabel(physicalVolumesMultiLabel),
  _protonBunchTimeLabel(protonBunchTimeLabel)
  {}

void ContentSelector::firstLoop()  //This is useful for now, but may be changed later
{
  TGLBEntry *entry;

  entry=_hitBox->FindEntry("StrawHit:makeSH:");
  if(entry==nullptr) entry=_hitBox->FindEntry("StrawHit:makeDcH:");
  if(entry!=nullptr) _hitBox->Select(entry->EntryId());

  entry=_caloHitBox->FindEntry("CaloHit:CaloHitsMaker:");
  if(entry!=nullptr) _caloHitBox->Select(entry->EntryId());

  entry=_trackBox->FindEntry("SimParticle:g4run:");
  if(entry!=nullptr) _trackBox->Select(entry->EntryId());

  _selectedHitFlagEntry="Nothing Selected";
  _selectedHitPositionEntry="Nothing Selected";
}

template<class CollectionType>
void ContentSelector::createNewEntries(std::vector<art::Handle<CollectionType> > &dataVector,
                                       const art::Event &event, const std::string &className,
                                       std::vector<entryStruct> &newEntries, int classID)
{
  dataVector=event.getMany<CollectionType>();
  typedef std::vector<art::Handle<CollectionType> > CollectionVector;
  typedef typename CollectionVector::const_iterator itertype;
  itertype iter;
  int vectorPos=0;
  for(iter=dataVector.begin(); iter!=dataVector.end(); iter++, vectorPos++)
  {
    entryStruct e;
    e.entryID=newEntries.size();
    e.classID=classID;
    e.vectorPos=vectorPos;
    e.entryText=className+":"+iter->provenance()->moduleLabel()
               +":"+iter->provenance()->productInstanceName();
    e.className=iter->provenance()->producedClassName();
    e.moduleLabel=iter->provenance()->moduleLabel();
    e.productInstanceName=iter->provenance()->productInstanceName();
    newEntries.push_back(e);
  }
}

void ContentSelector::setAvailableCollections(const art::Event& event)
{
  std::vector<entryStruct> newEntries;

  entryStruct nothingSelected;
  nothingSelected.entryID=0;
  nothingSelected.entryText="Nothing Selected";

//Hit Selection
  newEntries.push_back(nothingSelected);
  createNewEntries<mu2e::StepPointMCCollection>(_stepPointMCVector, event, "StepPointMC", newEntries, 1);
  createNewEntries<mu2e::StrawHitCollection>(_strawHitVector, event, "StrawHit", newEntries, 2);

  if(newEntries!=_hitEntries)
  {
    _hitEntries.clear();
    TGTextLBEntry *selectedEntryTmp=dynamic_cast<TGTextLBEntry*>(_hitBox->GetSelectedEntry());
    std::string selectedEntry="Nothing Selected";
    if(selectedEntryTmp) selectedEntry=selectedEntryTmp->GetText()->GetString();
    _hitBox->RemoveAll();
    for(unsigned int i=0; i<newEntries.size(); i++)
    {
      _hitEntries.push_back(newEntries[i]);
      _hitBox->AddEntry(newEntries[i].entryText.c_str(), newEntries[i].entryID);
      if(newEntries[i].entryText.compare(selectedEntry)==0) _hitBox->Select(newEntries[i].entryID);
    }
    _hitBox->GetListBox()->GetEntry(0)->SetBackgroundColor(0x00FF00);
  }

//Hit Flag Selection
  _hitFlagEntries.clear();
  _hitFlagEntries.push_back(nothingSelected);
  createNewEntries<mu2e::StrawHitFlagCollection>(_strawHitFlagVector, event, "StrawHitFlag", _hitFlagEntries, 1);

//Hit Position Selection
  _hitPositionEntries.clear();
  _hitPositionEntries.push_back(nothingSelected);
  createNewEntries<mu2e::StrawHitPositionCollection>(_strawHitPositionVector, event, "StrawHitPosition", _hitPositionEntries, 1);

//Calo Hit Selection
  newEntries.clear();
  newEntries.push_back(nothingSelected);
  createNewEntries<mu2e::StepPointMCCollection>(_caloStepPointMCVector, event, "StepPointMC", newEntries, 1);
  createNewEntries<mu2e::CaloHitCollection>(_caloHitVector, event, "CaloHit", newEntries, 3);

  if(newEntries!=_caloHitEntries)
  {
    _caloHitEntries.clear();
    TGTextLBEntry *selectedEntryTmp=dynamic_cast<TGTextLBEntry*>(_caloHitBox->GetSelectedEntry());
    std::string selectedEntry="Nothing Selected";
    if(selectedEntryTmp) selectedEntry=selectedEntryTmp->GetText()->GetString();
    _caloHitBox->RemoveAll();
    for(unsigned int i=0; i<newEntries.size(); i++)
    {
      _caloHitEntries.push_back(newEntries[i]);
      _caloHitBox->AddEntry(newEntries[i].entryText.c_str(), newEntries[i].entryID);
      if(newEntries[i].entryText.compare(selectedEntry)==0) _caloHitBox->Select(newEntries[i].entryID);
    }
    _caloHitBox->GetListBox()->GetEntry(0)->SetBackgroundColor(0x00FF00);
  }

//CRV Hit Selection
  newEntries.clear();
  newEntries.push_back(nothingSelected);
  createNewEntries<mu2e::CrvRecoPulseCollection>(_crvRecoPulseVector, event, "CrvRecoPulse", newEntries, 1);

  if(newEntries!=_crvHitEntries)
  {
    _crvHitEntries.clear();
    TGTextLBEntry *selectedEntryTmp=dynamic_cast<TGTextLBEntry*>(_crvHitBox->GetSelectedEntry());
    std::string selectedEntry="Nothing Selected";
    if(selectedEntryTmp) selectedEntry=selectedEntryTmp->GetText()->GetString();
    _crvHitBox->RemoveAll();
    for(unsigned int i=0; i<newEntries.size(); i++)
    {
      _crvHitEntries.push_back(newEntries[i]);
      _crvHitBox->AddEntry(newEntries[i].entryText.c_str(), newEntries[i].entryID);
      if(newEntries[i].entryText.compare(selectedEntry)==0) _crvHitBox->Select(newEntries[i].entryID);
    }
    _crvHitBox->GetListBox()->GetEntry(0)->SetBackgroundColor(0x00FF00);
  }

//CRV Waveforms (not inside a menu)
  _crvDigisVector=event.getMany<mu2e::CrvDigiCollection>();

//Track Selection
  newEntries.clear();
  createNewEntries<mu2e::SimParticleCollection>(_simParticleVector, event, "SimParticle", newEntries, 1);
  createNewEntries<mu2e::KalSeedCollection>(_kalSeedTrkVector, event, "KalSeed", newEntries, 4);

  if(newEntries!=_trackEntries)
  {
    _trackEntries.clear();
    TList selections;
    _trackBox->GetSelectedEntries(&selections);
    std::vector<std::string> oldSelections;
    for(int i=0; i<selections.GetSize(); i++)
    {
      std::string selectedEntry=(dynamic_cast<TGTextLBEntry*>(selections.At(i)))->GetText()->GetString();
      oldSelections.push_back(selectedEntry);
    }
    selections.RemoveAll();
    _trackBox->RemoveAll();
    for(unsigned int i=0; i<newEntries.size(); i++)
    {
      _trackEntries.push_back(newEntries[i]);
      _trackBox->AddEntry(newEntries[i].entryText.c_str(), newEntries[i].entryID);
      if(find(oldSelections.begin(),oldSelections.end(),newEntries[i].entryText)!=oldSelections.end())
           _trackBox->Select(newEntries[i].entryID);
    }
  }

//MCTrajectories
  _mcTrajectoryVector=event.getMany<mu2e::MCTrajectoryCollection>();

//Other
  _hasPhysicalVolumesMulti=event.getSubRun().getByLabel(_physicalVolumesMultiLabel, _physicalVolumesMulti);
  event.getByLabel(_protonBunchTimeLabel, _protonBunchTime);
}

bool ContentSelector::getSelectedHitsName(std::string &className,
                                          std::string &moduleLabel,
                                          std::string &productInstanceName) const
{
  int i=_hitBox->GetSelected();
  std::vector<entryStruct>::const_iterator iter;
  for(iter=_hitEntries.begin(); iter!=_hitEntries.end(); iter++)
  {
    if(iter->entryID==i)
    {
      className          =iter->className;
      moduleLabel        =iter->moduleLabel;
      productInstanceName=iter->productInstanceName;
      return(true);
    }
  }
  return(false);
}

std::vector<ContentSelector::trackInfoStruct> ContentSelector::getSelectedTrackNames() const
{
  std::vector<trackInfoStruct> to_return;
  TList selections;
  _trackBox->GetSelectedEntries(&selections);
  for(int i=0; i<selections.GetSize(); i++)
  {
    TGTextLBEntry *entry=dynamic_cast<TGTextLBEntry*>(selections.At(i));
    if(entry==nullptr) continue;
    trackInfoStruct t;
    std::vector<entryStruct>::const_iterator iter;
    for(iter=_trackEntries.begin(); iter!=_trackEntries.end(); iter++)
    {
      if(iter->entryID==entry->EntryId())
      {
        t.classID=iter->classID;
        t.index=iter->vectorPos;
        t.entryText=iter->entryText;
        t.moduleLabel=iter->moduleLabel;
        t.productInstanceName=iter->productInstanceName;
        if(t.classID>0) to_return.push_back(t);
        break;
      }
    }
  }
  return(to_return);
}

template<typename CollectionType>
const CollectionType* ContentSelector::getSelectedHitCollection() const
{
  int i=_hitBox->GetSelected();
  int classID=0;
  int index=0;
  std::vector<entryStruct>::const_iterator iter;
  for(iter=_hitEntries.begin(); iter!=_hitEntries.end(); iter++)
  {
    if(iter->entryID==i)
    {
      classID=iter->classID;
      index=iter->vectorPos;
      break;
    }
  }
  switch(classID)
  {
    case 1 : if(typeid(CollectionType)!=typeid(mu2e::StepPointMCCollection)) return(nullptr);
             if(index>=static_cast<int>(_stepPointMCVector.size())) return(nullptr);
             return(reinterpret_cast<const CollectionType*>(_stepPointMCVector[index].product()));
    case 2 : if(typeid(CollectionType)!=typeid(mu2e::StrawHitCollection)) return(nullptr);
             if(index>=static_cast<int>(_strawHitVector.size())) return(nullptr);
             return(reinterpret_cast<const CollectionType*>(_strawHitVector[index].product()));
    case 4 : if(typeid(CollectionType)!=typeid(mu2e::KalSeedCollection)) return(nullptr);
             if(index>=static_cast<int>(_kalSeedHitVector.size())) return(nullptr);
             return(reinterpret_cast<const CollectionType*>(_kalSeedHitVector[index].product()));
//Note about the use of reinterpret_cast: While it is generally unsafe to use it, in this case it is Ok.
//the typeid check makes that the program advances to the line with the reinterpret_cast ONLY if the
//type of the vector element and the CollectionType are identical. The compiler doesn't see this, the compiler
//knows that this template function (and therefore CollectionType) gets implemented with different types (see below),
//which created an incompatibility between the return type of the function and the object which is returned.
//However, the compiler doesn't see that the object gets returned only if both types match.
//In order to satisfy the compiler, an reinterpret_cast is used, while during run time the argument type and the
//return type of the reinterpret_cast will always be the same, i.e. no "reinterpretation" will happen
  };
  return(nullptr);
}
template const mu2e::StepPointMCCollection* ContentSelector::getSelectedHitCollection<mu2e::StepPointMCCollection>() const;
template const mu2e::StrawHitCollection*    ContentSelector::getSelectedHitCollection<mu2e::StrawHitCollection>() const;
template const mu2e::KalSeedCollection*  ContentSelector::getSelectedHitCollection<mu2e::KalSeedCollection>() const;

template<typename CollectionType>
const CollectionType* ContentSelector::getSelectedCaloHitCollection() const
{
  int i=_caloHitBox->GetSelected();
  int classID=0;
  int index=0;
  std::vector<entryStruct>::const_iterator iter;
  for(iter=_caloHitEntries.begin(); iter!=_caloHitEntries.end(); iter++)
  {
    if(iter->entryID==i)
    {
      classID=iter->classID;
      index=iter->vectorPos;
      break;
    }
  }
  switch(classID)
  {
    case 1 : if(typeid(CollectionType)!=typeid(mu2e::StepPointMCCollection)) return(nullptr);
             if(index>=static_cast<int>(_caloStepPointMCVector.size())) return(nullptr);
             return(reinterpret_cast<const CollectionType*>(_caloStepPointMCVector[index].product()));
    case 3 : if(typeid(CollectionType)!=typeid(mu2e::CaloHitCollection)) return(nullptr);
             if(index>=static_cast<int>(_caloHitVector.size())) return(nullptr);
             return(reinterpret_cast<const CollectionType*>(_caloHitVector[index].product()));
  };
  return(nullptr);
}
template const mu2e::StepPointMCCollection* ContentSelector::getSelectedCaloHitCollection<mu2e::StepPointMCCollection>() const;
template const mu2e::CaloHitCollection*    ContentSelector::getSelectedCaloHitCollection<mu2e::CaloHitCollection>() const;


template<typename CollectionType>
const CollectionType* ContentSelector::getSelectedCrvHitCollection() const
{
  int i=_crvHitBox->GetSelected();
  int classID=0;
  int index=0;
  std::vector<entryStruct>::const_iterator iter;
  for(iter=_crvHitEntries.begin(); iter!=_crvHitEntries.end(); iter++)
  {
    if(iter->entryID==i)
    {
      classID=iter->classID;
      index=iter->vectorPos;
      break;
    }
  }
  switch(classID)
  {
    case 1 : if(typeid(CollectionType)!=typeid(mu2e::CrvRecoPulseCollection)) return(nullptr);
             if(index>=static_cast<int>(_crvRecoPulseVector.size())) return(nullptr);
             return(reinterpret_cast<const CollectionType*>(_crvRecoPulseVector[index].product()));
  };
  return(nullptr);
}
template const mu2e::CrvRecoPulseCollection* ContentSelector::getSelectedCrvHitCollection<mu2e::CrvRecoPulseCollection>() const;

const std::vector<art::Handle<mu2e::CrvDigiCollection> > &ContentSelector::getSelectedCrvDigiCollection() const
{
  return _crvDigisVector;
}

template<typename CollectionType>
std::vector<const CollectionType*> ContentSelector::getSelectedTrackCollection(std::vector<trackInfoStruct> &v) const
{
  std::vector<const CollectionType*> to_return;

  TList selections;
  _trackBox->GetSelectedEntries(&selections);
  for(int i=0; i<selections.GetSize(); i++)
  {
    TGTextLBEntry *entry=dynamic_cast<TGTextLBEntry*>(selections.At(i));
    if(entry==nullptr) continue;
    trackInfoStruct t;
    int classID=0;
    int index=0;
    std::vector<entryStruct>::const_iterator iter;
    for(iter=_trackEntries.begin(); iter!=_trackEntries.end(); iter++)
    {
      if(iter->entryID==entry->EntryId())
      {
        classID=iter->classID;
        index=iter->vectorPos;
        t.classID=classID;
        t.index=index;
        t.entryText=iter->entryText;
        t.moduleLabel=iter->moduleLabel;
        t.productInstanceName=iter->productInstanceName;
        break;
      }
    }
    switch(classID)
    {
      case 1 : if(typeid(CollectionType)!=typeid(mu2e::SimParticleCollection)) break;
               if(index>=static_cast<int>(_simParticleVector.size())) break;
               to_return.push_back(reinterpret_cast<const CollectionType*>(_simParticleVector[index].product()));
               t.productId=_simParticleVector[index].id();
               v.push_back(t);
               break;
      case 4 : if(typeid(CollectionType)!=typeid(mu2e::KalSeedCollection)) break;
               if(index>=static_cast<int>(_kalSeedTrkVector.size())) break;
               to_return.push_back(reinterpret_cast<const CollectionType*>(_kalSeedTrkVector[index].product()));
               t.productId=_kalSeedTrkVector[index].id();
               v.push_back(t);
               break;
    };
  }
  return(to_return);
}
template std::vector<const mu2e::SimParticleCollection*> ContentSelector::getSelectedTrackCollection<mu2e::SimParticleCollection>(std::vector<trackInfoStruct> &v) const;
template std::vector<const mu2e::KalSeedCollection*> ContentSelector::getSelectedTrackCollection<mu2e::KalSeedCollection>(std::vector<trackInfoStruct> &v) const;

const mu2e::PhysicalVolumeInfoMultiCollection* ContentSelector::getPhysicalVolumeInfoMultiCollection() const
{
  if(_hasPhysicalVolumesMulti) return(_physicalVolumesMulti.product());
  else return(nullptr);
}

const double ContentSelector::getTDC0time() const
{
  if(!_protonBunchTime.isValid())
  {
    std::cout<<"ProtonBunchTime not found! Setting it to 0."<<std::endl;
    return 0;
  }
  double TDC0time = -_protonBunchTime->pbtime_;
  return TDC0time;
}

const mu2e::MCTrajectoryCollection* ContentSelector::getMCTrajectoryCollection(const trackInfoStruct &t) const
{
  typedef std::vector<art::Handle<mu2e::MCTrajectoryCollection> > CollectionVector;
  typedef typename CollectionVector::const_iterator itertype;
  itertype iter;
  for(iter=_mcTrajectoryVector.begin(); iter!=_mcTrajectoryVector.end(); iter++)
  {
//    if(t.moduleLabel==iter->provenance()->moduleLabel() && t.productInstanceName==iter->provenance()->productInstanceName())
    if(t.moduleLabel==iter->provenance()->moduleLabel()) //TODO: why no check for product instance name?
      return(iter->product());
  }
  return(nullptr);
}

const mu2e::StrawHitFlagCollection* ContentSelector::getStrawHitFlagCollection() const
{
  int classID=0;
  int index=0;
  std::vector<entryStruct>::const_iterator iter;
  for(iter=_hitFlagEntries.begin(); iter!=_hitFlagEntries.end(); iter++)
  {
    if(iter->entryText==_selectedHitFlagEntry)
    {
      classID=iter->classID;
      index=iter->vectorPos;
      break;
    }
  }
  if(classID==1 && index<static_cast<int>(_strawHitFlagVector.size()))
  {
    return(_strawHitFlagVector[index].product());
  }
  return(nullptr);
}

const mu2e::StrawHitPositionCollection* ContentSelector::getStrawHitPositionCollection() const
{
  return(nullptr);
}

}
