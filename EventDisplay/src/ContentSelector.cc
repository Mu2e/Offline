#include "ContentSelector.h"

namespace mu2e_eventdisplay
{

ContentSelector::ContentSelector(TGComboBox *hitBox, TGComboBox *caloHitBox, TGListBox *trackBox, std::string const &g4ModuleLabel)
  :
  _hasPhysicalVolumes(false),
  _hasPointTrajectories(false),
  _hitBox(hitBox),
  _caloHitBox(caloHitBox),
  _trackBox(trackBox),
  _g4ModuleLabel(g4ModuleLabel),
  _hitsAreSelected(false)
  {}

void ContentSelector::firstLoop()  //This is useful for now, but may be changed later
{
  TGLBEntry *entry;

  entry=_hitBox->FindEntry("TrkRecoTrk:KalFitTest:");
  if(entry!=NULL) _hitBox->Select(entry->EntryId());
  else
  {
    entry=_hitBox->FindEntry("StrawHit:makeSH:");
    if(entry!=NULL) _hitBox->Select(entry->EntryId());
  }

  entry=_caloHitBox->FindEntry("CaloHit:CaloReadoutHitsMaker:");
  if(entry!=NULL) _caloHitBox->Select(entry->EntryId());

  entry=_trackBox->FindEntry("TrkRecoTrk:KalFitTest:");
  if(entry!=NULL) _trackBox->Select(entry->EntryId());
  else
  {
    entry=_trackBox->FindEntry("SimParticle:g4run:");
    if(entry!=NULL) _trackBox->Select(entry->EntryId());
  }
}

bool ContentSelector::compareLists(const std::vector<entryStruct> &newEntries, const TGListBox *boxContent) const
{
  if(static_cast<int>(newEntries.size())!=boxContent->GetNumberOfEntries()) return(false);

  for(unsigned int i=0; i<newEntries.size(); i++)
  {
    int entryID=newEntries[i].entryID;
    TGTextLBEntry *entryTmp=dynamic_cast<TGTextLBEntry*>(boxContent->GetEntry(entryID));
    if(entryTmp==NULL) return(false);
    const char *entry=entryTmp->GetText()->GetString();
    if(newEntries[i].entryText.compare(entry)!=0) return(false);
  }
  return(true);
}

template<class CollectionType>
void ContentSelector::createNewEntries(std::vector<art::Handle<CollectionType> > &dataVector,
                                       const art::Event &event, const std::string &className,
                                       std::vector<entryStruct> &newEntries, int entryIDStart)
{
  int entryID;
  event.getManyByType(dataVector);
  typedef std::vector<art::Handle<CollectionType> > CollectionVector;
  typedef typename CollectionVector::const_iterator itertype;
  itertype iter;
  for(iter=dataVector.begin(), entryID=entryIDStart;
      iter!=dataVector.end();
      iter++, entryID++)
  {
    entryStruct e;
    e.entryID=entryID;
    e.entryText=className+":"+iter->provenance()->moduleLabel()
               +":"+iter->provenance()->productInstanceName();
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
  createNewEntries<mu2e::StepPointMCCollection>(_stepPointMCVector, event, "StepPointMC", newEntries, 1000);
  createNewEntries<mu2e::StrawHitCollection>(_strawHitVector, event, "StrawHit", newEntries, 2000);
#ifdef BABARINSTALLED
  createNewEntries<mu2e::TrkRecoTrkCollection>(_hitOnTrackVector, event, "TrkRecoTrk", newEntries, 3000);
#endif
  newEntries.push_back(nothingSelected);

  if(compareLists(newEntries,_hitBox->GetListBox())==false)
  {
    TGTextLBEntry *selectedEntryTmp=dynamic_cast<TGTextLBEntry*>(_hitBox->GetSelectedEntry());
    std::string selectedEntry="Nothing Selected";
    if(selectedEntryTmp) selectedEntry=selectedEntryTmp->GetText()->GetString();
    _hitBox->RemoveAll();
    for(unsigned int i=0; i<newEntries.size(); i++)
    {
      _hitBox->AddEntry(newEntries[i].entryText.c_str(), newEntries[i].entryID);
      if(newEntries[i].entryText.compare(selectedEntry)==0) _hitBox->Select(newEntries[i].entryID);
    }
    _hitBox->GetListBox()->GetEntry(0)->SetBackgroundColor(0x00FF00);
  }

  newEntries.clear();

//Calo Hit Selection
  createNewEntries<mu2e::CaloCrystalHitCollection>(_caloCrystalHitVector, event, "CaloCrystalHit", newEntries, 1000);
  createNewEntries<mu2e::CaloHitCollection>(_caloHitVector, event, "CaloHit", newEntries, 2000);
  newEntries.push_back(nothingSelected);

  if(compareLists(newEntries,_caloHitBox->GetListBox())==false)
  {
    TGTextLBEntry *selectedEntryTmp=dynamic_cast<TGTextLBEntry*>(_caloHitBox->GetSelectedEntry());
    std::string selectedEntry="Nothing Selected";
    if(selectedEntryTmp) selectedEntry=selectedEntryTmp->GetText()->GetString();
    _caloHitBox->RemoveAll();
    for(unsigned int i=0; i<newEntries.size(); i++)
    {
      _caloHitBox->AddEntry(newEntries[i].entryText.c_str(), newEntries[i].entryID);
      if(newEntries[i].entryText.compare(selectedEntry)==0) _caloHitBox->Select(newEntries[i].entryID);
    }
    _caloHitBox->GetListBox()->GetEntry(0)->SetBackgroundColor(0x00FF00);
  }

  newEntries.clear();

//Track Selection
  createNewEntries<mu2e::SimParticleCollection>(_simParticleVector, event, "SimParticle", newEntries, 1000);
#ifdef BABARINSTALLED
  createNewEntries<mu2e::TrkRecoTrkCollection>(_trkRecoTrkVector, event, "TrkRecoTrk", newEntries, 2000);
#endif

  if(compareLists(newEntries,_trackBox)==false)
  {
    TList selections;
    _trackBox->GetSelectedEntries(&selections);
    std::vector<std::string> oldSelections;
    for(int i=0; i<selections.GetSize(); i++)
    {
      std::string selectedEntry=(dynamic_cast<TGTextLBEntry*>(selections.At(i)))->GetText()->GetString();
      oldSelections.push_back(selectedEntry);
    }
    _trackBox->RemoveAll();
    for(unsigned int i=0; i<newEntries.size(); i++)
    {
      _trackBox->AddEntry(newEntries[i].entryText.c_str(), newEntries[i].entryID);
      if(find(oldSelections.begin(),oldSelections.end(),newEntries[i].entryText)!=oldSelections.end())
           _trackBox->Select(newEntries[i].entryID);
    }
  }


//Other
  _hasPhysicalVolumes=event.getRun().getByLabel(_g4ModuleLabel, _physicalVolumes);
  _hasPointTrajectories=event.getByLabel(_g4ModuleLabel, _pointTrajectories);
}

void ContentSelector::setSelectedHitsName()
{
  _hitsAreSelected=false;
  int i=_hitBox->GetSelected();
  int classtype=i/1000;
  int index=i%1000;
  art::Provenance const *provenance=NULL;
  switch(classtype)
  {
    case 1 : if(index>=static_cast<int>(_stepPointMCVector.size())) return;
             _hitsAreSelected=true;
             provenance=_stepPointMCVector[index].provenance();
             break;
    case 2 : if(index>=static_cast<int>(_strawHitVector.size())) return;
             _hitsAreSelected=true;
             provenance=_strawHitVector[index].provenance();
             break;
#ifdef BABARINSTALLED
    case 3 : if(index>=static_cast<int>(_hitOnTrackVector.size())) return;
             _hitsAreSelected=true;
             provenance=_hitOnTrackVector[index].provenance();
             break;
#endif
    default: return;
  };
  _hitsClassName          =provenance->className();
  _hitsModuleLabel        =provenance->moduleLabel();
  _hitsProductInstanceName=provenance->productInstanceName();
}

bool ContentSelector::getSelectedHitsName(std::string &className,
                                          std::string &moduleLabel,
                                          std::string &productInstanceName) const
{
  className          =_hitsClassName;
  moduleLabel        =_hitsModuleLabel;
  productInstanceName=_hitsProductInstanceName;
  return(_hitsAreSelected);
}

template<typename CollectionType>
const CollectionType* ContentSelector::getSelectedHitCollection() const
{
  int i=_hitBox->GetSelected();
  int classtype=i/1000;
  int index=i%1000;
  switch(classtype)
  {
    case 1 : if(typeid(CollectionType)!=typeid(mu2e::StepPointMCCollection)) return(NULL);
             if(index>=static_cast<int>(_stepPointMCVector.size())) return(NULL);
             return(reinterpret_cast<const CollectionType*>(_stepPointMCVector[index].product()));
    case 2 : if(typeid(CollectionType)!=typeid(mu2e::StrawHitCollection)) return(NULL);
             if(index>=static_cast<int>(_strawHitVector.size())) return(NULL);
             return(reinterpret_cast<const CollectionType*>(_strawHitVector[index].product()));
#ifdef BABARINSTALLED
    case 3 : if(typeid(CollectionType)!=typeid(mu2e::TrkRecoTrkCollection)) return(NULL);
             if(index>=static_cast<int>(_hitOnTrackVector.size())) return(NULL);
             return(reinterpret_cast<const CollectionType*>(_hitOnTrackVector[index].product()));
#endif
//Note about the use of reinterpret_cast: While it is generally unsafe to use it, in this case it is Ok.
//the typeid check makes that the program advances to the line with the reinterpret_cast ONLY if the
//type of the vector element and the CollectionType are identical. The compiler doesn't see this, the compiler
//knows that this template function (and therefore CollectionType) gets implemented with different types (see below),
//which created an incompatibility between the return type of the function and the object which is returned.
//However, the compiler doesn't see that the object gets returned only if both types match.
//In order to satisfy the compiler, an reinterpret_cast is used, while during run time the argument type and the
//return type of the reinterpret_cast will always be the same, i.e. no "reinterpretation" will happen
  };
  return(NULL);
}
template const mu2e::StepPointMCCollection* ContentSelector::getSelectedHitCollection<mu2e::StepPointMCCollection>() const;
template const mu2e::StrawHitCollection*    ContentSelector::getSelectedHitCollection<mu2e::StrawHitCollection>() const;
#ifdef BABARINSTALLED
template const mu2e::TrkRecoTrkCollection*  ContentSelector::getSelectedHitCollection<mu2e::TrkRecoTrkCollection>() const;
#endif

template<typename CollectionType>
const CollectionType* ContentSelector::getSelectedCaloHitCollection() const
{
  int i=_caloHitBox->GetSelected();
  int classtype=i/1000;
  int index=i%1000;
  switch(classtype)
  {
    case 1 : if(typeid(CollectionType)!=typeid(mu2e::CaloCrystalHitCollection)) return(NULL);
             if(index>=static_cast<int>(_caloCrystalHitVector.size())) return(NULL);
             return(reinterpret_cast<const CollectionType*>(_caloCrystalHitVector[index].product()));
    case 2 : if(typeid(CollectionType)!=typeid(mu2e::CaloHitCollection)) return(NULL);
             if(index>=static_cast<int>(_caloHitVector.size())) return(NULL);
             return(reinterpret_cast<const CollectionType*>(_caloHitVector[index].product()));
  };
  return(NULL);
}
template const mu2e::CaloCrystalHitCollection* ContentSelector::getSelectedCaloHitCollection<mu2e::CaloCrystalHitCollection>() const;
template const mu2e::CaloHitCollection*    ContentSelector::getSelectedCaloHitCollection<mu2e::CaloHitCollection>() const;

template<typename CollectionType>
std::vector<const CollectionType*> ContentSelector::getSelectedTrackCollection() const
{
  std::vector<const CollectionType*> to_return;

  TList selections;
  _trackBox->GetSelectedEntries(&selections);
  for(int i=0; i<selections.GetSize(); i++)
  {
    TGTextLBEntry *entry=dynamic_cast<TGTextLBEntry*>(selections.At(i));
    if(entry==NULL) continue;
    int id=entry->EntryId();
    int classtype=id/1000;
    int index=id%1000;
    switch(classtype)
    {
      case 1 : if(typeid(CollectionType)!=typeid(mu2e::SimParticleCollection)) break;
               if(index>=static_cast<int>(_simParticleVector.size())) break;
               to_return.push_back(reinterpret_cast<const CollectionType*>(_simParticleVector[index].product())); 
               break;
#ifdef BABARINSTALLED
      case 2 : if(typeid(CollectionType)!=typeid(mu2e::TrkRecoTrkCollection)) break;
               if(index>=static_cast<int>(_trkRecoTrkVector.size())) break;
               to_return.push_back(reinterpret_cast<const CollectionType*>(_trkRecoTrkVector[index].product()));
               break;
#endif
    };
  }
  return(to_return);
}
template std::vector<const mu2e::SimParticleCollection*> ContentSelector::getSelectedTrackCollection<mu2e::SimParticleCollection>() const;
#ifdef BABARINSTALLED
template std::vector<const mu2e::TrkRecoTrkCollection*> ContentSelector::getSelectedTrackCollection<mu2e::TrkRecoTrkCollection>() const;
#endif

const mu2e::PhysicalVolumeInfoCollection* ContentSelector::getPhysicalVolumeInfoCollection() const
{
  if(_hasPhysicalVolumes) return(_physicalVolumes.product());
  else return(NULL);
}

const mu2e::PointTrajectoryCollection* ContentSelector::getPointTrajectoryCollection() const
{
  if(_hasPointTrajectories) return(_pointTrajectories.product());
  else return(NULL);
}

}
