#include "ContentSelector.h"

namespace mu2e_eventdisplay
{

ContentSelector::ContentSelector(TGComboBox *hitBox, TGComboBox *caloHitBox, TGListBox *trackBox)
{
  _hitBox=hitBox;
  _caloHitBox=caloHitBox;
  _trackBox=trackBox;
  _hitsAreSelected=false;
  _hasPhysicalVolumes=false;
  _hasPointTrajectories=false;
}

void ContentSelector::firstLoop()  //This is useful for now, but may be changed later on
{
  for(int i=0; i<_hitBox->GetNumberOfEntries(); i++)
  {
    const char *entry=(dynamic_cast<TGTextLBEntry*>(_hitBox->GetListBox()->GetEntry(i)))->GetText()->GetString();
    if(strcmp(entry,"StepPointMC:g4run:tracker")==0) _hitBox->Select(i);
  }
  for(int i=0; i<_caloHitBox->GetNumberOfEntries(); i++)
  {
    const char *entry=(dynamic_cast<TGTextLBEntry*>(_caloHitBox->GetListBox()->GetEntry(i)))->GetText()->GetString();
    if(strcmp(entry,"CaloCrystalHit:CaloCrystalHitsMaker:")==0) _caloHitBox->Select(i);
  }
  for(int i=0; i<_trackBox->GetNumberOfEntries(); i++)
  {
    const char *entry=(dynamic_cast<TGTextLBEntry*>(_trackBox->GetEntry(i)))->GetText()->GetString();
    if(strcmp(entry,"SimParticle:g4run:")==0) _trackBox->Select(i);
  }
}

bool ContentSelector::compareLists(const std::vector<std::string> &newContent, const TGListBox *boxContent) const
{
  for(unsigned int i=0; i<newContent.size(); i++)
  {
    TGTextLBEntry *entryTmp=dynamic_cast<TGTextLBEntry*>(boxContent->GetEntry(i));
    if(entryTmp==NULL) return(false);
    const char *entry=entryTmp->GetText()->GetString();
    if(newContent[i].compare(entry)!=0) return(false);
  }
  return(true);
}

void ContentSelector::setAvailableCollections(const edm::Event& event)
{
  std::vector<std::string> newContent;

//Hit Selection
  event.getManyByType(_hitsVector);
  std::vector<edm::Handle<mu2e::StepPointMCCollection> >::const_iterator hitsIter;
  for(hitsIter=_hitsVector.begin(); hitsIter!=_hitsVector.end(); hitsIter++)
  {
    std::string s="StepPointMC:"+hitsIter->provenance()->moduleLabel()
                 +":"+hitsIter->provenance()->productInstanceName();
    newContent.push_back(s);
  }
  newContent.push_back("Nothing Selected");

  if(compareLists(newContent,_hitBox->GetListBox())==false)
  {
    TGTextLBEntry *selectedEntryTmp=dynamic_cast<TGTextLBEntry*>(_hitBox->GetSelectedEntry());
    std::string selectedEntry="Nothing Selected";
    if(selectedEntryTmp) selectedEntry=selectedEntryTmp->GetText()->GetString();
    _hitBox->RemoveAll();
    for(unsigned int i=0; i<newContent.size(); i++)
    {
      _hitBox->AddEntry(newContent[i].c_str(), i);
      if(newContent[i].compare(selectedEntry)==0) _hitBox->Select(i);
    }
    _hitBox->GetListBox()->GetEntry(newContent.size()-1)->SetBackgroundColor(0x00FF00);
  }

  newContent.clear();

//Calo Hit Selection
  event.getManyByType(_caloHitsVector);
  std::vector<edm::Handle<mu2e::CaloCrystalHitCollection> >::const_iterator caloHitsIter;
  for(caloHitsIter=_caloHitsVector.begin(); caloHitsIter!=_caloHitsVector.end(); caloHitsIter++)
  {
    std::string s="CaloCrystalHit:"+caloHitsIter->provenance()->moduleLabel()
                 +":"+caloHitsIter->provenance()->productInstanceName();
    newContent.push_back(s);
  }
  newContent.push_back("Nothing Selected");

  if(compareLists(newContent,_caloHitBox->GetListBox())==false)
  {
    TGTextLBEntry *selectedEntryTmp=dynamic_cast<TGTextLBEntry*>(_caloHitBox->GetSelectedEntry());
    std::string selectedEntry="Nothing Selected";
    if(selectedEntryTmp) selectedEntry=selectedEntryTmp->GetText()->GetString();
    _caloHitBox->RemoveAll();
    for(unsigned int i=0; i<newContent.size(); i++)
    {
      _caloHitBox->AddEntry(newContent[i].c_str(), i);
      if(newContent[i].compare(selectedEntry)==0) _caloHitBox->Select(i);
    }
    _caloHitBox->GetListBox()->GetEntry(newContent.size()-1)->SetBackgroundColor(0x00FF00);
  }

  newContent.clear();

//Track Selection
  event.getManyByType(_simParticlesVector);
  std::vector<edm::Handle<mu2e::SimParticleCollection> >::const_iterator simParticlesIter;
  for(simParticlesIter=_simParticlesVector.begin(); simParticlesIter!=_simParticlesVector.end(); simParticlesIter++)
  {
    std::string s="SimParticle:"+simParticlesIter->provenance()->moduleLabel()
                 +":"+simParticlesIter->provenance()->productInstanceName();
    newContent.push_back(s);
  }
  if(compareLists(newContent,_trackBox)==false)
  {
    TList *selections = new TList;
    _trackBox->GetSelectedEntries(selections);
    std::vector<std::string> oldSelections;
    for(int i=0; i<selections->GetSize(); i++)
    {
      std::string selectedEntry=(dynamic_cast<TGTextLBEntry*>(selections->At(i)))->GetText()->GetString();
      oldSelections.push_back(selectedEntry);
    }
    _trackBox->RemoveAll();
    for(unsigned int i=0; i<newContent.size(); i++)
    {
      _trackBox->AddEntry(newContent[i].c_str(), i);
      if(find(oldSelections.begin(),oldSelections.end(),newContent[i])!=oldSelections.end()) _trackBox->Select(i);
    }
  }


//Other
  _hasPhysicalVolumes=event.getRun().getByType(_physicalVolumes);
  _hasPointTrajectories=event.getByType(_pointTrajectories);
}

void ContentSelector::setSelectedHitsName()
{
  _hitsAreSelected        =false;
  int i=_hitBox->GetSelected();
  if(i<0 || i>=static_cast<int>(_hitsVector.size())) return;
  _hitsAreSelected        =true;
  _hitsClassName          =_hitsVector[i].provenance()->className();
  _hitsModuleLabel        =_hitsVector[i].provenance()->moduleLabel();
  _hitsProductInstanceName=_hitsVector[i].provenance()->productInstanceName();
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

const mu2e::StepPointMCCollection* ContentSelector::getSelectedHitCollection() const
{
  unsigned int i=_hitBox->GetSelected();
  if(i>=_hitsVector.size()) return(NULL);
  else return(_hitsVector[i].product());
}

const mu2e::CaloCrystalHitCollection* ContentSelector::getSelectedCaloHitCollection() const
{
  unsigned int i=_caloHitBox->GetSelected();
  if(i>=_caloHitsVector.size()) return(NULL);
  else return(_caloHitsVector[i].product());
}

std::vector<const mu2e::SimParticleCollection*> ContentSelector::getSelectedTrackCollection() const
{
  std::vector<const mu2e::SimParticleCollection*> to_return;

  TList *selections = new TList;
  _trackBox->GetSelectedEntries(selections);
  for(int i=0; i<selections->GetSize(); i++)
  {
    const char *entry=(dynamic_cast<TGTextLBEntry*>(selections->At(i)))->GetText()->GetString();
    std::vector<edm::Handle<mu2e::SimParticleCollection> >::const_iterator simParticlesIter;
    for(simParticlesIter=_simParticlesVector.begin(); 
        simParticlesIter!=_simParticlesVector.end(); 
        simParticlesIter++)
    {
      std::string s="SimParticle:"+simParticlesIter->provenance()->moduleLabel()
                   +":"+simParticlesIter->provenance()->productInstanceName();
      if(s.compare(entry)==0) to_return.push_back(simParticlesIter->product());
    }
  }

  return(to_return);
}

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
