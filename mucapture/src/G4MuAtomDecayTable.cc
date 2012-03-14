// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, cloned from G4DecayTable.cc
//      20 February 2010 Kevin Lynch
// ------------------------------------------------------------

#include "globals.hh"
#include "Randomize.hh"
#include <algorithm>

#include "G4MuAtomDecayTable.hh"
#include "G4MuAtom.hh"

G4MuAtomDecayTable::G4MuAtomDecayTable():parent(0)
{
  channels =  new G4VMuAtomDecayChannelVector;
}

G4MuAtomDecayTable::~G4MuAtomDecayTable()
{
  // remove and delete all contents  
  G4VMuAtomDecayChannelVector::iterator i;
  for (i = channels->begin(); i!= channels->end(); ++i) {
    delete (*i);
  }
  channels->clear();
  delete  channels;
  channels = 0;
}    

G4MuAtomDecayTable::G4MuAtomDecayTable(const G4MuAtomDecayTable & rhs):
  parent(rhs.parent), channels(new G4VMuAtomDecayChannelVector)
{
  G4VMuAtomDecayChannelVector::iterator i;
  for( i = rhs.channels->begin(); i != rhs.channels->end(); ++i){
    G4VMuAtomDecayChannel *chan = (*i)->Clone();
    Insert(chan);
  }
}

G4MuAtomDecayTable* G4MuAtomDecayTable::Clone(G4MuAtom const* muatom) const {
  G4MuAtomDecayTable* table = new G4MuAtomDecayTable(*this);
  G4VMuAtomDecayChannelVector *v = table->channels;
  G4VMuAtomDecayChannelVector::iterator i;
  for( i = v->begin(); i != v->end(); ++i )
    (*i)->SetParent(muatom->GetParticleName());
  return table;
}



G4MuAtomDecayTable& 
G4MuAtomDecayTable::operator=(const G4MuAtomDecayTable & rhs){
  G4MuAtomDecayTable lhs(rhs);
  swap(lhs);
  return *this;
}

void G4MuAtomDecayTable::Insert( G4VMuAtomDecayChannel * aChannel){
  if (parent == 0) { parent = (G4ParticleDefinition*)(aChannel->GetParent()); }
  if (parent != aChannel->GetParent()) {
#ifdef G4VERBOSE
    G4cout << " G4MuAtomDecayTable::Insert :: bad   G4VMuAtomDecayChannel (mismatch parent) ";
    G4cout << "       " << parent->GetParticleName();
    G4cout << " input:" << aChannel->GetParent()->GetParticleName() << G4endl;
#endif
  } else {
    G4double r = aChannel->GetBR();
    G4VMuAtomDecayChannelVector::iterator i;
    for (i = channels->begin(); i!= channels->end(); ++i) {
      if (r > (*i)->GetBR()) {
	channels->insert(i,aChannel);
	return;
      }
    }
    channels->push_back(aChannel);
  }
}

G4VMuAtomDecayChannel *G4MuAtomDecayTable::SelectADecayChannel()
{
  // check if contents exist
  if (channels->size()<1) return 0;

  while (1) {
    G4double sumBR = 0.0;
    G4double r= G4UniformRand();
    // select decay channel
    G4VMuAtomDecayChannelVector::iterator i;
    for (i = channels->begin(); i!= channels->end(); ++i) {
      sumBR += (*i)->GetBR();
      if (r < sumBR) {
	return (*i);
      }
    }
  }
  return 0;
}

void G4MuAtomDecayTable::DumpInfo() const
{
  G4cout << "G4MuAtomDecayTable:  " << parent->GetParticleName() << G4endl;
  G4int index =0;
  G4VMuAtomDecayChannelVector::iterator i;
  for (i = channels->begin(); i!= channels->end(); ++i) {
    G4cout << index << ": ";
    (*i)->DumpInfo();
    index +=1;
  }
  G4cout << G4endl;
}

void G4MuAtomDecayTable::swap(G4MuAtomDecayTable & rhs){
  std::swap(parent, rhs.parent);
  std::swap(channels, rhs.channels);
}









