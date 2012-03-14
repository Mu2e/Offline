#ifndef G4MuAtomDecayTable_h
#define G4MuAtomDecayTable_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, cloned from G4DecayTable.hh
//      20 February 2010  Kevin Lynch
// ------------------------------------------------------------

#include "G4ios.hh"
#include <vector>
#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4VMuAtomDecayChannel.hh"

class G4MuAtom;

class G4MuAtomDecayTable 
{
 // Class Description
 //   G4MuAtomDecayTable is the table of pointer to G4VMuAtomDecayChannel.
 //   Decay channels inside is sorted by using decay branching ratio
 //
  
public:
  typedef std::vector<G4VMuAtomDecayChannel*> G4VMuAtomDecayChannelVector;
  
  //constructors
public:
  G4MuAtomDecayTable();
  ~G4MuAtomDecayTable();
  
public:
  // FIXME ... reexamine the behavior of these member functions in
  // light of the Clone fixes to G4TCaptureKineticsChannel.  Do I even
  // need this functionality anymore?  If the answer is NO, then I
  // need to modify the DIO channel, and I need to jetison the new
  // decay table ... all of which would be nice, as it would eliminate
  // a big cut and paste headache
  G4MuAtomDecayTable* Clone(G4MuAtom const*) const;
  G4MuAtomDecayTable(const G4MuAtomDecayTable &);
  G4MuAtomDecayTable & operator=(const G4MuAtomDecayTable & rhs);
  
public:
  // equality operators
  G4int operator==(const G4MuAtomDecayTable &right) const {return (this == &right);};
  G4int operator!=(const G4MuAtomDecayTable &right) const {return (this != &right);};
  
public: // With Description
  void  Insert( G4VMuAtomDecayChannel* aChannel);
  // Insert a decay channel at proper position 
  // (i.e. sorted by using branching ratio ) 
  
  G4int entries() const;
  // Returns number of decay channels inside
  
public: // With Description
  G4VMuAtomDecayChannel* SelectADecayChannel();
  // A decay channel is selected at random according to the branching ratio 
  
  G4VMuAtomDecayChannel* GetDecayChannel(G4int index) const;
  G4VMuAtomDecayChannel* operator[](G4int index);
  // Get index-th Decay channel
  
  void DumpInfo() const;
  
private:
  G4ParticleDefinition       *parent;
  G4VMuAtomDecayChannelVector       *channels;

  void swap(G4MuAtomDecayTable &);
};

inline     
G4int G4MuAtomDecayTable::entries() const
{
  return channels->size();
}

inline     
G4VMuAtomDecayChannel* G4MuAtomDecayTable::operator[](G4int index)
{
  return (*channels)[index];
}


inline     
G4VMuAtomDecayChannel* G4MuAtomDecayTable::GetDecayChannel(G4int index) const
{
  G4VMuAtomDecayChannel* selectedChannel = 0;
  if ( (index>=0) && (index<G4int(channels->size())) ){
    selectedChannel = (*channels)[index];
  }
  return selectedChannel;
}


#endif
