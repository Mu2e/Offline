#ifndef G4TCaptureKineticsTable_h
#define G4TCaptureKineticsTable_h

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, Kevin Lynch, February 12 2010
// ----------------------------------------------------------------

#include "globals.hh"

#include "G4TCaptureKineticsChannel.hh"

#include <vector>

// FIXME add verboseLevel

// FIXME ... templatize so it works for MuAtom and MuMolecule, with
// typedefs in the appropriate headers.  also, update the clone, and
// add an IsApplicable which is called from Insert

template<class T> class G4TCaptureKineticsTable {
public:

  G4TCaptureKineticsTable(T const*);
  ~G4TCaptureKineticsTable();

  G4TCaptureKineticsTable* Clone(T const*) const;

  void Insert(G4TCaptureKineticsChannel<T>* aChannel);

  G4int entries() const { return channels.size(); }
  
  G4TCaptureKineticsChannel<T>* SelectAChannel() const;
  
  G4TCaptureKineticsChannel<T>* GetChannel(G4int index) const { return channels[index]; }
  G4TCaptureKineticsChannel<T>* operator[](G4int index) { return channels[index]; }

  G4double GetTotalRate() const { return totalRate; }
  
  // IMPLEMENT_ME
  void DumpInfo() const {}

private:
  G4TCaptureKineticsTable(G4TCaptureKineticsTable const&);
  G4TCaptureKineticsTable& operator=(G4TCaptureKineticsTable const&);

  T const* particle;

  typedef std::vector<G4TCaptureKineticsChannel<T>*> G4TCaptureKineticsChannelVector;
  G4TCaptureKineticsChannelVector channels;
  typedef typename G4TCaptureKineticsChannelVector::const_iterator const_iterator;
  std::vector<G4double> BRs;

  G4double totalRate;

  void swap(G4TCaptureKineticsTable &);
};

#define G4TCAPTUREKINETICSTABLE_ICC
#include "G4TCaptureKineticsTable.icc"
#undef G4TCAPTUREKINETICSTABLE_ICC

#endif

