#ifndef G4MuAtomTable_h
#define G4MuAtomTable_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, Kevin Lynch, February 12 2010
// ----------------------------------------------------------------

#include "globals.hh"
#include "G4MuAtom.hh"
#include "G4MuonMinusAtomicCaptureCascadeModel.hh"
#include "G4MuonMinusAtomicCaptureStateModel.hh"
#include "G4MuonMinusAtomicCaptureChargeModel.hh"
#include "G4MuAtomLifetimeModel.hh"
#include "G4MuAtomCaptureRateModel.hh"

#include <map>

class G4MuAtomTable {

  // Class Description:

  // This is a singleton class that collaborates with the G4MuAtom
  // class to track the current G4MuAtom instances

public:

  static G4MuAtomTable* GetInstance();

  // does not create
  G4MuAtom* FindMuAtom(G4int Z, G4int A, G4int iSpin = 0);
  // does create
  G4MuAtom* GetMuAtom(G4int Z, G4int A, G4int iSpin = 0);

  G4MuonMinusAtomicCaptureCascadeModel*
  CascadeModel(G4int Z, G4int A, G4MuonMinusAtomicCaptureCascadeModel*);
  G4MuonMinusAtomicCaptureCascadeModel*
  CascadeModel(G4int Z, G4int A) const;

  G4MuonMinusAtomicCaptureStateModel*
  StateModel(G4int Z, G4int A, G4MuonMinusAtomicCaptureStateModel*);
  G4MuonMinusAtomicCaptureStateModel*
  StateModel(G4int Z, G4int A) const;

  G4MuonMinusAtomicCaptureChargeModel*
  ChargeModel(G4int Z, G4int A, G4MuonMinusAtomicCaptureChargeModel*);
  G4MuonMinusAtomicCaptureChargeModel*
  ChargeModel(G4int Z, G4int A) const;

  G4MuAtomLifetimeModel*
  LifetimeModel(G4int Z, G4int A, G4MuAtomLifetimeModel*);
  G4MuAtomLifetimeModel*
  LifetimeModel(G4int Z, G4int A) const;

  G4MuAtomCaptureRateModel*
  CaptureRateModel(G4int Z, G4int A, G4MuAtomCaptureRateModel*);
  G4MuAtomCaptureRateModel*
  CaptureRateModel(G4int Z, G4int A) const;

protected:
  G4MuAtom* CreateMuAtom(G4int Z, G4int A, G4int iSpin);

private:
  G4MuAtomTable();
  ~G4MuAtomTable();

  // FIXME ... there are ownership issues here ...
  typedef std::multimap<G4int, G4MuAtom*> G4MuAtomList;
  typedef std::multimap<G4int, G4MuAtom*>::iterator iterator;
  typedef std::multimap<G4int, G4MuAtom*>::const_iterator const_iterator;
  G4MuAtomList fMuAtomList;

  // FIXME ... there are ownership issues here ... 
  typedef std::map<G4int, G4MuonMinusAtomicCaptureCascadeModel*> G4CascadeModelList;
  typedef std::map<G4int, G4MuonMinusAtomicCaptureCascadeModel*>::iterator cascade_iterator;
  typedef std::map<G4int, G4MuonMinusAtomicCaptureCascadeModel*>::const_iterator 
  cascade_const_iterator;
  G4CascadeModelList fMuCascadeModelList;
  G4MuonMinusAtomicCaptureCascadeModel *defaultCascadeModel;

  // FIXME ... there are ownership issues here ...
  typedef std::map<G4int, G4MuonMinusAtomicCaptureStateModel*> G4StateModelList;
  typedef std::map<G4int, G4MuonMinusAtomicCaptureStateModel*>::iterator state_iterator;
  typedef std::map<G4int, G4MuonMinusAtomicCaptureStateModel*>::const_iterator 
  state_const_iterator;
  G4StateModelList fMuStateModelList;
  G4MuonMinusAtomicCaptureStateModel *defaultStateModel;

  // FIXME ... there are ownership issues here ...
  typedef std::map<G4int, G4MuonMinusAtomicCaptureChargeModel*> G4ChargeModelList;
  typedef std::map<G4int, G4MuonMinusAtomicCaptureChargeModel*>::iterator charge_iterator;
  typedef std::map<G4int, G4MuonMinusAtomicCaptureChargeModel*>::const_iterator 
  charge_const_iterator;
  G4ChargeModelList fMuChargeModelList;
  G4MuonMinusAtomicCaptureChargeModel *defaultChargeModel;

  // FIXME ... there are ownership issues here ...
  typedef std::map<G4int, G4MuAtomLifetimeModel*> G4LifetimeModelList;
  typedef std::map<G4int, G4MuAtomLifetimeModel*>::iterator lifetime_iterator;
  typedef std::map<G4int, G4MuAtomLifetimeModel*>::const_iterator 
  lifetime_const_iterator;
  G4LifetimeModelList fMuLifetimeModelList;
  G4MuAtomLifetimeModel *defaultLifetimeModel;

  // FIXME ... there are ownership issues here ...
  typedef std::map<G4int, G4MuAtomCaptureRateModel*> G4CaptureRateModelList;
  typedef std::map<G4int, G4MuAtomCaptureRateModel*>::iterator capture_iterator;
  typedef std::map<G4int, G4MuAtomCaptureRateModel*>::const_iterator 
  capture_const_iterator;
  G4CaptureRateModelList fMuCaptureRateModelList;
  G4MuAtomCaptureRateModel *defaultCaptureRateModel;

  // don't implement!
  G4MuAtomTable(G4MuAtomTable const&);
  G4MuAtomTable& operator=(G4MuAtomTable const&);

  friend class G4MuAtom;

  void Insert(G4MuAtom* muatom);
  void Insert(G4int Z, G4int A, G4MuonMinusAtomicCaptureCascadeModel* cascademodel);
  void Insert(G4int Z, G4int A, G4MuonMinusAtomicCaptureStateModel* statemodel);
  void Insert(G4int Z, G4int A, G4MuonMinusAtomicCaptureChargeModel* chargemodel);
  void Insert(G4int Z, G4int A, G4MuAtomLifetimeModel* lifemodel);
  void Insert(G4int Z, G4int A, G4MuAtomCaptureRateModel* capturemodel);
};

#endif
