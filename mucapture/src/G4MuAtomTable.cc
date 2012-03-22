// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Kevin Lynch, February 12 2010
// ----------------------------------------------------------------


#include "G4MuAtomTable.hh"
#include "G4IonTable.hh"

#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4RunManager.hh"
#include "G4VUserPhysicsList.hh"
#include "G4MuAtomDecayTable.hh"
#include "G4HadronicProcessType.hh"

#include <utility>

G4MuAtomTable* G4MuAtomTable::GetInstance(){
  static G4MuAtomTable table;
  return &table;
}

G4MuAtomTable::G4MuAtomTable() :
  defaultCascadeModel(new G4MuonMinusAtomicCaptureCascadeModel()),
  defaultStateModel(new G4MuonMinusAtomicCaptureStateModel()),
  defaultChargeModel(new G4MuonMinusAtomicCaptureChargeModel()),
  defaultLifetimeModel(new G4MuAtomLifetimeModel()),
  defaultCaptureRateModel(new G4MuAtomCaptureRateModel()) {
}

G4MuAtomTable::~G4MuAtomTable(){
  delete defaultCaptureRateModel;
  delete defaultLifetimeModel;
  delete defaultChargeModel;
  delete defaultStateModel;
  delete defaultCascadeModel;
}

G4MuAtom* G4MuAtomTable::FindMuAtom(G4int Z, G4int A, G4int iSpin){
  G4MuAtom* muatom = 0;
  G4int const encoding = G4MuAtom::MakeEncoding(Z,A);

  const_iterator l = fMuAtomList.lower_bound(encoding);
  const_iterator u = fMuAtomList.upper_bound(encoding);
  //    end = fMuAtomList.end();

  while( l!=u ){
    if( l->second->GetPDGiSpin() == iSpin ){
      muatom = l->second;
      break;
    }
    ++l;
  }

  return muatom;
}

G4MuAtom* G4MuAtomTable::GetMuAtom(G4int Z, G4int A, G4int iSpin){
  G4MuAtom* muatom = FindMuAtom(Z,A,iSpin);
  if( muatom == 0 )
    muatom = CreateMuAtom(Z,A,iSpin);
  return muatom;
}

G4MuAtom* G4MuAtomTable::CreateMuAtom(G4int Z, G4int A, G4int iSpin){
  G4MuAtom* muatom = G4MuAtom::MuAtomDefinition(Z,A,iSpin);

  // There's a spurious, hardcoded warning in the G4ParticleDefinition
  // constructor about creating particles outside the PreInit state,
  // even though I appear to be doing everything right from here
  // down.  The design has a few flaws here...
  G4ParticleTable *table = G4ParticleTable::GetParticleTable();

  G4MuAtom* gmuatom = 
    static_cast<G4MuAtom*>(table->FindParticle("GenericMuAtom"));
  // FIXME ... this whole process is less clean than it should/could
  // be.  Also, one should be somewhat careful about setting BR/Rates
  // in clones.  Figure that out and write it down this time
  if( gmuatom != 0 ){
    // If G4GenericMuAtom exists, we copy the tables and process manager
    // entries from there...
    G4ProcessManager *proc = gmuatom->GetProcessManager();
    if( proc ){
      // Copy the process manager
      G4ProcessManager *newproc = new G4ProcessManager(*proc);
      newproc->SetParticleType(muatom);
      muatom->SetProcessManager(newproc);
      // "Copy" the Decay Table
      G4MuAtomDecayTable *dt = 
	gmuatom->GetMuAtomDecayTable()->Clone(muatom);
      muatom->SetMuAtomDecayTable(dt);
      // "Copy" the Capture Kinetics Table
      G4MuAtomCaptureKineticsTable* gckt = gmuatom->CaptureKineticsTable();
      G4MuAtomCaptureKineticsTable *cckt = 0;
      if( gckt ){
	cckt = gckt->Clone(muatom);
      }
      muatom->CaptureKineticsTable(cckt);

      // FIXME remove this crosscheck
      G4ProcessManager* pmanager = muatom->GetProcessManager();
      G4ProcessVector const* pVector = pmanager->GetProcessList();

      G4cout << "G4MuAtomTable::GetMuAtom " << muatom->GetParticleName() 
             << " processes: " << G4endl;
      for( G4int j=0; j<pmanager->GetProcessListLength(); ++j ) {
        G4VProcess* proc = (*pVector)[j];
        G4String const & name  = proc->GetProcessName();
        G4cout << "G4MuAtomTable::GetMuAtom " << name << G4endl;
      }

    } else {
      G4cout << "GenericMuAtom needs a process manager!\n";
    }
  } else {
    // If not, we complain
    G4cout << "GenericMuAtom needs to be instantiated for implicit MuAtom creation!\n"; 
  }

  //G4MuonMinus::Definition()
  //SetProcessSubType(fCapture)

  // G4MuonMinusCapture::GetMuonMinusAtomicCaptureProcess()

//   G4VUserPhysicsList const* cphysics = G4RunManager::GetRunManager()->GetUserPhysicsList(); 
//   G4VUserPhysicsList *physics = const_cast<G4VUserPhysicsList*>(cphysics);
//   physics->PreparePhysicsTable(muatom);
//   physics->BuildPhysicsTable(muatom);

//  G4HadronicProcess* theProcess = FindProcess(G4MuonMinus::Definition(), fCapture);
//  theProcess->PreparePhysicsTable(muatom);
//  theProcess->BuildPhysicsTable(muatom);

  // well those functions are not equivalent... 
  // the one in the G4VUserPhysicsList does the equivalent of
  
  // void G4VUserPhysicsList::PreparePhysicsTable(G4ParticleDefinition* particle){
  //    G4ProcessVector* pVector = pManager->GetProcessList();
  //  for (G4int j=0; j < pVector->size(); ++j) {
  //    (*pVector)[j]->PreparePhysicsTable(*particle);
  //  }
  // }

  G4ProcessManager* pManager = muatom->GetProcessManager();
  G4ProcessVector* pVector = pManager->GetProcessList();

  if (!pVector) {
    G4cout << "G4VUserPhysicsList::PreparePhysicsTable  "
           << ": No Process Vector for " 
           << muatom->GetParticleName() <<G4endl;
    G4Exception("G4VUserPhysicsList::PreparePhysicsTable",
                "Run0274", FatalException,
                "No process Vector");
    return muatom;
  }

  for (G4int j=0; j < pVector->size(); ++j) {
    (*pVector)[j]->PreparePhysicsTable(*muatom);
    (*pVector)[j]->BuildPhysicsTable(*muatom);
  }

  return muatom;

}


void G4MuAtomTable::Insert(G4MuAtom* muatom){
  G4int const Z = muatom->GetAtomicNumber();
  G4int const A = muatom->GetAtomicMass();  
  G4int const encoding = G4MuAtom::MakeEncoding(Z,A);

  fMuAtomList.insert( std::make_pair(encoding, muatom) );
}

void G4MuAtomTable::Insert(G4int Z, G4int A,
			   G4MuonMinusAtomicCaptureCascadeModel* cascademodel){
  G4int const encoding = G4MuAtom::MakeEncoding(Z,A);
  fMuCascadeModelList.insert( std::make_pair(encoding, cascademodel) );
}

void G4MuAtomTable::Insert(G4int Z, G4int A,
			   G4MuonMinusAtomicCaptureStateModel* statemodel){
  G4int const encoding = G4MuAtom::MakeEncoding(Z,A);
  fMuStateModelList.insert( std::make_pair(encoding, statemodel) );
}

void G4MuAtomTable::Insert(G4int Z, G4int A,
			   G4MuonMinusAtomicCaptureChargeModel* chargemodel){
  G4int const encoding = G4MuAtom::MakeEncoding(Z,A);
  fMuChargeModelList.insert( std::make_pair(encoding, chargemodel) );
}

void G4MuAtomTable::Insert(G4int Z, G4int A,
			   G4MuAtomLifetimeModel* lifemodel){
  G4int const encoding = G4MuAtom::MakeEncoding(Z,A);
  fMuLifetimeModelList.insert( std::make_pair(encoding, lifemodel) );
}

void G4MuAtomTable::Insert(G4int Z, G4int A,
			   G4MuAtomCaptureRateModel* lifemodel){
  G4int const encoding = G4MuAtom::MakeEncoding(Z,A);
  fMuCaptureRateModelList.insert( std::make_pair(encoding, lifemodel) );
}

// FIXME ... I coudl really clean these implementations up with a
// template function...
G4MuonMinusAtomicCaptureCascadeModel*
G4MuAtomTable::CascadeModel(G4int Z, G4int A, 
			  G4MuonMinusAtomicCaptureCascadeModel* model){
  G4int const encoding = G4MuAtom::MakeEncoding(Z,A);
  cascade_iterator i = fMuCascadeModelList.find(encoding);
  G4MuonMinusAtomicCaptureCascadeModel *old = defaultCascadeModel;
  if( i != fMuCascadeModelList.end() ){
    old = i->second;
  } else {
    Insert(Z, A, model);
  }
  return old;
}

G4MuonMinusAtomicCaptureCascadeModel*
G4MuAtomTable::CascadeModel(G4int Z, G4int A) const {
  G4int const encoding = G4MuAtom::MakeEncoding(Z,A);
  cascade_const_iterator i = fMuCascadeModelList.find(encoding);
  G4MuonMinusAtomicCaptureCascadeModel *old = defaultCascadeModel;
  if( i != fMuCascadeModelList.end() )
    old = i->second;
  return old;
}

G4MuonMinusAtomicCaptureStateModel*
G4MuAtomTable::StateModel(G4int Z, G4int A, 
			  G4MuonMinusAtomicCaptureStateModel* model){
  G4int const encoding = G4MuAtom::MakeEncoding(Z,A);
  state_iterator i = fMuStateModelList.find(encoding);
  G4MuonMinusAtomicCaptureStateModel *old = defaultStateModel;
  if( i != fMuStateModelList.end() ){
    old = i->second;
  } else {
    Insert(Z, A, model);
  }
  return old;
}

G4MuonMinusAtomicCaptureStateModel*
G4MuAtomTable::StateModel(G4int Z, G4int A) const {
  G4int const encoding = G4MuAtom::MakeEncoding(Z,A);
  state_const_iterator i = fMuStateModelList.find(encoding);
  G4MuonMinusAtomicCaptureStateModel *old = defaultStateModel;
  if( i != fMuStateModelList.end() )
    old = i->second;
  return old;
}

G4MuonMinusAtomicCaptureChargeModel*
G4MuAtomTable::ChargeModel(G4int Z, G4int A, 
			  G4MuonMinusAtomicCaptureChargeModel* model){
  G4int const encoding = G4MuAtom::MakeEncoding(Z,A);
  charge_iterator i = fMuChargeModelList.find(encoding);
  G4MuonMinusAtomicCaptureChargeModel *old = defaultChargeModel;
  if( i != fMuChargeModelList.end() ){
    old = i->second;
  } else {
    Insert(Z, A, model);
  }
  return old;
}

G4MuonMinusAtomicCaptureChargeModel*
G4MuAtomTable::ChargeModel(G4int Z, G4int A) const {
  G4int const encoding = G4MuAtom::MakeEncoding(Z,A);
  charge_const_iterator i = fMuChargeModelList.find(encoding);
  G4MuonMinusAtomicCaptureChargeModel *old = defaultChargeModel;
  if( i != fMuChargeModelList.end() )
    old = i->second;
  return old;
}

G4MuAtomLifetimeModel*
G4MuAtomTable::LifetimeModel(G4int Z, G4int A, 
			  G4MuAtomLifetimeModel* model){
  G4int const encoding = G4MuAtom::MakeEncoding(Z,A);
  lifetime_iterator i = fMuLifetimeModelList.find(encoding);
  G4MuAtomLifetimeModel *old = defaultLifetimeModel;
  if( i != fMuLifetimeModelList.end() ){
    old = i->second;
  } else {
    Insert(Z, A, model);
  }
  return old;
}

G4MuAtomLifetimeModel*
G4MuAtomTable::LifetimeModel(G4int Z, G4int A) const {
  G4int const encoding = G4MuAtom::MakeEncoding(Z,A);
  lifetime_const_iterator i = fMuLifetimeModelList.find(encoding);
  G4MuAtomLifetimeModel *old = defaultLifetimeModel;
  if( i != fMuLifetimeModelList.end() )
    old = i->second;
  return old;
}


G4MuAtomCaptureRateModel*
G4MuAtomTable::CaptureRateModel(G4int Z, G4int A,
			  G4MuAtomCaptureRateModel* model){
  G4int const encoding = G4MuAtom::MakeEncoding(Z,A);
  capture_iterator i = fMuCaptureRateModelList.find(encoding);
  G4MuAtomCaptureRateModel *old = defaultCaptureRateModel;
  if( i != fMuCaptureRateModelList.end() ){
    old = i->second;
  } else {
    Insert(Z, A, model);
  }
  return old;
}

G4MuAtomCaptureRateModel*
G4MuAtomTable::CaptureRateModel(G4int Z, G4int A) const {
  G4int const encoding = G4MuAtom::MakeEncoding(Z,A);
  capture_const_iterator i = fMuCaptureRateModelList.find(encoding);
  G4MuAtomCaptureRateModel *old = defaultCaptureRateModel;
  if( i != fMuCaptureRateModelList.end() )
    old = i->second;
  return old;
}

