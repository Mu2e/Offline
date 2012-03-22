// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Kevin Lynch, February 12 2010
// ----------------------------------------------------------------

#include "G4MuAtom.hh"
#include "G4MuAtomTable.hh"

#include "G4NucleiProperties.hh"
#include "G4IonTable.hh"

#include <sstream>
#include <string>

G4MuAtom::~G4MuAtom(){
}

G4MuAtom* G4MuAtom::MuAtom(G4int Z, G4int A, G4int iSpin){
  return MuAtomDefinition(Z,A,iSpin);
}

G4MuAtom* G4MuAtom::Definition(G4int Z, G4int A, G4int iSpin){
  return MuAtomDefinition(Z,A,iSpin);
}

G4MuAtom* G4MuAtom::MuAtomDefinition(G4int Z, G4int A, G4int iSpin){
  // FIXME Scrub the Z, A, iSpin values...

  G4MuAtomTable* table = G4MuAtomTable::GetInstance();
  G4MuAtom* muatom = table->FindMuAtom(Z,A,iSpin);
  if( muatom != 0 ){
    //    G4Exception("Attempt to instantiate G4MuAtom already in G4MuAtomTable!");
    return muatom;
  }

  G4String const name = G4MuAtom::MakeName(Z,A,iSpin);

  // FIXME ... ought to be able to get this from somewhere!
  G4double const muon_mass_c2 = 0.1056584*GeV;
  // FIXME ... subtract the binding energy!
  G4double const mass = G4NucleiProperties::GetNuclearMass(A,Z) + muon_mass_c2;

  G4int const encoding = MakeEncoding(Z,A);

  muatom = new G4MuAtom(name, mass, Z, A, iSpin, encoding); //FIXME
							    //leak
							    //... table
							    //takes
							    //ownership? 

  table->Insert(muatom);

  return muatom;
}

G4MuAtom::G4MuAtom(const G4String&  aName,  
		   G4double         mass,     
		   G4int            Z,
		   G4int            A,
		   G4int            iSpin,
		   G4int            encoding)
  //               name             mass          width         charge
  //             2*spin           parity  C-conjugation
  //          2*Isospin       2*Isospin3       G-parity
  //               type    lepton number  baryon number   PDG encoding
  //             stable         lifetime    decay table
  //             shortlived      subType    anti_encoding
  : G4ParticleDefinition( aName,mass,0,(Z-1)*eplus,iSpin,0,0,0,0,0,
			  "MuAtom",-1,A,encoding,false,0.,0,
			  false,"",0),
    fMuAtomDecayTable(0), fCaptureTable(0)
{  
  SetAtomicNumber(Z);
  SetAtomicMass(A);
  G4MuAtomTable *table = G4MuAtomTable::GetInstance();
  SetPDGLifeTime(table->LifetimeModel(Z,A)->GetLifetime(Z,A));
}


G4MuAtomCaptureKineticsTable* 
G4MuAtom::CaptureKineticsTable(G4MuAtomCaptureKineticsTable* ckt){
  G4MuAtomCaptureKineticsTable* old = fCaptureTable;
  fCaptureTable = ckt;
  return old;
}


G4double G4MuAtom::GetCaptureRate() const {
  if( !fCaptureTable )
    return 0;
  return fCaptureTable->GetTotalRate();
}


// Cascade Model
G4MuonMinusAtomicCaptureCascadeModel* 
G4MuAtom::CascadeModel(G4int Z, G4int A, G4MuonMinusAtomicCaptureCascadeModel* m){
  G4MuAtomTable *table = G4MuAtomTable::GetInstance();
  return table->CascadeModel(Z,A,m);
}

G4MuonMinusAtomicCaptureCascadeModel* 
G4MuAtom::CascadeModel(G4int Z, G4int A){
  G4MuAtomTable *table = G4MuAtomTable::GetInstance();
  return table->CascadeModel(Z,A);
}

G4MuonMinusAtomicCaptureCascadeModel* 
G4MuAtom::CascadeModel(G4MuonMinusAtomicCaptureCascadeModel* m){
  G4MuAtomTable *table = G4MuAtomTable::GetInstance();
  return table->CascadeModel(GetAtomicNumber(),GetAtomicMass(),m);
}

G4MuonMinusAtomicCaptureCascadeModel* G4MuAtom::CascadeModel() const{
  G4MuAtomTable *table = G4MuAtomTable::GetInstance();
  return table->CascadeModel(GetAtomicNumber(),GetAtomicMass());
}
  

// State Model
G4MuonMinusAtomicCaptureStateModel* 
G4MuAtom::StateModel(G4int Z, G4int A, G4MuonMinusAtomicCaptureStateModel* m){
  G4MuAtomTable *table = G4MuAtomTable::GetInstance();
  return table->StateModel(Z,A,m);
}

G4MuonMinusAtomicCaptureStateModel* 
G4MuAtom::StateModel(G4int Z, G4int A){
  G4MuAtomTable *table = G4MuAtomTable::GetInstance();
  return table->StateModel(Z,A);
}

G4MuonMinusAtomicCaptureStateModel* 
G4MuAtom::StateModel(G4MuonMinusAtomicCaptureStateModel* m){
  G4MuAtomTable *table = G4MuAtomTable::GetInstance();
  return table->StateModel(GetAtomicNumber(),GetAtomicMass(),m);
}

G4MuonMinusAtomicCaptureStateModel* G4MuAtom::StateModel() const{
  G4MuAtomTable *table = G4MuAtomTable::GetInstance();
  return table->StateModel(GetAtomicNumber(),GetAtomicMass());
}
  
// Charge Model
G4MuonMinusAtomicCaptureChargeModel* 
G4MuAtom::ChargeModel(G4int Z, G4int A, G4MuonMinusAtomicCaptureChargeModel* m){
  G4MuAtomTable *table = G4MuAtomTable::GetInstance();
  return table->ChargeModel(Z,A,m);
}

G4MuonMinusAtomicCaptureChargeModel* 
G4MuAtom::ChargeModel(G4int Z, G4int A){
  G4MuAtomTable *table = G4MuAtomTable::GetInstance();
  return table->ChargeModel(Z,A);
}

G4MuonMinusAtomicCaptureChargeModel* 
G4MuAtom::ChargeModel(G4MuonMinusAtomicCaptureChargeModel* m){
  G4MuAtomTable *table = G4MuAtomTable::GetInstance();
  return table->ChargeModel(GetAtomicNumber(),GetAtomicMass(),m);
}

G4MuonMinusAtomicCaptureChargeModel* G4MuAtom::ChargeModel() const{
  G4MuAtomTable *table = G4MuAtomTable::GetInstance();
  return table->ChargeModel(GetAtomicNumber(),GetAtomicMass());
}
  
// Lifetime Model
G4MuAtomLifetimeModel* 
G4MuAtom::LifetimeModel(G4int Z, G4int A, G4MuAtomLifetimeModel* m){
  G4MuAtomTable *table = G4MuAtomTable::GetInstance();
  if( G4MuAtom * ma = table->FindMuAtom(Z,A) )
    ma->SetPDGLifeTime( m->GetLifetime(Z,A) );
  return table->LifetimeModel(Z,A,m);
}

G4MuAtomLifetimeModel* G4MuAtom::LifetimeModel(G4int Z, G4int A){
  G4MuAtomTable *table = G4MuAtomTable::GetInstance();
  return table->LifetimeModel(Z,A);
}

G4MuAtomLifetimeModel* G4MuAtom::LifetimeModel(G4MuAtomLifetimeModel* m){
  G4MuAtomTable *table = G4MuAtomTable::GetInstance();
  SetPDGLifeTime( m->GetLifetime( GetAtomicNumber(), GetAtomicMass() ) );
  return table->LifetimeModel( GetAtomicNumber(), GetAtomicMass(), m);
}

G4MuAtomLifetimeModel* G4MuAtom::LifetimeModel() const{
  G4MuAtomTable *table = G4MuAtomTable::GetInstance();
  return table->LifetimeModel(GetAtomicNumber(), GetAtomicMass());
}


// From G4MuMinusCaptureCascade
G4double G4MuAtom::GetKShellEnergy() const { 
  // Calculate the Energy of K Mesoatom Level for this Element using
  // the Energy of Hydrogen Atom taken into account finite size of the
  // nucleus (V.Ivanchenko)
  const G4int ListK = 28;
  static G4double ListZK[ListK] = {
      1., 2.,  4.,  6.,  8., 11., 14., 17., 18., 21., 24.,
     26., 29., 32., 38., 40., 41., 44., 49., 53., 55.,
     60., 65., 70., 75., 81., 85., 92.};
  static G4double ListKEnergy[ListK] = {
     0.00275, 0.011, 0.043, 0.098, 0.173, 0.326,
     0.524, 0.765, 0.853, 1.146, 1.472,
     1.708, 2.081, 2.475, 3.323, 3.627, 
     3.779, 4.237, 5.016, 5.647, 5.966,
     6.793, 7.602, 8.421, 9.249, 10.222,
    10.923,11.984};

  // Energy with finit size corrections
  G4double KEnergy = GetLinApprox(ListK,ListZK,ListKEnergy,GetAtomicNumber());

  return KEnergy;
}


G4String G4MuAtom::MakeName(G4int Z, G4int A, G4int iSpin){
  // mu + ion table name - "[0.0]"
  G4String name;
  std::ostringstream o;
  if( Z==1 ){
    if( A == 1 )
      name = "mu_P";
    else if ( A == 2 )
      name = "mu_D";
    else if ( A == 3 )
      name = "mu_T";
  } else {
    o << "mu_";
    o << G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonName(Z,A,0.0);
    name = o.str();
    G4int pos = name.find('[');
    // string [nnnn]
    name = name.substr(0,pos);
  }
  if( iSpin!=0 ){
    o.str("");
    o << name;
    if( iSpin%2 == 0 )
      o << "_" << iSpin/2;
    else
      o<< "_" << iSpin << "/2";
    name = o.str();
  }
  
  //  G4cout << "G4MuAtom::MakeName" << " making " << name << G4endl;
  return name;
}

G4int G4MuAtom::MakeEncoding(G4int Z, G4int A){
  G4int encoding = 
    G4ParticleTable::GetParticleTable()->GetIonTable()->GetNucleusEncoding(Z, A, 0);
  encoding += 1000000000;
  //  G4cout  << "G4MuAtom::MakeEncoding" << " Z " << Z << " A " << A << " encoding " << encoding << G4endl;

  return encoding;
}


// From G4MuMinusCaptureCascade
G4double G4MuAtom::GetLinApprox(G4int N, 
				const G4double* X, 
				const G4double* Y, 
				G4double Xuser) const {
  G4double Yuser;
  if(Xuser <= X[0])        Yuser = Y[0];
  else if(Xuser >= X[N-1]) Yuser = Y[N-1];
  else {
    G4int i;
    for (i=1; i<N; i++){
      if(Xuser <= X[i]) break; 
    }    

    if(Xuser == X[i]) Yuser = Y[i];
    else Yuser = Y[i-1] + (Y[i] - Y[i-1])*(Xuser - X[i-1])/(X[i] - X[i-1]);
  }
  return Yuser;
}

