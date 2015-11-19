#include "WLSMaterials.hh"

#include "G4SystemOfUnits.hh"

WLSMaterials::WLSMaterials()
{
  nistMan = G4NistManager::Instance();

  nistMan->SetVerbose(2);

  CreateMaterials();
}

WLSMaterials::~WLSMaterials()
{
  delete    PMMA;
  delete    Pethylene;
  delete    FPethylene;
  delete    Polystyrene;
}

WLSMaterials* WLSMaterials::instance = NULL;

WLSMaterials* WLSMaterials::GetInstance()
{
  if(instance == 0) instance = new WLSMaterials();
  return instance;
}

G4Material* WLSMaterials::GetMaterial(const G4String material)
{
  G4Material* mat =  nistMan->FindOrBuildMaterial(material);

  if (!mat) mat = G4Material::GetMaterial(material);
  if (!mat) 
  {
     std::ostringstream o;
     o << "Material " << material << " not found!";
     G4Exception("WLSMaterials::GetMaterial","", FatalException,o.str().c_str());
  }

  return mat;
}

void WLSMaterials::CreateMaterials()
{
  G4double density;
  G4int ncomponents;
  G4double fractionmass;
  std::vector<G4int> natoms;
  std::vector<G4double> fractionMass;
  std::vector<G4String> elements;

  // Materials Definitions
  // =====================

  //--------------------------------------------------
  // Vacuum
  //--------------------------------------------------

  nistMan->FindOrBuildMaterial("G4_Galactic");

  //--------------------------------------------------
  // Air
  //--------------------------------------------------

  Air = nistMan->FindOrBuildMaterial("G4_AIR");

  //--------------------------------------------------
  // PVC
  //--------------------------------------------------

  PVC = nistMan->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");

  //--------------------------------------------------
  // WLSfiber PMMA
  //--------------------------------------------------

  elements.push_back("C");     natoms.push_back(5);
  elements.push_back("H");     natoms.push_back(8);
  elements.push_back("O");     natoms.push_back(2);

  density = 1.190*g/cm3;

  PMMA = nistMan->
          ConstructNewMaterial("PMMA", elements, natoms, density);

  elements.clear();
  natoms.clear();

  //--------------------------------------------------
  // Cladding (polyethylene)
  //--------------------------------------------------

  elements.push_back("C");     natoms.push_back(2);
  elements.push_back("H");     natoms.push_back(4);

  density = 1.200*g/cm3;

  Pethylene = nistMan->
          ConstructNewMaterial("Pethylene", elements, natoms, density);

  elements.clear();
  natoms.clear();

  //--------------------------------------------------
  // Double Cladding (fluorinated polyethylene)
  //--------------------------------------------------

  elements.push_back("C");     natoms.push_back(2);
  elements.push_back("H");     natoms.push_back(4);

  density = 1.400*g/cm3;

  FPethylene = nistMan->
          ConstructNewMaterial("FPethylene", elements, natoms, density);

  elements.clear();
  natoms.clear();

  //--------------------------------------------------
  // Polystyrene
  //--------------------------------------------------
 
  elements.push_back("C");     natoms.push_back(8);
  elements.push_back("H");     natoms.push_back(8);

  density = 1.050*g/cm3;

  Polystyrene = nistMan->
          ConstructNewMaterial("Polystyrene", elements, natoms, density);

  elements.clear();
  natoms.clear();

  //--------------------------------------------------
  // Aluminium
  //--------------------------------------------------

  nistMan->FindOrBuildMaterial("G4_Al");

  //--------------------------------------------------
  // TiO2
  //--------------------------------------------------

  elements.push_back("Ti");     natoms.push_back(1);
  elements.push_back("O");      natoms.push_back(2);

  density     = 4.26*g/cm3;

  G4Material* TiO2 = nistMan->
          ConstructNewMaterial("TiO2", elements, natoms, density);

  elements.clear();
  natoms.clear();

  //--------------------------------------------------
  // Scintillator Coating - 15% TiO2 and 85% polystyrene by weight.
  //--------------------------------------------------

  density = 1.52*g/cm3;

  Coating =
          new G4Material("Coating", density, ncomponents=2);

  Coating->AddMaterial(TiO2,        fractionmass = 15*perCent);
  Coating->AddMaterial(Polystyrene, fractionmass = 85*perCent);

  //
  // ------------ Generate & Add Material Properties Table ------------
  //

  const G4int nEntries = 60;

  G4double PhotonEnergy[nEntries] =
  {2.00*eV,2.03*eV,2.06*eV,2.09*eV,2.12*eV,
   2.15*eV,2.18*eV,2.21*eV,2.24*eV,2.27*eV,
   2.30*eV,2.33*eV,2.36*eV,2.39*eV,2.42*eV,
   2.45*eV,2.48*eV,2.51*eV,2.54*eV,2.57*eV,
   2.60*eV,2.63*eV,2.66*eV,2.69*eV,2.72*eV,
   2.75*eV,2.78*eV,2.81*eV,2.84*eV,2.87*eV,
   2.90*eV,2.93*eV,2.96*eV,2.99*eV,3.02*eV,
   3.05*eV,3.08*eV,3.11*eV,3.14*eV,3.17*eV,
   3.20*eV,3.23*eV,3.26*eV,3.29*eV,3.32*eV,
   3.35*eV,3.38*eV,3.41*eV,3.44*eV,3.47*eV,
   3.50*eV,3.53*eV,3.56*eV,3.59*eV,3.62*eV,
   3.65*eV,3.68*eV,3.71*eV,3.74*eV,3.77*eV};

  //--------------------------------------------------
  // Air
  //--------------------------------------------------

  G4double RefractiveIndex[nEntries] =
  { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00};

  G4MaterialPropertiesTable* MPT = new G4MaterialPropertiesTable();
  MPT->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex, nEntries);

  Air->SetMaterialPropertiesTable(MPT);

  //--------------------------------------------------
  //  PMMA for WLSfibers
  //--------------------------------------------------

  G4double RefractiveIndexWLSfiber[nEntries] =
  { 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59,
    1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59,
    1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59,
    1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59,
    1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59,
    1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59};


  G4double AbsWLSfiber[nEntries] =
  {24.0*m, 7.0*m, 7.0*m,20.0*m,25.0*m,25.0*m,22.0*m,21.5*m,21.0*m,19.0*m,  //these photons are lost, i.e. no wavelength shifting
   16.0*m,12.0*m,10.5*m, 9.5*m, 6.5*m, 4.5*m, 3.0*m, 0.5*m, 0.1*m, 1.*km,  //...
    1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km,  //these photons are not lost - they will be wavelength shifted (see below)
    1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km,
    1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km,
    1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km};
  for(int i=0; i<nEntries; i++) AbsWLSfiber[i]*=1.4;    //override the original values

  G4double AbsWLSfiberWLS[nEntries] =
  { 1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km,  //no wavelength shifting in these photons - they will be lost (see above)
    1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*km, 1.*mm,  //...
    1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm,  //these photons will be wavelength shifted
    1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm,  //...
    1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm,
    1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm};

  G4double EmissionFib[nEntries] =
  { 0.00, 0.50, 0.75, 0.75, 0.75, 1.00, 1.50, 2.00, 3.00, 4.00,
    5.00, 6.00, 7.50,10.50,13.00,15.00,14.00,13.00,13.50,16.50,
   18.50,16.00,10.50, 3.50, 1.50, 0.50, 0.00, 0.00, 0.00, 0.00, 
    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};

  G4MaterialPropertiesTable* MPTWLSfiber = new G4MaterialPropertiesTable();
  MPTWLSfiber->AddProperty("RINDEX",PhotonEnergy,RefractiveIndexWLSfiber,nEntries);
  MPTWLSfiber->AddProperty("ABSLENGTH",PhotonEnergy,AbsWLSfiber,nEntries);
  MPTWLSfiber->AddProperty("WLSABSLENGTH",PhotonEnergy,AbsWLSfiberWLS,nEntries);
  MPTWLSfiber->AddProperty("WLSCOMPONENT",PhotonEnergy,EmissionFib,nEntries);
  MPTWLSfiber->AddConstProperty("WLSTIMECONSTANT", 7.4*ns);

  PMMA->SetMaterialPropertiesTable(MPTWLSfiber);

  //--------------------------------------------------
  //  Polyethylene
  //--------------------------------------------------

  G4double RefractiveIndexClad1[nEntries] =
  { 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49};

  G4double AbsClad[nEntries] =
  {20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,
   20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,
   20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,
   20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,
   20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,
   20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m,20.0*m};

  G4MaterialPropertiesTable* MPTClad1 = new G4MaterialPropertiesTable();
  MPTClad1->AddProperty("RINDEX",PhotonEnergy,RefractiveIndexClad1,nEntries);
  MPTClad1->AddProperty("ABSLENGTH",PhotonEnergy,AbsClad,nEntries);

  Pethylene->SetMaterialPropertiesTable(MPTClad1);

  //--------------------------------------------------
  // Fluorinated Polyethylene
  //--------------------------------------------------

   G4double RefractiveIndexClad2[nEntries] =
   { 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
     1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
     1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
     1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
     1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
     1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42};

  // Add entries into properties table
  G4MaterialPropertiesTable* MPTClad2 = new G4MaterialPropertiesTable();
  MPTClad2->AddProperty("RINDEX",PhotonEnergy,RefractiveIndexClad2,nEntries);
  MPTClad2->AddProperty("ABSLENGTH",PhotonEnergy,AbsClad,nEntries);

  FPethylene->SetMaterialPropertiesTable(MPTClad2);

  //--------------------------------------------------
  //  Polystyrene
  //--------------------------------------------------

  G4double RefractiveIndexPS[nEntries] =
  { 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59,
    1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59,
    1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59,
    1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59,
    1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59,
    1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59, 1.59};

  G4double AbsPS[nEntries] =
  {50.*cm,50.*cm,50.*cm,50.*cm,50.*cm,50.*cm,50.*cm,50.*cm,50.*cm,50.*cm,
   50.*cm,40.*cm,33.*cm,33.*cm,33.*cm,33.*cm,33.*cm,33.*cm,33.*cm,28.*cm,
   28.*cm,28.*cm,25.*cm,25.*cm,22.*cm,22.*cm,20.*cm,20.*cm,18.*cm,16.*cm,
   15.*cm,12.*cm,8.6*cm,5.4*cm,2.8*cm,1.3*cm,0.3*cm,0.2*cm,0.2*cm,0.2*cm,
   0.2*cm,0.2*cm,0.2*cm,0.2*cm,0.2*cm,0.2*cm,0.2*cm,0.2*cm,0.2*cm,0.2*cm,
   0.2*cm,0.2*cm,0.2*cm,0.2*cm,0.2*cm,0.2*cm,0.2*cm,0.2*cm,0.2*cm,0.2*cm};
  for(int i=36; i<nEntries; i++) AbsPS[i]=0.;  //(the original values lead to probabilities which do not agree with the test beam data)

  G4double ScintilFast[nEntries] =
  {   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
      0,    0,    0,   20,   60,  100,  170,  250,  370,  500,
    620,  750,  900, 1120, 1470, 1950, 2300, 2490, 2530, 2700,
   3200, 4000, 4600, 4600, 4250, 3950, 4100, 4750, 5050, 4600,
   4100, 3850, 3700, 3500, 3400, 3500, 3950, 4200, 3900, 3300,
   2600, 2250, 2300, 2200, 1600,  700,  300,   50,    0,    0};

  G4MaterialPropertiesTable* MPTPolystyrene = new G4MaterialPropertiesTable();
  MPTPolystyrene->AddProperty("RINDEX",PhotonEnergy,RefractiveIndexPS,nEntries);
  MPTPolystyrene->AddProperty("ABSLENGTH",PhotonEnergy,AbsPS,nEntries);
//no need to to simulate the entire emission spectrum, if all photons with i>=36 are absorbed
  MPTPolystyrene->AddProperty("FASTCOMPONENT",PhotonEnergy, ScintilFast,36);
  MPTPolystyrene->AddProperty("SLOWCOMPONENT",PhotonEnergy, ScintilFast,36);  //assumed to be delayed flourescence with same spectrum, not used
  MPTPolystyrene->AddConstProperty("SCINTILLATIONYIELD",4000./MeV);
  MPTPolystyrene->AddConstProperty("RESOLUTIONSCALE",1.0);
  MPTPolystyrene->AddConstProperty("FASTTIMECONSTANT", 3.*ns);    //includes WLS components in the scintillator
  MPTPolystyrene->AddConstProperty("SLOWTIMECONSTANT", 100.*ns);  //unknown, not used
  MPTPolystyrene->AddConstProperty("YIELDRATIO", 1.0); 

  Polystyrene->SetMaterialPropertiesTable(MPTPolystyrene);

  Polystyrene->GetIonisation()->SetBirksConstant(0.126*mm/MeV);
}
