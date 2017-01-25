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
  // Silicon
  //--------------------------------------------------

  nistMan->FindOrBuildMaterial("G4_Si");

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

  const G4int nEntries2  = 2;
  const G4int nEntries4  = 4;
  const G4int nEntries55 = 55;
  const G4int nEntries60 = 60;

  G4double PhotonEnergy2[nEntries2] = {2.00*eV, 3.77*eV};
  G4double PhotonEnergy4Fiber[nEntries4] = {2.00*eV, 2.55*eV, 2.56*eV, 3.77*eV};
  G4double PhotonEnergy4Scintillator[nEntries4] = {2.00*eV, 3.00*eV, 3.01*eV, 3.77*eV};

  G4double PhotonEnergy55[nEntries55] =
  {2.00*eV,2.03*eV,2.06*eV,2.09*eV,2.12*eV,
   2.13*eV,2.14*eV,2.15*eV,2.16*eV,2.17*eV,
   2.18*eV,2.19*eV,2.20*eV,2.21*eV,2.22*eV,
   2.23*eV,2.24*eV,2.25*eV,2.26*eV,2.27*eV,
   2.28*eV,2.29*eV,2.30*eV,2.31*eV,2.32*eV,
   2.33*eV,2.34*eV,2.35*eV,2.36*eV,2.37*eV,
   2.38*eV,2.39*eV,2.40*eV,2.41*eV,2.42*eV,
   2.43*eV,2.44*eV,2.45*eV,2.46*eV,2.47*eV,
   2.48*eV,2.49*eV,2.50*eV,2.51*eV,2.52*eV,
   2.53*eV,2.56*eV,2.59*eV,2.62*eV,2.65*eV,
   2.68*eV,2.71*eV,2.74*eV,2.77*eV,3.77*eV};

  G4double PhotonEnergy60[nEntries60] =
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

  G4double RefractiveIndex[nEntries2] = {1.00, 1.00};

  G4MaterialPropertiesTable* MPT = new G4MaterialPropertiesTable();
  MPT->AddProperty("RINDEX", PhotonEnergy2, RefractiveIndex, nEntries2);

  Air->SetMaterialPropertiesTable(MPT);

  //--------------------------------------------------
  //  PMMA for WLSfibers
  //--------------------------------------------------

  G4double RefractiveIndexWLSfiber[nEntries2] = {1.59, 1.59};

  G4double AbsWLSfiber[nEntries55] =                                       //as measured by Yuri
  {23.5*m, 8.0*m,12.0*m,22.0*m,26.3*m,26.5*m,26.2*m,25.0*m,24.0*m,23.3*m,  //these photons are lost, i.e. no wavelength shifting
   23.1*m,23.0*m,22.9*m,22.8*m,22.6*m,22.3*m,21.7*m,21.0*m,20.1*m,19.1*m,  //...
   17.9*m,16.2*m,14.3*m,12.4*m,11.1*m,10.7*m,10.6*m,10.4*m,10.1*m, 9.6*m,  //...
    8.9*m, 8.0*m, 7.0*m, 6.0*m, 5.0*m, 4.0*m, 3.2*m, 2.4*m, 1.8*m, 1.3*m,  //...
    .95*m, .75*m, .55*m, 0.4*m, 0.3*m, 0.2*m, 0.1*m, 0.1*m, 0.1*m, 0.1*m,  //these photons are not lost - they will be wavelength shifted (see below)
    0.1*m, 0.1*m, 0.1*m, 0.1*m, 0.1*m};                                    //...
  for(int i=0; i<nEntries55; i++) AbsWLSfiber[i]*=1.8;                     //after tuning it to the spectrometer scan results

  G4double AbsWLSfiberWLS[nEntries4] = {1.0e6*km, 1.0e6*km,   //no wavelength shifting in these photons - they will be lost (see above)
                                          0.4*mm,   0.4*mm};  //these photons will be wavelength shifted
                                                              //0.4 for 175 ppm K27 (from NOvA TDR)
                                                              //photon yield is proportional to surface --> short WLS absorption length

  G4double EmissionFib[nEntries55] =
  { 0.04, 0.06, 0.09, 0.14, 0.17, 0.18, 0.19, 0.20, 0.22, 0.23,
    0.23, 0.24, 0.24, 0.26, 0.30, 0.32, 0.35, 0.38, 0.41, 0.45,
    0.48, 0.56, 0.63, 0.63, 0.67, 0.63, 0.64, 0.63, 0.62, 0.61,
    0.61, 0.62, 0.64, 0.70, 0.68, 0.72, 0.80, 0.82, 0.91, 0.91,
    0.93, 0.84, 0.77, 0.67, 0.68, 0.71, 0.20, 0.98, 0.93, 0.70,
    0.28, 0.11, 0.04, 0.00, 0.00};                            //after tuning it to the spectrometer scan results

  G4MaterialPropertiesTable* MPTWLSfiber = new G4MaterialPropertiesTable();
  MPTWLSfiber->AddProperty("RINDEX",PhotonEnergy2,RefractiveIndexWLSfiber,nEntries2);
  MPTWLSfiber->AddProperty("ABSLENGTH",PhotonEnergy55,AbsWLSfiber,nEntries55);
  MPTWLSfiber->AddProperty("WLSABSLENGTH",PhotonEnergy4Fiber,AbsWLSfiberWLS,nEntries4);
  MPTWLSfiber->AddProperty("WLSCOMPONENT",PhotonEnergy55,EmissionFib,nEntries55);
  MPTWLSfiber->AddConstProperty("WLSTIMECONSTANT", 7.4*ns);

  PMMA->SetMaterialPropertiesTable(MPTWLSfiber);

  //--------------------------------------------------
  //  Polyethylene

  G4double AbsClad[nEntries2] = {40.0*m, 40.0*m};

  G4double RefractiveIndexClad1[nEntries2] = {1.49, 1.49};

  G4MaterialPropertiesTable* MPTClad1 = new G4MaterialPropertiesTable();
  MPTClad1->AddProperty("RINDEX",PhotonEnergy2,RefractiveIndexClad1,nEntries2);
  MPTClad1->AddProperty("ABSLENGTH",PhotonEnergy2,AbsClad,nEntries2);  

  Pethylene->SetMaterialPropertiesTable(MPTClad1);

  //--------------------------------------------------
  // Fluorinated Polyethylene
  //--------------------------------------------------

   G4double RefractiveIndexClad2[nEntries2] = {1.42, 1.42};

  // Add entries into properties table
  G4MaterialPropertiesTable* MPTClad2 = new G4MaterialPropertiesTable();
  MPTClad2->AddProperty("RINDEX",PhotonEnergy2,RefractiveIndexClad2,nEntries2);
  MPTClad2->AddProperty("ABSLENGTH",PhotonEnergy2,AbsClad,nEntries2); 

  FPethylene->SetMaterialPropertiesTable(MPTClad2);

  //--------------------------------------------------
  //  Polystyrene
  //--------------------------------------------------

  G4double RefractiveIndexPS[nEntries2] = {1.59, 1.59};

  G4double AbsPS[nEntries4] = {300*cm, 300*cm, 0.01*mm, 0.01*mm};

  G4double ScintilFast[nEntries60] =
  {   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
      0,    0,    0,   20,   60,  100,  170,  250,  370,  500,
    620,  750,  900, 1120, 1470, 1950, 2300, 2490, 2530, 2700,
   3200, 4000, 4600, 4600, 4250, 3950, 4100, 4750, 5050, 4600,
   4100, 3850, 3700, 3500, 3400, 3500, 3950, 4200, 3900, 3300,
   2600, 2250, 2300, 2200, 1600,  700,  300,   50,    0,    0};

  G4MaterialPropertiesTable* MPTPolystyrene = new G4MaterialPropertiesTable();
  MPTPolystyrene->AddProperty("RINDEX",PhotonEnergy2,RefractiveIndexPS,nEntries2);
  MPTPolystyrene->AddProperty("ABSLENGTH",PhotonEnergy4Scintillator,AbsPS,nEntries4);
//  MPTPolystyrene->AddConstProperty("SCINTILLATIONYIELD",15000./MeV); 
//  MPTPolystyrene->AddProperty("FASTCOMPONENT",PhotonEnergy60, ScintilFast,nEntries60); 
//  MPTPolystyrene->AddProperty("SLOWCOMPONENT",PhotonEnergy60, ScintilFast,nEntries60);  //assumed to be delayed flourescence, not used
  MPTPolystyrene->AddConstProperty("SCINTILLATIONYIELD",8000./MeV); 
  MPTPolystyrene->AddProperty("FASTCOMPONENT",PhotonEnergy60, ScintilFast,36); 
  MPTPolystyrene->AddProperty("SLOWCOMPONENT",PhotonEnergy60, ScintilFast,36);  //assumed to be delayed flourescence, not used
  MPTPolystyrene->AddConstProperty("RESOLUTIONSCALE",1.0);
  MPTPolystyrene->AddConstProperty("FASTTIMECONSTANT", 3.*ns);    //includes WLS components in the scintillator
  MPTPolystyrene->AddConstProperty("SLOWTIMECONSTANT", 100.*ns);  //unknown, not used
  MPTPolystyrene->AddConstProperty("YIELDRATIO", 1.0); 

  Polystyrene->SetMaterialPropertiesTable(MPTPolystyrene);

  Polystyrene->GetIonisation()->SetBirksConstant(0.126*mm/MeV);
}
