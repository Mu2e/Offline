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
  // Scintillator Coating - 30% TiO2 and 70% polystyrene by weight.
  //--------------------------------------------------

  density = 2.01*g/cm3;

  Coating =
          new G4Material("Coating", density, ncomponents=2);

  Coating->AddMaterial(TiO2,        fractionmass = 30*perCent);
  Coating->AddMaterial(Polystyrene, fractionmass = 70*perCent);

  //
  // ------------ Generate & Add Material Properties Table ------------
  //

  const G4int nEntries2  = 2;
  const G4int nEntries279= 279;
  const G4int nEntries65 = 65;
  const G4int nEntries43 = 43;
  const G4int nEntries17 = 17;
  const G4int nEntries91 = 91;

  G4double PhotonEnergy2[nEntries2] = {2.00*eV, 3.80*eV};

  G4double PhotonEnergy279[nEntries279] =
  {
 2.000,  2.002,  2.003,  2.005,  2.006,  2.008,  2.009,  2.011,  2.013,  2.014,   //as measured with fiber tester by Yuri
 2.016,  2.017,  2.019,  2.020,  2.022,  2.024,  2.025,  2.027,  2.028,  2.030,   //...
 2.031,  2.033,  2.035,  2.036,  2.038,  2.039,  2.041,  2.043,  2.044,  2.046, 
 2.047,  2.049,  2.051,  2.052,  2.054,  2.056,  2.057,  2.059,  2.060,  2.062, 
 2.064,  2.065,  2.067,  2.069,  2.070,  2.072,  2.074,  2.075,  2.077,  2.079, 
 2.080,  2.082,  2.084,  2.085,  2.087,  2.089,  2.090,  2.092,  2.094,  2.095, 
 2.097,  2.099,  2.100,  2.102,  2.104,  2.105,  2.107,  2.109,  2.110,  2.112, 
 2.114,  2.116,  2.117,  2.119,  2.121,  2.123,  2.124,  2.126,  2.128,  2.129, 
 2.131,  2.133,  2.135,  2.136,  2.138,  2.140,  2.142,  2.143,  2.145,  2.147, 
 2.149,  2.150,  2.152,  2.154,  2.156,  2.158,  2.159,  2.161,  2.163,  2.165, 
 2.167,  2.168,  2.170,  2.172,  2.174,  2.176,  2.177,  2.179,  2.181,  2.183, 
 2.185,  2.186,  2.188,  2.190,  2.192,  2.194,  2.196,  2.197,  2.199,  2.201, 
 2.203,  2.205,  2.207,  2.209,  2.210,  2.212,  2.214,  2.216,  2.218,  2.220, 
 2.222,  2.224,  2.225,  2.227,  2.229,  2.231,  2.233,  2.235,  2.237,  2.239, 
 2.241,  2.243,  2.244,  2.246,  2.248,  2.250,  2.252,  2.254,  2.256,  2.258, 
 2.260,  2.262,  2.264,  2.266,  2.268,  2.270,  2.272,  2.274,  2.276,  2.278, 
 2.279,  2.281,  2.283,  2.285,  2.287,  2.289,  2.291,  2.293,  2.295,  2.297, 
 2.299,  2.301,  2.303,  2.305,  2.307,  2.310,  2.312,  2.314,  2.316,  2.318, 
 2.320,  2.322,  2.324,  2.326,  2.328,  2.330,  2.332,  2.334,  2.336,  2.338, 
 2.340,  2.342,  2.344,  2.347,  2.349,  2.351,  2.353,  2.355,  2.357,  2.359, 
 2.361,  2.363,  2.365,  2.368,  2.370,  2.372,  2.374,  2.376,  2.378,  2.380, 
 2.383,  2.385,  2.387,  2.389,  2.391,  2.393,  2.396,  2.398,  2.400,  2.402, 
 2.404,  2.406,  2.409,  2.411,  2.413,  2.415,  2.417,  2.420,  2.422,  2.424, 
 2.426,  2.429,  2.431,  2.433,  2.435,  2.437,  2.440,  2.442,  2.444,  2.446, 
 2.449,  2.451,  2.453,  2.456,  2.458,  2.460,  2.462,  2.465,  2.467,  2.469, 
 2.472,  2.474,  2.476,  2.479,  2.481,  2.483,  2.486,  2.488,  2.490,  2.493, 
 2.495,  2.497,  2.500,  2.502,  2.504,  2.507,  2.509,  2.511,  2.514,  2.516, 
 2.519,  2.521,  2.523,  2.526,  2.528,  2.531,                                    //...
 2.55999,   2.56,  3.800};                                                       //additional values
  for(int i=0; i<nEntries279; i++) PhotonEnergy279[i]*=eV;

  G4double PhotonEnergy65[nEntries65] =
  {
  2.00, 2.55999,                                               //additional values
  2.56, 2.58, 2.60, 2.62, 2.64, 2.66, 2.68, 2.70, 2.72, 2.74,  //as given by Kuraray
  2.76, 2.78, 2.80, 2.82, 2.84, 2.86, 2.88, 2.90, 2.92, 2.94,  //...
  2.96, 2.98, 3.00, 3.02, 3.04, 3.06, 3.08, 3.10, 3.12, 3.14, 
  3.16, 3.18, 3.20, 3.22, 3.24, 3.26, 3.28, 3.30, 3.32, 3.34, 
  3.36, 3.38, 3.40, 3.42, 3.44, 3.46, 3.48, 3.50, 3.52, 3.54, 
  3.56, 3.58, 3.60, 3.62, 3.64, 3.66, 3.68, 3.70, 3.72, 3.74, 
  3.76, 3.78, 3.80};
  for(int i=0; i<nEntries65; i++) PhotonEnergy65[i]*=eV;

  G4double PhotonEnergy43[nEntries43] =
  {
  2.00, 2.02, 2.04, 2.06, 2.08, 2.10, 2.12, 2.14, 2.16, 2.18,  //as given by Kuraray
  2.20, 2.22, 2.24, 2.26, 2.28, 2.30, 2.32, 2.34, 2.36, 2.38,  //...
  2.40, 2.42, 2.44, 2.46, 2.48, 2.50, 2.52, 2.54, 2.56, 2.58, 
  2.60, 2.62, 2.64, 2.66, 2.68, 2.70, 2.72, 2.74, 2.76, 2.78, 
  2.80, 
  2.82, 3.80};                                                 //additional values
  for(int i=0; i<nEntries43; i++) PhotonEnergy43[i]*=eV;


  G4double PhotonEnergy17[nEntries17] =
  {
  2.00, 2.25, 2.40, 2.60, 2.70, 2.80, 2.90, 2.95, 3.00, 3.02,
  3.04, 3.06, 3.08, 3.10, 3.12, 3.40, 3.80};
  for(int i=0; i<nEntries17; i++) PhotonEnergy17[i]*=eV;

  G4double PhotonEnergy91[nEntries91] =
  {
  2.00, 2.02, 2.04, 2.06, 2.08, 2.10, 2.12, 2.14, 2.16, 2.18,
  2.20, 2.22, 2.24, 2.26, 2.28, 2.30, 2.32, 2.34, 2.36, 2.38,
  2.40, 2.42, 2.44, 2.46, 2.48, 2.50, 2.52, 2.54, 2.56, 2.58,
  2.60, 2.62, 2.64, 2.66, 2.68, 2.70, 2.72, 2.74, 2.76, 2.78,
  2.80, 2.82, 2.84, 2.86, 2.88, 2.90, 2.92, 2.94, 2.96, 2.98,
  3.00, 3.02, 3.04, 3.06, 3.08, 3.10, 3.12, 3.14, 3.16, 3.18,
  3.20, 3.22, 3.24, 3.26, 3.28, 3.30, 3.32, 3.34, 3.36, 3.38,
  3.40, 3.42, 3.44, 3.46, 3.48, 3.50, 3.52, 3.54, 3.56, 3.58,
  3.60, 3.62, 3.64, 3.66, 3.68, 3.70, 3.72, 3.74, 3.76, 3.78,
  3.80};
  for(int i=0; i<nEntries91; i++) PhotonEnergy91[i]*=eV;

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

  G4double AbsWLSfiber[nEntries279] = 
  {
 22.3,  21.6,  20.2,  20.0,  19.2,  19.4,  17.8,  17.0,  15.9,  15.2,    //as measured by Yuri 
 14.6,  13.6,  12.7,  11.5,  10.5,  10.3,   9.4,   8.7,   8.5,   7.8,    //these photons are lost, i.e. no wavelength shifting
  7.4,   7.5,   5.7,   6.8,   6.8,   6.8,   7.3,   7.0,   7.1,   7.6,    //...
  7.4,   7.9,   8.4,   8.7,   9.0,   9.5,  10.5,  10.3,  11.0,  11.5, 
 12.1,  12.9,  12.8,  13.4,  13.7,  14.5,  15.1,  15.3,  15.9,  16.9, 
 17.0,  17.3,  18.1,  18.8,  19.6,  19.7,  19.6,  20.3,  20.5,  21.2, 
 21.3,  21.5,  22.1,  22.1,  22.4,  22.7,  22.8,  23.9,  23.2,  23.5, 
 23.5,  23.4,  23.6,  23.7,  23.9,  24.1,  24.0,  24.3,  24.5,  24.4, 
 24.2,  24.4,  24.1,  24.4,  24.4,  24.7,  24.4,  24.1,  24.1,  23.6, 
 24.1,  23.6,  23.3,  23.3,  23.1,  22.9,  23.0,  22.5,  22.4,  22.2, 
 22.2,  21.8,  21.8,  21.8,  21.5,  21.5,  21.6,  21.3,  21.4,  21.7, 
 22.0,  21.3,  21.2,  21.2,  21.3,  21.3,  21.2,  21.1,  21.2,  21.6, 
 21.3,  21.2,  21.1,  21.2,  21.2,  21.3,  21.0,  21.0,  20.8,  20.8, 
 21.1,  20.8,  20.6,  20.6,  20.6,  20.5,  20.4,  20.3,  20.3,  20.1, 
 20.1,  20.0,  19.8,  19.8,  19.7,  19.5,  19.4,  19.4,  19.1,  19.0, 
 18.8,  19.0,  18.4,  18.3,  18.3,  18.1,  17.9,  17.7,  17.5,  17.3, 
 17.1,  16.9,  16.6,  16.5,  16.2,  16.0,  15.8,  15.5,  15.2,  14.9, 
 14.6,  14.2,  13.9,  13.5,  13.1,  12.8,  12.4,  12.0,  11.7,  11.4, 
 11.2,  10.9,  10.7,  10.5,  10.4,  10.3,  10.2,  10.1,  10.1,  10.1, 
 10.1,  10.0,  10.0,  10.0,  10.0,  10.0,   9.9,   9.9,   9.9,   9.8, 
  9.8,   9.7,   9.6,   9.5,   9.4,   9.4,   9.3,   9.1,   9.0,   8.9, 
  8.7,   8.6,   8.4,   8.3,   8.1,   8.0,   7.8,   7.7,   7.4,   7.3, 
  7.1,   6.9,   6.7,   6.5,   6.3,   6.1,   5.9,   5.7,   5.5,   5.3, 
  5.1,   4.9,   4.7,   4.5,   4.3,   4.1,   3.9,   3.7,   3.5,   3.3, 
  3.1,   3.0,   2.8,   2.6,   2.5,   2.3,   2.2,   2.0,   1.9,   1.8, 
  1.7,   1.6,   1.5,   1.4,   1.3,   1.2,   1.1,   1.0,   1.0,   0.9, 
  0.8,   0.8,   0.7,   0.7,   0.6,   0.6,   0.6,   0.5,   0.5,   0.5, 
  0.4,   0.4,   0.4,   0.4,   0.4,   0.3,                                 //...
  0.1,   1e3,   1e3};                                                     //these photons are not lost - they will be wavelength shifted (see below)

  for(int i=0; i<nEntries279; i++) AbsWLSfiber[i]*=0.84*m;                //fudge factor 
  AbsWLSfiber[nEntries279-3]=111.19*mm;                                   //to match the value in AbsWLSfiberWLS

  G4double AbsWLSfiberWLS[nEntries65] =
  {
  1e6,    1e6,                                                            //no wavelength shifting in these photons - they will be lost (see above)
  111.19, 36.71, 13.44,  5.67,  2.64,  1.57,  0.99,  0.70,  0.61,  0.64,  //these photons will be wavelength shifted
    0.72,  0.80,  0.80,  0.66,  0.52,  0.43,  0.39,  0.48,  0.64,  0.80,
    0.96,  1.12,  1.21,  1.24,  1.37,  1.49,  1.66,  1.86,  2.16,  2.46,
    2.83,  3.30,  3.93,  4.78,  5.67,  6.91,  8.06,  9.63, 11.88, 14.37,
   16.66, 17.49, 19.07, 20.96, 21.82, 24.30, 26.07, 29.67, 31.39, 34.39,
   38.00, 42.44, 48.05, 55.33, 65.19, 79.27,101.03,139.11,222.90,558.02,
  1e3,    1e3,   1e3};
  for(int i=0; i<nEntries65; i++) AbsWLSfiberWLS[i]*=1.00*mm;

  G4double EmissionFib[nEntries43] =
  {
     1,   2,   4,   5,   5,   6,   8,  10,  13,  17,
    23,  28,  34,  39,  44,  50,  54,  60,  75,  94,
   108, 123, 129, 126, 117, 110, 112, 122, 134, 147,
   148, 131, 100,  71,  42,  23,  14,   7,   2,   1,
     1,
     0,   0};

  G4MaterialPropertiesTable* MPTWLSfiber = new G4MaterialPropertiesTable();
  MPTWLSfiber->AddProperty("RINDEX",PhotonEnergy2,RefractiveIndexWLSfiber,nEntries2);
  MPTWLSfiber->AddProperty("ABSLENGTH",PhotonEnergy279,AbsWLSfiber,nEntries279);
  MPTWLSfiber->AddProperty("WLSABSLENGTH",PhotonEnergy65,AbsWLSfiberWLS,nEntries65);
  MPTWLSfiber->AddProperty("WLSCOMPONENT",PhotonEnergy43,EmissionFib,nEntries43);
  MPTWLSfiber->AddConstProperty("WLSTIMECONSTANT", 8.7*ns);

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

  G4double AbsPS[nEntries17] = 
  {49.5, 49.5, 32.8, 28.1, 24.5, 19.5, 16.2, 9.5, 4.5, 2.8,
    2.0,  0.8,  0.4,  0.3,  0.2,  0.2, 0.0};
  for(int i=0; i<nEntries17; i++) AbsPS[i]*=3.0*cm;  //seems to give the best result

  G4double ScintilFast[nEntries91] =
  {
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,  10,  10,  10,   9,   9,   9,   9,  14,
    19,  24,  27,  31,  36,  43,  52,  62,  75,  87,
    99, 111, 122, 130, 146, 171, 210, 246, 276, 304,
   313, 314, 316, 318, 344, 392, 446, 509, 531, 529,
   500, 466, 433, 422, 443, 479, 503, 507, 481, 440,
   400, 381, 362, 348, 335, 321, 308, 314, 328, 350,
   363, 351, 330, 292, 248, 212, 186, 180, 185, 179,
   162, 114,  68,  37,  22,   7,   0,   0,   0,   0,
   0};

  for(int i=55; i<nEntries91; i++) ScintilFast[i]=0;  //Absorption lengths for these energies are very short, but unknown

  G4MaterialPropertiesTable* MPTPolystyrene = new G4MaterialPropertiesTable();
  MPTPolystyrene->AddProperty("RINDEX",PhotonEnergy2,RefractiveIndexPS,nEntries2);
  MPTPolystyrene->AddProperty("ABSLENGTH",PhotonEnergy17,AbsPS,nEntries17);
  MPTPolystyrene->AddProperty("FASTCOMPONENT",PhotonEnergy91, ScintilFast,nEntries91); 
  MPTPolystyrene->AddProperty("SLOWCOMPONENT",PhotonEnergy91, ScintilFast,nEntries91);  //assumed to be delayed flourescence, not used
  MPTPolystyrene->AddConstProperty("SCINTILLATIONYIELD",47000./MeV); //to match the testbeam number of PEs
  MPTPolystyrene->AddConstProperty("RESOLUTIONSCALE",1.0);
  MPTPolystyrene->AddConstProperty("FASTTIMECONSTANT", 3.*ns);    //includes WLS components in the scintillator
  MPTPolystyrene->AddConstProperty("SLOWTIMECONSTANT", 100.*ns);  //unknown, not used
  MPTPolystyrene->AddConstProperty("YIELDRATIO", 1.0); 
  //MPTPolystyrene->AddConstProperty("ISOTHERMAL_COMPRESSIBILITY", 1.0e-19*m3/MeV);  //for Rayleigh scattering 
                                                                                   //mean free path length is inverse proportional to this value
                                                                                   //value of 1.0e-19*m3/MeV and n=1.59 gives a 
                                                                                   //mean free path of 14.3cm [1.1cm] for 2.0eV [3.8eV] photons

  Polystyrene->SetMaterialPropertiesTable(MPTPolystyrene);

  Polystyrene->GetIonisation()->SetBirksConstant(0.126*mm/MeV);
}
