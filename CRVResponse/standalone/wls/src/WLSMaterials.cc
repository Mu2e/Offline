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
  delete    FPethylene;
  delete    PolystyreneFiber;
  delete    PolystyreneScint;
  delete    Epoxy;
  delete    Coating;
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
  // PMMA (Fiber Inner Cladding)
  //--------------------------------------------------

  elements.push_back("C");     natoms.push_back(5);
  elements.push_back("H");     natoms.push_back(8);
  elements.push_back("O");     natoms.push_back(2);

  density = 1.190*g/cm3;

  PMMA = nistMan->ConstructNewMaterial("PMMA", elements, natoms, density);

  elements.clear();
  natoms.clear();

  //--------------------------------------------------
  //Fluorinated Polyethylene (Fiber Outer Cladding)
  //--------------------------------------------------

  elements.push_back("C");     natoms.push_back(2);
  elements.push_back("H");     natoms.push_back(4);
//TODO: where is the flour?

  density = 1.430*g/cm3;

  FPethylene = nistMan->ConstructNewMaterial("FPethylene", elements, natoms, density);

  elements.clear();
  natoms.clear();

  //--------------------------------------------------
  // Polystyrene: 
  //--------------------------------------------------
 
  elements.push_back("C");     natoms.push_back(8);
  elements.push_back("H");     natoms.push_back(8);

  density = 1.050*g/cm3;

  PolystyreneFiber = nistMan->ConstructNewMaterial("PolystyreneFiber", elements, natoms, density);
  PolystyreneScint = nistMan->ConstructNewMaterial("PolystyreneScint", elements, natoms, density);

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
  // Epoxy 
  //--------------------------------------------------
 
  elements.push_back("C");     natoms.push_back(21);
  elements.push_back("H");     natoms.push_back(25);
  elements.push_back("Cl");    natoms.push_back(1);
  elements.push_back("O");     natoms.push_back(5);

  density = 1.25*g/cm3;

  Epoxy = nistMan->ConstructNewMaterial("Epoxy", elements, natoms, density);

  elements.clear();
  natoms.clear();

  //--------------------------------------------------
  // TiO2
  //--------------------------------------------------

  elements.push_back("Ti");     natoms.push_back(1);
  elements.push_back("O");      natoms.push_back(2);

  density = 4.26*g/cm3;

  G4Material* TiO2 = nistMan->ConstructNewMaterial("TiO2", elements, natoms, density);

  elements.clear();
  natoms.clear();

  //--------------------------------------------------
  // Scintillator Coating - 30% TiO2 and 70% polystyrene by weight.
  //--------------------------------------------------

  density = 2.01*g/cm3;

  Coating = new G4Material("Coating", density, ncomponents=2);

  Coating->AddMaterial(TiO2,             fractionmass = 30*perCent);
  Coating->AddMaterial(PolystyreneScint, fractionmass = 70*perCent);

  //
  // ------------ Generate & Add Material Properties Table ------------
  //

  //--------------------------------------------------
  // Air
  //--------------------------------------------------

  const G4int nEntriesRefractiveIndexAir  = 2;
  G4double PhotonEnergyRefractiveIndexAir[nEntriesRefractiveIndexAir] = {2.00*eV, 15.75*eV};
  G4double RefractiveIndexAir[nEntriesRefractiveIndexAir] = {1.00, 1.00};

  G4MaterialPropertiesTable* MPTAir = new G4MaterialPropertiesTable();
  MPTAir->AddProperty("RINDEX", PhotonEnergyRefractiveIndexAir, RefractiveIndexAir, nEntriesRefractiveIndexAir);

  Air->SetMaterialPropertiesTable(MPTAir);

  //--------------------------------------------------
  // Fibers Core (polystyrene)
  //--------------------------------------------------

  //refractive index polystyrene 
  //up to 4eV: https://refractiveindex.info/?shelf=organic&book=poly(methyl_methacrylate)&page=Sultanova
  //above 4eV: (http://zeus.phys.uconn.edu/halld/tagger/fp-prototype/polystyrene_abs.pdf
  //assumption: refractive index is the same for the PS+Y11 and PS+PPO+POPOP mixtures
  const G4int nEntriesRefractiveIndexPS = 61;
  G4double PhotonEnergyRefractiveIndexPS[nEntriesRefractiveIndexPS] = 
  {
     2.00,  2.05,  2.10,  2.15,  2.20,  2.25,  2.30,  2.35,  2.40,  2.45,
     2.50,  2.55,  2.60,  2.65,  2.70,  2.75,  2.80,  2.85,  2.90,  2.95,
     3.00,  3.05,  3.10,  3.15,  3.20,  3.25,  3.30,  3.40,  3.50,  3.60,
     3.70,  3.80,  3.90,  4.00,
     4.50,  5.00,  5.25,  5.50,  5.75,  6.00,  6.25,  6.50,  7.00,  7.50, 
     8.00,  8.50,  9.00,  9.50, 10.00, 10.50, 11.00, 11.50, 12.00, 12.50, 
    13.00, 13.50, 14.00, 14.50, 15.00, 15.50, 15.75
  };
  for(int i=0; i<nEntriesRefractiveIndexPS; i++) PhotonEnergyRefractiveIndexPS[i]*=eV;
  G4double RefractiveIndexPS[nEntriesRefractiveIndexPS] = 
  {
    1.589, 1.590, 1.591, 1.593, 1.594, 1.596, 1.597, 1.599, 1.601, 1.602,
    1.604, 1.606, 1.608, 1.609, 1.611, 1.613, 1.615, 1.617, 1.620, 1.622,
    1.624, 1.626, 1.629, 1.631, 1.633, 1.636, 1.639, 1.644, 1.649, 1.655,
    1.661, 1.668, 1.675, 1.682,
    1.72, 1.81, 1.95, 2.10, 2.10, 2.20, 1.70, 1.35, 1.24, 1.40, 
    1.55, 1.62, 1.64, 1.62, 1.60, 1.57, 1.54, 1.50, 1.44, 1.38, 
    1.32, 1.26, 1.20, 1.14, 1.10, 1.04, 1.00
  };


  //Absorption Y11
  //from 2.0eV to 2.531eV -as measured with fiber tester by Yuri
  //                      -values in m
  //                      -scaled by 0.84 to match test beam resuls
  //from 2.531eV to 2.56eV -linear interpolation
  //from 2.56eV to 3.74eV -from Kuraray
  //                      -in mm

  const G4int nEntriesAbsorptionY11 = 325;
  G4double PhotonEnergyAbsorptionY11[nEntriesAbsorptionY11] =
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
  2.519,  2.521,  2.523,                                                           //... 
  2.56, 2.58, 2.60, 2.62, 2.64, 2.66, 2.68, 2.70, 2.72, 2.74,  //as given by Kuraray
  2.76, 2.78, 2.80, 2.82, 2.84, 2.86, 2.88, 2.90, 2.92, 2.94,  //...
  2.96, 2.98, 3.00, 3.02, 3.04, 3.06, 3.08, 3.10, 3.12, 3.14, 
  3.16, 3.18, 3.20, 3.22, 3.24, 3.26, 3.28, 3.30, 3.32, 3.34, 
  3.36, 3.38, 3.40, 3.42, 3.44, 3.46, 3.48, 3.50, 3.52, 3.54, 
  3.5401, 15.75
  };
  for(int i=0; i<nEntriesAbsorptionY11; i++) PhotonEnergyAbsorptionY11[i]*=eV;
  G4double AbsorptionY11[nEntriesAbsorptionY11] = 
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
   3.13,  2.96,  2.77,  2.61,  2.46,  2.31,  2.16,  2.03,  1.90,  1.78, 
   1.67,  1.56,  1.45,  1.36,  1.27,  1.19,  1.11,  1.03,  0.97,  0.90, 
   0.84,  0.78,  0.73,  0.684, 0.644, 0.596, 0.561, 0.526, 0.494, 0.467, 
   0.440, 0.411, 0.399,                                                   //... 
   38.90, 12.97,  4.86,  2.16,  1.11,  0.75,  0.56,  0.47,  0.44,  0.45,  //these photons will be wavelength shifted
    0.47,  0.50,  0.50,  0.46,  0.42,  0.40,  0.39,  0.41,  0.45,  0.50,
    0.55,  0.60,  0.63,  0.64,  0.68,  0.72,  0.78,  0.85,  0.95,  1.05,
    1.18,  1.34,  1.56,  1.85,  2.16,  2.59,  2.99,  3.54,  4.32,  5.19,
    5.98,  6.27,  6.82,  7.48,  7.78,  8.64,  9.26, 10.51, 11.11, 12.16,
    0.0,   0.0
  };
  for(int i=0; i<273; i++) AbsorptionY11[i]*=0.85*m;
  for(int i=273; i<nEntriesAbsorptionY11; i++) AbsorptionY11[i]*=1.00*mm;

  //quantum yield for Y11
  //(http://hallaweb.jlab.org/experiment/PVDIS/SoLID/EC/meetings/kuraray/Absorption&Emission%20of%20Y7811-1.pdf)
  //note: no quantum yields were available, so the quantum yield was set to 1
  //      and a cut was made at the low energy end where the absorption curve ends in the Kuraray plot
  const G4int nEntriesQuantumYieldY11 = 6;
  G4double PhotonEnergyQuantumYieldY11[nEntriesQuantumYieldY11] = 
  {
     2.00,  2.55999, 2.56, 3.54, 3.5401, 15.75
  };
  for(int i=0; i<nEntriesQuantumYieldY11; i++) PhotonEnergyQuantumYieldY11[i]*=eV;
  G4double QuantumYieldY11[nEntriesQuantumYieldY11] =
  {
     0, 0, 1, 1, 0, 0
  };

  //Emission Y11 (from Kuraray)
  //(http://hallaweb.jlab.org/experiment/PVDIS/SoLID/EC/meetings/kuraray/Absorption&Emission%20of%20Y7811-1.pdf)
  const G4int nEntriesEmissionY11 = 43;
  G4double PhotonEnergyEmissionY11[nEntriesEmissionY11] =
  {
  2.00, 2.02, 2.04, 2.06, 2.08, 2.10, 2.12, 2.14, 2.16, 2.18,
  2.20, 2.22, 2.24, 2.26, 2.28, 2.30, 2.32, 2.34, 2.36, 2.38,
  2.40, 2.42, 2.44, 2.46, 2.48, 2.50, 2.52, 2.54, 2.56, 2.58,
  2.60, 2.62, 2.64, 2.66, 2.68, 2.70, 2.72, 2.74, 2.76, 2.78,
  2.80, 
  2.82, 15.75
  };
  for(int i=0; i<nEntriesEmissionY11; i++) PhotonEnergyEmissionY11[i]*=eV;
  G4double EmissionY11[nEntriesEmissionY11] =
  {
     1,   2,   4,   5,   5,   6,   8,  10,  13,  17,
    23,  28,  34,  39,  44,  50,  54,  60,  75,  94,
   108, 123, 129, 126, 117, 110, 112, 122, 134, 147,
   148, 131, 100,  71,  42,  23,  14,   7,   2,   1,
     1,
     0,   0
  };

  G4MaterialPropertiesTable* MPTWLSfiber = new G4MaterialPropertiesTable();
  MPTWLSfiber->AddProperty("RINDEX",PhotonEnergyRefractiveIndexPS,RefractiveIndexPS,nEntriesRefractiveIndexPS);
  MPTWLSfiber->AddProperty("WLSY11ABSLENGTH",PhotonEnergyAbsorptionY11,AbsorptionY11,nEntriesAbsorptionY11);
  MPTWLSfiber->AddProperty("WLSY11QUANTUMYIELD",PhotonEnergyQuantumYieldY11,QuantumYieldY11,nEntriesQuantumYieldY11);
  MPTWLSfiber->AddProperty("WLSY11COMPONENT",PhotonEnergyEmissionY11,EmissionY11,nEntriesEmissionY11);
  MPTWLSfiber->AddConstProperty("WLSY11TIMECONSTANT", 10.0*ns);

  PolystyreneFiber->SetMaterialPropertiesTable(MPTWLSfiber);

  //--------------------------------------------------
  // Fiber Inner Cladding (PMMA)
  //--------------------------------------------------

  //refractive index PMMA 
  //up to 4eV: https://refractiveindex.info/?shelf=organic&book=poly%28methyl_methacrylate%29&page=Sultanova
  const G4int nEntriesRefractiveIndexClad1 = 48;
  G4double PhotonEnergyRefractiveIndexClad1[nEntriesRefractiveIndexClad1] = 
  {
     2.00,  2.10,  2.20,  2.30,  2.40,  2.50,  2.60,  2.70,  2.80,  2.90,
     3.00,  3.10,  3.20,  3.30,  3.40,  3.50,  3.60,  3.70,  3.80,  3.90,
     4.00,
     4.50,  5.00,  5.25,  5.50,  5.75,  6.00,  6.25,  6.50,  7.00,  7.50, 
     8.00,  8.50,  9.00,  9.50, 10.00, 10.50, 11.00, 11.50, 12.00, 12.50, 
    13.00, 13.50, 14.00, 14.50, 15.00, 15.50, 15.75
  };
  for(int i=0; i<nEntriesRefractiveIndexClad1; i++) PhotonEnergyRefractiveIndexClad1[i]*=eV;
  G4double RefractiveIndexClad1[nEntriesRefractiveIndexClad1] = 
  {
    1.489, 1.490, 1.492, 1.493, 1.495, 1.496, 1.498, 1.500, 1.501, 1.503,
    1.505, 1.507, 1.509, 1.512, 1.514, 1.516, 1.519, 1.521, 1.524, 1.527,
    1.530,
    1.72, 1.81, 1.95, 2.10, 2.10, 2.20, 1.70, 1.35, 1.24, 1.40,  //from polystyrene
    1.55, 1.62, 1.64, 1.62, 1.60, 1.57, 1.54, 1.50, 1.44, 1.38, 
    1.32, 1.26, 1.20, 1.14, 1.10, 1.04, 1.00
  };
  //0.1 less than polystyrene starting at 4.5eV (due to lack of data)
  for(int i=21; i<nEntriesRefractiveIndexClad1; i++) RefractiveIndexClad1[i]-=0.15;

  G4MaterialPropertiesTable* MPTClad1 = new G4MaterialPropertiesTable();
  MPTClad1->AddProperty("RINDEX",PhotonEnergyRefractiveIndexClad1,RefractiveIndexClad1,nEntriesRefractiveIndexClad1);

  PMMA->SetMaterialPropertiesTable(MPTClad1);

  //--------------------------------------------------
  // Fiber Outer Cladding (Fluorinated Polyethylene)
  //--------------------------------------------------

  //refractive index outer cladding
  //0.07 less than PMMA (due to the lack of data) 
  G4double RefractiveIndexClad2[nEntriesRefractiveIndexClad1];
  for(int i=0; i<nEntriesRefractiveIndexClad1; i++) RefractiveIndexClad2[i]=RefractiveIndexClad1[i]-0.07;

  G4MaterialPropertiesTable* MPTClad2 = new G4MaterialPropertiesTable();
  MPTClad2->AddProperty("RINDEX",PhotonEnergyRefractiveIndexClad1,RefractiveIndexClad2,nEntriesRefractiveIndexClad1);

  FPethylene->SetMaterialPropertiesTable(MPTClad2);

  //--------------------------------------------------
  // Epoxy (SiPM window)
  //--------------------------------------------------

  //refractive index of Epoxy SiPM window is 1.55 according to Hamamatsu - assumed to be at 589nm. 
  //since no other information is given, the fiber polystyrene index of refraction is used -0.04.
  G4double RefractiveIndexEpoxy[nEntriesRefractiveIndexPS];
  for(int i=0; i<nEntriesRefractiveIndexPS; i++) RefractiveIndexEpoxy[i]=RefractiveIndexPS[i]-0.04;

  G4MaterialPropertiesTable* MPTEpoxy = new G4MaterialPropertiesTable();
  MPTEpoxy->AddProperty("RINDEX", PhotonEnergyRefractiveIndexPS, RefractiveIndexEpoxy, nEntriesRefractiveIndexPS);

  Epoxy->SetMaterialPropertiesTable(MPTEpoxy);

  //--------------------------------------------------
  //  Polystyrene
  //--------------------------------------------------

  //refractive index polystyrene same as for fiber (see above) 


  //emission spectrum for PPO (used for WLS [for PS+PPO] and for scintillation [for PS+PPO])
  //(https://pubs.acs.org/doi/suppl/10.1021/ac062160k/suppl_file/ac062160ksi20061218_105400.pdf)
  //notes: -most emissions (due to WLS and scintillation) happen via an non-radiative energy transfer from PS to PPO
  //        followed by a PPO emission, so that the PS emissions can be neglected.
  //       -for WLS, no distingtion is made between photons being absorped by PS or PPO.
  const G4int nEntriesEmissionPPO = 70;
  G4double PhotonEnergyEmissionPPO[nEntriesEmissionPPO] = 
  {
    2.00, 
    2.48, 2.50, 2.52, 2.54, 2.56, 2.58, 2.60, 2.62, 2.64, 2.66, 
    2.68, 2.70, 2.72, 2.74, 2.76, 2.78, 2.80, 2.82, 2.84, 2.86, 
    2.88, 2.90, 2.92, 2.94, 2.96, 2.98, 3.00, 3.02, 3.04, 3.06, 
    3.08, 3.10, 3.12, 3.14, 3.16, 3.18, 3.20, 3.22, 3.24, 3.26, 
    3.28, 3.30, 3.32, 3.34, 3.36, 3.38, 3.40, 3.42, 3.44, 3.46, 
    3.48, 3.50, 3.52, 3.54, 3.56, 3.58, 3.60, 3.62, 3.64, 3.66, 
    3.68, 3.70, 3.72, 3.74, 3.76, 3.78, 3.80, 
    3.82, 15.75
  };
  for(int i=0; i<nEntriesEmissionPPO; i++) PhotonEnergyEmissionPPO[i]*=eV;
  G4double EmissionPPO[nEntriesEmissionPPO] = 
  {   0, 
      0,  0,  0,  0,  0,  0,  1,  1,  1,  1,
      1,  1,  1,  1,  2,  2,  2,  3,  3,  4,
      4,  5,  6,  7,  8,  9, 10, 12, 13, 16,
     19, 21, 24, 27, 29, 31, 33, 36, 42, 48,
     55, 59, 63, 64, 64, 59, 62, 70, 80, 92,
    100, 96, 89, 77, 66, 55, 51, 60, 73, 85,
     82, 62, 40, 20, 11,  6,  5, 
      0,  0
  };

  //quantum yield for PS+PPO (used for WLS)
  //(https://aip.scitation.org/doi/pdf/10.1063/1.1840616)  
  //(https://www.researchgate.net/publication/232974179_Photodynamics_of_OLED_triplet_emitters_lrppy3_and_PtOEP)
  //(https://pubs.acs.org/doi/suppl/10.1021/ac062160k/suppl_file/ac062160ksi20061218_105400.pdf)
  //notes: -quantum efficiciencies for PS+PPO not given for energies less than 4.35eV, 
  //       -use a value of 0.94 for energies below 4.35eV (from Boens et al., Supporting Information for 
  //        Fluorescence Lifetime Standards for Time and Frequency Domain Fluorescence Spectroscopy), 
  //        which is the quantum yield of PPO (can be justified, because at these energies, almost all photons
  //        will be absorbed by PPO and not by PS due to the shorter absorption length of PPO)
  //       -quantum efficiciencies for PS+PPO for energies less than 3.4eV are set to zero
  //        (absorption length for PPO is much greater than for PS for energies less than 3.4eV,
  //        absorption at these transition energies is dominated by POPOP, so that this sudden jump doesn't matter)
  //       -no attenuation lengths (ABSLENGTH) will be specified for PS+PPO, the attenuation is handled by the quantum yield
  const G4int nEntriesQuantumYieldPSPPO = 27;
  G4double PhotonEnergyQuantumYieldPSPPO[nEntriesQuantumYieldPSPPO] = 
  {
     2.00,  3.39,  3.40,  4.34,  4.35,  4.43,  4.59,  4.77,  4.96,  5.17,
     5.39,  5.64,  5.90,  6.20,  6.53,  6.89,  7.29,  7.75,  8.27,  8.86,
     9.54, 10.33, 11.27, 12.40, 13.78, 15.50, 17.71
  };
  for(int i=0; i<nEntriesQuantumYieldPSPPO; i++) PhotonEnergyQuantumYieldPSPPO[i]*=eV;
  G4double QuantumYieldPSPPO[nEntriesQuantumYieldPSPPO] =
  {
     0.00, 0.00, 0.94, 0.94, 0.78, 0.74, 0.67, 0.63, 0.60, 0.65,
     0.68, 0.57, 0.48, 0.39, 0.38, 0.48, 0.53, 0.58, 0.56, 0.49,
     0.43, 0.40, 0.38, 0.35, 0.35, 0.38, 0.43
  };

  //absorption length for PS+PPO (used for attenuation [via quantum yield] and WLS)
  //(http://zeus.phys.uconn.edu/halld/tagger/fp-prototype/polystyrene_abs.pdf)
  //(https://arxiv.org/pdf/1112.5941.pdf)
  //(https://www.sciencedirect.com/science/article/pii/S0168900206011946)
  //notes: -for energies greater/equal than 3.72eV, the absorption lengths for PS and PPO are so small
  //        that a value of 1e-3cm was used for all energies greater than 3.72eV
  //       -for energies less than 3.72eV, the absorption of PS+PPO is dominated by PPO,
  //        so that the PPO absorption spectrum is used
  //       -for energies less than 3.40eV, the absorption lengths for PS+PPO are dominated by PS
  //        (the absorption lengths for PS+PPO are negligible compared to POPOP in around this transition energy,
  //         so that this sudden jump doesn't matter)
  const G4int nEntriesAbsorptionPSPPO = 45;
  G4double PhotonEnergyAbsorptionPSPPO[nEntriesAbsorptionPSPPO] = 
  { 
    2.00, 2.04, 2.08, 2.12, 2.16, 2.20, 2.24, 2.28, 2.32, 2.36, //PS
    2.40, 2.44, 2.48, 2.52, 2.56, 2.60, 2.64, 2.68, 2.72, 2.76, //PS
    2.80, 2.84, 2.88, 2.92, 2.96, 3.00, 3.04, 3.08, 3.12, 3.16, //PS
    3.20, 3.24, 3.28, 2.32, 3.36, 3.40, 3.44, 3.48, 3.52, 3.56, //PS/PPO
    3.60, 3.64, 3.68, 3.72, 15.75                               //PPO/PS+PPO
  };
  for(int i=0; i<nEntriesAbsorptionPSPPO; i++) PhotonEnergyAbsorptionPSPPO[i]*=eV;
  G4double AbsorptionPSPPO[nEntriesAbsorptionPSPPO] = 
  {
    449,  361,  449,  449,  439,  430,  421,  421,  396,  404,  //PS 
    396,  389,  381,  367,  355,  342,  331,  321,  311,  297,  //PS
    285,  273,  253,  241,  230,  220,  202,  180,  135,   92,  //PS
     72,   56,   39,   29,   23, 8.43, 2.05, .637, .191, .110,  //PS/PPO
2.87e-2, 3.82e-3, 1.06e-3, 1e-3, 1e-3                           //PPO/PS+PPO
  };
  for(int i=0; i<nEntriesAbsorptionPSPPO; i++) AbsorptionPSPPO[i]*=cm;


  //emission spectrum for POPOP (used for WLS)
  //(https://pubs.acs.org/doi/suppl/10.1021/ac062160k/suppl_file/ac062160ksi20061218_105400.pdf)

  const G4int nEntriesEmissionPOPOP = 63;
  G4double PhotonEnergyEmissionPOPOP[nEntriesEmissionPOPOP] = 
  {
    2.00, 2.24,
    2.26, 2.28, 2.30, 2.32, 2.34, 2.36, 2.38, 2.40, 2.42, 2.44, 
    2.46, 2.48, 2.50, 2.52, 2.54, 2.56, 2.58, 2.60, 2.62, 2.64, 
    2.66, 2.68, 2.70, 2.72, 2.74, 2.76, 2.78, 2.80, 2.82, 2.84, 
    2.86, 2.88, 2.90, 2.92, 2.94, 2.96, 2.98, 3.00, 3.02, 3.04, 
    3.06, 3.08, 3.10, 3.12, 3.14, 3.16, 3.18, 3.20, 3.22, 3.24, 
    3.26, 3.28, 3.30, 3.32, 3.34, 3.36, 3.38, 3.40, 3.42, 3.44, 
   15.75
  };
  for(int i=0; i<nEntriesEmissionPOPOP; i++) PhotonEnergyEmissionPOPOP[i]*=eV;
  G4double EmissionPOPOP[nEntriesEmissionPOPOP] = 
  {
      0,  0,
      1,  1,  1,  1,  1,  2,  2,  3,  3,  4, 
      5,  6,  6,  7,  7,  9, 10, 13, 15, 19, 
     21, 23, 24, 23, 26, 28, 35, 43, 52, 56, 
     58, 58, 54, 49, 50, 57, 74, 89,100, 97, 
     85, 73, 56, 46, 44, 56, 78, 94, 88, 60, 
     36, 18,  9,  4,  2,  1,  1,  1,  1,  0, 
      0
  };

  //quantum yield for POPOP (used for WLS)
  //(https://pubs.acs.org/doi/suppl/10.1021/ac062160k/suppl_file/ac062160ksi20061218_105400.pdf)
  //note: only one quantum yield is available, so it will be used for all wavelength
  //      even though it is probably not constant

  const G4int nEntriesQuantumYieldPOPOP = 2;
  G4double PhotonEnergyQuantumYieldPOPOP[nEntriesQuantumYieldPOPOP] = 
  {
     2.00,  15.75
  };
  for(int i=0; i<nEntriesQuantumYieldPOPOP; i++) PhotonEnergyQuantumYieldPOPOP[i]*=eV;
  G4double QuantumYieldPOPOP[nEntriesQuantumYieldPOPOP] =
  {
     0.97, 0.97    //POPOP in cyclohexane (Boens et al.)
  };

  //absorption length for POPOP (used for attenuation [via quantum yield] and WLS)
  //(https://www.sciencedirect.com/science/article/pii/S0168900206011946)
  //notes: -absorption length for energies less than 2.6eV are uncertain/unknown
  //        therefore, they are set to 20m
  //       -absorption length for energies greater than 3.52eV are unknown
  //        therefore, they are set to 0.030cm (=last known value) all the way to the highest energy
  //        (above 3.7eV, the POPOP absorption is negligible compared to PS+PPO)

  const G4int nEntriesAbsorptionPOPOP = 26;
  G4double PhotonEnergyAbsorptionPOPOP[nEntriesAbsorptionPOPOP] = 
  { 
    2.00, 2.60,
    2.64, 2.68, 2.72, 2.76, 2.80, 2.84, 2.88, 2.92, 2.96, 3.00, 
    3.04, 3.08, 3.12, 3.16, 3.20, 3.24, 3.28, 3.32, 3.36, 3.40,
    3.44, 3.48, 3.52,
   15.75
  };
  for(int i=0; i<nEntriesAbsorptionPOPOP; i++) PhotonEnergyAbsorptionPOPOP[i]*=eV;
  G4double AbsorptionPOPOP[nEntriesAbsorptionPOPOP] = 
  {
     2000,  2000,
     1658,  1324,   946,   880,   661,   529,   358,   149,  42.6,  14.7,
     4.95,  1.79, 0.517, 0.192, 0.084, 0.056, 0.045, 0.038, 0.034, 0.031,
    0.030, 0.029, 0.030,
    0.030
  };
  for(int i=0; i<nEntriesAbsorptionPOPOP; i++) AbsorptionPOPOP[i]*=cm;


  G4MaterialPropertiesTable* MPTPolystyrene = new G4MaterialPropertiesTable();

  //for scintillation emission of PS+PPO
  MPTPolystyrene->AddProperty("FASTCOMPONENT",PhotonEnergyEmissionPPO,EmissionPPO,nEntriesEmissionPPO);
  MPTPolystyrene->AddConstProperty("FASTTIMECONSTANT", 1.36*ns);  //PPO in cyclohexane (Boens et al.)
  MPTPolystyrene->AddConstProperty("FASTSCINTILLATIONRISETIME", 1.0*ns);  //exact value not know, but the rise time of PS based scintillators tend to be around 1ns
  MPTPolystyrene->AddConstProperty("SCINTILLATIONYIELD",35900./MeV); //to match the testbeam number of PEs
  MPTPolystyrene->AddConstProperty("RESOLUTIONSCALE",1.0);
  MPTPolystyrene->AddConstProperty("YIELDRATIO", 1.0);  //100% fast component 

  //for cerenkov emission of PS, and light refraction at the PS borders
  MPTPolystyrene->AddProperty("RINDEX",PhotonEnergyRefractiveIndexPS,RefractiveIndexPS,nEntriesRefractiveIndexPS);

  //for WLS absorption/emission of PS+PPO
  MPTPolystyrene->AddProperty("WLSPSPPOCOMPONENT",PhotonEnergyEmissionPPO,EmissionPPO,nEntriesEmissionPPO);
  MPTPolystyrene->AddProperty("WLSPSPPOQUANTUMYIELD",PhotonEnergyQuantumYieldPSPPO,QuantumYieldPSPPO,nEntriesQuantumYieldPSPPO);
  MPTPolystyrene->AddProperty("WLSPSPPOABSLENGTH",PhotonEnergyAbsorptionPSPPO,AbsorptionPSPPO,nEntriesAbsorptionPSPPO);
  MPTPolystyrene->AddConstProperty("WLSPSPPOTIMECONSTANT", 1.36*ns);  //PPO in cyclohexane (Boens et al.)

  //for WLS absorption/emission of POPOP
  MPTPolystyrene->AddProperty("WLSPOPOPCOMPONENT",PhotonEnergyEmissionPOPOP,EmissionPOPOP,nEntriesEmissionPOPOP);
  MPTPolystyrene->AddProperty("WLSPOPOPQUANTUMYIELD",PhotonEnergyQuantumYieldPOPOP,QuantumYieldPOPOP,nEntriesQuantumYieldPOPOP);
  MPTPolystyrene->AddProperty("WLSPOPOPABSLENGTH",PhotonEnergyAbsorptionPOPOP,AbsorptionPOPOP,nEntriesAbsorptionPOPOP);
  MPTPolystyrene->AddConstProperty("WLSPOPOPTIMECONSTANT", 1.12*ns);  //POPOP in cyclohexane (Boesn et al.)

  PolystyreneScint->SetMaterialPropertiesTable(MPTPolystyrene);

  PolystyreneScint->GetIonisation()->SetBirksConstant(0.126*mm/MeV); //https://arxiv.org/pdf/1106.5649v2.pdf
}
