// This is a simplified copy of the ATLAS SiliconProperties class.
//
// Adapted for Mu2e by Andrei Gaponenko.

#include "ExtinctionMonitorFNAL/Digitization/inc/SiliconProperties.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include <cmath>

namespace mu2e {
  namespace ExtMonFNAL {

    namespace {
      using namespace CLHEP;

      // Constants used in the formula
      const double elecHallFactZero  = 1.13;
      const double elecHallFact_drdt = 8e-4;
      const double elecV_sat_0       = 1.53e9 * cm/s;
      const double elecV_sat_exp     = -0.87;
      const double elecE_crit_0      = 1.01   * volt/cm;
      const double elecE_crit_exp    = 1.55;
      const double elecBeta_0        = 2.57e-2;
      const double elecBeta_exp      = 0.66;

      const double holeHallFactZero  = 0.72;
      const double holeHallFact_drdt = -5e-4;
      const double holeV_sat_0       = 1.62e8 * cm/s;
      const double holeV_sat_exp     = -0.52;
      const double holeE_crit_0      = 1.24   * volt/cm;
      const double holeE_crit_exp    = 1.68;
      const double holeBeta_0        = 0.46;
      const double holeBeta_exp      = 0.17;

      const double temperatureZero   = 273.15 * kelvin;

    } // anonymous namespace


    // This value for the number of eh pairs per deposited energy is fairly standard I think.
    // In reality there is some temperture dependence but for the temperature ranges we deal with
    // I don't think variations are too signifcant.
    const double SiliconProperties::s_ehPairsPerEnergyDefault = 1. / (3.62 * eV); // 1 eh pair per 3.62 eV.

    const double SiliconProperties::s_fanoFactorDefault = 0.1;

    SiliconProperties::SiliconProperties()
      :  m_temperature(0),
         m_electricField(0),
         m_electronDriftMobility(0),
         m_holeDriftMobility(0),
         m_electronHallMobility(0),
         m_holeHallMobility(0),
         m_electronDiffusionConstant(0),
         m_holeDiffusionConstant(0),
         m_ehPairsPerEnergy(s_ehPairsPerEnergyDefault),
         m_fanoFactor(s_fanoFactorDefault),
         m_override(false)
    {}


    SiliconProperties::SiliconProperties(double temperature, double electricField)
      :  m_temperature(temperature),
         m_electricField(electricField),
         m_electronDriftMobility(0),
         m_holeDriftMobility(0),
         m_electronHallMobility(0),
         m_holeHallMobility(0),
         m_electronDiffusionConstant(0),
         m_holeDiffusionConstant(0),
         m_ehPairsPerEnergy(s_ehPairsPerEnergyDefault),
         m_override(false)
    {
      setConditions(temperature, electricField);
    }

    void
    SiliconProperties::setConditions(double temperature, double electricField)
    {
      if (!m_override) {
        m_temperature = temperature;
        m_electricField = electricField;
        m_electronDriftMobility = calcElectronDriftMobility(temperature, electricField);
        m_holeDriftMobility = calcHoleDriftMobility(temperature, electricField);
        m_electronHallMobility = calcElectronHallFactor(temperature) * m_electronDriftMobility;
        m_holeHallMobility = calcHoleHallFactor(temperature) * m_holeDriftMobility;
        m_electronDiffusionConstant = calcDiffusionConstant(temperature, m_electronDriftMobility);
        m_holeDiffusionConstant = calcDiffusionConstant(temperature, m_holeDriftMobility);
      }
    }

    double
    SiliconProperties::calcElectronHallFactor(double temperature) const
    {
      // Equation from ATL-INDET-2001-004
      return elecHallFactZero + elecHallFact_drdt * (temperature - temperatureZero);
    }

    double
    SiliconProperties::calcHoleHallFactor(double temperature) const
    {
      // Equation from ATL-INDET-2001-004
      return holeHallFactZero + holeHallFact_drdt * (temperature - temperatureZero);
    }

    // driftMobility
    double
    SiliconProperties::calcDriftMobility(double electricField, double electricField_critical,
                                         double saturationVelocity, double beta) const
    {
      // Equation from ATL-INDET-2001-004
      return saturationVelocity / electricField_critical /
        pow(1. + pow(electricField/electricField_critical, beta), 1./beta);
    }

    double
    SiliconProperties::calcElectronDriftMobility(double temperature, double electricField) const
    {
      // Equations from ATL-INDET-2001-004
      double saturationVelocity = elecV_sat_0 * pow(temperature, elecV_sat_exp);
      double electricField_critical = elecE_crit_0 * pow(temperature, elecE_crit_exp);
      double beta = elecBeta_0 * pow(temperature, elecBeta_exp);

      return calcDriftMobility(electricField, electricField_critical, saturationVelocity, beta);
    }

    double
    SiliconProperties::calcHoleDriftMobility(double temperature, double electricField) const
    {
      // Equations from ATL-INDET-2001-004
      double saturationVelocity = holeV_sat_0 * pow(temperature, holeV_sat_exp);
      double electricField_critical = holeE_crit_0 * pow(temperature, holeE_crit_exp);
      double beta = holeBeta_0 * pow(temperature, holeBeta_exp);

      return calcDriftMobility(electricField, electricField_critical, saturationVelocity, beta);
    }


    double
    SiliconProperties::calcDiffusionConstant(double temperature, double mobility) const
    {
      // Einstein's relationship (in many text books)
      return -k_Boltzmann * temperature / electron_charge * mobility; // k_Boltzmann and electron_charge
      // are defined in CLHEP/PhysicalConstants.h
    }

    double
    SiliconProperties::electronDriftMobility() const
    {
      return m_electronDriftMobility;
    }

    double
    SiliconProperties::holeDriftMobility() const
    {
      return m_holeDriftMobility;
    }

    double
    SiliconProperties::electronHallMobility() const
    {
      return m_electronHallMobility;
    }

    double
    SiliconProperties::holeHallMobility() const
    {
      return m_holeHallMobility;
    }

    double
    SiliconProperties::electronDiffusionConstant() const
    {
      return m_electronDiffusionConstant;
    }

    double
    SiliconProperties::holeDiffusionConstant() const
    {
      return m_holeDiffusionConstant;
    }

    void
    SiliconProperties::setElectronDriftMobility(double mobility)
    {
      m_override = true;
      m_electronDriftMobility = mobility;
    }


    void
    SiliconProperties::setHoleDriftMobility(double mobility)
    {
      m_override = true;
      m_holeDriftMobility = mobility;
    }

    void
    SiliconProperties::setElectronHallMobility(double mobility)
    {
      m_override = true;
      m_electronHallMobility = mobility;
    }

    void
    SiliconProperties::setHoleHallMobility(double mobility)
    {
      m_override = true;
      m_holeHallMobility = mobility;
    }


    void
    SiliconProperties::setElectronDiffusionConstant(double diffusionConstant)
    {
      m_override = true;
      m_electronDiffusionConstant = diffusionConstant;
    }

    void
    SiliconProperties::setHoleDiffusionConstant(double diffusionConstant)
    {
      m_override = true;
      m_holeDiffusionConstant = diffusionConstant;
    }

    void
    SiliconProperties::setElectronHolePairsPerEnergy(double ehPairsPerEnergy)
    {
      m_ehPairsPerEnergy = ehPairsPerEnergy;
    }


    double SiliconProperties::electronHolePairsPerEnergy() const
    {
      return m_ehPairsPerEnergy;
    }

  } // namespace ExtMonFNAL
} // namespace mu2e
