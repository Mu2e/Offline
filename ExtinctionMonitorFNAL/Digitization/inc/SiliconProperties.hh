// This is a simplified copy of the ATLAS SiliconProperties class.
//
// Adapted for Mu2e by Andrei Gaponenko.

#ifndef SiliconProperties_hh
#define SiliconProperties_hh

namespace mu2e {

  namespace ExtMonFNAL {

    class SiliconProperties {
    public:

      SiliconProperties();
      SiliconProperties(double temperature, double electricField);

      void setConditions(double temperature, double electricField);
      double temperature() const { return m_temperature; }
      double electricField() const { return m_electricField; }


      double electronDriftMobility() const;
      double holeDriftMobility() const;
      double electronHallMobility() const;
      double holeHallMobility() const;
      double electronDiffusionConstant() const;
      double holeDiffusionConstant() const;

      double electronHolePairsPerEnergy() const;
      double fanoFactor() const { return m_fanoFactor; }

      // These are mainly for use internally but are left public
      double calcElectronHallFactor(double temperature) const;
      double calcHoleHallFactor(double temperature) const;
      double calcDriftMobility(double electricField, double electricField_critical,
                               double saturationVelocity, double beta) const;
      double calcElectronDriftMobility(double temperature, double electricField) const;
      double calcHoleDriftMobility(double temperature, double electricField) const;
      double calcElectronHallMobility(double temperature, double mobility) const;
      double calcHoleHallMobility(double temperature, double mobility) const;
      double calcDiffusionConstant(double temperature, double mobility) const;

      // Allow overriding calculated quantities.
      // Setting any one (other than setElectronHolePairsPerEnergy) will disable recalculations
      // ie setConditions(temperature, electricField) will have no further effect.
      void setElectronDriftMobility(double mobility);
      void setHoleDriftMobility(double mobility);
      void setElectronHallMobility(double mobility);
      void setHoleHallMobility(double mobility);
      void setElectronDiffusionConstant(double diffusionConstant);
      void setHoleDiffusionConstant(double diffusionConstant);
      void setElectronHolePairsPerEnergy(double ehPairsPerEnergy);

    private:
      double m_temperature;
      double m_electricField;

      double m_electronDriftMobility;
      double m_holeDriftMobility;
      double m_electronHallMobility;
      double m_holeHallMobility;
      double m_electronDiffusionConstant;
      double m_holeDiffusionConstant;
      double m_ehPairsPerEnergy;
      double m_fanoFactor;

      bool   m_override; // signal that quantities are overriden and recalculation are disabled.

      static const double s_ehPairsPerEnergyDefault;
      static const double s_fanoFactorDefault;
    };

  } // namespace ExtMonFNAL
} // namespace mu2e

#endif/*SiliconProperties_hh*/
