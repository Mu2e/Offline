// Original author: Andrei Gaponenko, 2012

#include "Offline/Sources/inc/ExtMonFNALMARSUtils.hh"

#include <map>
#include <cmath>

#include "cetlib_except/exception.h"

#include "CLHEP/Vector/LorentzVector.h"

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    std::ostream& operator<<(std::ostream& os, const MARSParticle& mp) {
      return os<<"MARSParticle("<<mp.protonNumber
               <<", pid="<<mp.pid
               <<", kE="<<mp.kineticEnergy
               <<", w="<<mp.weight
               <<", tof="<<mp.tof
               <<")";
    }

    //================================================================
    bool readMARSLine(std::istream& file, MARSParticle& res) {
      std::string line;
      std::getline(file, line);

      if(line.empty()) {
        res.protonNumber = -1U;
        return false;
      }
      else {
        std::istringstream is(line);

        if(!(is>>res.protonNumber
             >>res.pid>>res.kineticEnergy>>res.weight
             >>res.x>>res.y>>res.z
             >>res.dcx>>res.dcy>>res.dcz
             >>res.tof)
           ) {

          throw cet::exception("BADINPUT")<<" parseMARSLine(): error processing line: "<<line<<"\n";
        }
        return true;
      }
    }

    //================================================================
    GenParticle MARSMu2eConverter::marsToMu2eParticle(const MARSParticle& mp) {
      PDGCode::type pdgId = marsToMu2eParticleCode(mp.pid);

      const double mass = pdt_->particle(pdgId).mass();

      const double energy = mass + marsToMu2eEnergy(mp.kineticEnergy);
      const double p3mag = sqrt((energy-mass)*(energy+mass));

      const CLHEP::HepLorentzVector p4(mp.dcx * p3mag,
                                       mp.dcy * p3mag,
                                       mp.dcz * p3mag,
                                       energy
                                       );

      return GenParticle(pdgId,
                         GenId::MARS,
                         marsToMu2ePosition(mp.x, mp.y, mp.z),
                         p4,
                         marsToMu2eTime(mp.tof)
                         );
    }

    //================================================================
    // returns pdgId
    PDGCode::type MARSMu2eConverter::marsToMu2eParticleCode(int marsPID) {
      static std::map<int,int> table;
      if(table.empty()) {
        table[1] = PDGCode::proton; // proton
        table[2] = PDGCode::n0; // neutron
        table[3] = PDGCode::pi_plus; // pi+
        table[4] = PDGCode::pi_minus; // pi-
        table[5] = PDGCode::K_plus; // K+
        table[6] = PDGCode::K_minus; // K-
        table[7] = PDGCode::mu_plus; // mu+
        table[8] = PDGCode::mu_minus; // mu-
        table[9] = PDGCode::gamma; // gamma
        table[10] = PDGCode::e_minus; // e-
        table[11] = PDGCode::e_plus; // e+
        table[12] = PDGCode::anti_proton; // p- (antiproton)
        table[13] = PDGCode::pi0; // pi0
        table[14] = PDGCode::deuteron; // deutron
        table[15] = PDGCode::tritium; // H3
        table[16] = PDGCode::He3; // He3
        table[17] = PDGCode::alpha; // He4
        table[18] = PDGCode::nu_mu; // numu
        table[19] = PDGCode::anti_nu_mu; // numubar
        table[20] = PDGCode::nu_e; // nue
        table[21] = PDGCode::anti_nu_e; // nuebar
        table[22] = PDGCode::K_L0; // K0L
        table[23] = PDGCode::K_S0; // K0S
        //table[24] = ; // K0
        //table[25] = ; // K0bar
        table[26] = PDGCode::Lambda0; // Lambda
        table[27] = PDGCode::anti_Lambda0; // Lambdabar
        table[28] = PDGCode::Sigma_plus; // Sigma+
        table[29] = PDGCode::Sigma0; // Sigma0
        table[30] = PDGCode::Sigma_minus; // Sigma-
        table[31] = PDGCode::anti_n0; // nbar
        table[32] = PDGCode::Xi0; // Xi0
        table[33] = PDGCode::Xi_minus; // Xi-
        table[34] = PDGCode::Omega_minus; // Omega-
        table[35] = PDGCode::anti_Sigma_plus; // anti-Sigma-
        table[36] = PDGCode::anti_Sigma0; // anti-Sigma0
        table[37] = PDGCode::anti_Sigma_minus; // anti-Sigma+
        table[38] = PDGCode::anti_Xi0; // anti-Xi0
        table[39] = PDGCode::anti_Xi_minus; // anti-Xi+
        table[40] = PDGCode::anti_Omega_minus; // anti-Omega+
      }
      std::map<int,int>::const_iterator ip = table.find(marsPID);
      if(ip == table.end()) {
        throw cet::exception("BADINPUT")<<" marsToMu2eParticleCode(): unknonw MARS particle code "<<marsPID<<"\n";
      }
      return PDGCode::type(ip->second);
    }
  } // namespace ExtMonFNAL
} // namespace mu2e
