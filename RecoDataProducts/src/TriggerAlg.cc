#include "RecoDataProducts/inc/TriggerAlg.hh"
#include <stdexcept>
#include <iostream>
#include <stdio.h>

namespace mu2e {

  std::string const& TriggerAlgDetail::typeName()
  {
     static std::string type("TriggerAlg");
     return type;
  }

  std::map<std::string,TriggerAlgDetail::mask_type> const& TriggerAlgDetail::bitNames()
  {
     static std::map<std::string,mask_type> bitnames;
     if (bitnames.size()==0)
     {
         bitnames[std::string("Trig0")]   = bit_to_mask(trig0);
         bitnames[std::string("Trig1")]   = bit_to_mask(trig1);
         bitnames[std::string("Trig2")]   = bit_to_mask(trig2);
         bitnames[std::string("Trig3")]   = bit_to_mask(trig3);
         bitnames[std::string("Trig4")]   = bit_to_mask(trig4);
         bitnames[std::string("Trig5")]   = bit_to_mask(trig5);
         bitnames[std::string("Trig6")]   = bit_to_mask(trig6);
         bitnames[std::string("Trig7")]   = bit_to_mask(trig7);
         bitnames[std::string("Trig8")]   = bit_to_mask(trig8);
         bitnames[std::string("Trig9")]   = bit_to_mask(trig9);
	 bitnames[std::string("Trig10")]  = bit_to_mask(trig10);
         bitnames[std::string("Trig11")]  = bit_to_mask(trig11);
         bitnames[std::string("Trig12")]  = bit_to_mask(trig12);
         bitnames[std::string("Trig13")]  = bit_to_mask(trig13);
         bitnames[std::string("Trig14")]  = bit_to_mask(trig14);
         bitnames[std::string("Trig15")]  = bit_to_mask(trig15);
         bitnames[std::string("Trig16")]  = bit_to_mask(trig16);
         bitnames[std::string("Trig17")]  = bit_to_mask(trig17);
         bitnames[std::string("Trig18")]  = bit_to_mask(trig18);
         bitnames[std::string("Trig19")]  = bit_to_mask(trig19);
	 bitnames[std::string("Trig20")]  = bit_to_mask(trig20);
         bitnames[std::string("Trig21")]  = bit_to_mask(trig21);
         bitnames[std::string("Trig22")]  = bit_to_mask(trig22);
         bitnames[std::string("Trig23")]  = bit_to_mask(trig23);
         bitnames[std::string("Trig24")]  = bit_to_mask(trig24);
         bitnames[std::string("Trig25")]  = bit_to_mask(trig25);
         bitnames[std::string("Trig26")]  = bit_to_mask(trig26);
         bitnames[std::string("Trig27")]  = bit_to_mask(trig27);
         bitnames[std::string("Trig28")]  = bit_to_mask(trig28);
         bitnames[std::string("Trig30")]  = bit_to_mask(trig30);
         bitnames[std::string("Trig31")]  = bit_to_mask(trig31);
       }
     return bitnames;
  }

}
