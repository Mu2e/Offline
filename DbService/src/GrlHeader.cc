#include "Offline/DbService/inc/GrlHeader.hh"
#include "Offline/GeneralUtilities/inc/TimeUtility.hh"
#include <iostream>
#include <iomanip>
using namespace mu2e;

//**************************************************
std::string GrlHeader::formatted() const {
  std::ostringstream ss;
  ss << std::setw(20) << std::left << name() << std::right;
  if(locked()) {
    ss << "   locked";
  } else {
    ss << " unlocked";
  }
  ss << std::setw(28) << TimeUtility::reformat1(createTime());
  ss << std::setw(12) << createUser();

  return ss.str();

}
