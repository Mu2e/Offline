#include "DataProducts/inc/STMChannel.hh"

void testSTMChannel() {
  mu2e::STMChannel ch_default;
  std::cout << "Default constructor should be \"unknown\": " << ch_default << " : ";
  if (ch_default.name().find("unknown") != std::string::npos) {
    std::cout << "OK" << std::endl;
  }
  else {
    std::cout << "FAILED" << std::endl;
    return -1;
  }
  std::cout << "and isValid() should be false: isValid() = " << std::boolalpha << ch_default.isValid() << " : ";
  if (ch_default.isValid() == false) {
    std::cout << "OK" << std::endl;
  }
  else {
    std::cout << "FAILED" << std::endl;
    return -1;
  }

  std::cout << std::endl;

  mu2e::STMChannel ch_HPGe(mu2e::STMChannel::HPGe);
  std::cout << "This should be \"HPGe\": " << ch_HPGe << " : ";
  if (ch_HPGe.name().find("HPGe") != std::string::npos) {
    std::cout << "OK" << std::endl;
  }
  else {
    std::cout << "FAILED" << std::endl;
    return -1;
  }
  std::cout << "and isValid() should be true: isValid() = " << std::boolalpha << ch_HPGe.isValid() << " : ";
  if (ch_HPGe.isValid() == true) {
    std::cout << "OK" << std::endl;
  }
  else {
    std::cout << "FAILED" << std::endl;
    return -1;
  }

  std::cout << std::endl;

  mu2e::STMChannel ch_LaBr(mu2e::STMChannel::LaBr);
  std::cout << "This should be \"LaBr\": " << ch_LaBr << " : ";
  if (ch_LaBr.name().find("LaBr") != std::string::npos) {
    std::cout << "OK" << std::endl;
  }
  else {
    std::cout << "FAILED" << std::endl;
    return -1;
  }
  std::cout << "and isValid() should be true: isValid() = " << std::boolalpha << ch_LaBr.isValid() << " : ";
  if (ch_LaBr.isValid() == true) {
    std::cout << "OK" << std::endl;
  }
  else {
    std::cout << "FAILED" << std::endl;
    return -1;
  }

  std::cout << std::endl;

  mu2e::STMChannel ch_unknown(mu2e::STMChannel::unknown);
  std::cout << "This should be \"unknown\": " << ch_unknown << " : ";
  if (ch_unknown.name().find("unknown") != std::string::npos) {
    std::cout << "OK" << std::endl;
  }
  else {
    std::cout << "FAILED" << std::endl;
    return -1;
  }
  std::cout << "and isValid() should be false: isValid() = " << std::boolalpha << ch_unknown.isValid() << " : ";
  if (ch_unknown.isValid() == false) {
    std::cout << "OK" << std::endl;
  }
  else {
    std::cout << "FAILED" << std::endl;
    return -1;
  }

  std::cout << std::endl;

  std::cout << "All tests passed" << std::endl;
}
