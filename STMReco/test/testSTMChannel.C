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

  mu2e::STMChannel ch_HPGe(mu2e::STMChannel::HPGe);
  std::cout << "This should be \"HPGe\": " << ch_HPGe << " : ";
  if (ch_HPGe.name().find("HPGe") != std::string::npos) {
    std::cout << "OK" << std::endl;
  }
  else {
    std::cout << "FAILED" << std::endl;
    return -1;
  }

  mu2e::STMChannel ch_LaBr(mu2e::STMChannel::LaBr);
  std::cout << "This should be \"LaBr\": " << ch_LaBr << " : ";
  if (ch_LaBr.name().find("LaBr") != std::string::npos) {
    std::cout << "OK" << std::endl;
  }
  else {
    std::cout << "FAILED" << std::endl;
    return -1;
  }

  mu2e::STMChannel ch_unknown(mu2e::STMChannel::unknown);
  std::cout << "This should be \"unknown\": " << ch_unknown << " : ";
  if (ch_unknown.name().find("unknown") != std::string::npos) {
    std::cout << "OK" << std::endl;
  }
  else {
    std::cout << "FAILED" << std::endl;
    return -1;
  }
}
