#include "DataProducts/inc/STMChannel.hh"
#include "canvas/Utilities/InputTag.h"

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

  // Test findByName()
  const int n_channels = 4;
  std::string names[n_channels] = {"HPGe", "LaBr", "unknown", "a channel that does not exist"};
  mu2e::STMChannel expectedChannels[n_channels] = {mu2e::STMChannel::HPGe, mu2e::STMChannel::LaBr, mu2e::STMChannel::unknown, mu2e::STMChannel::unknown};
  bool expectedValids[n_channels] = {true, true, false, false};
  for (int i_channel = 0; i_channel < n_channels; ++i_channel) {
    std::string name = names[i_channel];
    mu2e::STMChannel expectedChannel = expectedChannels[i_channel];
    bool expectedValid = expectedValids[i_channel];

    std::cout << "Looking for \"" << name << "\"... " << std::endl;
    mu2e::STMChannel ch = mu2e::STMChannel::findByName(name);
    std::cout << "Expecting to find " << expectedChannel << " and found: " << ch << " : ";
    if (ch == expectedChannel) {
      std::cout << "OK" << std::endl;
    }
    else {
      std::cout << "FAILED" << std::endl;
      return -1;
    }

    std::cout << "isValid() should be " << std::boolalpha << expectedValid << ": isValid() = " << ch.isValid() << " : ";
    if (ch.isValid() == expectedValid) {
      std::cout << "OK" << std::endl;
    }
    else {
      std::cout << "FAILED" << std::endl;
      return -1;
    }

    std::cout << std::endl;
  }

  std::string expected_all_printed = "List of STM channels: \n 0 HPGe\n 1 LaBr\n 2 unknown\n";
  std::stringstream all_printed;
  mu2e::STMChannel::printAll(all_printed);
  std::cout << "Output of printAll(): " << std::endl;
  std::cout << all_printed.str();
  std::cout << "Is it what we expected? : ";
  if (all_printed.str().compare(expected_all_printed) == 0) {
    std::cout << "YES" << std::endl;
  }
  else {
    std::cout << "NO" << std::endl;
    return -1;
  }
  std::cout << std::endl;

  art::InputTag hpgeInstance("", "HPGe");
  std::cout << "Testing STMChannel::getChannel() with art::InputTag " << hpgeInstance << ": ";
  if (mu2e::STMChannel::getChannel(hpgeInstance) == mu2e::STMChannel::HPGe) {
    std::cout << "PASSED" << std::endl;
  }
  else {
    std::cout << "FAILED" << std::endl;
    return -1;
  }
  art::InputTag labrInstance("", "LaBr");
  std::cout << "Testing STMChannel::getChannel() with art::InputTag " << labrInstance << ": ";
  if (mu2e::STMChannel::getChannel(labrInstance) == mu2e::STMChannel::LaBr) {
    std::cout << "PASSED" << std::endl;
  }
  else {
    std::cout << "FAILED" << std::endl;
    return -1;
  }
  art::InputTag hpgeLabel("labelledHPGe", "");
  std::cout << "Testing STMChannel::getChannel() with art::InputTag " << hpgeLabel << ": ";
  if (mu2e::STMChannel::getChannel(hpgeLabel) == mu2e::STMChannel::HPGe) {
    std::cout << "PASSED" << std::endl;
  }
  else {
    std::cout << "FAILED" << std::endl;
    return -1;
  }
  art::InputTag labrLabel("labelledLaBr", "");
  std::cout << "Testing STMChannel::getChannel() with art::InputTag " << labrLabel << ": ";
  if (mu2e::STMChannel::getChannel(labrLabel) == mu2e::STMChannel::LaBr) {
    std::cout << "PASSED" << std::endl;
  }
  else {
    std::cout << "FAILED" << std::endl;
    return -1;
  }

  std::cout << "All tests passed" << std::endl;
}
