//
// The volume names in the CSV can contain memory address (i.e. 0xAAAAAAA).
// For comparing geometries these need to be removed
//

std::string clean_volume_names(std::string name) {
  bool all_found = false;
  while (!all_found) {
    auto pos_0x = name.find("0x");
    if (pos_0x == std::string::npos) {
      all_found = true;
    }
    else {
      name = name.erase(pos_0x, 9); // always 9 chars
    }
  }
  return name;
}
