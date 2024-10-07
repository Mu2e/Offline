//
// The volume names in the CSV can contain memory address (i.e. 0xAAAAAAA).
// For comparing geometries these need to be removed
//

std::string clean_volume_names(std::string name) {
  auto pos_0x = name.find("0x");
  if (pos_0x != std::string::npos) {
    name = name.erase(pos_0x, name.size()-pos_0x);
  }

  return name;
}

std::string clean_path_names(std::string name) {

  std::string clean_path;

  std::vector<std::string> tokens;
  std::stringstream ssname(name);
  std::string temp;
  while (std::getline(ssname, temp, '/')) {
    clean_path += clean_volume_names(temp) + "/";
  }

  return clean_path;
}
