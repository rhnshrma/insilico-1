/*
 core/configuration.hpp - insilico Configuration header and source

 Copyright (C) 2014 Collins Assisi, Collins Assisi Lab, IISER, Pune
 Copyright (C) 2014-2015 Pranav Kulkarni, Collins Assisi Lab, IISER, Pune <pranavcode@gmail.com>

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef INCLUDED_INSILICO_INCLUDE_CORE_CONFIGURATION_SERIAL_HPP
#define INCLUDED_INSILICO_INCLUDE_CORE_CONFIGURATION_SERIAL_HPP

#include "core/engine.hpp"
#include "core/helper/file.hpp"
#include "core/helper/string.hpp"

#include "core/configuration/parser.hpp"
#include "core/injector.hpp"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace insilico { namespace configuration {

std::ofstream outstream;

// initialization, check and handle commandline arguments
void initialize(int argc, char **argv) {
  bool repeat[4] = {false, false, false, false};
  std::string sargv, output_file, neuron_file, synapse_file, external_current_file;
  std::string error_msg = "[insilico::configuration::initialize] Supply of file more than once not allowed.\n"; 
  std::string usage_error_msg = "[insilico::configuration::initialize] USAGE: ";
  usage_error_msg += argv[0];
  usage_error_msg += " -o<output_file.csv> -n<neuron_file.isf> -s<synapse_file.isf> -e<external_file.isfc>\n\n \
    Options:\n\t-o   Output file\n\t-n   Neuron configuration file\n\t-s   Synapse configuration file (optional)\n \
    \t-e   External current configuration file (optional)\n             \
    \nNo space allowed between option and its value.\n";

  for(int i=1; i < argc; ++i) { sargv += argv[i]; sargv += ' '; }
  std::vector<std::string> cmds = split(sargv, ' ');
  if(cmds.size() < 2 || cmds.size() > 4) {
    std::cerr << usage_error_msg; exit(1);
  }

  std::cerr << "[insilico::configuration::initialize] SUCCESS: Initializing with following parameters:\n";
  for(std::string cmd : cmds) {
    switch(cmd.at(1)) {
      case 'o':
        if(repeat[0]) { std::cerr << error_msg; exit(1); }
        output_file.assign(cmd.begin()+2, cmd.end());
        std::cerr << "Output file: " << output_file << '\n';
        repeat[0] = true;
        break;
      case 'n':
        if(repeat[1]) { std::cerr << error_msg; exit(1); }
        neuron_file.assign(cmd.begin()+2, cmd.end());
        std::cerr << "Neuron file: " << neuron_file << '\n';
        repeat[1] = true;
        break;
      case 's':
        if(repeat[2]) { std::cerr << error_msg; exit(1); }
        synapse_file.assign(cmd.begin()+2, cmd.end());
        std::cerr << "Synapse file: " << synapse_file << '\n';
        repeat[2] = true;
        break;
      case 'e':
        if(repeat[3]) { std::cerr << error_msg; exit(1); }
        external_current_file.assign(cmd.begin()+2, cmd.end());
        std::cerr << "External current file: " << external_current_file << '\n';
        repeat[3] = true;
        break;
      default:
        std::cout << "Default error: " << cmd.at(1) << std::endl;
        std::cerr << usage_error_msg; exit(1);
    }
  }

  if(repeat[0]) { outstream.open(output_file, std::ios::out); }
  else {
    std::cerr << "[insilico::configuration::initialize] Output file is required.\n"; 
    exit(1);
  }
  
  if(repeat[2]) { read(neuron_file, synapse_file); }
  else if(repeat[1]) { read(neuron_file); }
  else {
    std::cerr << "[insilico::configuration::initialize] Neuron configuration file is required.\n";
    exit(1);
  }
 
  if(repeat[3]) { injector::read(external_current_file); }

  std::cerr << "[insilico::configuration::read] SUCCESS: Input file read complete." << '\n';
}
                 
// close all output streams
void finalize() {
  outstream.close();
  std::cerr << "[insilico::configuration::finalize] SUCCESS: Releasing resources."<<'\n';
}

} } // namespace insilico::configuration

#endif
