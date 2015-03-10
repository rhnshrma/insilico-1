/*
  core/engine.hpp - insilico's Simulation Engine API header and source

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

#ifndef INCLUDED_INSILICO_CORE_GLOBAL_TYPE_HPP
#define INCLUDED_INSILICO_CORE_GLOBAL_TYPE_HPP

#include <iostream>
#include <vector>

namespace insilico {

using state_type = std::vector<double>;

} // namespace insilico

#endif