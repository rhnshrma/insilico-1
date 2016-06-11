/*
 neuron/Izh.hpp Izhikevich model(2003)
 http://www.izhikevich.org/publications/spikes.htm
 look up above website for more details on equation

 Copyright (C) 2015 Pratyush Ramakrishna, Suhita Nadkarni Lab,
                    IISER, Pune <pratyush.ramakrishna@students.iiserpune.ac.in>

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

#include "insilico/core/engine.hpp"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

namespace insilico
{

class Izh : public Neuron
{

public:

    void ode_set(state_type &variables, state_type &dxdt,
                 const double t, unsigned index) {

        // Getting parameter values from file. Pick file of the
        // neuron type you want to use while calling the neuron file
        // (while using -n<filename>) and also pick injected current
        // file if needed. Template codes for injector current is
        // written in examples

        // v is voltage
        unsigned v_index = engine::neuron_index(index, "v");
        double v = variables[v_index];

        // u is recovery variable
        unsigned u_index = engine::neuron_index(index, "u");
        double u = variables[u_index];

        double a = engine::neuron_value(index, "a");
        double b = engine::neuron_value(index, "b");
        double c = engine::neuron_value(index, "c");
        double d = engine::neuron_value(index, "d");
        double I = engine::neuron_value(index, "I");

        if (v >= 30) {
            variables[v_index] = c;
            variables[u_index] = u + d;
        }

        v = variables[v_index];
        u = variables[u_index];

        // To inject current for short time periods or to inject
        // variable current, below code gets current values saved
        // in a text file. This text file must be included while
        // executing the program by including a  -e (external current)
        // and adding file name. Code for creating the text file
        // which has variable external current values should be written
        // by user. Sample codes for different types have been given.

        std::vector<double> I_inj = injector::external_current(index, t);

        double current = 0;
        if (!I_inj.empty()) {
            current = I_inj[0];
        }

        dxdt[u_index] = a * ((b * v) - u);
        dxdt[v_index] = (0.04 * v * v) + (5.0 * v) + 140 - u + I + current;

        // in above equation, constant external current and variable
        // injected current have been shown in two different variables
        // (I and current).

        // There are many neuron types. Each neuron behaves in required
        // way only if you put in the right parameter values for the
        // injected current.

    } // function ode_set

}; // class Izh

} // insilico
