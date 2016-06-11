/*
 neuron/Morris_Lecar.hpp - Morris-Lecar model for Spiking Neurons
                                 (Morris-Lecar, 1981)

 Copyright (C) 2015 Pranav Kulkarni, Collins Assisi Lab,
                    IISER, Pune <pranavcode@gmail.com>
 Copyright (C) 2015 Nishant Singh, Suhita Nadkarni Lab,
                    IISER, Pune <nishant.singh@students.iiserpune.ac.in>

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

#ifndef INCLUDED_N_Morris_Lecar_HPP
#define INCLUDED_N_Morris_Lecar_HPP

#include "insilico/core/engine.hpp"

namespace insilico {

class Morris_Lecar : public Neuron {
 public:
  void ode_set(state_type &variables, state_type &dxdt,
               const double t, unsigned index) {
    

    //Initializing Values        
    double g_L = engine::neuron_value(index,"g_L"), v_L = engine::neuron_value(index,"v_L"), g_Ca = engine::neuron_value(index,"g_Ca") , v_Ca =engine::neuron_value(index,"v_Ca"),g_K = engine::neuron_value(index,"g_K"), v_K = engine::neuron_value(index,"v_K");
    double v1 = engine::neuron_value(index,"v1"), v2 = engine::neuron_value(index,"v2");
    double v3 = engine::neuron_value(index,"v3"), v4 = engine::neuron_value(index,"v4"), lamda_n_bar = engine::neuron_value(index,"lamda_n_bar");

    unsigned v_index = engine::neuron_index(index, "v");
    double v = variables[v_index];

	unsigned n_index = engine::neuron_index(index, "n");
    double n = variables[n_index];
	
	double m_ss = (1.0 + tanh((v-v1)/v2))/2;
	double n_ss = (1.0 + tanh((v-v3)/v4))/2;
    double lamda_n = lamda_n_bar*cosh((v-v3)/(2*v4));
    
    // ODE set
    double I_L = g_L*(v - v_L);
    double I_Ca = g_Ca*m_ss*(v - v_Ca);
    double I_K = g_K*n*(v - v_K);
    double I_Ext = engine::neuron_value(index, "iext");
    
    // incoming synaptic currents - commented out because we are not running network simulations
/*    double I_Syn = 0;
    std::vector<unsigned> g1_indices = engine::get_pre_neuron_indices(index, "g1");
    state_type esyn_values = engine::get_pre_neuron_values(index, "esyn");
    
    double g_Syn=0.13;

    for(unsigned iterator = 0; iterator < g1_indices.size(); ++iterator) {
      I_Syn += g_Syn*variables[g1_indices[iterator]] * (v - esyn_values[iterator]);
    }
*/

	dxdt[n_index] = lamda_n*(n_ss - n);
	dxdt[v_index] = I_Ext - I_L - I_Ca - I_K;
//  dxdt[v_index] = I_Ext - I_L - I_Ca - I_K + I_Syn;
    
  } // function ode_set
}; // class N_Morris_Lecar

} // insilico

#endif
