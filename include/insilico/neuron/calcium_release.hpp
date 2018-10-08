/*
 neuron/N_FA.hpp - KCNQ regulation with M1 activation (FalkenburgerJGP2010)

 Copyright (C) 2014 Pranav Kulkarni, Collins Assisi Lab, IISER, Pune <pranavcode@gmail.com>

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

#ifndef INCLUDED_N_FA_HPP
#define INCLUDED_N_FA_HPP

#include "insilico/core/engine.hpp"



namespace insilico {

        


class calcium_release: public Neuron{
 public:
  void ode_set(state_type &variables, state_type &dxdt, const double t,const unsigned index) {


    unsigned D28kB_high_Cytosol_index = engine::neuron_index(index,"D28kB_high_Cytosol");
    unsigned IP3_D_Cytosol_index = engine::neuron_index(index,"IP3_D_Cytosol");
    unsigned h_ERM_index = engine::neuron_index(index,"h_ERM");
    unsigned CGB_Cytosol_index = engine::neuron_index(index,"CGB_Cytosol");
    unsigned h_D_ERM_index = engine::neuron_index(index,"h_D_ERM");
    unsigned Ca_D_Cytosol_index = engine::neuron_index(index,"Ca_D_Cytosol");
    unsigned D28kB_Cytosol_index = engine::neuron_index(index,"D28kB_Cytosol");
    unsigned PABMg_D_Cytosol_index = engine::neuron_index(index,"PABMg_D_Cytosol");
    unsigned IP3_Cytosol_index = engine::neuron_index(index,"IP3_Cytosol");
    unsigned PA_D_Cytosol_index = engine::neuron_index(index,"PA_D_Cytosol");
    unsigned D28kB_D_Cytosol_index = engine::neuron_index(index,"D28kB_D_Cytosol");
    unsigned CG_D_Cytosol_index = engine::neuron_index(index,"CG_D_Cytosol");
    unsigned D28kB_high_D_Cytosol_index = engine::neuron_index(index,"D28kB_high_D_Cytosol");
    unsigned PABCa_Cytosol_index = engine::neuron_index(index,"PABCa_Cytosol");
    unsigned D28k_high_Cytosol_index = engine::neuron_index(index,"D28k_high_Cytosol");
    unsigned D28k_high_D_Cytosol_index = engine::neuron_index(index,"D28k_high_D_Cytosol");
    unsigned CG_Cytosol_index = engine::neuron_index(index,"CG_Cytosol");
    unsigned PABMg_Cytosol_index = engine::neuron_index(index,"PABMg_Cytosol");
    unsigned CGB_D_Cytosol_index = engine::neuron_index(index,"CGB_D_Cytosol");
    unsigned D28k_Cytosol_index = engine::neuron_index(index,"D28k_Cytosol");
    unsigned PABCa_D_Cytosol_index = engine::neuron_index(index,"PABCa_D_Cytosol");
    unsigned Ca_Cytosol_index = engine::neuron_index(index,"Ca_Cytosol");
    unsigned PA_Cytosol_index = engine::neuron_index(index,"PA_Cytosol");
    unsigned D28k_D_Cytosol_index = engine::neuron_index(index,"D28k_D_Cytosol");


    double D28kB_high_Cytosol = variables[D28kB_high_Cytosol_index];
    double IP3_D_Cytosol = variables[IP3_D_Cytosol_index];
    double h_ERM = variables[h_ERM_index];
    double CGB_Cytosol = variables[CGB_Cytosol_index];
    double h_D_ERM = variables[h_D_ERM_index];
    double Ca_D_Cytosol = variables[Ca_D_Cytosol_index];
    double D28kB_Cytosol = variables[D28kB_Cytosol_index];
    double PABMg_D_Cytosol = variables[PABMg_D_Cytosol_index];
    double IP3_Cytosol = variables[IP3_Cytosol_index];
    double PA_D_Cytosol = variables[PA_D_Cytosol_index];
    double D28kB_D_Cytosol = variables[D28kB_D_Cytosol_index];
    double CG_D_Cytosol = variables[CG_D_Cytosol_index];
    double D28kB_high_D_Cytosol = variables[D28kB_high_D_Cytosol_index];
    double PABCa_Cytosol = variables[PABCa_Cytosol_index];
    double D28k_high_Cytosol = variables[D28k_high_Cytosol_index];
    double D28k_high_D_Cytosol = variables[D28k_high_D_Cytosol_index];
    double CG_Cytosol = variables[CG_Cytosol_index];
    double PABMg_Cytosol = variables[PABMg_Cytosol_index];
    double CGB_D_Cytosol = variables[CGB_D_Cytosol_index];
    double D28k_Cytosol = variables[D28k_Cytosol_index];
    double PABCa_D_Cytosol = variables[PABCa_D_Cytosol_index];
    double Ca_Cytosol = variables[Ca_Cytosol_index];
    double PA_Cytosol = variables[PA_Cytosol_index];
    double D28k_D_Cytosol = variables[D28k_D_Cytosol_index];

    // Load Parameters
    double l_star_IP3deg = engine::neuron_value(index,"l_star_IP3deg");
    double PABMg_F = engine::neuron_value(index,"PABMg_F");
    double r_n_PA_deg = engine::neuron_value(index,"r_n_PA_deg");
    double r_neck_Ca_d = engine::neuron_value(index,"r_neck_Ca_d");
    double mlabfix_F_ = engine::neuron_value(index,"mlabfix_F_");
    double r_spine_PABMg_d = engine::neuron_value(index,"r_spine_PABMg_d");
    double delta = engine::neuron_value(index,"delta");
    double Mg_Cytosol_init_uM = engine::neuron_value(index,"Mg_Cytosol_init_uM");
    double IP3_D_Cytosol_init_uM = engine::neuron_value(index,"IP3_D_Cytosol_init_uM");
    double D_PA_deg = engine::neuron_value(index,"D_PA_deg");
    double lc_IP3deg = engine::neuron_value(index,"lc_IP3deg");
    double h_ERM_init_molecules_um_2 = engine::neuron_value(index,"h_ERM_init_molecules_um_2");
    double CGB_F = engine::neuron_value(index,"CGB_F");
    double l_n_CG_deg = engine::neuron_value(index,"l_n_CG_deg");
    double l_star_PA_deg = engine::neuron_value(index,"l_star_PA_deg");
    double Kf_PA_Ca = engine::neuron_value(index,"Kf_PA_Ca");
    double D_D28k_d = engine::neuron_value(index,"D_D28k_d");
    double r_neck_D28k_high_d = engine::neuron_value(index,"r_neck_D28k_high_d");
    double l_star_PABCa_deg = engine::neuron_value(index,"l_star_PABCa_deg");
    double r_spine_D28kB_high_d = engine::neuron_value(index,"r_spine_D28kB_high_d");
    double Kr_CD28k_high = engine::neuron_value(index,"Kr_CD28k_high");
    double D_D28kB_deg = engine::neuron_value(index,"D_D28kB_deg");
    double r_neck_PABCa_d = engine::neuron_value(index,"r_neck_PABCa_d");
    double Mg_D_Cytosol_init_uM = engine::neuron_value(index,"Mg_D_Cytosol_init_uM");
    double Kr_D28kBDbinding = engine::neuron_value(index,"Kr_D28kBDbinding");
    double r_d_PABMg_deg = engine::neuron_value(index,"r_d_PABMg_deg");
    double lc_PA_deg = engine::neuron_value(index,"lc_PA_deg");
    double PABMg_Cytosol_init_uM = engine::neuron_value(index,"PABMg_Cytosol_init_uM");
    double l_star_D28k_deg = engine::neuron_value(index,"l_star_D28k_deg");
    double D28kB_Cytosol_init_uM = engine::neuron_value(index,"D28kB_Cytosol_init_uM");
    double D_PA_d = engine::neuron_value(index,"D_PA_d");
    double r_D_D28kB_high_deg = engine::neuron_value(index,"r_D_D28kB_high_deg");
    double netValence_flux1 = engine::neuron_value(index,"netValence_flux1");
    double netValence_flux0 = engine::neuron_value(index,"netValence_flux0");
    double l_n_D28k_high_deg = engine::neuron_value(index,"l_n_D28k_high_deg");
    double r_neck_CGB_d = engine::neuron_value(index,"r_neck_CGB_d");
    double Ca_D_ER_init_uM = engine::neuron_value(index,"Ca_D_ER_init_uM");
    double l_CG_d = engine::neuron_value(index,"l_CG_d");
    double D_PABMg_deg = engine::neuron_value(index,"D_PABMg_deg");
    double r_d_CG_deg = engine::neuron_value(index,"r_d_CG_deg");
    double Kf_CGbinding = engine::neuron_value(index,"Kf_CGbinding");
    double l_D28k_d = engine::neuron_value(index,"l_D28k_d");
    double lc_D28kB_high_deg = engine::neuron_value(index,"lc_D28kB_high_deg");
    double l_IP3_d = engine::neuron_value(index,"l_IP3_d");
    double l_n_D28kB_deg = engine::neuron_value(index,"l_n_D28kB_deg");
    double r_neck_PABMg_d = engine::neuron_value(index,"r_neck_PABMg_d");
    double Jch = engine::neuron_value(index,"Jch");
    double r_spine_CG_d = engine::neuron_value(index,"r_spine_CG_d");
    double r_n_D28k_deg = engine::neuron_value(index,"r_n_D28k_deg");
    double Kf_CD28k_high = engine::neuron_value(index,"Kf_CD28k_high");
    double l_star_D28kB_deg = engine::neuron_value(index,"l_star_D28kB_deg");
    double r_neck_PA_d = engine::neuron_value(index,"r_neck_PA_d");
    double Ca_F = engine::neuron_value(index,"Ca_F");
    double PABCa_D_Cytosol_init_uM = engine::neuron_value(index,"PABCa_D_Cytosol_init_uM");
    double Kr_PA_Dbinding = engine::neuron_value(index,"Kr_PA_Dbinding");
    double Ca_Extracellular_init_uM = engine::neuron_value(index,"Ca_Extracellular_init_uM");
    double D_CGB_deg = engine::neuron_value(index,"D_CGB_deg");
    double Size_Extracellular = engine::neuron_value(index,"Size_Extracellular");
    double Size_ERM = engine::neuron_value(index,"Size_ERM");
    double CGB_Cytosol_init_uM = engine::neuron_value(index,"CGB_Cytosol_init_uM");
    double netValence_reaction1 = engine::neuron_value(index,"netValence_reaction1");
    double netValence_reaction0 = engine::neuron_value(index,"netValence_reaction0");
    double r_n_D28kB_high_deg = engine::neuron_value(index,"r_n_D28kB_high_deg");
    double l_star_PABMg_deg = engine::neuron_value(index,"l_star_PABMg_deg");
    double D_D28kB_d = engine::neuron_value(index,"D_D28kB_d");
    double netValence_SERCA_flux = engine::neuron_value(index,"netValence_SERCA_flux");
    double Kr_CG_Dbinding = engine::neuron_value(index,"Kr_CG_Dbinding");
    double D28kB_D_Cytosol_init_uM = engine::neuron_value(index,"D28kB_D_Cytosol_init_uM");
    double Jmax2_IP3R_flux = engine::neuron_value(index,"Jmax2_IP3R_flux");
    double l_n_Ca_deg = engine::neuron_value(index,"l_n_Ca_deg");
    double lc_PABCa_deg = engine::neuron_value(index,"lc_PABCa_deg");
    double lc_D28k_deg = engine::neuron_value(index,"lc_D28k_deg");
    double r_n_PABCa_deg = engine::neuron_value(index,"r_n_PABCa_deg");
    double D28k_Cytosol_init_uM = engine::neuron_value(index,"D28k_Cytosol_init_uM");
    double Kr_PA_Mg = engine::neuron_value(index,"Kr_PA_Mg");
    double K_millivolts_per_volt = engine::neuron_value(index,"K_millivolts_per_volt");
    double mlabfix_K_GHK_ = engine::neuron_value(index,"mlabfix_K_GHK_");
    double Size_PM = engine::neuron_value(index,"Size_PM");
    double netValence_IP3R_fluxD = engine::neuron_value(index,"netValence_IP3R_fluxD");
    double D_D28kB_high_deg = engine::neuron_value(index,"D_D28kB_high_deg");
    double D28kB_F = engine::neuron_value(index,"D28kB_F");
    double Ca_Cytosol_init_uM = engine::neuron_value(index,"Ca_Cytosol_init_uM");
    double r_neck_D28k_d = engine::neuron_value(index,"r_neck_D28k_d");
    double l_CGB_d = engine::neuron_value(index,"l_CGB_d");
    double r_spine_D28kB_d = engine::neuron_value(index,"r_spine_D28kB_d");
    double D_D28k_deg = engine::neuron_value(index,"D_D28k_deg");
    double l_D28k_high_d = engine::neuron_value(index,"l_D28k_high_d");
    double Kdegr_IP3_degr = engine::neuron_value(index,"Kdegr_IP3_degr");
    double lc_D28k_high_deg = engine::neuron_value(index,"lc_D28k_high_deg");
    double D_D28k_high_deg = engine::neuron_value(index,"D_D28k_high_deg");
    double mlabfix_F_nmol_ = engine::neuron_value(index,"mlabfix_F_nmol_");
    double PA_F = engine::neuron_value(index,"PA_F");
    double Voltage_ER = engine::neuron_value(index,"Voltage_ER");
    double Kinh_reaction1 = engine::neuron_value(index,"Kinh_reaction1");
    double Kinh_reaction0 = engine::neuron_value(index,"Kinh_reaction0");
    double lc_D28kB_deg = engine::neuron_value(index,"lc_D28kB_deg");
    double netValence_SERCA_fluxD = engine::neuron_value(index,"netValence_SERCA_fluxD");
    double PABCa_Cytosol_init_uM = engine::neuron_value(index,"PABCa_Cytosol_init_uM");
    double r_n_D28kB_deg = engine::neuron_value(index,"r_n_D28kB_deg");
    double PA_Cytosol_init_uM = engine::neuron_value(index,"PA_Cytosol_init_uM");
    double CGB_D_Cytosol_init_uM = engine::neuron_value(index,"CGB_D_Cytosol_init_uM");
    double D_D28kB_high_d = engine::neuron_value(index,"D_D28kB_high_d");
    double p9 = engine::neuron_value(index,"p9");
    double IP3_Cytosol_init_uM = engine::neuron_value(index,"IP3_Cytosol_init_uM");
    double p8 = engine::neuron_value(index,"p8");
    double p7 = engine::neuron_value(index,"p7");
    double l_n_D28kB_high_deg = engine::neuron_value(index,"l_n_D28kB_high_deg");
    double p6 = engine::neuron_value(index,"p6");
    double p5 = engine::neuron_value(index,"p5");
    double vL_ER_leak_flux = engine::neuron_value(index,"vL_ER_leak_flux");
    double p4 = engine::neuron_value(index,"p4");
    double p3 = engine::neuron_value(index,"p3");
    double p2 = engine::neuron_value(index,"p2");
    double l_n_CGB_deg = engine::neuron_value(index,"l_n_CGB_deg");
    double p1 = engine::neuron_value(index,"p1");
    double p0 = engine::neuron_value(index,"p0");
    double netValence_IP3R_flux = engine::neuron_value(index,"netValence_IP3R_flux");
    double vL_ER_leak_fluxD = engine::neuron_value(index,"vL_ER_leak_fluxD");
    double lc_PABMg_deg = engine::neuron_value(index,"lc_PABMg_deg");
    double Kact_IP3R_flux = engine::neuron_value(index,"Kact_IP3R_flux");
    double IP3_CytosolS = engine::neuron_value(index,"IP3_CytosolS");
    double r_n_PABMg_deg = engine::neuron_value(index,"r_n_PABMg_deg");
    double SVR = engine::neuron_value(index,"SVR");
    double D28kB_high_Cytosol_init_uM = engine::neuron_value(index,"D28kB_high_Cytosol_init_uM");
    double IP3_CytosolD = engine::neuron_value(index,"IP3_CytosolD");
    double JchD = engine::neuron_value(index,"JchD");
    double r_n_CG_deg = engine::neuron_value(index,"r_n_CG_deg");
    double netValence_ER_leak_flux = engine::neuron_value(index,"netValence_ER_leak_flux");
    double r_neck_D28kB_d = engine::neuron_value(index,"r_neck_D28kB_d");
    double mlabfix_T_ = engine::neuron_value(index,"mlabfix_T_");
    double vP_SERCA_fluxD = engine::neuron_value(index,"vP_SERCA_fluxD");
    double CG_Cytosol_init_uM = engine::neuron_value(index,"CG_Cytosol_init_uM");
    double Kdegr_IP3_degr1 = engine::neuron_value(index,"Kdegr_IP3_degr1");
    double D_CG_deg = engine::neuron_value(index,"D_CG_deg");
    double Kf_CG_Dbinding = engine::neuron_value(index,"Kf_CG_Dbinding");
    double l_star_CG_deg = engine::neuron_value(index,"l_star_CG_deg");
    double l_star_D28kB_high_deg = engine::neuron_value(index,"l_star_D28kB_high_deg");
    double r_D_D28k_high_deg = engine::neuron_value(index,"r_D_D28k_high_deg");
    double D_CG_d = engine::neuron_value(index,"D_CG_d");
    double D_IP3_d = engine::neuron_value(index,"D_IP3_d");
    double lc_CG_deg = engine::neuron_value(index,"lc_CG_deg");
    double D28kB_high_F = engine::neuron_value(index,"D28kB_high_F");
    double PABMg_D_Cytosol_init_uM = engine::neuron_value(index,"PABMg_D_Cytosol_init_uM");
    double D28kB_high_D_Cytosol_init_uM = engine::neuron_value(index,"D28kB_high_D_Cytosol_init_uM");
    double kP_SERCA_fluxD = engine::neuron_value(index,"kP_SERCA_fluxD");
    double mlabfix_R_ = engine::neuron_value(index,"mlabfix_R_");
    double Rs = engine::neuron_value(index,"Rs");
    double Kf_PA_MgD = engine::neuron_value(index,"Kf_PA_MgD");
    double r_D_Ca_deg = engine::neuron_value(index,"r_D_Ca_deg");
    double l_Ca_d = engine::neuron_value(index,"l_Ca_d");
    double r_D_D28k_deg = engine::neuron_value(index,"r_D_D28k_deg");
    double r_neck_D28kB_high_d = engine::neuron_value(index,"r_neck_D28kB_high_d");
    double l_n_IP3deg = engine::neuron_value(index,"l_n_IP3deg");
    double r_spine_IP3_d = engine::neuron_value(index,"r_spine_IP3_d");
    double r_spine_Ca_d = engine::neuron_value(index,"r_spine_Ca_d");
    double Ca_ER_init_uM = engine::neuron_value(index,"Ca_ER_init_uM");
    double kP_SERCA_flux = engine::neuron_value(index,"kP_SERCA_flux");
    double tau = engine::neuron_value(index,"tau");
    double h_D_ERM_init_molecules_um_2 = engine::neuron_value(index,"h_D_ERM_init_molecules_um_2");
    double r_neck_CG_d = engine::neuron_value(index,"r_neck_CG_d");
    double l_star_D28k_high_deg = engine::neuron_value(index,"l_star_D28k_high_deg");
    double l_n_PA_deg = engine::neuron_value(index,"l_n_PA_deg");
    double ERDensity_D_ERM_init_molecules_um_2 = engine::neuron_value(index,"ERDensity_D_ERM_init_molecules_um_2");
    double D28k_D_Cytosol_init_uM = engine::neuron_value(index,"D28k_D_Cytosol_init_uM");
    double p11 = engine::neuron_value(index,"p11");
    double D28k_high_D_Cytosol_init_uM = engine::neuron_value(index,"D28k_high_D_Cytosol_init_uM");
    double p10 = engine::neuron_value(index,"p10");
    double l_star_CGB_deg = engine::neuron_value(index,"l_star_CGB_deg");
    double Kf_PA_Mg = engine::neuron_value(index,"Kf_PA_Mg");
    double r_n_Ca_deg = engine::neuron_value(index,"r_n_Ca_deg");
    double CG_D_Cytosol_init_uM = engine::neuron_value(index,"CG_D_Cytosol_init_uM");
    double D_Ca_deg = engine::neuron_value(index,"D_Ca_deg");
    double l_star_Ca_deg = engine::neuron_value(index,"l_star_Ca_deg");
    double Kr_CaD28k_med = engine::neuron_value(index,"Kr_CaD28k_med");
    double Kr_PA_Ca = engine::neuron_value(index,"Kr_PA_Ca");
    double D_CGB_d = engine::neuron_value(index,"D_CGB_d");
    double l_PA_d = engine::neuron_value(index,"l_PA_d");
    double PA_D_Cytosol_init_uM = engine::neuron_value(index,"PA_D_Cytosol_init_uM");
    double l_n_D28k_deg = engine::neuron_value(index,"l_n_D28k_deg");
    double Kf_CD28k_highDbinding = engine::neuron_value(index,"Kf_CD28k_highDbinding");
    double lc_Ca_deg = engine::neuron_value(index,"lc_Ca_deg");
    double Kr_CD28k_highDbinding = engine::neuron_value(index,"Kr_CD28k_highDbinding");
    double Kon_reaction1 = engine::neuron_value(index,"Kon_reaction1");
    double Kon_reaction0 = engine::neuron_value(index,"Kon_reaction0");
    double Size_ER = engine::neuron_value(index,"Size_ER");
    double l_n_PABCa_deg = engine::neuron_value(index,"l_n_PABCa_deg");
    double r_spine_PA_d = engine::neuron_value(index,"r_spine_PA_d");
    double Kf_PA_Dbinding = engine::neuron_value(index,"Kf_PA_Dbinding");
    double D28k_high_F = engine::neuron_value(index,"D28k_high_F");
    double Ca_D_Cytosol_init_uM = engine::neuron_value(index,"Ca_D_Cytosol_init_uM");
    double r_n_CGB_deg = engine::neuron_value(index,"r_n_CGB_deg");
    double r_spine_CGB_d = engine::neuron_value(index,"r_spine_CGB_d");
    double CG_F = engine::neuron_value(index,"CG_F");
    double r_n_D28k_high_deg = engine::neuron_value(index,"r_n_D28k_high_deg");
    double ERDensity_ERM_init_molecules_um_2 = engine::neuron_value(index,"ERDensity_ERM_init_molecules_um_2");
    double vP_SERCA_flux = engine::neuron_value(index,"vP_SERCA_flux");
    double D28k_F = engine::neuron_value(index,"D28k_F");
    double Kf_CaD28k_med = engine::neuron_value(index,"Kf_CaD28k_med");
    double D28k_high_Cytosol_init_uM = engine::neuron_value(index,"D28k_high_Cytosol_init_uM");
    double mlabfix_PI_ = engine::neuron_value(index,"mlabfix_PI_");
    double r_spine_D28k_d = engine::neuron_value(index,"r_spine_D28k_d");
    double l_PABCa_d = engine::neuron_value(index,"l_PABCa_d");
    double Kr_CGbinding = engine::neuron_value(index,"Kr_CGbinding");
    double Ks = engine::neuron_value(index,"Ks");
    double l_D28kB_high_d = engine::neuron_value(index,"l_D28kB_high_d");
    double IP3_F = engine::neuron_value(index,"IP3_F");
    double D_PABCa_d = engine::neuron_value(index,"D_PABCa_d");
    double Js = engine::neuron_value(index,"Js");
    double Size_Cytosol = engine::neuron_value(index,"Size_Cytosol");
    double r_spine_D28k_high_d = engine::neuron_value(index,"r_spine_D28k_high_d");
    double mlabfix_N_pmol_ = engine::neuron_value(index,"mlabfix_N_pmol_");
    double Kf_D28kBDbinding = engine::neuron_value(index,"Kf_D28kBDbinding");
    double l_PABMg_d = engine::neuron_value(index,"l_PABMg_d");
    double l_n_PABMg_deg = engine::neuron_value(index,"l_n_PABMg_deg");
    double l_D28kB_d = engine::neuron_value(index,"l_D28kB_d");
    double t2_flux1 = engine::neuron_value(index,"t2_flux1");
    double t2_flux0 = engine::neuron_value(index,"t2_flux0");
    double r_d_PABCa_deg = engine::neuron_value(index,"r_d_PABCa_deg");
    double Kact_IP3R_fluxD = engine::neuron_value(index,"Kact_IP3R_fluxD");
    double r_neck_IP3_d = engine::neuron_value(index,"r_neck_IP3_d");
    double KMOLE = engine::neuron_value(index,"KMOLE");
    double D_PABMg_d = engine::neuron_value(index,"D_PABMg_d");
    double Voltage_PM = engine::neuron_value(index,"Voltage_PM");
    double dI_IP3R_flux = engine::neuron_value(index,"dI_IP3R_flux");
    double D_Ca_d = engine::neuron_value(index,"D_Ca_d");
    double D_D28k_high_d = engine::neuron_value(index,"D_D28k_high_d");
    double PABCa_F = engine::neuron_value(index,"PABCa_F");
    double Ca_D_Extracellular_init_uM = engine::neuron_value(index,"Ca_D_Extracellular_init_uM");
    double r_d_CGB_deg = engine::neuron_value(index,"r_d_CGB_deg");
    double r_D_D28kB_deg = engine::neuron_value(index,"r_D_D28kB_deg");
    double r_spine_PABCa_d = engine::neuron_value(index,"r_spine_PABCa_d");
    double r_d_IP3deg = engine::neuron_value(index,"r_d_IP3deg");
    double D_PABCa_deg = engine::neuron_value(index,"D_PABCa_deg");
    double r_n_IP3deg = engine::neuron_value(index,"r_n_IP3deg");
    double Jmax2_IP3R_fluxD = engine::neuron_value(index,"Jmax2_IP3R_fluxD");
    double Kr_PA_MgD = engine::neuron_value(index,"Kr_PA_MgD");
    double netValence_ER_leak_fluxD = engine::neuron_value(index,"netValence_ER_leak_fluxD");
    double r_d_PA_deg = engine::neuron_value(index,"r_d_PA_deg");
    double D_IP3deg = engine::neuron_value(index,"D_IP3deg");
    double t1_flux1 = engine::neuron_value(index,"t1_flux1");
    double dI_IP3R_fluxD = engine::neuron_value(index,"dI_IP3R_fluxD");
    double t1_flux0 = engine::neuron_value(index,"t1_flux0");
    double lc_CGB_deg = engine::neuron_value(index,"lc_CGB_deg");

    double J_IP3=0;
    for(int idx=0; idx<=11; idx++){

        double p = engine::neuron_value(index,"p"+std::to_string(idx));
        short flg = 0;
        if (t > idx*tau + delta){
            flg = 1;
        }
        else{
            flg = 0;
        }
        double J_pulse = p*flg*exp(-1*((t - (idx*tau +delta))*Ks)) ;
        
        J_IP3 = J_IP3 + J_pulse;
    }


    // Functions for calculations

    //double J_pulses = (SVR * Js * ( (p0 * (t > delta) * exp( - ((t - delta) * Ks))) + (p1 * (t > (tau + delta)) * exp( - ((t - (tau + delta)) * Ks))) + (p2 * (t > ((2.0 * tau) + delta)) * exp( - (Ks * (t - ((2.0 * tau) + delta))))) + (p3 * (t > ((3.0 * tau) + delta)) * exp( - (Ks * (t - ((3.0 * tau) + delta))))) + (p4 * (t > ((4.0 * tau) + delta)) * exp( - (Ks * (t - ((4.0 * tau) + delta))))) + (p5 * (t > ((5.0 * tau) + delta)) * exp( - (Ks * (t - ((5.0 * tau) + delta))))) + (p6 * (t > ((6.0 * tau) + delta)) * exp( - (Ks * (t - ((6.0 * tau) + delta))))) + (p7 * (t > ((7.0 * tau) + delta)) * exp( - (Ks * (t - ((7.0 * tau) + delta))))) + (p8 * (t > ((8.0 * tau) + delta)) * exp( - (Ks * (t - ((8.0 * tau) + delta))))) + (p9 * (t > ((9.0 * tau) + delta)) * exp( - (Ks * (t - ((9.0 * tau) + delta))))) + (p10 * (t > ((10.0 * tau) + delta)) * exp( - (Ks * (t - ((10.0 * tau) + delta))))) + (p11 * (t > ((11.0 * tau) + delta)) * exp( - (Ks * (t - ((11.0 * tau) + delta)))))) / Rs);
    double J_pulses = SVR*Js*J_IP3;
    double ERDensity_ERM = ERDensity_ERM_init_molecules_um_2;
    double J_SERCA_flux = (ERDensity_ERM * vP_SERCA_flux * Ca_Cytosol * Ca_Cytosol / ((kP_SERCA_flux * kP_SERCA_flux) + (Ca_Cytosol * Ca_Cytosol)));
    double J_CGB_deg = ((D_CGB_deg * pow(r_n_CGB_deg,2.0) * (CGB_D_Cytosol - CGB_Cytosol) / l_n_CGB_deg / pow(r_d_CGB_deg,2.0) / l_star_CGB_deg) + (D_CGB_deg * (CGB_D_Cytosol - CGB_F) / l_star_CGB_deg / lc_CGB_deg));
    double KFlux_PM_Cytosol = (Size_PM / Size_Cytosol);
    double J_D28kB_high_d = (0.75 * D_D28kB_high_d * (D28kB_high_Cytosol - D28kB_high_D_Cytosol) * pow(r_neck_D28kB_high_d,2.0) / l_D28kB_high_d / pow(r_spine_D28kB_high_d,3.0));
    double J_D28kB_d = (0.75 * D_D28kB_d * (D28kB_Cytosol - D28kB_Cytosol) * pow(r_neck_D28kB_d,2.0) / l_D28kB_d / pow(r_spine_D28kB_d,3.0));
    double J_reaction1 =  - ((Kinh_reaction1 - ((Ca_D_Cytosol + Kinh_reaction1) * h_D_ERM)) * Kon_reaction1);
    double J_reaction0 =  - ((Kinh_reaction0 - ((Ca_Cytosol + Kinh_reaction0) * h_ERM)) * Kon_reaction0);
    double ERDensity_D_ERM = ERDensity_D_ERM_init_molecules_um_2;
    double Ca_D_ER = Ca_D_ER_init_uM;
    double J_IP3R_fluxD =  - (ERDensity_D_ERM * Jmax2_IP3R_fluxD * (1.0 - (Ca_D_Cytosol / Ca_D_ER)) * pow((h_D_ERM * IP3_D_Cytosol * Ca_D_Cytosol / (IP3_D_Cytosol + dI_IP3R_fluxD) / (Ca_D_Cytosol + Kact_IP3R_fluxD)),3.0));
    double Mg_Cytosol = Mg_Cytosol_init_uM;
    double J_IP3degrade_IP3_degr = (Kdegr_IP3_degr * (IP3_Cytosol - IP3_CytosolS));
    double KFlux_ERM_Cytosol = (Size_ERM / Size_Cytosol);
    double J_D28k_high_d = (0.75 * D_D28k_high_d * (D28k_high_Cytosol - D28k_high_D_Cytosol) * pow(r_neck_D28k_high_d,2.0) / l_D28k_high_d / pow(r_spine_D28k_high_d,3.0));
    double Ca_ER = Ca_ER_init_uM;
    double J_IP3R_flux =  - (ERDensity_ERM * Jmax2_IP3R_flux * (1.0 - (Ca_Cytosol / Ca_ER)) * (pow(h_ERM * IP3_Cytosol * Ca_Cytosol / (IP3_Cytosol + dI_IP3R_flux) / (Ca_Cytosol + Kact_IP3R_flux),3.0)));
    double J_PABCa_deg = ((D_PABCa_deg * pow(r_n_PABCa_deg,2.0) * (PABCa_D_Cytosol - PABCa_Cytosol) / l_n_PABCa_deg / pow(r_d_PABCa_deg,2.0) / l_star_PABCa_deg) + (D_PABCa_deg * (PABCa_D_Cytosol - PABCa_F) / l_star_PABCa_deg / lc_PABCa_deg));
    double J_CG_deg = ((D_CG_deg * pow(r_n_CG_deg,2.0) * (CG_D_Cytosol - CG_Cytosol) / l_n_CG_deg / pow(r_d_CG_deg,2.0) / l_star_CG_deg) + (D_CG_deg * (CG_D_Cytosol - CG_F) / l_star_CG_deg / lc_CG_deg));
    double J_Ca_d = (0.75 * D_Ca_d * (Ca_Cytosol - Ca_D_Cytosol) * pow(r_neck_Ca_d,2.0) / l_Ca_d / pow(r_spine_Ca_d,3.0));
    double J_D28kB_deg = ((D_D28kB_deg * pow(r_n_D28kB_deg,2.0) * (D28kB_D_Cytosol - D28kB_Cytosol) / l_n_D28kB_deg / pow(r_D_D28kB_deg,2.0) / l_star_D28kB_deg) + (D_D28kB_deg * (D28kB_D_Cytosol - D28kB_F) / l_star_D28kB_deg / lc_D28kB_deg));
    double J_D28k_deg = ((D_D28k_deg * pow(r_n_D28k_deg,2.0) * (D28k_D_Cytosol - D28k_Cytosol) / l_n_D28k_deg / pow(r_D_D28k_deg,2.0) / l_star_D28k_deg) + (D_D28k_deg * (D28k_D_Cytosol - D28k_F) / l_star_D28k_deg / lc_D28k_deg));
    double J_ER_leak_fluxD =  - (ERDensity_D_ERM * vL_ER_leak_fluxD * (1.0 - (Ca_D_Cytosol / Ca_D_ER)));
    double J_D28k_high_deg = ((D_D28k_high_deg * pow(r_n_D28k_high_deg,2.0) * (D28k_high_D_Cytosol - D28k_high_Cytosol) / l_n_D28k_high_deg / pow(r_D_D28k_high_deg,2.0) / l_star_D28k_high_deg) + (D_D28k_high_deg * (D28k_high_D_Cytosol - D28k_high_F) / l_star_D28k_high_deg / lc_D28k_high_deg));
    double J_D28kB_high_deg = ((D_D28kB_high_deg * pow(r_n_D28kB_high_deg,2.0) * (D28kB_high_D_Cytosol - D28kB_high_Cytosol) / l_n_D28kB_high_deg / pow(r_D_D28kB_high_deg,2.0) / l_star_D28kB_high_deg) + (D_D28kB_high_deg * (D28kB_high_D_Cytosol - D28kB_high_F) / l_star_D28kB_high_deg / lc_D28kB_high_deg));
    double J_PABMg_deg = ((D_PABMg_deg * pow(r_n_PABMg_deg,2.0) * (PABMg_D_Cytosol - PABMg_Cytosol) / l_n_PABMg_deg / pow(r_d_PABMg_deg,2.0) / l_star_PABMg_deg) + (D_PABMg_deg * (PABMg_D_Cytosol - PABMg_F) / l_star_PABMg_deg / lc_PABMg_deg));
    double J_D28kBDbinding = (((Kf_D28kBDbinding * D28k_D_Cytosol) * Ca_D_Cytosol) - (Kr_D28kBDbinding * D28kB_D_Cytosol));
    double Ca_D_Extracellular = Ca_D_Extracellular_init_uM;
    double J_flux1 = (JchD * (t > t1_flux1) * (t < t2_flux1) * (Ca_D_Extracellular - Ca_D_Cytosol));
    double Ca_Extracellular = Ca_Extracellular_init_uM;
    double J_flux0 = (Jch * (t > t1_flux0) * (t < t2_flux0) * (Ca_Extracellular - Ca_Cytosol));
    double J_PA_d = (0.75 * D_PA_d * (PA_Cytosol - PA_D_Cytosol) * pow(r_neck_PA_d,2.0) / l_PA_d / pow(r_spine_PA_d,3.0));
    double J_PA_Dbinding = (((Kf_PA_Dbinding * PA_D_Cytosol) * Ca_D_Cytosol) - (Kr_PA_Dbinding * PABCa_D_Cytosol));
    double J_Ca_deg = ((D_Ca_deg * pow(r_n_Ca_deg,2.0) * (Ca_D_Cytosol - Ca_Cytosol) / l_n_Ca_deg / pow(r_D_Ca_deg,2.0) / l_star_Ca_deg) + (D_Ca_deg * (Ca_D_Cytosol - Ca_F) / l_star_Ca_deg / lc_Ca_deg));
    double J_CG_Dbinding = (((Kf_CG_Dbinding * Ca_D_Cytosol) * CG_D_Cytosol) - (Kr_CG_Dbinding * CGB_D_Cytosol));
    double J_CD28k_highDbinding = (((Kf_CD28k_highDbinding * Ca_D_Cytosol) * D28k_high_D_Cytosol) - (Kr_CD28k_highDbinding * D28kB_high_D_Cytosol));
    double J_CaD28k_med = (((Kf_CaD28k_med * D28k_Cytosol) * Ca_Cytosol) - (Kr_CaD28k_med * D28kB_Cytosol));
    double Mg_D_Cytosol = Mg_D_Cytosol_init_uM;
    double J_PACa = (((Kf_PA_Ca * PA_Cytosol) * Ca_Cytosol) - (Kr_PA_Ca * PABCa_Cytosol));
    double J_PABCa_d = (0.75 * D_PABCa_d * (PABCa_Cytosol - PABCa_D_Cytosol) * pow(r_neck_PABCa_d,2.0) / l_PABCa_d / pow(r_spine_PABCa_d,3.0));
    double J_IP3_d = (0.75 * D_IP3_d * (IP3_Cytosol - IP3_D_Cytosol) * pow(r_neck_IP3_d,2.0) / l_IP3_d / pow(r_spine_IP3_d,3.0));
    double J_SERCA_fluxD = (ERDensity_D_ERM * vP_SERCA_fluxD * Ca_D_Cytosol * Ca_D_Cytosol / ((kP_SERCA_fluxD * kP_SERCA_fluxD) + (Ca_D_Cytosol * Ca_D_Cytosol)));
    double J_ER_leak_flux =  - (ERDensity_ERM * vL_ER_leak_flux * (1.0 - (Ca_Cytosol / Ca_ER)));
    double J_PABMg_d = (0.75 * D_PABMg_d * (PABMg_Cytosol - PABMg_D_Cytosol) * pow(r_neck_PABMg_d,2.0) / l_PABMg_d / pow(r_spine_PABMg_d,3.0));
    double J_IP3degrade_IP3_degr1 = (Kdegr_IP3_degr1 * (IP3_D_Cytosol - IP3_CytosolD));
    double J_PAMg_PA_MgD = (((Kf_PA_MgD * Mg_D_Cytosol) * PA_D_Cytosol) - (Kr_PA_MgD * PABMg_D_Cytosol));
    double J_IP3deg = ((D_IP3deg * pow(r_n_IP3deg,2.0) * (IP3_D_Cytosol - IP3_Cytosol) / l_n_IP3deg / pow(r_d_IP3deg,2.0) / l_star_IP3deg) + (D_IP3deg * (IP3_D_Cytosol - IP3_F) / l_star_IP3deg / lc_IP3deg));
    double J_CGB_d = (0.75 * D_CGB_d * (CGB_Cytosol - CGB_D_Cytosol) * pow(r_neck_CGB_d,2.0) / l_CGB_d / pow(r_spine_CGB_d,3.0));
    double J_PAMg_PA_Mg = (((Kf_PA_Mg * PA_Cytosol) * Mg_Cytosol) - (Kr_PA_Mg * PABMg_Cytosol));
    double J_PA_deg = ((D_PA_deg * pow(r_n_PA_deg,2.0) * (PA_D_Cytosol - PA_Cytosol) / l_n_PA_deg / pow(r_d_PA_deg,2.0) / l_star_PA_deg) + (D_PA_deg * (PA_D_Cytosol - PA_F) / l_star_PA_deg / lc_PA_deg));
    double J_D28k_d = (0.75 * D_D28k_d * (D28k_Cytosol - D28k_D_Cytosol) * pow(r_neck_D28k_d,2.0) / l_D28k_d / pow(r_spine_D28k_d,3.0));
    double J_CGbinding = (((Kf_CGbinding * Ca_Cytosol) * CG_Cytosol) - (Kr_CGbinding * CGB_Cytosol));
    double J_CG_d = (0.75 * D_CG_d * (CG_Cytosol - CG_D_Cytosol) * pow(r_neck_CG_d,2.0) / l_CG_d / pow(r_spine_CG_d,3.0));
    double J_CD28k_high = (((Kf_CD28k_high * Ca_Cytosol) * D28k_high_Cytosol) - (Kr_CD28k_high * D28kB_high_Cytosol));

    

    //Storing all functional values
    engine::neuron_value(index,"J_pulses",J_pulses);
    engine::neuron_value(index,"ERDensity_ERM",ERDensity_ERM);
    engine::neuron_value(index,"J_SERCA_flux",J_SERCA_flux);
    engine::neuron_value(index,"J_CGB_deg",J_CGB_deg);
    engine::neuron_value(index,"KFlux_PM_Cytosol",KFlux_PM_Cytosol);
    engine::neuron_value(index,"J_D28kB_high_d",J_D28kB_high_d);
    engine::neuron_value(index,"J_D28kB_d",J_D28kB_d);
    engine::neuron_value(index,"J_reaction1",J_reaction1);
    engine::neuron_value(index,"J_reaction0",J_reaction0);
    engine::neuron_value(index,"ERDensity_D_ERM",ERDensity_D_ERM);
    engine::neuron_value(index,"Ca_D_ER",Ca_D_ER);
    engine::neuron_value(index,"J_IP3R_fluxD",J_IP3R_fluxD);
    engine::neuron_value(index,"Mg_Cytosol",Mg_Cytosol);
    engine::neuron_value(index,"J_IP3degrade_IP3_degr",J_IP3degrade_IP3_degr);
    engine::neuron_value(index,"KFlux_ERM_Cytosol",KFlux_ERM_Cytosol);
    engine::neuron_value(index,"J_D28k_high_d",J_D28k_high_d);
    engine::neuron_value(index,"Ca_ER",Ca_ER);
    engine::neuron_value(index,"J_IP3R_flux",J_IP3R_flux);
    engine::neuron_value(index,"J_PABCa_deg",J_PABCa_deg);
    engine::neuron_value(index,"J_CG_deg",J_CG_deg);
    engine::neuron_value(index,"J_Ca_d",J_Ca_d);
    engine::neuron_value(index,"J_D28kB_deg",J_D28kB_deg);
    engine::neuron_value(index,"J_D28k_deg",J_D28k_deg);
    engine::neuron_value(index,"J_ER_leak_fluxD",J_ER_leak_fluxD);
    engine::neuron_value(index,"J_D28k_high_deg",J_D28k_high_deg);
    engine::neuron_value(index,"J_D28kB_high_deg",J_D28kB_high_deg);
    engine::neuron_value(index,"J_PABMg_deg",J_PABMg_deg);
    engine::neuron_value(index,"J_D28kBDbinding",J_D28kBDbinding);
    engine::neuron_value(index,"Ca_D_Extracellular",Ca_D_Extracellular);
    engine::neuron_value(index,"J_flux1",J_flux1);
    engine::neuron_value(index,"Ca_Extracellular",Ca_Extracellular);
    engine::neuron_value(index,"J_flux0",J_flux0);
    engine::neuron_value(index,"J_PA_d",J_PA_d);
    engine::neuron_value(index,"J_PA_Dbinding",J_PA_Dbinding);
    engine::neuron_value(index,"J_Ca_deg",J_Ca_deg);
    engine::neuron_value(index,"J_CG_Dbinding",J_CG_Dbinding);
    engine::neuron_value(index,"J_CD28k_highDbinding",J_CD28k_highDbinding);
    engine::neuron_value(index,"J_CaD28k_med",J_CaD28k_med);
    engine::neuron_value(index,"Mg_D_Cytosol",Mg_D_Cytosol);
    engine::neuron_value(index,"J_PACa",J_PACa);
    engine::neuron_value(index,"J_PABCa_d",J_PABCa_d);
    engine::neuron_value(index,"J_IP3_d",J_IP3_d);
    engine::neuron_value(index,"J_SERCA_fluxD",J_SERCA_fluxD);
    engine::neuron_value(index,"J_ER_leak_flux",J_ER_leak_flux);
    engine::neuron_value(index,"J_PABMg_d",J_PABMg_d);
    engine::neuron_value(index,"J_IP3degrade_IP3_degr1",J_IP3degrade_IP3_degr1);
    engine::neuron_value(index,"J_PAMg_PA_MgD",J_PAMg_PA_MgD);
    engine::neuron_value(index,"J_IP3deg",J_IP3deg);
    engine::neuron_value(index,"J_CGB_d",J_CGB_d);
    engine::neuron_value(index,"J_PAMg_PA_Mg",J_PAMg_PA_Mg);
    engine::neuron_value(index,"J_PA_deg",J_PA_deg);
    engine::neuron_value(index,"J_D28k_d",J_D28k_d);
    engine::neuron_value(index,"J_CGbinding",J_CGbinding);
    engine::neuron_value(index,"J_CG_d",J_CG_d);
    engine::neuron_value(index,"J_CD28k_high",J_CD28k_high);
    



    // ODE
    dxdt[D28kB_high_Cytosol_index] =        (J_CD28k_high - J_D28kB_high_d);    // rate for D28kB_high_Cytosol
    dxdt[IP3_D_Cytosol_index] =         ( - J_IP3deg - J_IP3degrade_IP3_degr1);    // rate for IP3_D_Cytosol
    dxdt[h_ERM_index] =          - J_reaction0;    // rate for h_ERM
    dxdt[CGB_Cytosol_index] =       (J_CGbinding - J_CGB_d);    // rate for CGB_Cytosol
    dxdt[h_D_ERM_index] =        - J_reaction1;    // rate for h_D_ERM
    dxdt[Ca_D_Cytosol_index] =      ( - J_CD28k_highDbinding - J_PA_Dbinding - J_D28kBDbinding - J_CG_Dbinding - J_Ca_deg - (KFlux_ERM_Cytosol * J_SERCA_fluxD) - (KFlux_ERM_Cytosol * J_IP3R_fluxD) + (KFlux_PM_Cytosol * J_flux1) - (KFlux_ERM_Cytosol * J_ER_leak_fluxD));    // rate for Ca_D_Cytosol
    dxdt[D28kB_Cytosol_index] =         ( - J_D28kB_d + J_CaD28k_med);    // rate for D28kB_Cytosol
    dxdt[PABMg_D_Cytosol_index] =       (J_PAMg_PA_MgD - J_PABMg_deg);    // rate for PABMg_D_Cytosol
    dxdt[IP3_Cytosol_index] =       ( - J_IP3_d - J_IP3degrade_IP3_degr + J_pulses);    // rate for IP3_Cytosol
    dxdt[PA_D_Cytosol_index] =      ( - J_PA_Dbinding - J_PAMg_PA_MgD - J_PA_deg);    // rate for PA_D_Cytosol
    dxdt[D28kB_D_Cytosol_index] =       ( - J_D28kB_deg + J_D28kBDbinding);    // rate for D28kB_D_Cytosol
    dxdt[CG_D_Cytosol_index] =      ( - J_CG_deg - J_CG_Dbinding);    // rate for CG_D_Cytosol
    dxdt[D28kB_high_D_Cytosol_index] =      (J_CD28k_highDbinding - J_D28kB_high_deg);    // rate for D28kB_high_D_Cytosol
    dxdt[PABCa_Cytosol_index] =         (J_PACa - J_PABCa_d);    // rate for PABCa_Cytosol
    dxdt[D28k_high_Cytosol_index] =         ( - J_D28k_high_d - J_CD28k_high);    // rate for D28k_high_Cytosol
    dxdt[D28k_high_D_Cytosol_index] =       ( - J_CD28k_highDbinding - J_D28k_high_deg);    // rate for D28k_high_D_Cytosol
    dxdt[CG_Cytosol_index] =        ( - J_CG_d - J_CGbinding);    // rate for CG_Cytosol
    dxdt[PABMg_Cytosol_index] =         ( - J_PABMg_d + J_PAMg_PA_Mg);    // rate for PABMg_Cytosol
    dxdt[CGB_D_Cytosol_index] =         ( - J_CGB_deg + J_CG_Dbinding);    // rate for CGB_D_Cytosol
    dxdt[D28k_Cytosol_index] =      ( - J_CaD28k_med - J_D28k_d);    // rate for D28k_Cytosol
    dxdt[PABCa_D_Cytosol_index] =       (J_PA_Dbinding - J_PABCa_deg);    // rate for PABCa_D_Cytosol
    dxdt[Ca_Cytosol_index] =        ( - J_PACa - J_Ca_d - J_CGbinding - J_CaD28k_med - J_CD28k_high - (KFlux_ERM_Cytosol * J_ER_leak_flux) - (KFlux_ERM_Cytosol * J_SERCA_flux) + (KFlux_PM_Cytosol * J_flux0) - (KFlux_ERM_Cytosol * J_IP3R_flux));    // rate for Ca_Cytosol
    dxdt[PA_Cytosol_index] =        ( - J_PACa - J_PA_d - J_PAMg_PA_Mg);    // rate for PA_Cytosol
    dxdt[D28k_D_Cytosol_index] =        ( - J_D28k_deg - J_D28kBDbinding);    // rate for D28k_D_Cytosol
  }
};

} // insilico

#endif
