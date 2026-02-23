/**
 * @file mctruth.h
 * @brief Definitions of analysis variables which can extract information from
 * the SRTrueInteraction object.
 * @details This file contains definitions of analysis variables which can be
 * used to extract information from the SRTrueInteraction object. Each variable
 * is implemented as a function which takes an SRTrueInteraction object as an
 * argument and returns a double. The association of an SRInteractionTruthDLP
 * object to an SRTrueInteraction object is handled upstream in the SpineVar
 * functions.
 * @author mueller@fnal.gov
 */
#ifndef MCTRUTH_H
#define MCTRUTH_H
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"

#include "framework.h"

/**
 * @namespace mctruth
 * @brief Namespace for organizing variables which act on true interactions.
 * @details This namespace is intended to be used for organizing variables
 * which act on true interactions. Each variable is implemented as a function
 * which takes an SRTrueInteraction object as an argument and returns a double.
 */
namespace mctruth
{
    /**
     * @brief Variable for the true neutrino energy.
     * @details This variable is intended to provide the true energy of the
     * parent neutrino that produced the interaction.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true neutrino energy.
     */
    template<typename T>
        double neutrino_energy(const T & obj) { return obj.E; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, neutrino_energy, neutrino_energy);

    /**
     * @brief Variable for the true neutrino baseline.
     * @details This variable is intended to provide the true baseline of the
     * parent neutrino that produced the interaction.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true neutrino baseline.
     */
    template<typename T>
        double baseline(const T & obj) { return obj.baseline; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, baseline, baseline);

    /**
     * @brief Variable for the true neutrino PDG code.
     * @details This variable is intended to provide the true PDG code of the
     * parent neutrino that produced the interaction.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true neutrino PDG code.
     */
    template<typename T>
        double pdg(const T & obj) { return obj.pdg; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, pdg, pdg);

    /**
     * @brief Variable for the PDG code of the parent of the neutrino.
     * @details This variable is intended to provide the PDG code of the
     * parent of the neutrino that produced the interaction.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the PDG code of the parent of the neutrino.
     */
    template<typename T>
        double parent_pdg(const T & obj) { return obj.parent_pdg; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, parent_pdg, parent_pdg);

    /**
     * @brief Variable for the true neutrino current value.
     * @details This variable is intended to provide the true current value of
     * the parent neutrino that produced the interaction.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the true neutrino current value.
     */
    template<typename T>
        double cc(const T & obj) { return obj.iscc; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, cc, cc);

    /**
     * @brief Variable for the interaction mode of the interaction.
     * @details This variable is intended to provide the interaction mode of the
     * interaction. This is based on the GENIE interaction mode enumeration 
     * defined in the LArSoft MCNeutrino class.
     * @tparam T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the interaction mode.
     * Here for numeration: https://code-doc.larsoft.org/docs/latest/html/MCNeutrino_8h_source.html
     */
    template<typename T>
        double interaction_mode(const T & obj) { return obj.genie_mode; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, interaction_mode, interaction_mode);

    /**
     * @brief Variable for the interaction type of the interaction.
     * @details This variable is intended to provide the interaction type of the
     * interaction. This is based on the GENIE interaction type enumeration 
     * defined in the LArSoft MCNeutrino class.
     * @param T the type of the object to apply the variable on.
     * @param obj the SRTrueInteraction to apply the variable on.
     * @return the interaction type.
     */
    template<typename T>
        double interaction_type(const T & obj) { return obj.genie_inttype; }
    REGISTER_VAR_SCOPE(RegistrationScope::MCTruth, interaction_type, interaction_type);

    /**
     * @brief Cut for exactly one true final state muon above threshold.
     * @param obj the SRTrueInteraction to apply the cut on.
     * @param params KE threshold in MeV, defaults to 143.425 MeV.
     * @return true if exactly one muon above threshold.
     */
    template<typename T>
    bool single_muon_srtruth(const T & obj, std::vector<double> params={143.425,})
    {
        int num_muons(0);
        for(const auto & p : obj.prim)
        {
            if(abs(p.pdg) == 13)
            {
                double ke = 1000. * (p.genE - (MUON_MASS/1000.));
                if(ke >= params[0])
                    num_muons++;
            }
        }
        return num_muons == 1;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, single_muon_srtruth, single_muon_srtruth);

    /**
     * @brief Cut for at least one true final state proton above threshold.
     * @param obj the SRTrueInteraction to apply the cut on.
     * @param params KE threshold in MeV, defaults to 50 MeV.
     * @return true if at least one proton above threshold.
     */
    template<typename T>
    bool has_proton_srtruth(const T & obj, std::vector<double> params={50.0,})
    {
        for(const auto & p : obj.prim)
        {
            if(p.pdg == 2212)
            {
                double ke = 1000. * (p.genE - (PROTON_MASS/1000.));
                if(ke >= params[0])
                    return true;
            }
        }
        return false;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, has_proton_srtruth, has_proton_srtruth);

    /**
     * @brief Cut for zero true final state charged pions above threshold.
     * @param obj the SRTrueInteraction to apply the cut on.
     * @param params KE threshold in MeV, defaults to 25 MeV.
     * @return true if no charged pions above threshold.
     */
    template<typename T>
    bool no_charged_pions_srtruth(const T & obj, std::vector<double> params={25.0,})
    {
        for(const auto & p : obj.prim)
        {
            if(abs(p.pdg) == 211)
            {
                double ke = 1000. * (p.genE - (PION_MASS/1000.));
                if(ke >= params[0])
                    return false;
            }
        }
        return true;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, no_charged_pions_srtruth, no_charged_pions_srtruth);

    /**
     * @brief Cut for zero true final state neutral pions.
     * @param obj the SRTrueInteraction to apply the cut on.
     * @param params KE threshold in MeV, defaults to 0 MeV.
     * @return true if no neutral pions above threshold.
     */
    template<typename T>
    bool no_neutral_pions_srtruth(const T & obj, std::vector<double> params={0.0,})
    {
        for(const auto & p : obj.prim)
        {
            if(p.pdg == 111)
            {
                double ke = 1000. * p.genE - 134.977;
                if(ke >= params[0])
                    return false;
            }
        }
        return true;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, no_neutral_pions_srtruth, no_neutral_pions_srtruth);

    /**
     * @brief Cut for zero true final state photons above threshold.
     * @param obj the SRTrueInteraction to apply the cut on.
     * @param params energy threshold in MeV, defaults to 25 MeV.
     * @return true if no photons above threshold.
     */
    template<typename T>
    bool no_photons_srtruth(const T & obj, std::vector<double> params={25.0,})
    {
        for(const auto & p : obj.prim)
        {
            if(p.pdg == 22)
            {
                double energy = 1000. * p.genE;
                if(energy >= params[0])
                    return false;
            }
        }
        return true;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, no_photons_srtruth, no_photons_srtruth);

    /**
     * @brief Cut for zero true final state electrons above threshold.
     * @param obj the SRTrueInteraction to apply the cut on.
     * @param params KE threshold in MeV, defaults to 25 MeV.
     * @return true if no electrons above threshold.
     */
    template<typename T>
    bool no_electrons_srtruth(const T & obj, std::vector<double> params={25.0,})
    {
        for(const auto & p : obj.prim)
        {
            if(abs(p.pdg) == 11)
            {
                double ke = 1000. * (p.genE - (ELECTRON_MASS/1000.));
                if(ke >= params[0])
                    return false;
            }
        }
        return true;
    }
    REGISTER_CUT_SCOPE(RegistrationScope::MCTruth, no_electrons_srtruth, no_electrons_srtruth);

} // namespace mctruth
#endif