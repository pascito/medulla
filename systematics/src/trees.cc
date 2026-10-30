/**
 * @file trees.cc
 * @brief Implementation file for the trees namespace.
 * @details This file contains the implementation of the functions that read
 * and interface with the TTrees produced by the CAFAna analysis framework.
 * Different "copying" actions can be performed on the TTrees, such as a simple
 * copy or adding systematics to the output file based on the selected signal
 * candidates and the configured systematics.
 * @author mueller@fnal.gov
 */
#include <iostream>
#include <stdexcept>
#include <cstdint>
#include <cmath>
#include "TMD5.h"

#include "trees.h"
#include "detsys.h"
#include "utilities.h"
#include "configuration.h"
#include "systematic.h"
#include "weight_reader.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"

// Copy the input TTree to the output TTree.
void sys::trees::copy_tree(cfg::ConfigurationTable & table, TFile * output, TFile * input)
{
    /**
     * @brief Create the output subdirectory following the nesting outlined
     * in the configuration file.
     */
    TDirectory * directory = (TDirectory *) output;
    directory = create_directory(directory, table.get_string_field("destination").c_str());
    directory->cd();

    /**
     * @brief Check if the exposure information ("POT", "Livetime") has
     * alread been copied and saved. If not, copy the exposure information
     * to the output TTree.
     */
    if(!directory->GetListOfKeys()->Contains("POT"))
    {
        std::cout << "Copying POT and Livetime histograms." << std::endl;
        TDirectory * parent = (TDirectory *) input;
        parent = get_parent_directory(parent, table.get_string_field("origin").c_str());
        TH1D * pot = (TH1D *) parent->Get("POT");
        TH1D * livetime = (TH1D *) parent->Get("Livetime");
        directory->WriteObject(pot, "POT");
        directory->WriteObject(livetime, "Livetime");
    }
    
    /**
     * @brief Create the output TTree with the name specified in the
     * configuration file.
     */
    TTree * output_tree = new TTree(table.get_string_field("name").c_str(), table.get_string_field("name").c_str());

    /**
     * @brief Connect to the input TTree and associated branches.
     * @details Three are N+3 branches in the input TTree, where N is the
     * number of branches in the input TTree of type double. The three
     * other branches (Run, Subrun, Evt) are of type int. This block uses
     * a single array to store the values of the double branches and three
     * separate variables to store the values of the int branches.
     */
    TTree * input_tree = (TTree *) input->Get(table.get_string_field("origin").c_str());
    int run, subrun, event;
    double br[input_tree->GetNbranches()-3];
    for (int i = 0; i < input_tree->GetNbranches()-3; i++)
        input_tree->SetBranchAddress(input_tree->GetListOfBranches()->At(i)->GetName(), br+i);
    input_tree->SetBranchAddress("Run", &run);
    input_tree->SetBranchAddress("Subrun", &subrun);
    input_tree->SetBranchAddress("Evt", &event);

    /**
     * @brief Create the branches in the output TTree following the same
     * structure as the input TTree.
     * @details The same array (type double) and three variables (type int)
     * as used above with the input TTree are used to create the branches
     * in the output TTree. The branches are created in the same order as
     * the input TTree. This process streamlines the copying of the values
     * from the input TTree to the output TTree.
     */
    for (int i = 0; i < input_tree->GetNbranches()-3; i++)
        output_tree->Branch(input_tree->GetListOfBranches()->At(i)->GetName(), br+i);
    output_tree->Branch("Run", &run);
    output_tree->Branch("Subrun", &subrun);
    output_tree->Branch("Evt", &event);

    /**
     * @brief Loop over the input TTree and copy the values to the output
     * TTree.
     */
    for(int i(0); i < input_tree->GetEntries(); ++i)
    {
        input_tree->GetEntry(i);
        output_tree->Fill();
    }

    /**
     * @brief Write the output TTree to the output ROOT file.
     */
    directory->WriteObject(output_tree, table.get_string_field("name").c_str());
}

// Add reweightable systematics to the output TTree.
void sys::trees::copy_with_weight_systematics(cfg::ConfigurationTable & config, cfg::ConfigurationTable & table, TFile * output, TFile * input, sys::detsys::DetsysCalculator & calc)
{
    /**
     * @brief Create the output subdirectory following the nesting outlined
     * in the configuration file.
     */
    TDirectory * directory = (TDirectory *) output;
    directory = create_directory(directory, table.get_string_field("destination").c_str());
    directory->cd();
    
    /**
     * @brief Check if the exposure information ("POT", "Livetime") has
     * alread been copied and saved. If not, copy the exposure information
     * to the output TTree.
     */
    if(!directory->GetListOfKeys()->Contains("POT"))
    {
        std::cout << "Copying POT and Livetime histograms." << std::endl;
        TDirectory * parent = (TDirectory *) input;
        parent = get_parent_directory(parent, table.get_string_field("origin").c_str());
        TH1D * pot = (TH1D *) parent->Get("POT");
        TH1D * livetime = (TH1D *) parent->Get("Livetime");
        directory->WriteObject(pot, "POT");
        directory->WriteObject(livetime, "Livetime");
    }

    /**
     * @brief Connect to the input TTree and associated branches.
     * @details Three are N+3 branches in the input TTree, where N is the
     * number of branches in the input TTree of type double. The three
     * other branches (Run, Subrun, Evt) are of type int. This block uses
     * a single array to store the values of the double branches and three
     * separate variables to store the values of the int branches. There is
     * one quirk, however, as we would also like to have access to the
     * "true_neutrino_id" branch in the input TTree directly.
     */
    TTree * input_tree = (TTree *) input->Get(table.get_string_field("origin").c_str());
    std::map<std::string, double> brs;
    double & nu_id = brs["true_neutrino_id"];
    Int_t run, subrun, event;
    for(int i(0); i < input_tree->GetNbranches()-3; ++i)
    {
        std::string brname = input_tree->GetListOfBranches()->At(i)->GetName();
        brs[brname] = 0;
        input_tree->SetBranchAddress(brname.c_str(), &brs[brname]);
    }
    input_tree->SetBranchAddress("Run", &run);
    input_tree->SetBranchAddress("Subrun", &subrun);
    input_tree->SetBranchAddress("Evt", &event);


    /**
     * @brief Create the output TTree with the name specified in the
     * configuration file.
     * @details The output TTree is created with the same branches as the
     * input TTree, plus the Run, Subrun, and Evt branches. The same array
     * (type double) and three variables (type int) as used above with the
     * input TTree are used to create the branches in the output TTree. The
     * branches are created in the same order as the input TTree. This
     * process streamlines the copying of the values from the input TTree to
     * the output TTree.
     */
    TTree * output_tree = new TTree(table.get_string_field("name").c_str(), table.get_string_field("name").c_str());
    for(auto & br : brs)
        output_tree->Branch(br.first.c_str(), &br.second);
    output_tree->Branch("Run", &run);
    output_tree->Branch("Subrun", &subrun);
    output_tree->Branch("Evt", &event);

    /**
     * @brief Create an output TTree for non-matched signal candidates.
     * @details This differs from the output TTree above in that it is meant to
     * capture signal candidates that do not have a matching neutrino in the
     * input CAF files (e.g., cosmics or failed truth matching). This TTree is
     * only created if the "output.save_nonmatched" field in the configuration
     * file is set to true.
     */
    TTree * nonmatched_tree = nullptr;
    if(config.get_bool_field("output.nonmatched", true))
    {
        nonmatched_tree = new TTree((table.get_string_field("name") + "_nonmatched").c_str(), (table.get_string_field("name") + "_nonmatched").c_str());
        for(auto & br : brs)
            nonmatched_tree->Branch(br.first.c_str(), &br.second);
        nonmatched_tree->Branch("Run", &run);
        nonmatched_tree->Branch("Subrun", &subrun);
        nonmatched_tree->Branch("Evt", &event);
    }
    
    /**
     * @brief Create the map of selected signal candidates.
     * @details This block creates a map of selected signal candidates. The
     * map is built by looping over the input TTree and storing an index of the
     * run, subrun, event, nu_id, and nu_energy branches as the key. The
     * value is the index of the entry in the input TTree.
     */
    std::map<index_t, std::vector<Long64_t>> candidates;
    bool use_additional_hash = config.get_bool_field("input.use_additional_hash", false);
    const bool use_vertex_matching = config.get_bool_field("input.use_vertex_matching", false);
    const double vertex_tolerance = config.has_field("input.vertex_tolerance_cm")
        ? config.get_double_field("input.vertex_tolerance_cm") : 0.001;
    const char * vertex_branches[] = {"true_vertex_x", "true_vertex_y", "true_vertex_z"};
    if(use_vertex_matching)
    {
        if(!std::isfinite(vertex_tolerance) || vertex_tolerance <= 0)
            throw std::runtime_error("input.vertex_tolerance_cm must be finite and positive");
        for(const char * branch : vertex_branches)
            if(!input_tree->GetBranch(branch))
                throw std::runtime_error(std::string("Vertex matching requires input branch ") + branch);
        std::cout << "[MATCHING] vertex tolerance=" << vertex_tolerance << " cm" << std::endl;
    }
    for(Long64_t i = 0; i < input_tree->GetEntries(); ++i)
    {
        if(input_tree->GetEntry(i) <= 0)
            throw std::runtime_error("Failed to read input candidate " + std::to_string(i));
        // Preserve negative sentinel IDs through a defined signed-to-unsigned cast.
        index_t key = std::make_tuple(run, subrun, event,
            static_cast<uint64_t>(static_cast<int64_t>(nu_id)),
            use_additional_hash && nu_id >= 0 ? brs["true_neutrino_energy"] : 0.0);
        candidates[key].push_back(i);
    }

    /**
     * @brief Configure the weight-based systematics.
     * @details This block configures the weight-based systematics. The
     * systematics are split (by type) into separate TTrees, which is
     * enforced by the "type" field in the configuration block for each
     * systematic. Because we do not wish to loop over the selected signal
     * candidates multiple times, we must store the systematic information
     * in such a way that we can easily accomodate this scheme. The variable
     * "systematics" is a map of Systematic objects keyed by the name of the
     * systematic parameter. Each Systematic object contains metadata about
     * the systematic parameter (name, index, type, etc.), some
     * configuration information, and a pointer to the output TTree,
     * weights vector, and zscores vector.
     */
    std::map<std::string, Systematic *> systematics;
    std::map<std::string, TTree *> systrees;

    /**
     * @brief Create histograms for storing the systematic results as a 
     * function of a collection of variables.
     * @details This block creates histograms for storing the systematic
     * weights / selected ratios as a function of a the variables specified
     * in the configuration file. The histograms are stored in a map with
     * the key being a pair of the variable name and the systematic name.
     * The 1D histogram contains a single entry per universe with a fill
     * value corresponding to the ratio of the selected signal candidates
     * with the universe weight to the nominal count. The 2D histogram
     * contains a 2D histogram with the variable on the x-axis and the
     * universe index on the y-axis. The fill value is the universe weight.
     * The 1D histograms can be easily inspected to see the one-bin effect
     * (uncertainty) of the systematic on the selected signal candidates.
     * The 2D histograms contain similar information, but can additionally
     * be used to inspect the effect of the systematic as a function of the
     * variable or calculate a covariance matrix.
     */
    std::vector<SysVariable> sysvariables;
    std::map<syst_t, TH2D *> results2d;
    std::map<syst_t, TH1D *> results1d;
    for(cfg::ConfigurationTable & t : config.get_subtables("sysvar"))
    {
        sysvariables.push_back(SysVariable(t));
        calc.add_variable(sysvariables.back());
    }

    /**
     * @brief Loop over the systematic types in the configuration file.
     * @details This block loops over the systematic types in the
     * configuration file. The "table" field in the configuration file
     * specifies the name of the table lists in the configuration file that
     * contain the exact definition of the systematics. Principally, this
     * loop is used to load and configure the systematics of each type in
     * sequential order.
     */
    std::vector<std::string> table_types = table.get_string_vector("table_types");
    for(const std::string & s : table_types)
    {
        std::string tname = table.get_string_field("name") + '_' + s;
        systrees[tname] = new TTree(
            (tname + "Tree").c_str(),
            (tname + "Tree").c_str());
        systrees[tname]->SetDirectory(nullptr);
        systrees[tname]->Branch("Run", &run);
        systrees[tname]->Branch("Subrun", &subrun);
        systrees[tname]->Branch("Evt", &event);
        systrees[tname]->SetDirectory(directory);
        systrees[tname]->SetAutoFlush(1000);
    }

    for(cfg::ConfigurationTable & t : config.get_subtables("sys"))
    {
        std::string tname = table.get_string_field("name") + '_' + t.get_string_field("type");
        systematics.insert(std::make_pair<std::string, Systematic *>(t.get_string_field("name"), new Systematic(t, systrees[tname])));
        Systematic * tmp = systematics[t.get_string_field("name")];
        tmp->get_tree()->Branch(t.get_string_field("name").c_str(), &systematics[t.get_string_field("name")]->get_weights());
        if(tmp->get_nsigma()->size() > 0)
        {
            tmp->get_tree()->Branch((t.get_string_field("name") + "_sigma").c_str(), &systematics[t.get_string_field("name")]->get_nsigma());
        }
    }

    sys::WeightReader reader(config.get_string_field("input.weights"));
    // Keep a digest rather than every universe vector for every matched key.
    // This is a consistency fingerprint, not a security boundary.
    // Track each candidate separately: one coarse event/energy key can
    // contain distinct interactions with different vertices and weights.
    std::map<Long64_t, std::pair<std::string, Long64_t>> matched_signatures;
    size_t duplicate_records = 0;
    Long64_t matched_rows = 0, nonmatched_rows = 0;
    double nominal_count(0);
    while(reader.next())
    {
        for(size_t idn = 0; idn < reader.get_nnu(); ++idn)
        {
            index_t index = std::make_tuple(reader.get_run(), reader.get_subrun(),
                reader.get_event(), idn,
                use_additional_hash ? static_cast<double>(reader.get_energy(idn)) : 0.0);
            auto candidate = candidates.find(index);
            if(candidate == candidates.end()) continue;

            std::vector<Long64_t> matching_entries;
            if(use_vertex_matching)
            {
                const auto position = reader.get_vertex(idn);
                for(double coordinate : position)
                    if(!std::isfinite(coordinate))
                        throw std::runtime_error("Nonfinite CAF vertex in vertex matching");
                for(Long64_t input_entry : candidate->second)
                {
                    if(input_tree->GetEntry(input_entry) <= 0)
                        throw std::runtime_error("Failed to read candidate for vertex matching");
                    bool matches = true;
                    for(size_t axis = 0; axis < position.size(); ++axis)
                    {
                        const double coordinate = brs.at(vertex_branches[axis]);
                        if(!std::isfinite(coordinate))
                            throw std::runtime_error("Nonfinite candidate vertex in vertex matching");
                        if(std::abs(coordinate - position[axis]) > vertex_tolerance)
                            matches = false;
                    }
                    if(matches) matching_entries.push_back(input_entry);
                }
            }
            else matching_entries = candidate->second;
            if(matching_entries.empty()) continue;

            TMD5 digest;
            for(auto & [key, value] : systematics)
            {
                value->get_weights()->clear();
                if(value->get_type() == Type::kVARIATION) continue;
                reader.set(value->get_index());
                const uint64_t group = value->get_index();
                const uint64_t count = reader.get_nuniv(idn);
                if(count == 0)
                    throw std::runtime_error("Empty weight vector for " + key);
                if(value->get_type() == Type::kMULTISIGMA &&
                   count != value->get_nsigma()->size())
                    throw std::runtime_error("Weight/sigma length mismatch for " + key +
                        ": weights=" + std::to_string(count) +
                        ", sigma=" + std::to_string(value->get_nsigma()->size()));
                digest.Update(reinterpret_cast<const UChar_t*>(&group), sizeof(group));
                digest.Update(reinterpret_cast<const UChar_t*>(&count), sizeof(count));
                for(size_t u = 0; u < count; ++u)
                {
                    const float weight = reader.get_weight(idn, u);
                    value->get_weights()->push_back(weight);
                    digest.Update(reinterpret_cast<const UChar_t*>(&weight), sizeof(weight));
                }
            }
            digest.Final();
            const std::string signature = digest.AsString();
            bool duplicate_record = false;
            // Preserve every candidate that matches this interaction's vertex.
            for(Long64_t input_entry : matching_entries)
            {
                auto previous = matched_signatures.find(input_entry);
                if(previous != matched_signatures.end())
                {
                    if(previous->second.first != signature)
                        throw std::runtime_error("Conflicting duplicate weights for run=" +
                            std::to_string(reader.get_run()) + " subrun=" +
                            std::to_string(reader.get_subrun()) + " event=" +
                            std::to_string(reader.get_event()) + " neutrino=" + std::to_string(idn) +
                            " candidate=" + std::to_string(input_entry) +
                            " first_CAF_entry=" + std::to_string(previous->second.second) +
                            " second_CAF_entry=" + std::to_string(reader.get_entry()));
                    duplicate_record = true;
                    continue;
                }
                matched_signatures.emplace(input_entry, std::make_pair(signature, reader.get_entry()));
                if(input_tree->GetEntry(input_entry) <= 0)
                    throw std::runtime_error("Failed to reread matched input candidate");
                if(calc.is_initialized()) calc.increment_nominal_count(1.0);
                nominal_count += 1.0;
                for(auto & [key, value] : systematics)
                {
                    if(value->get_type() != Type::kVARIATION)
                    {
                        const auto & weights = *value->get_weights();
                        for(SysVariable & sv : sysvariables)
                        {
                            syst_t syskey = std::make_pair(sv.name, value->get_index());
                            if(results1d.find(syskey) == results1d.end())
                            {
                                results1d[syskey] = new TH1D((sv.name + "_" + key + "_1d").c_str(), (sv.name + "_" + key + "_1d").c_str(), 1000, -0.25, 0.25);
                                results1d[syskey]->SetDirectory(nullptr);
                                results2d[syskey] = new TH2D((sv.name + "_" + key + "_2d").c_str(), (sv.name + "_" + key + "_2d").c_str(), sv.nbins, sv.min, sv.max, weights.size(), 0, weights.size());
                                results2d[syskey]->SetDirectory(nullptr);
                            }
                            for(size_t u = 0; u < weights.size(); ++u)
                                results2d[syskey]->Fill(brs[sv.name], u, weights[u]);
                        }
                    }
                    else
                    {
                        value->get_weights()->clear();
                        for(double z : calc.get_zscores(key))
                            value->get_weights()->push_back(calc.get_weight(key, brs[calc.get_variable()], z));
                        for(SysVariable & sv : sysvariables)
                            calc.add_value(sv.name, brs[sv.name], key, brs[calc.get_variable()]);
                    }
                }
                output_tree->Fill();
                for(auto & [key, value] : systrees) value->Fill();
                ++matched_rows;
            }
            if(duplicate_record) ++duplicate_records;
        }
    }
    if(nonmatched_tree)
    {
        for(auto & [key, entries] : candidates)
        {
            for(Long64_t input_entry : entries)
            {
                if(matched_signatures.count(input_entry)) continue;
                if(input_tree->GetEntry(input_entry) <= 0)
                    throw std::runtime_error("Failed to reread nonmatched input candidate");
                nonmatched_tree->Fill();
                ++nonmatched_rows;
            }
        }
        if(matched_rows + nonmatched_rows != input_tree->GetEntries())
            throw std::runtime_error("Candidate conservation check failed");
        directory->WriteObject(nonmatched_tree, nonmatched_tree->GetName());
        delete nonmatched_tree;
    }
    std::cout << "[MATCHING] input=" << input_tree->GetEntries()
              << " matched=" << matched_rows << " nonmatched=" << nonmatched_rows
              << " identical_duplicate_weight_records=" << duplicate_records << std::endl;

    // Write the output TTree to the output file.
    directory->WriteObject(output_tree, table.get_string_field("name").c_str());
    for(auto & [key, value] : systrees)
        directory->WriteObject(value, (key+"Tree").c_str());
    
    // Write the systematic histograms to the output file.
    std::string destination = config.get_string_field("output.histogram_destination", "");
    if(destination != "")
    {
        TDirectory * histogram_directory = create_directory(output, destination.c_str());
        for(auto & [key, value] : results2d)
        {
            std::string name = value->GetName();
            histogram_directory->WriteObject(value, name.c_str());
            for(int i(0); i < value->GetNbinsY(); ++i)
            {
                double sum(0);
                for(int j(0); j < value->GetNbinsX(); ++j)
                    sum += value->GetBinContent(j+1, i+1);
                results1d[key]->Fill((sum - nominal_count) / nominal_count);
            }
            delete value;
        }
        for(auto & [key, value] : results1d)
        {
            std::string name = value->GetName();
            histogram_directory->WriteObject(value, name.c_str());
            delete value;
        }
    }

    // Write detector systematic histograms to the output file.
    if(calc.is_initialized())
        calc.write_results();
}
