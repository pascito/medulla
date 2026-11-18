/**
 * @file main.cc
 * @brief Main function for the code that adds TTrees with systematics to the
 * output ROOT file from CAFAna analyses.
 * @details This code is designed to add TTrees with systematics to the output
 * ROOT file. The code reads in the input ROOT file, which contains the TTrees
 * produced by the CAFAna analysis framework with the sBruce tree format. The 
 * code then matches the selected interactions in the sBruce trees with the 
 * corresponding universe weights from the original CAF files. The final output
 * is a ROOT file with TTrees that contain the selected interactions and some
 * additional TTrees containing the universe weights for the configured
 * systematics.
 * @author mueller@fnal.gov
 */
#include <iostream>
#include <string>
#include <memory>
#include <filesystem>
#include <stdexcept>

#include "configuration.h"
#include "trees.h"
#include "detsys.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

namespace
{
    // Check required selection objects before creating or replacing output.
    void validate_input(cfg::ConfigurationTable & config, TFile * input)
    {
        if(!input || input->IsZombie())
            throw std::runtime_error("Cannot open selection input ROOT file: " +
                config.get_string_field("input.path"));
        if(input->TestBit(TFile::kRecovered))
            throw std::runtime_error("Selection input ROOT file required recovery: " +
                config.get_string_field("input.path"));
        if(config.has_field("variations"))
        {
            const std::string origin = config.get_string_field("variations.origin");
            const std::string tree = config.get_string_field("variations.tree");
            const std::string variable = config.get_string_field("variations.variable");
            for(const std::string & key : config.get_string_vector("variations.keys"))
            {
                const std::string name = origin + key + "/" + tree;
                auto * selected = input->Get<TTree>(name.c_str());
                if(!selected)
                    throw std::runtime_error("Missing detector-variation TTree: " + name +
                        " in " + input->GetName() + ". input.path must point to the merged selection file.");
                if(!selected->GetBranch(variable.c_str()))
                    throw std::runtime_error("Missing detector-variation branch " + variable + " in " + name);
                const std::string exposure_name = name + "_exposure";
                auto * exposure_object = input->Get(exposure_name.c_str());
                if(exposure_object)
                {
                    auto * exposure = dynamic_cast<TTree *>(exposure_object);
                    if(!exposure || !exposure->GetBranch("pot"))
                        throw std::runtime_error("Invalid POT exposure TTree: " + exposure_name);
                }
                else
                {
                    const std::string pot_name = origin + key + "/POT";
                    auto * pot = dynamic_cast<TH1 *>(input->Get(pot_name.c_str()));
                    if(!pot || pot->GetDimension() != 1 || pot->GetNbinsX() < 1)
                        throw std::runtime_error("Missing POT exposure TTree or histogram for " + name);
                }
            }
        }
        for(auto & table : config.get_subtables("tree"))
        {
            if(table.get_string_field("action") != "add_weights") continue;
            const std::string origin = table.get_string_field("origin");
            if(!input->Get<TTree>(origin.c_str()))
                throw std::runtime_error("Missing selection TTree for weights: " + origin);
        }
    }

    bool same_file(const std::string & input, const std::string & output)
    {
        if(input == output) return true;
        if(input.find("://") != std::string::npos || output.find("://") != std::string::npos)
            return false;
        std::error_code error;
        if(std::filesystem::equivalent(input, output, error) && !error) return true;
        return std::filesystem::weakly_canonical(input) == std::filesystem::weakly_canonical(output);
    }
}

int main(int argc, char * argv[])
{
    gErrorIgnoreLevel = kError;
    const bool check_only = argc == 3 && std::string(argv[1]) == "--check-input";
    if(argc != 2 && !check_only)
    {
        std::cerr << "Usage: " << argv[0] << " [--check-input] <configuration.toml>" << std::endl;
        return 1;
    }
    try
    {
        cfg::ConfigurationTable config;
        config.set_config(argv[argc - 1]);
        const std::string input_path = config.get_string_field("input.path");
        const std::string output_path = config.get_string_field("output.path");
        if(same_file(input_path, output_path))
            throw std::runtime_error("Input and output ROOT paths refer to the same file; refusing to overwrite selection input.");
        std::unique_ptr<TFile> input(TFile::Open(input_path.c_str(), "READ"));
        validate_input(config, input.get());
        std::cout << "[INPUT] Selection input validated: " << input_path << std::endl;
        if(check_only)
        {
            std::cout << "[INPUT] Check complete; no output file created or modified." << std::endl;
            return 0;
        }
        std::unique_ptr<TFile> output(TFile::Open(output_path.c_str(), "RECREATE"));
        if(!output || output->IsZombie() || !output->IsWritable())
            throw std::runtime_error("Cannot create output ROOT file: " + output_path);
        sys::detsys::DetsysCalculator calc;
        if(config.has_field("variations"))
        {
            calc = sys::detsys::DetsysCalculator(config, output.get(), input.get());
            calc.write();
        }
        for(auto & table : config.get_subtables("tree"))
        {
            const std::string origin = table.get_string_field("origin");
            std::cout << "Processing tree: " << origin << std::endl;
            if(!input->Get(origin.c_str()))
            {
                std::cerr << "Info: Tree " << origin << " not found in input file. Skipping..." << std::endl;
                continue;
            }
            const std::string type = table.get_string_field("action");
            if(type == "copy")
                sys::trees::copy_tree(table, output.get(), input.get());
            else if(type == "add_weights")
                sys::trees::copy_with_weight_systematics(config, table, output.get(), input.get(), calc);
        }
        output->Close();
        input->Close();
    }
    catch(const std::exception & error)
    {
        std::cerr << "[FATAL] " << error.what() << std::endl;
        return 2;
    }
    return 0;
}
