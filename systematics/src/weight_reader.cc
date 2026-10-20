/**
 * @file weight_reader.cc
 * @brief A simple weight reader for CAF / flat CAF files.
 *
 * @details
 * This class reads weight information from ROOT TTrees and provides
 * accessor methods for retrieving metadata and systematic weights.
 *
 * Debugging additions:
 *   1. ROOT file registration diagnostics
 *   2. TChain entry count verification
 *   3. Flat/structured CAF identification
 *   4. TTreeReader initialization diagnostics
 *   5. Protection against division by zero in progress bar
 *
 * @author mueller@fnal.gov
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <sstream>
#include <stdexcept>
#include <string>

#include "weight_reader.h"

#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include "sbnanaobj/StandardRecord/SRTrueInteraction.h"


// ============================================================
// 1. Constructor
// ============================================================

sys::WeightReader::WeightReader(const std::string & input)
    : chain("recTree"),
      entry(0),
      idx(0),
      progress_started(false)
{
    std::cout << "\n========================================\n";
    std::cout << "[WeightReader] Initializing\n";
    std::cout << "Input: " << input << std::endl;
    std::cout << "========================================\n";


    // --------------------------------------------------------
    // Case 1: Wildcard input
    // --------------------------------------------------------

    if(input.find("*") != std::string::npos)
    {
        std::cout << "[DEBUG] Input type: wildcard ROOT files"
                  << std::endl;

        isflat = input.find("flat") != std::string::npos;

        int added = chain.Add(input.c_str());

        std::cout << "[DEBUG] Files added: "
                  << added << std::endl;
    }


    // --------------------------------------------------------
    // Case 2: Text file containing ROOT file paths
    // --------------------------------------------------------

    else if(input.find(".txt") != std::string::npos)
    {
        std::cout << "[DEBUG] Input type: ROOT file list"
                  << std::endl;

        std::ifstream infile(input);

        if(!infile.is_open())
        {
            std::cerr
                << "[ERROR] Cannot open weights list: "
                << input
                << std::endl;

            throw std::runtime_error(
                "WeightReader: Cannot open weights list"
            );
        }

        std::string line;
        std::string last_file;

        int file_count = 0;
        int total_added = 0;


        // ----------------------------------------------------
        // Read each ROOT file from the list
        // ----------------------------------------------------

        while(std::getline(infile, line))
        {
            // Remove possible carriage return
            if(!line.empty() && line.back() == '\r')
                line.pop_back();

            // Skip empty or whitespace-only lines
            if(line.find_first_not_of(" \t") == std::string::npos)
                continue;

            std::cout
                << "[DEBUG] Adding file "
                << file_count
                << ": "
                << line
                << std::endl;

            int added = chain.Add(line.c_str());

            std::cout
                << "[DEBUG] chain.Add returned: "
                << added
                << std::endl;

            if(added <= 0)
            {
                std::cerr
                    << "[WARNING] Failed to register file: "
                    << line
                    << std::endl;
            }

            total_added += added;

            last_file = line;

            file_count++;
        }


        // ----------------------------------------------------
        // Check for empty file list
        // ----------------------------------------------------

        if(file_count == 0)
        {
            throw std::runtime_error(
                "WeightReader: Input weights list is empty"
            );
        }


        // ----------------------------------------------------
        // Determine whether files are flat CAFs
        // ----------------------------------------------------

        isflat = last_file.find("flat") != std::string::npos;


        // ----------------------------------------------------
        // Print ROOT file diagnostics
        // ----------------------------------------------------

        std::cout
            << "\n========== WEIGHT READER DEBUG =========="
            << std::endl;

        std::cout
            << "Input file list: "
            << input
            << std::endl;

        std::cout
            << "Files in list: "
            << file_count
            << std::endl;

        std::cout
            << "Total files registered: "
            << total_added
            << std::endl;

        std::cout
            << "TChain number of trees: "
            << chain.GetNtrees()
            << std::endl;

        std::cout
            << "[DEBUG] Calculating TChain entries..."
            << std::endl;

        Long64_t total_entries = chain.GetEntries();

        std::cout
            << "TChain total entries: "
            << total_entries
            << std::endl;

        std::cout
            << "Is flat CAF: "
            << isflat
            << std::endl;

        if(total_entries == 0)
        {
            std::cerr
                << "[ERROR] TChain contains zero entries!"
                << std::endl;
        }

        std::cout
            << "=========================================\n"
            << std::endl;
    }


    // --------------------------------------------------------
    // Case 3: Single ROOT file
    // --------------------------------------------------------

    else
    {
        std::cout
            << "[DEBUG] Input type: Single ROOT file"
            << std::endl;

        isflat = input.find("flat") != std::string::npos;

        int added = chain.Add(input.c_str());

        std::cout
            << "[DEBUG] chain.Add returned: "
            << added
            << std::endl;
    }


    // ========================================================
    // 2. Create TTreeReader
    // ========================================================

    std::cout
        << "[DEBUG] Creating TTreeReader..."
        << std::endl;

    reader = std::make_unique<TTreeReader>(&chain);


    // --------------------------------------------------------
    // Metadata branches
    // --------------------------------------------------------

    run = std::make_unique<TTreeReaderValue<uint32_t>>(
        *reader,
        "rec.hdr.run"
    );

    subrun = std::make_unique<TTreeReaderValue<uint32_t>>(
        *reader,
        "rec.hdr.subrun"
    );

    event = std::make_unique<TTreeReaderValue<uint32_t>>(
        *reader,
        "rec.hdr.evt"
    );


    // ========================================================
    // 3. Flat CAF configuration
    // ========================================================

    if(isflat)
    {
        std::cout
            << "[DEBUG] Configuring flat CAF branches..."
            << std::endl;


        // ----------------------------------------------------
        // Event-level indexing
        // ----------------------------------------------------

        chain.SetBranchAddress(
            "rec.mc.nu..length",
            &nnu
        );


        // ----------------------------------------------------
        // Neutrino-level indexing
        // ----------------------------------------------------

        chain.SetBranchAddress(
            "rec.mc.nu.wgt..length",
            nwgt
        );

        chain.SetBranchAddress(
            "rec.mc.nu.wgt..idx",
            &iwgt
        );

        chain.SetBranchAddress(
            "rec.mc.nu.E",
            &nu_energy
        );


        // ----------------------------------------------------
        // Systematic-level indexing
        // ----------------------------------------------------

        chain.SetBranchAddress(
            "rec.mc.nu.wgt.univ..length",
            &nuniv
        );

        chain.SetBranchAddress(
            "rec.mc.nu.wgt.univ..idx",
            &iuniv
        );

        chain.SetBranchAddress(
            "rec.mc.nu.wgt.univ",
            &wgts
        );


        // Load first entry
        std::cout
            << "[DEBUG] Loading first flat CAF entry..."
            << std::endl;

        Long64_t bytes = chain.GetEntry(0);

        std::cout
            << "[DEBUG] chain.GetEntry(0) returned: "
            << bytes
            << std::endl;
    }


    // ========================================================
    // 4. Structured CAF configuration
    // ========================================================

    else
    {
        std::cout
            << "[DEBUG] Configuring structured CAF branches..."
            << std::endl;


        // Corrected ROOT branch type:
        // unsigned long long instead of uint64_t

        nnu_structured =
            std::make_unique<
                TTreeReaderValue<unsigned long long>
            >(
                *reader,
                "rec.mc.nnu"
            );


        mc =
            std::make_unique<
                TTreeReaderArray<caf::SRTrueInteraction>
            >(
                *reader,
                "rec.mc.nu"
            );


        nu_energy_structured =
            std::make_unique<
                TTreeReaderArray<Float_t>
            >(
                *reader,
                "rec.mc.nu.E"
            );
    }


    // ========================================================
    // 5. Initialize reader
    // ========================================================

    std::cout
        << "\n[DEBUG] Before reader->Next()"
        << std::endl;

    bool success = reader->Next();

    std::cout
        << "[DEBUG] reader->Next() returned: "
        << success
        << std::endl;

    if(!success)
    {
        std::cerr
            << "[WARNING] TTreeReader failed to initialize "
            << "or read the first entry."
            << std::endl;
    }

    std::cout
        << "[DEBUG] WeightReader constructor finished."
        << std::endl;
}


// ============================================================
// 6. Advance to next entry
// ============================================================

bool sys::WeightReader::next()
{
    Long64_t total_entries = chain.GetEntries();

    this->progress_bar(
        entry + 1,
        total_entries
    );


    // Check reader and tree
    if(!chain.GetTree() || !reader)
        return false;


    // Check entry range
    if(entry >= static_cast<size_t>(total_entries))
        return false;


    // Advance reader
    if(!reader->Next())
        return false;


    // Load next entry
    chain.GetEntry(++entry);

    return true;
}


// ============================================================
// 7. Set weight group index
// ============================================================

void sys::WeightReader::set(size_t index)
{
    idx = index;
}


// ============================================================
// 8. Get number of neutrinos
// ============================================================

uint32_t sys::WeightReader::get_nnu() const
{
    return isflat ? nnu : **nnu_structured;
}


// ============================================================
// 9. Get number of weight groups
// ============================================================

uint32_t sys::WeightReader::get_nwgt(Int_t i) const
{
    if(i < 0 || i >= static_cast<Int_t>(this->get_nnu()))
    {
        throw std::out_of_range(
            "WeightReader: Index out of range in get_nwgt()"
        );
    }

    return isflat
        ? nwgt[i]
        : (*mc)[i].wgt.size();
}


// ============================================================
// 10. Get number of universes
// ============================================================

uint32_t sys::WeightReader::get_nuniv(size_t idn) const
{
    if(idn >= get_nnu() || idx >= get_nwgt(idn))
    {
        throw std::out_of_range(
            "WeightReader: Index out of range in get_nuniv()"
        );
    }

    return isflat
        ? nuniv[iwgt[idn] + idx]
        : (*mc)[idn].wgt[idx].univ.size();
}


// ============================================================
// 11. Get systematic weight
// ============================================================

float sys::WeightReader::get_weight(
    size_t idn,
    size_t idu
) const
{
    if(isflat)
    {
        size_t n = iwgt[idn] + idx;

        size_t univ_offset = iuniv[n];

        return wgts[univ_offset + idu];
    }

    return (*mc)[idn].wgt[idx].univ[idu];
}


// ============================================================
// 12. Get neutrino energy
// ============================================================

float sys::WeightReader::get_energy(size_t idn) const
{
    return isflat
        ? nu_energy[idn]
        : (*nu_energy_structured)[idn];
}


// ============================================================
// 13. Progress bar
// ============================================================

void sys::WeightReader::progress_bar(
    size_t entry,
    size_t total
) const
{
    // --------------------------------------------------------
    // Prevent division by zero
    // --------------------------------------------------------

    if(total == 0)
    {
        std::cerr
            << "\n[ERROR] WeightReader: "
            << "TChain contains zero entries!"
            << std::endl;

        return;
    }


    // --------------------------------------------------------
    // Start timer
    // --------------------------------------------------------

    if(!progress_started)
    {
        progress_start_time =
            std::chrono::steady_clock::now();

        progress_started = true;
    }


    // --------------------------------------------------------
    // Calculate progress
    // --------------------------------------------------------

    float percent =
        static_cast<float>(entry) / total;

    int percent_int =
        static_cast<int>(percent * 1000.0);


    if(percent_int == last_printed_percent &&
       entry != total)
    {
        return;
    }

    last_printed_percent = percent_int;


    // --------------------------------------------------------
    // Print progress bar
    // --------------------------------------------------------

    std::cout << "\r\033[K[";

    int pos =
        static_cast<int>(50 * percent);

    for(int i = 0; i < 50; ++i)
    {
        if(i < pos)
            std::cout << "=";

        else if(i == pos)
            std::cout << ">";

        else
            std::cout << " ";
    }

    std::cout
        << "] "
        << std::fixed
        << std::setprecision(2)
        << percent * 100.0
        << "%  ";


    // --------------------------------------------------------
    // Calculate elapsed time and ETA
    // --------------------------------------------------------

    auto now =
        std::chrono::steady_clock::now();

    double elapsed =
        std::chrono::duration_cast<
            std::chrono::duration<double>
        >(
            now - progress_start_time
        ).count();


    double eta =
        percent > 0.0
        ? elapsed / percent - elapsed
        : 0.0;


    // --------------------------------------------------------
    // Format time
    // --------------------------------------------------------

    auto format_time =
        [](double seconds) -> std::string
    {
        int h =
            static_cast<int>(seconds) / 3600;

        int m =
            (static_cast<int>(seconds) % 3600) / 60;

        double s =
            seconds - h * 3600 - m * 60;

        std::ostringstream oss;

        oss << std::setfill('0')
            << std::setw(2) << h << ":"
            << std::setw(2) << m << ":"
            << std::setw(4)
            << std::fixed
            << std::setprecision(1)
            << s;

        return oss.str();
    };


    // --------------------------------------------------------
    // Print time information
    // --------------------------------------------------------

    std::cout
        << "Elapsed: "
        << format_time(elapsed)
        << ", ETA: "
        << format_time(eta)
        << std::flush;


    // --------------------------------------------------------
    // Finish progress bar
    // --------------------------------------------------------

    if(entry == total)
    {
        std::cout << std::endl;

        progress_started = false;
    }
}