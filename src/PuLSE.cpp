/*
PuLSE version 1.2
Copyright(c) 2017 Steven Shave

Distributed under the MIT license

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files(the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and / or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions :

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "PuLSE-HTMLWriter.hpp"
#include "PuLSE-RunData.hpp"
#include <array>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// Version number
constexpr char PuLSEVersion[] = "1.2";

// 2 different usages, depending on platform (linux has the ability to read gz
// compressed input files)
#ifdef __linux__
constexpr char usage[] =
    "./pulse inFile.fastq(.gz) libraryDefinition[triplet "
    "residue][...]\nwhere...\ninFile.fastq(.gz) --> the file containing "
    "library sequence reads (can be gz.\nlibrary Definition --> denotes the "
    "sequence preceeding the randomized stretch of bases.\n[triplet residue] "
    "--> optional argument to change the standard codon triplet to AA residue "
    "mapping.  Useful if using a page display system which incorporates "
    "nonsense surpression\n";
#else
constexpr char usage[] =
    "./pulse inFile.fastq libraryDefinition[triplet "
    "residue][...]\nwhere...\ninFile.fastq --> the file containing library "
    "sequence reads.\nlibrary Definition --> denotes the sequence preceeding "
    "the randomized stretch of bases.\n[triplet residue] --> optional argument "
    "to change the standard codon triplet to AA residue mapping.  Useful if "
    "using a page display system which incorporates nonsense surpression\n";

#endif

// Display version number, usage, and then exit.
void ShowUsageAndExit() {
  std::cerr << "PuLSE v";
  std::cerr << PuLSEVersion;
  std::cerr << "\n";
  std::cerr << usage;
  exit(-1);
}

int main(int argc, char **argv) {
  // Check arguments
  if (argc < 3 || ((argc - 3) % 2 != 0.0f)) {
    if (argc == 2) {
      if (std::strcmp(argv[1], "-v") == 0 || std::strcmp(argv[1], "--v") == 0) {
        std::cerr << "PuLSE v";
        std::cerr << PuLSEVersion;
        std::cerr << "\n";
        exit(0);
      }
    }
    ShowUsageAndExit();
  }

  // Make a RunData object, then pass to a PuLSEHTMLWriter object
  auto rundata =
      std::make_shared<RunData>(std::string(argv[1]), std::string(argv[2]));
  PuLSEHTMLWriter htmlwriter(rundata);

  // Output some essentails to the termial
  std::cerr << "Data file:\t\t\t" << rundata->fastaFilename << "\n";
  std::cerr << "Library definition:\t\t" << rundata->libraryDefinition << "\n";
  std::cerr << "Upstream marker:\t" << rundata->sequenceBeginMarker << "\n";
  std::cerr << "Downstream marker:\t" << rundata->sequenceEndMarker << "\n";
  std::cerr << "Randomized DNA positions:\t" << rundata->dnalength << "\n";

  // In systems with nonsense surpression, we can supply trailing parameters
  // such as: "UAG Q".  Only if present on the command line are they added
  rundata->populateNonStandardTriplets(argc, argv);

  // Read DNA sequences.
  rundata->ReadSequences();

  // Write some simple parameters to the HTML report.
  htmlwriter.WriteRunInfo();

  // Output basic statistics
  rundata->PopulateBasicStats();
  htmlwriter.WriteBasicStats();

  // Output cumulative occurances
  rundata->PopulateCumulativeOccurances();
  htmlwriter.WriteCumulativeCounts();

  // Output common occurances
  rundata->PopulateCommonOccurances();
  htmlwriter.WriteCommonOccurances();

  // Output positional counts and heatmaps.
  rundata->PopulateHeatmaps();
  htmlwriter.WriteHeatMaps();

  return 0;
}
