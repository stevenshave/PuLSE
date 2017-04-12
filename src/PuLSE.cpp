// Steven R Shave (stevenshave@gmail.com) - 16/10/2016

#include "PuLSE-RunData.hpp"
#include "PuLSE-HTMLWriter.hpp"
#include <array>
#include <cmath>
#include <fstream>
#include <utility>
#include <iostream>

#include <memory>
#include <set>
#include <cstdlib>
#include <string>
#include <tuple>
#include <vector>



void ShowUsageAndExit() {
  std::cerr << "./profilePhageLibrary inFile.fastq libraryDefinition"
               "[triplet residue][...]\n";
  std::cerr << "where...\n";
  std::cerr << "inFile.fastq --> the file containing library sequence reads.\n";
  std::cerr << "library Definition --> denotes the sequence preceeding the randomized stretch of bases.\n";

  std::cerr << "[triplet residue] --> optional argument to change the standard codon triplet to AA residue mappin\n";
  std::cerr << "                      useful if using a page display system which incorporates nonsense surpression\n";
  exit(-1);
}

int main(int argc, char **argv) {

  if (argc < 3 || ((argc - 3) % 2 != 0.0f)) {
    ShowUsageAndExit();
  }

  // Make a RunData object, then pass to a PuLSEHTMLWriter object
  auto rundata =
      std::make_shared<RunData>(std::string(argv[1]), std::string(argv[2]));
  PuLSEHTMLWriter htmlwriter(rundata);

  // Output some essentails to the termial
  std::cerr << "Data file:\t\t\t" << rundata->fastaFilename << "\n";
  std::cerr << "Library definition:\t\t" << rundata->libraryDefinition<< "\n";
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

  //Output positional counts and heatmaps.
  rundata->PopulateHeatmaps();
  htmlwriter.WriteHeatMaps();
  
  return 0;
}
