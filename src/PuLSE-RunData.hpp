#pragma once
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <cmath>

#include <vector>
class RunData {
public:
  std::string residues = "FLIMVSPTAYHQENKDCWRG*";
  std::string bases = "ATCG";
  std::string fastaFilename;
  std::string libraryDefinition;
  std::string nonstandardCodesUsed = "";
  unsigned dnalength;
  unsigned numReads;
  std::string sequenceBeginMarker, sequenceEndMarker;
  std::unordered_map<std::string, unsigned> map_proteinSequences;
  std::unordered_map<std::string, unsigned> map_DNASequences;
  std::vector<std::pair<unsigned int, std::string>> commondna, commonprot;
  std::vector<std::vector<std::string>> protheatmap, dnaheatmap;
  int nProtMax;
  int nProtNotFound;
  int nDNAMax;
  int nDNANotFound;
  unsigned dnamaxoccurances = 0;
  unsigned protmaxoccurances = 0;
  using seenCount = unsigned int;
  using seenCountCumulative = unsigned int;
  std::vector<unsigned int> nTimesProtSeen, nTimesDNASeen,
      nTimesProtSeenCumulative, nTimesDNASeenCumulative;
  std::unordered_map<std::string, char> tripletToAAMap{
      std::make_pair("UUU", 'F'), std::make_pair("UUC", 'F'),
      std::make_pair("UUA", 'L'), std::make_pair("UUG", 'L'),
      std::make_pair("CUU", 'L'), std::make_pair("CUC", 'L'),
      std::make_pair("CUA", 'L'), std::make_pair("CUG", 'L'),
      std::make_pair("AUU", 'I'), std::make_pair("AUC", 'I'),
      std::make_pair("AUA", 'I'), std::make_pair("AUG", 'M'),
      std::make_pair("GUU", 'V'), std::make_pair("GUC", 'V'),
      std::make_pair("GUA", 'V'), std::make_pair("GUG", 'V'),
      std::make_pair("UCU", 'S'), std::make_pair("UCC", 'S'),
      std::make_pair("UCA", 'S'), std::make_pair("UCG", 'S'),
      std::make_pair("AGU", 'S'), std::make_pair("AGC", 'S'),
      std::make_pair("CCU", 'P'), std::make_pair("CCC", 'P'),
      std::make_pair("CCA", 'P'), std::make_pair("CCG", 'P'),
      std::make_pair("ACU", 'T'), std::make_pair("ACC", 'T'),
      std::make_pair("ACA", 'T'), std::make_pair("ACG", 'T'),
      std::make_pair("GCU", 'A'), std::make_pair("GCC", 'A'),
      std::make_pair("GCA", 'A'), std::make_pair("GCG", 'A'),
      std::make_pair("UAU", 'Y'), std::make_pair("UAC", 'Y'),
      std::make_pair("UAA", '*'), std::make_pair("UAG", '*'),
      std::make_pair("UGA", '*'), std::make_pair("CAU", 'H'),
      std::make_pair("CAC", 'H'), std::make_pair("CAA", 'Q'),
      std::make_pair("CAG", 'Q'), std::make_pair("GAA", 'E'),
      std::make_pair("GAG", 'E'), std::make_pair("AAU", 'N'),
      std::make_pair("AAC", 'N'), std::make_pair("AAA", 'K'),
      std::make_pair("AAG", 'K'), std::make_pair("GAU", 'D'),
      std::make_pair("GAC", 'D'), std::make_pair("UGU", 'C'),
      std::make_pair("UGC", 'C'), std::make_pair("UGG", 'W'),
      std::make_pair("CGU", 'R'), std::make_pair("CGC", 'R'),
      std::make_pair("CGA", 'R'), std::make_pair("CGG", 'R'),
      std::make_pair("AGA", 'R'), std::make_pair("AGG", 'R'),
      std::make_pair("GGU", 'G'), std::make_pair("GGC", 'G'),
      std::make_pair("GGA", 'G'), std::make_pair("GGG", 'G')};

  // Constructor
  explicit RunData(std::string &&ffn, std::string &&libDef) {
    fastaFilename = std::move(ffn);
    libraryDefinition = std::move(libDef);
    std::transform(libraryDefinition.begin(), libraryDefinition.end(),
                   libraryDefinition.begin(), toupper);
    dnalength =
        std::count(libraryDefinition.begin(), libraryDefinition.end(), 'X');
    sequenceBeginMarker = std::string(
        libraryDefinition.begin(),
        std::find(libraryDefinition.begin(), libraryDefinition.end(), 'X'));
    sequenceEndMarker = std::string(libraryDefinition.begin() +
                                        libraryDefinition.find_last_of("X") + 1,
                                    libraryDefinition.end());
  };

  // Alter the tripletToAAMap, which converts DNA triplets into protein AA
  // residue characters.  Using this, expression systems with nonsense
  // surpression may be analysed; for example by specifying "UAG Q" on the
  // command line.
  void populateNonStandardTriplets(const int argc, char **argv) {
    for (int i = 3; i < argc; i += 2) {
      nonstandardCodesUsed +=
          std::string(argv[i]) + "-->" + std::string(argv[i + 1]) + ", ";
      tripletToAAMap[argv[i]] = argv[i + 1][0];
      std::cerr << "Added customized tripplet to DNA mapping: " << argv[i]
                << " --> " << tripletToAAMap[argv[i]] << "\n";
    }
    if (nonstandardCodesUsed.length() == 0) {
      nonstandardCodesUsed = "None";
    } else {
      nonstandardCodesUsed = std::string(nonstandardCodesUsed.begin(),
                                         nonstandardCodesUsed.end() - 2);
    }
  };

  //Start reading sequences
  void ReadSequences() {
    std::string protsequence;
    unsigned ntimesfound;
    auto insertIntoProtMap = [&]() {
      auto index = map_proteinSequences.find(protsequence);
      if (index == map_proteinSequences.end()) {
        map_proteinSequences[protsequence] = ntimesfound;
      } else {
        index->second += ntimesfound;
      }
    };

    map_DNASequences =
        ReadPhageLibraryQC(fastaFilename, sequenceBeginMarker,
                           sequenceEndMarker, dnalength, numReads);

    for (auto &&i : map_DNASequences) {
      protsequence = translateDNA(replaceChar(i.first, 'T', 'U'));
      ntimesfound = i.second;
      insertIntoProtMap();
    }
    dnamaxoccurances =
        std::max_element(map_DNASequences.begin(), map_DNASequences.end(),
                         [](const auto &p1, const auto &p2) {
                           return p1.second < p2.second;
                         })
            ->second;
    protmaxoccurances = std::max_element(map_proteinSequences.begin(),
                                         map_proteinSequences.end(),
                                         [](const auto &p1, const auto &p2) {
                                           return p1.second < p2.second;
                                         })
                            ->second;
  };

  //Populate simple count statistics
  void PopulateBasicStats() {
    nProtMax = static_cast<int>((std::pow(20, (dnalength / 3))));
    nProtNotFound = nProtMax - map_proteinSequences.size();
    nDNAMax = static_cast<int>(std::pow(4, dnalength));
    nDNANotFound = nDNAMax - map_DNASequences.size();
  };

  //Populate cumulative counts
  void PopulateCumulativeOccurances() {
    nTimesDNASeen = std::vector<unsigned int>(dnamaxoccurances + 1, 0);
    nTimesDNASeenCumulative = nTimesDNASeen;
    nTimesProtSeen = std::vector<unsigned int>(protmaxoccurances + 1, 0);
    nTimesProtSeenCumulative = nTimesProtSeen;
    nTimesDNASeen[0] = nDNANotFound;
    nTimesProtSeen[0] = nProtNotFound;

    for (auto &&i : map_DNASequences) {
      ++nTimesDNASeen[i.second];
    }
    for (auto &&i : map_proteinSequences) {
      ++nTimesProtSeen[i.second];
    }
    nTimesProtSeenCumulative[nTimesProtSeenCumulative.size() - 1] =
        nTimesProtSeen[nTimesProtSeen.size() - 1];
    nTimesDNASeenCumulative[nTimesDNASeenCumulative.size() - 1] =
        nTimesDNASeen[nTimesDNASeen.size() - 1];
    for (int i = nTimesProtSeenCumulative.size() - 2; i >= 0; --i) {
      nTimesProtSeenCumulative[i] =
          nTimesProtSeenCumulative[i + 1] + nTimesProtSeen[i];
    }
    for (int i = nTimesDNASeenCumulative.size() - 2; i >= 0; --i) {
      nTimesDNASeenCumulative[i] =
          nTimesDNASeenCumulative[i + 1] + nTimesDNASeen[i];
    }
  };

  //Populate the most common protein and DNA sequences
  void PopulateCommonOccurances() {
    for (auto &&i : map_DNASequences) {
      if (i.second == 1) {
        continue;
      }
      commondna.push_back(std::make_pair(i.second, i.first));
    }
    std::sort(commondna.rbegin(), commondna.rend());
    if (commondna.size() > 100) {
      commondna.resize(100);
    }
    for (auto &&i : map_proteinSequences) {
      if (i.second == 1) {
        continue;
      }
      commonprot.push_back(std::make_pair(i.second, i.first));
    }
    std::sort(commonprot.rbegin(), commonprot.rend());
    if (commonprot.size() > 100) {
      commonprot.resize(100);
    }
  };

  //Make heatmaps
  void PopulateHeatmaps() {
    std::vector<std::unordered_map<char, unsigned int>> tmp_protmap(
        dnalength / 3, std::unordered_map<char, unsigned int>());
    std::vector<std::unordered_map<char, unsigned int>> tmp_dnamap(
        dnalength, std::unordered_map<char, unsigned int>());

    // Initialise temp maps to zero
    for (auto &&phm : tmp_protmap) {
      for (auto &&res : residues) {
        phm[res] = 0;
      }
    }
    for (auto &&dhm : tmp_dnamap) {
      for (auto &&base : bases) {
        dhm[base] = 0;
      }
    }

    // Count chars in appropriate positions
    for (auto &&i : map_proteinSequences) {
      for (unsigned aa = 0; aa < dnalength / 3; ++aa) {
        tmp_protmap[aa][i.first[aa]] += i.second;
      }
    }
    for (auto &&i : map_DNASequences) {
      for (unsigned aa = 0; aa < dnalength; ++aa) {
        tmp_dnamap[aa][i.first[aa]] += i.second;
      }
    }

    // Populate heatmap vectors
    protheatmap = std::vector<std::vector<std::string>>(
        residues.length(), std::vector<std::string>());
    for (unsigned i = 0; i < residues.length(); ++i) {
      protheatmap[i].push_back(std::string(1, residues[i]));
      for (unsigned j = 0; j < dnalength / 3; ++j) {
        protheatmap[i].push_back(
            std::to_string(tmp_protmap[j][residues.at(i)]));
      }
    }
    dnaheatmap = std::vector<std::vector<std::string>>(
        bases.length(), std::vector<std::string>());
    for (unsigned i = 0; i < bases.length(); ++i) {
      dnaheatmap[i].push_back(std::string(1, bases[i]));
      for (unsigned j = 0; j < dnalength; ++j) {
        dnaheatmap[i].push_back(std::to_string(tmp_dnamap[j][bases.at(i)]));
      }
    }
  };

private:
  // Read fastq library file containing reads and populate
  std::unordered_map<std::string, unsigned int>
  ReadPhageLibraryQC(const std::string &fastaFilename,
                     const std::string &beginMarker,
                     const std::string &endMarker, const unsigned &dnalength,
                     unsigned &numReads) {
    numReads = 0;

    std::string shortBeginMarker = beginMarker;

    if (shortBeginMarker.length() > 3) {
      shortBeginMarker =
          shortBeginMarker.substr(shortBeginMarker.size() - 3, 3);
    }

    std::string shortEndMarker = endMarker;
    if (shortEndMarker.length() > 3) {
      shortEndMarker = shortEndMarker.substr(0, 3);
    }

    std::unordered_map<std::string, unsigned int> sequences;
    std::string line;
    std::ifstream inputfile(fastaFilename);
    std::string sequence;

    auto insertSequenceIntoMap = [&]() {
      auto index = sequences.find(sequence);
      if (index == sequences.end()) {
        sequences[sequence] = 1;
      } else {
        ++(index->second);
      }
    };

    auto getForwardSequence = [&]() {
      auto loc = line.find(beginMarker);
      while (loc != std::string::npos) {
        if (line.size() >
            loc + beginMarker.size() + dnalength + shortEndMarker.size()) {
          if (line.substr(loc + beginMarker.size() + dnalength,
                          shortEndMarker.size())
                  .compare(shortEndMarker) == 0) {
            return line.substr(loc + beginMarker.size(), dnalength);
          }
        }
        loc = line.find(beginMarker, loc + 1);
      }
      return std::string();
    };

    auto getReverseSequence = [&]() {
      auto loc = line.find(shortBeginMarker);
      while (loc != std::string::npos) {
        if (line.size() >
            loc + shortBeginMarker.size() + dnalength + endMarker.size()) {
          if (line.substr(loc + shortBeginMarker.size() + dnalength,
                          endMarker.size())
                  .compare(endMarker) == 0) {
            return line.substr(loc + shortBeginMarker.size(), dnalength);
          }
        }
        loc = line.find(shortBeginMarker, loc + 1);
      }
      return std::string();
    };

    while (getline(inputfile, line)) {
      std::transform(line.begin(), line.end(), line.begin(), toupper);
      sequence = getForwardSequence();
      if (sequence.size() == dnalength) {
        ++numReads;
        insertSequenceIntoMap();
        continue;
      }
      std::reverse(line.begin(), line.end());
      if (swapBasesInPlace(line)) {
        sequence = getReverseSequence();
        if (sequence.size() == dnalength) {
          ++numReads;
          insertSequenceIntoMap();
          continue;
        }
      }
    }
    // Copy elision means that this is not an expensive copy
    return sequences;
  };

  // Change DNA sequence to complimentary strand in place
  bool swapBasesInPlace(std::string &in) {
	  for (auto &&c : in) {
		  if (c == 'A') {
			  c = 'T';
			  continue;
		  }
		  if (c == 'T') {
			  c = 'A';
			  continue;
		  }
		  if (c == 'C') {
			  c = 'G';
			  continue;
		  }
		  if (c == 'G') {
			  c = 'C';
			  continue;
		  }
		  return false;
	  }
	  return true;
  };

  // Replace a character in a string, returning the new string
  std::string replaceChar(std::string str, const char ch1, const char ch2) {
	  std::replace(str.begin(), str.end(), ch1, ch2);
	  return str;
  };

  // Convert a DNA triplet (string of length 3) into protein AA residue, such
  // as: "CAT" --> 'H'
  std::string translateDNATriplet(const std::string &in) {
	  auto location = tripletToAAMap.find(in);
	  if (location == tripletToAAMap.end()) {
		  return "-";
	  }
	  return std::string(1, location->second);
  };

  //Converts string of DNA into AA sequence;
  std::string translateDNA(std::string in) {
	  std::string out = "";
	  if (in.length() % 3 != 0) {
		  return out;
	  }
	  for (unsigned int i = 0; i < in.length(); i += 3) {
		  out += translateDNATriplet(in.substr(i, 3));
	  }
	  return out;
  };
};
