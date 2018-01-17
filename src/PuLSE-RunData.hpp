/*
PuLSE version 1.4
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

#pragma once
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#ifdef __linux__
#include "zstr/src/zstr.hpp"
#endif

#include <vector>
class RunData {
public:
  std::string residues = "FLIMVSPTAYHQENKDCWRG*";
  std::string bases = "ATCG";
  std::string fastaFilename;
  std::string libraryDefinition;
  std::string nonstandardCodesUsed = "";
  unsigned randomDNALength;
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

  std::vector<std::vector<char>> validBases;
  std::vector<std::unordered_map<char, float>> protExpectedRates;
  std::vector<std::unordered_map<char, float>> dnaExpectedRates;
  std::string randomStretchDefinition;

  std::vector<unsigned int> nTimesProtSeen, nTimesDNASeen,
      nTimesProtSeenCumulative, nTimesDNASeenCumulative;

  // DNA triplet to AA mapping.  Entries can be changed/remapped via command
  // line arguments
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
    // Ensure library definition is uppercase
    std::transform(libraryDefinition.begin(), libraryDefinition.end(),
                   libraryDefinition.begin(), toupper);

    auto rndBegin = libraryDefinition.begin();
    while (*rndBegin != 'N' && *rndBegin != 'K' && *rndBegin != 'X' &&
           *rndBegin != 'R' && *rndBegin != 'Y' && *rndBegin != 'S' &&
           *rndBegin != 'W' && *rndBegin != 'M' && *rndBegin != 'B' &&
           *rndBegin != 'D' && *rndBegin != 'H' && *rndBegin != 'V') {
      ++rndBegin;
    }

    auto rndEnd = libraryDefinition.end() - 1;
    while (*rndEnd != 'N' && *rndEnd != 'K' && *rndEnd != 'X' &&
           *rndEnd != 'R' && *rndEnd != 'Y' && *rndEnd != 'S' &&
           *rndEnd != 'W' && *rndEnd != 'M' && *rndEnd != 'B' &&
           *rndEnd != 'D' && *rndEnd != 'H' && *rndEnd != 'V') {
      --rndEnd;
    }
    ++rndEnd;

    sequenceBeginMarker = std::string(libraryDefinition.begin(), rndBegin);
    sequenceEndMarker = std::string(rndEnd, libraryDefinition.end());
    randomStretchDefinition = std::string(rndBegin, rndEnd);
    randomDNALength = randomStretchDefinition.length();

    for (const auto &i : randomStretchDefinition) {

      if (i == 'X') {
        validBases.push_back({'A', 'C', 'T', 'G'});
        dnaExpectedRates.push_back(std::unordered_map<char, float>{
            std::make_pair('A', 0.25), std::make_pair('T', 0.25),
            std::make_pair('C', 0.25), std::make_pair('G', 0.25)});
        continue;
      }
      if (i == 'N') {
        validBases.push_back({'A', 'C', 'T', 'G'});
        dnaExpectedRates.push_back(std::unordered_map<char, float>{
            std::make_pair('A', 0.25), std::make_pair('T', 0.25),
            std::make_pair('C', 0.25), std::make_pair('G', 0.25)});
        continue;
      }
      if (i == 'K') {
        validBases.push_back({'T', 'G'});
        dnaExpectedRates.push_back(std::unordered_map<char, float>{
            std::make_pair('A', 0), std::make_pair('T', 0.5),
            std::make_pair('C', 0), std::make_pair('G', 0.5)});
        continue;
      }
      if (i == 'A') {
        validBases.push_back({'A'});
        dnaExpectedRates.push_back(std::unordered_map<char, float>{
            std::make_pair('A', 1.0), std::make_pair('T', 0),
            std::make_pair('C', 0), std::make_pair('G', 0)});
        continue;
      }
      if (i == 'T') {
        validBases.push_back({'T'});
        dnaExpectedRates.push_back(std::unordered_map<char, float>{
            std::make_pair('A', 0), std::make_pair('T', 1),
            std::make_pair('C', 0), std::make_pair('G', 0)});
        continue;
      }
      if (i == 'C') {
        validBases.push_back({'C'});
        dnaExpectedRates.push_back(std::unordered_map<char, float>{
            std::make_pair('A', 0.0), std::make_pair('T', 0),
            std::make_pair('C', 1), std::make_pair('G', 0)});
        continue;
      }
      if (i == 'G') {
        validBases.push_back({'G'});
        dnaExpectedRates.push_back(std::unordered_map<char, float>{
            std::make_pair('A', 0), std::make_pair('T', 0),
            std::make_pair('C', 0), std::make_pair('G', 1)});
        continue;
      }

      if (i == 'R') {
        validBases.push_back({'G', 'A'});
        dnaExpectedRates.push_back(std::unordered_map<char, float>{
            std::make_pair('A', 0.5), std::make_pair('T', 0),
            std::make_pair('C', 0), std::make_pair('G', 0.5)});
        continue;
      }
      if (i == 'Y') {
        validBases.push_back({'C', 'T'});
        dnaExpectedRates.push_back(std::unordered_map<char, float>{
            std::make_pair('A', 0), std::make_pair('T', 0.5),
            std::make_pair('C', 0.5), std::make_pair('G', 0)});
        continue;
      }
      if (i == 'S') {
        validBases.push_back({'G', 'C'});
        dnaExpectedRates.push_back(std::unordered_map<char, float>{
            std::make_pair('A', 0), std::make_pair('T', 0),
            std::make_pair('C', 0.5), std::make_pair('G', 0.5)});
        continue;
      }
      if (i == 'W') {
        validBases.push_back({'A', 'T'});
        dnaExpectedRates.push_back(std::unordered_map<char, float>{
            std::make_pair('A', 0.5), std::make_pair('T', 0.5),
            std::make_pair('C', 0), std::make_pair('G', 0)});
        continue;
      }
      if (i == 'M') {
        validBases.push_back({'A', 'C'});
        dnaExpectedRates.push_back(std::unordered_map<char, float>{
            std::make_pair('A', 0.5), std::make_pair('T', 0),
            std::make_pair('C', 0.5), std::make_pair('G', 0)});
        continue;
      }
      if (i == 'B') {
        validBases.push_back({'C', 'G', 'T'});
        dnaExpectedRates.push_back(std::unordered_map<char, float>{
            std::make_pair('A', 0), std::make_pair('T', 0.333333),
            std::make_pair('C', 0.333333), std::make_pair('G', 0.333333)});
        continue;
      }
      if (i == 'D') {
        validBases.push_back({'A', 'T', 'G'});
        dnaExpectedRates.push_back(std::unordered_map<char, float>{
            std::make_pair('A', 0.333333), std::make_pair('T', 0.333333),
            std::make_pair('C', 0), std::make_pair('G', 0.333333)});
        continue;
      }
      if (i == 'H') {
        validBases.push_back({'A', 'C', 'T'});
        dnaExpectedRates.push_back(std::unordered_map<char, float>{
            std::make_pair('A', 0.333333), std::make_pair('T', 0.333333),
            std::make_pair('C', 0.333333), std::make_pair('G', 0)});
        continue;
      }
      if (i == 'V') {
        validBases.push_back({'A', 'C', 'G'});
        dnaExpectedRates.push_back(std::unordered_map<char, float>{
            std::make_pair('A', 0.333333), std::make_pair('T', 0),
            std::make_pair('C', 0.333333), std::make_pair('G', 0.333333)});
        continue;
      }
    }

    // Now populate the protExpectedRates
    {
      // Populate std::vector<std::unordered_map<std::string, float>>
      // protExpectedRates based on
      //  std::vector<std::vector<char>> validBases;
      // and std::unordered_map<std::string, char> tripletToAAMap
      std::string triplet = "";
      std::unordered_map<char, float> blankmap;
      for (auto &&i : residues) {
        blankmap[i] = 0.0;
      }
      for (unsigned i = 0; i < randomDNALength; i += 3) {
        protExpectedRates.push_back(blankmap);
        for (auto &&i1 : validBases[i]) {
          for (auto &&i2 : validBases[i + 1]) {
            for (auto &&i3 : validBases[i + 2]) {
              triplet =
                  std::string(1, i1) + std::string(1, i2) + std::string(1, i3);
              for (auto &&i : triplet) {
                if (i == 'T')
                  i = 'U';
              }

              protExpectedRates.back()[tripletToAAMap[triplet]] += 1;
            }
          }
        }
        // Now normalise validResidues.back;
        float total = 0;

        for (auto &&freq : protExpectedRates.back()) {
          total += freq.second;
        }
        for (auto &&freq : protExpectedRates.back()) {
          freq.second /= total;
        }
      }
    }
  }

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

  // Start reading sequences
  void ReadSequences() {
    std::string protsequence;
    unsigned ntimesfound;

    // Lambda to insert prot into map (if it does not already exist),
    // otherwise, increment the occurance count
    auto insertIntoProtMap = [&]() {
      auto index = map_proteinSequences.find(protsequence);
      if (index == map_proteinSequences.end()) {
        map_proteinSequences[protsequence] = ntimesfound;
      } else {
        index->second += ntimesfound;
      }
    };

    // Collect DNA reads in a map (DNA sequence --> count)
    map_DNASequences =
        ReadPhageLibraryQC(fastaFilename, sequenceBeginMarker,
                           sequenceEndMarker, randomDNALength, numReads);

    // Map DNA sequences to protein sequences map
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

  // Populate simple count statistics
  void PopulateBasicStats() {
    nProtMax = static_cast<int>((std::pow(20, (randomDNALength / 3))));
    nProtNotFound = nProtMax - map_proteinSequences.size();
    nDNAMax = static_cast<int>(std::pow(4, randomDNALength));
    nDNANotFound = nDNAMax - map_DNASequences.size();
  };

  // Populate cumulative counts
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

  // Populate the most common protein and DNA sequences
  void PopulateCommonOccurances() {
    for (auto &&i : map_DNASequences) {
      if (i.second == 1) {
        continue;
      }
      commondna.push_back(std::make_pair(i.second, i.first));
    }
    std::sort(commondna.rbegin(), commondna.rend());
    // Only keep the top 100 most frequently occuring DNA sequences.
    if (commondna.size() > 100) {
      commondna.resize(100);
    }
    for (auto &&i : map_proteinSequences) {
      if (i.second == 1) {
        continue;
      }
      commonprot.push_back(std::make_pair(i.second, i.first));
    }
    // Only keep the top 100 most commonly occuring protein sequences.
    std::sort(commonprot.rbegin(), commonprot.rend());
    if (commonprot.size() > 100) {
      commonprot.resize(100);
    }
  };

  // Make heatmaps
  void PopulateHeatmaps() {
    std::vector<std::unordered_map<char, unsigned int>> tmp_protmap(
        randomDNALength / 3, std::unordered_map<char, unsigned int>());
    std::vector<std::unordered_map<char, unsigned int>> tmp_dnamap(
        randomDNALength, std::unordered_map<char, unsigned int>());

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
      for (unsigned aa = 0; aa < randomDNALength / 3; ++aa) {
        tmp_protmap[aa][i.first[aa]] += i.second;
      }
    }
    for (auto &&i : map_DNASequences) {
      for (unsigned aa = 0; aa < randomDNALength; ++aa) {
        tmp_dnamap[aa][i.first[aa]] += i.second;
      }
    }

    /* for (unsigned i = 0; i < randomStretchDefinition.length(); ++i) {
       for (auto &&i : tmp_protmap[i]) {
         i.second *=
       }
     }*/

    // Populate heatmap vectors
    protheatmap = std::vector<std::vector<std::string>>(
        residues.length(), std::vector<std::string>());
    for (unsigned i = 0; i < residues.length(); ++i) {
      protheatmap[i].push_back(std::string(1, residues[i]));
      for (unsigned j = 0; j < randomDNALength / 3; ++j) {
        protheatmap[i].push_back(
            std::to_string(tmp_protmap[j][residues.at(i)]));
      }
    }
    dnaheatmap = std::vector<std::vector<std::string>>(
        bases.length(), std::vector<std::string>());
    for (unsigned i = 0; i < bases.length(); ++i) {
      dnaheatmap[i].push_back(std::string(1, bases[i]));
      for (unsigned j = 0; j < randomDNALength; ++j) {
        dnaheatmap[i].push_back(std::to_string(tmp_dnamap[j][bases.at(i)]));
      }
    }
  };

private:
  // Read fastq library file containing reads and populate a map containing
  // DNA sequences and occurances
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
// If on linux, then we can optionally read gz compressed input fastq files.
#ifdef __linux__
    zstr::ifstream inputfile(fastaFilename);
#else
    // If not, then just read a standard, plaintext fastq file.
    std::ifstream inputfile(fastaFilename);
#endif
    std::string sequence;

    // Lambda to insert DNA sequence into map (also containing occurance
    // counts)
    auto insertSequenceIntoMap = [&]() {
      auto index = sequences.find(sequence);
      if (index == sequences.end()) {
        sequences[sequence] = 1;
      } else {
        ++(index->second);
      }
    };

    // Lambda to extract the forward sequence
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

    // Lambda to extract the reverse sequence
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

    // Read the input fastq file
    while (getline(inputfile, line)) {
      // Capitalise
      std::transform(line.begin(), line.end(), line.begin(), toupper);
      // Read forwards
      sequence = getForwardSequence();
      if (sequence.size() == dnalength) {
        ++numReads;
        insertSequenceIntoMap();
        continue;
      }
      // Now try to read the reverse
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

  // Converts string of DNA into AA sequence;
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
