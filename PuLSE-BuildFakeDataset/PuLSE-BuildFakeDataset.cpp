/*
PuLSE-BuildFakeDataset version 1.3
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
#include "../src/PuLSE-RunData.hpp"
#include "../src/PuLSE-Writer.hpp"

#include <iostream>

#include <random>
#include <vector>

bool parseInput(const std::string &input,
                std::vector<std::vector<float>> &probabilities,
                std::string &upstreamMarker, std::string &downstreamMarker) {
  size_t numCommas = std::count(input.begin(), input.end(), ',') - 1;
  upstreamMarker = input.substr(0, input.find(",", 0));
  if (numCommas % 4 != 0)
    return false;
  std::vector<float> quad(4);
  auto pos1 = input.find(",", 0) + 1;
  auto pos2 = input.find(",", pos1);

  for (unsigned i = 0; i < numCommas; i += 4) {
    quad[0] = std::atof(input.substr(pos1, pos2 - pos1).c_str());
    pos1 = pos2 + 1;
    pos2 = input.find(",", pos1);

    quad[1] = std::atof(input.substr(pos1, pos2 - pos1).c_str());
    pos1 = pos2 + 1;
    pos2 = input.find(",", pos1);

    quad[2] = std::atof(input.substr(pos1, pos2 - pos1).c_str());
    pos1 = pos2 + 1;
    pos2 = input.find(",", pos1);

    quad[3] = std::atof(input.substr(pos1, pos2 - pos1).c_str());
    pos1 = pos2 + 1;
    pos2 = input.find(",", pos1);

    probabilities.push_back(quad);
  }
  downstreamMarker = input.substr(input.find_last_of(',') + 1);
  return true;
};

int main(int argc, char **argv) {

  std::random_device rd;
  std::mt19937 gen(rd());

  // Generate n fairly chosen bases and return as a string.
  std::discrete_distribution<int> fair({25, 25, 25, 25});
  auto genFairBases = [&gen, &fair](unsigned n = 1) {
    unsigned basenumber;
    std::string out;
    for (unsigned i = 0; i < n; ++i) {
      basenumber = fair(gen);
      if (basenumber == 0) {
        out += "A";
        continue;
      }
      if (basenumber == 1) {
        out += "T";
        continue;
      }
      if (basenumber == 2) {
        out += "C";
        continue;
      }
      if (basenumber == 3) {
        out += "G";
        continue;
      }
    }
    return out;
  };

  std::cerr
      << "Generate a fake dataset\n\n"
         "Pass arguments to as strings of comma separated likelihoods for the "
         "four ATCG bases appearing\nflanked by the up and downstream "
         "markers.\n"
         "Output (in the fastq file format) should be piped to a file.\n"
         "As an example the following input string represents a library with "
         "10,000 members, flanked by up and downstream markers,\n"
         "flanking a singular randomised amino acid.\n"
         "amino acids up and downstream respectively:\n"
         "CGTTGC,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,"
         "TGTGCT 10000\n"
         "Whereby CGTTGC is the upstream marker sequence\n"
         "The first four 0.25 values represent occurance rates for A, T, C and "
         "G bases respectively in the first base position.\n"
         "The next four 0.25 values represent occurance rates for A, T, C and "
         "G bases respectively in the second base position.\n"
         "The last four 0.25 values represent occurance rates for A, T, C and "
         "G bases respectively in the third base position.\n"
         "\nThe sequence is simply expanded to add more randomised amino acid "
         "positions, for example a perfectly random trimer could be specified "
         "as follows:\n"
         "CGTTGC,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0."
         "25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,"
         "0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25,TGTGCT\n";

  unsigned nSequences = std::atoi(argv[2]);
  std::vector<std::vector<float>> probabilities;
  std::string upstreamMarker = "CGTTGC";
  std::string downstreamMarker = "TGTGCT";

  // Parse big input string.

  parseInput(std::string(argv[1]), probabilities, upstreamMarker,
             downstreamMarker);

  std::vector<std::discrete_distribution<>> positions;
  for (auto i : probabilities) {
    positions.push_back(std::discrete_distribution<int>(i.begin(), i.end()));
  }
  auto genBases = [&gen, &positions]() {
    std::string out;
    for (auto &&pos : positions) {
      unsigned basenumber = pos(gen);
      if (basenumber == 0) {
        out += "A";
        continue;
      }
      if (basenumber == 1) {
        out += "T";
        continue;
      }
      if (basenumber == 2) {
        out += "C";
        continue;
      }
      if (basenumber == 3) {
        out += "G";
        continue;
      }
    }
    return out;
  };

  std::discrete_distribution<int> corruptDistribution({ 0.01,0.99 });
  std::normal_distribution<> normalDist(5, 3);
  std::string line;
  int numToDeleteUp = 0, numToDeleteDown = 0;
  for (unsigned sequenceNumber = 0; sequenceNumber < nSequences;
       sequenceNumber++) {
    line = "";
    line += genFairBases(15);
    line += upstreamMarker;
    line += genBases();
    line += downstreamMarker;
    line += genFairBases(15);
    numToDeleteUp = normalDist(gen);
    numToDeleteDown = normalDist(gen);
    if (numToDeleteUp < 0)
      numToDeleteUp = 0;
    if (numToDeleteDown < 0)
      numToDeleteDown = 0;

	line.erase(0, numToDeleteUp);
	line.erase(line.length() - numToDeleteDown);
	
	if (corruptDistribution(gen) == 0) {
		std::uniform_int_distribution<int> dist(0, line.length()-1);
		line.erase(dist(gen), 1);
	}
		
	std::cout << "@Simulated line " << sequenceNumber+1 << "\n"
		<< line<<"\n+\n#\n";
  }

  std::cout << "\n";

  return 0;
}
