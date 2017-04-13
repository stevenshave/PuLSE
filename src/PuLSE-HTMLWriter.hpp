#pragma once
#include "PuLSE-RunData.hpp"
#include <cctype>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

class PuLSEHTMLWriter {
private:
  std::ofstream file;
  std::shared_ptr<RunData> rundata;

public:
  explicit PuLSEHTMLWriter(std::shared_ptr<RunData> indata) {
    rundata = std::move(indata);
    file.open(
        rundata->fastaFilename.substr(0, rundata->fastaFilename.rfind(".")) +
        ".html");
    if (!file) {
      std::cerr << "Error opening HTML output file: "
                << rundata->fastaFilename.substr(
                       0, rundata->fastaFilename.rfind(".")) +
                       ".html"
                << ", exiting\n";
      exit(-1);
    }
    file << "<!DOCTYPE html>"
            "<html lang = \"en\">"
            "<head>"
            "<title>PuLSE output - Phage Library Sequence Evaluation </title>"
            "<meta charset = \"utf-8\">"
            "<meta name = \"viewport\" content = \"width=device-width, "
            "initial-scale=1\">"
            "<link rel = \"stylesheet\" href = "
            "\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/"
            "bootstrap.min.css\">"
            "<script src = "
            "\"https://ajax.googleapis.com/ajax/libs/jquery/3.2.0/"
            "jquery.min.js\"></script>"
            "<script src = "
            "\"https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/"
            "bootstrap.min.js\"></script>"
            "</head>"
            "<body>"
            "<div class = \"container\">"
            "<h1>PuLSE output</h1>"
            "</div>";
  };

  void WriteRunInfo() {
    file << "<div class = \"container\">"
            "<h3>Run information</h3>"
            "</div>"
            "<div class=\"container\"><table class=\"table\">"
            "<thread><tr><th>Parameter</th><th>Value</th></tr></thread>"
            "<tbody><tr>"

            "<tr><td>Data file</td><td>"
         << rundata->fastaFilename
         << "</td></tr>"
            "<tr><td>Library definition</td><td>"
         << rundata->libraryDefinition
         << "</td></tr>"
            "<tr><td>Forward upstream marker</td><td>"
         << rundata->sequenceBeginMarker
         << "</td></tr>"
            "<tr><td>Forward downstream  marker</td><td>"
         << rundata->sequenceEndMarker
         << "</td></tr>"
            "<tr><td>Randomized DNA positions</td><td>"
         << rundata->dnalength << "</td></tr>"
         << "<tr><td>Non-standard AA codes used</td><td>"
         << rundata->nonstandardCodesUsed << "</td></tr>"
         << "</tbody></table>"
            "</div>";
  };

  void CloseHTMLFile() {
    file << "</body>"
            "</html>";
    file.close();
  };

  void WriteBasicStats() {
    std::vector<std::vector<std::string>> data;
    insertTableRow(data, std::string("Number of valid reads"),
                   rundata->numReads);

    insertTableRow(data, "Found unique at dna level",
                   std::to_string(rundata->map_DNASequences.size()) + " (" +
                       std::to_string((rundata->map_DNASequences.size() /
                                       static_cast<float>(rundata->nDNAMax)) *
                                      100.0) +
                       "%)");
    insertTableRow(data, "Found unique at protein level",
                   std::to_string(rundata->map_proteinSequences.size()) + " (" +
                       std::to_string((rundata->map_proteinSequences.size() /
                                       static_cast<float>(rundata->nProtMax)) *
                                      100.0) +
                       "%)");
    insertTableRow(data, "DNA sequences not found",
                   std::to_string(rundata->nDNANotFound) + " (" +
                       std::to_string((rundata->nDNANotFound /
                                       static_cast<float>(rundata->nDNAMax)) *
                                      100) +
                       " %)");
    insertTableRow(data, "Protein sequences not found",
                   std::to_string(rundata->nProtNotFound) + " (" +
                       std::to_string((rundata->nProtNotFound /
                                       static_cast<float>(rundata->nProtMax)) *
                                      100) +
                       " %)");

    WriteHTMLTable(file, "Basic statistics",
                   "Properties of sequences found and not found",
                   std::vector<std::string>{"Property", "Value"}, data);
  };

  void WriteCumulativeCounts() {
    std::vector<std::vector<std::string>> data1, data2;
    unsigned int i = 0;
    for (; i < rundata->nTimesDNASeen.size() && i < 51; ++i) {
      insertTableRow(data1, i, rundata->nTimesDNASeen[i],
                     rundata->nTimesDNASeenCumulative[i],
                     (rundata->nTimesDNASeenCumulative[i] /
                      static_cast<float>(rundata->nDNAMax)) *
                         100);
    }
    for (; i < rundata->nTimesDNASeen.size(); ++i) {
      if (rundata->nTimesDNASeen[i] != 0) {
        insertTableRow(data1, i, rundata->nTimesDNASeen[i],
                       rundata->nTimesDNASeenCumulative[i],
                       (rundata->nTimesDNASeenCumulative[i] /
                        static_cast<float>(rundata->nDNAMax)) *
                           100);
      }
    }

    i = 0;

    for (; i < rundata->nTimesProtSeen.size() && i < 51; ++i) {
      insertTableRow(data2, i, rundata->nTimesProtSeen[i],
                     rundata->nTimesProtSeenCumulative[i],
                     (rundata->nTimesProtSeenCumulative[i] /
                      static_cast<float>(rundata->nProtMax)) *
                         100);
    }
    for (; i < rundata->nTimesProtSeen.size(); ++i) {
      if (rundata->nTimesProtSeen[i] != 0) {
        insertTableRow(data2, i, rundata->nTimesProtSeen[i],
                       rundata->nTimesProtSeenCumulative[i],
                       (rundata->nTimesProtSeenCumulative[i] /
                        static_cast<float>(rundata->nProtMax)) *
                           100);
      }
    }

    Write2HTMLTables(file, "Cumulative counts",
                     "The number of sequences (count) seen N times. Cumulative "
                     "count gives sequences seen at least N times.  N of 1 to "
                     "50 shown, then only when sequence found N times up to "
                     "the max occurance.",
                     {"DNA Occurances (N)", "Count", "Cumulative count",
                      "Cumulative count %"},
                     data1,
                     {"Prot Occurances (N)", "Count", "Cumulative count",
                      "Cumulative count %"},
                     data2);
  };

  void WriteCommonOccurances() {
    std::vector<std::vector<std::string>> data1, data2;
    for (unsigned i = 0; i < 100 && i < rundata->commondna.size(); ++i) {
      insertTableRow(data1, std::get<1>(rundata->commondna[i]),
                     std::get<0>(rundata->commondna[i]));
    }
    for (unsigned i = 0; i < 100 && i < rundata->commonprot.size(); ++i) {
      insertTableRow(data2, std::get<1>(rundata->commonprot[i]),
                     std::get<0>(rundata->commonprot[i]));
    }

    Write2HTMLTables(file, "Most common sequences",
                     "Counts of the top 100 occuring DNA and protein sequences",
                     {"DNA sequence", "Count"}, data1,
                     {"Protein sequence", "Count"}, data2);
  };

  void WriteHeatMaps() {

    std::vector<std::string> headings{"Residue"};
    for (unsigned i = 1; i <= rundata->dnalength / 3; ++i) {
      headings.push_back("Position #" + std::to_string(i));
    }
    WriteHTMLTable(
        file, "Protein residue counts",
        "Occurance counts of each AA resdue at every library position",
        headings, rundata->protheatmap);

    headings = {"Residue"};
    for (unsigned i = 1; i <= rundata->dnalength; ++i) {
      headings.push_back("Position #" + std::to_string(i));
    }
    WriteHTMLTable(
        file, "DNA base counts",
        "Occurance counts of each DNA base at every library position", headings,
        rundata->dnaheatmap);

    auto getResidueOccurance = [&rd = rundata](const char &c) {
      unsigned nTimesFound = 0;
      for (auto &&i : rd->tripletToAAMap) {
        if (i.second == c) {
          ++nTimesFound;
        }
      }
      return nTimesFound;
    };

    // Now write the normalised/enriched heatmaps
    std::unordered_map<char, unsigned int> residueExpectedOccurance;
    for (auto &&i : rundata->residues) {
      residueExpectedOccurance[i] = getResidueOccurance(i);
    }

    auto normalisedProtHeatmap = rundata->protheatmap;

    for (auto &&y : normalisedProtHeatmap) {
      for (unsigned x = 1; x < y.size(); ++x) {
        y[x] = std::to_string(
            std::stoi(y[x]) /
            (rundata->numReads * (residueExpectedOccurance[y[0][0]] / 64.0)));
      }
    }

    headings = {"Residue"};
    for (unsigned i = 1; i <= rundata->dnalength / 3; ++i) {
      headings.push_back("Position #" + std::to_string(i));
    }
    WriteHTMLTable(file, "Normalised protein residue heatmap",
                   "Enrichment over expected occurance of each AA resdue at "
                   "every library position",
                   headings, normalisedProtHeatmap);

    auto normalisedDNAHeatmap = rundata->dnaheatmap;

    for (auto &&y : normalisedDNAHeatmap) {
      for (unsigned x = 1; x < y.size(); ++x) {
        y[x] = std::to_string(std::stoi(y[x]) / (rundata->numReads * (0.25)));
      }
    }

    headings = {"Residue"};
    for (unsigned i = 1; i <= rundata->dnalength; ++i) {
      headings.push_back("Position #" + std::to_string(i));
    }
    WriteHTMLTable(file, "Normalised DNA base heatmap",
                   "Enrichment over expected occurance of each base at every "
                   "library position",
                   headings, normalisedDNAHeatmap);
  };

  ~PuLSEHTMLWriter() {
    file << "</body></html>\n";
    file.close();
  };

private:
  void insertTableRow(std::vector<std::vector<std::string>> &vec,
                      const std::string &in1, const std::string &in2) {
    vec.push_back(std::vector<std::string>{in1, in2});
  };

  template <typename T>
  void insertTableRow(std::vector<std::vector<std::string>> &vec,
                      const std::string &in1, const T &in2) {
    vec.push_back(std::vector<std::string>{in1, std::to_string(in2)});
  };

  template <typename T1, typename T2, typename T3, typename T4>
  void insertTableRow(std::vector<std::vector<std::string>> &vec, const T1 &in1,
                      const T2 &in2, const T3 &in3, const T4 &in4) {
    vec.push_back(
        std::vector<std::string>{std::to_string(in1), std::to_string(in2),
                                 std::to_string(in3), std::to_string(in4)});
  };

  template <typename T1, typename T2, typename T3>
  void insertTableRow(std::vector<std::vector<std::string>> &vec, const T1 &in1,
                      const T2 &in2, const T3 &in3, const std::string &in4) {
    vec.push_back(std::vector<std::string>{
        std::to_string(in1), std::to_string(in2), std::to_string(in3), in4});
  };

  void insertTableRow(std::vector<std::vector<std::string>> &vec,
                      const std::string &in1, const std::string &in2,
                      const std::string &in3, const std::string &in4) {
    vec.push_back(std::vector<std::string>{in1, in2, in3, in4});
  };

  void WriteHTMLTable(std::ofstream &out, const std::string &title,
                      const std::string &subtitle,
                      const std::vector<std::string> &headings,
                      const std::vector<std::vector<std::string>> &data) {
    out << "<div class = \"container\">"
           "<h2>"
        << title
        << "</h2>"
           "<p>"
        << subtitle
        << "</p>"
           "<table class = \"table table-bordered\">"
           "<thead><tr>";
    for (auto &&i : headings) {
      out << "<th>" << i << "</th>";
    }
    out << "</tr></thead><tbody>";
    for (auto &&i : data) {
      out << "<tr>";
      for (auto &&j : i) {
        out << "<td>" << j << "</td>";
      }
      out << "</tr>";
    }
    out << "</tbody></table></div>";
  };

  void
  WriteNormalisedHTMLTable(std::ofstream &out, const std::string &title,
                           const std::string &subtitle,
                           const std::vector<std::string> &headings,
                           const std::vector<std::vector<std::string>> &data) {

    auto getIntensities = [](float val, unsigned int &r, unsigned int &g,
                             unsigned int &b) {
      constexpr float upperbound = 2.0f;
      constexpr float lowerbound = 0.0f;
	  r = 100; g = 100; b = 100;
      unsigned dif;
      if (val < 1.0f) {
        val = (val < lowerbound ? 0.0f : val);
        dif = static_cast<int>(val * 100);
		r -= dif;
		g -= dif; return;
      } else {
        val = (val > upperbound ? 2.0f : val);
		dif = static_cast<int>((val - 1.0f) * 100);
		g -= 100;
		b -= 100;
		return;
      }
	  return;
    };

    auto getHexColor = [&getIntensities](const std::string &in) {
      unsigned int r, g, b;
      std::stringstream ss;
      std::string retval = "#000000";
      if (in.find(".") == std::string::npos) {
        return std::string("#FFFFFF");
      }
      float val = std::stof(in);
      getIntensities(val, r, g, b);
	  r = static_cast<float>((r / 100.0f)*255.0f);
	  g = static_cast<float>((g / 100.0f)*255.0f);
	  b = static_cast<float>((b / 100.0f)*255.0f);

        ss << std::hex << static_cast<int>(r);
		ss << std::hex << static_cast<int>(g);
		ss << std::hex << static_cast<int>(b);

        retval = ss.str();
      
      std::transform(retval.begin(), retval.end(), retval.begin(), toupper);
      return retval;
    };
    out << "<div class = \"container\">"
           "<h2>"
        << title
        << "</h2>"
           "<p>"
        << subtitle
        << "</p>"
           "<table class = \"table table-bordered\">"
           "<thead><tr>";
    for (auto &&i : headings) {
      out << "<th>" << i << "</th>";
    }
    out << "</tr></thead><tbody>";
    for (auto &&i : data) {
      out << "<tr>";
      for (auto &&j : i) {
        out << "<td bgcolor=\"" << getHexColor(j) << "\">" << j << "</td>";
      }
      out << "</tr>";
    }
    out << "</tbody></table></div>";
  };

  void Write2HTMLTables(std::ofstream &out, const std::string &title,
                        const std::string &subtitle,
                        const std::vector<std::string> &headings1,
                        const std::vector<std::vector<std::string>> &data1,
                        const std::vector<std::string> &headings2,
                        const std::vector<std::vector<std::string>> &data2) {

    out << "<div class = \"container\"><h2>" << title << "</h2><p>" << subtitle
        << "</p>"
           "<div class = \"col-xs-6\">"
           "<div class=\"table-responsive\"><table class = \"table "
           "table-bordered\">"
           "<thead><tr>";
    for (auto &&i : headings1) {
      out << "<th>" << i << "</th>";
    }
    out << "</tr></thead><tbody>";
    for (auto &&i : data1) {
      out << "<tr>";
      for (auto &&j : i) {
        out << "<td>" << j << "</td>";
      }
      out << "</tr>";
    }
    out << "</tbody></table></div></div>"

        << "<div class = \"col-xs-6\"><div class=\"table-responsive\"><table "
           "class = \"table table-bordered\">"
           "<thead><tr>";
    for (auto &&i : headings2) {
      out << "<th>" << i << "</th>";
    }
    out << "</tr></thead><tbody>";
    for (auto &&i : data2) {
      out << "<tr>";
      for (auto &&j : i) {
        out << "<td>" << j << "</td>";
      }
      out << "</tr>";
    }
    out << "</tbody></table></div></div>"
           "</div>";
  };
};
