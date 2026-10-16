#ifndef CANDIDATE_TREE_ENSEMBLE_H
#define CANDIDATE_TREE_ENSEMBLE_H

#include "Rtypes.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// Tree-ensemble BDT read by StCandidateFemtoMaker: one file per fold, written by common/macro/datadriven_bdt
// (Python + xgboost; tree_ensemble.py is the reference implementation of this file format and of the score).
//
// File (text; empty lines and lines starting with '#' are ignored; every record is required, in this order):
//   format candidateTreeEnsemble 1
//   learner <free text>
//   baseMargin <double>
//   inputs <n>
//   input <i> <expression>            input i, an expression over the candidate-tree variables (CandidateExpression)
//   trees <nTrees>
//   tree <t> <nNodes>
//   node <id> <input> <threshold> <yes> <no> <missing> <leaf>
//       split node: input >= 0; value < threshold (float comparison, as xgboost) -> yes, otherwise no, NaN -> missing
//       leaf node:  input = -1, yes = no = missing = -1, leaf value
//   end
// Score (the Maker contract): margin = sum of the leaves (double, tree order) + baseMargin,
// score = 2 / (1 + exp(-margin)) - 1, in [-1, 1].
class CandidateTreeEnsemble {
public:
  CandidateTreeEnsemble() : mBaseMargin(0.) {}

  bool Load(const std::string& file, std::string& error) {
    mInputs.clear();
    mTrees.clear();
    std::ifstream in(file.c_str());
    if (!in) { error = "cannot read " + file; return false; }
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
      const size_t b = line.find_first_not_of(" \t\r");
      if (b == std::string::npos || line[b] == '#') continue;
      const size_t e = line.find_last_not_of(" \t\r");
      lines.push_back(line.substr(b, e - b + 1));
    }
    size_t k = 0;
    std::string rest;
    if (k >= lines.size() || lines[k++] != "format candidateTreeEnsemble 1") { error = file + ": not a candidateTreeEnsemble 1 file"; return false; }
    if (!Expect(lines, k, "learner ", mLearner, error, file)) return false;
    if (!Expect(lines, k, "baseMargin ", rest, error, file)) return false;
    mBaseMargin = std::strtod(rest.c_str(), 0);
    if (!Expect(lines, k, "inputs ", rest, error, file)) return false;
    const int nIn = std::atoi(rest.c_str());
    for (int i = 0; i < nIn; ++i) {
      if (!Expect(lines, k, "input ", rest, error, file)) return false;
      const size_t sp = rest.find(' ');
      if (sp == std::string::npos || std::atoi(rest.substr(0, sp).c_str()) != i) { error = file + ": input records out of order"; return false; }
      mInputs.push_back(rest.substr(rest.find_first_not_of(' ', sp)));
    }
    if (!Expect(lines, k, "trees ", rest, error, file)) return false;
    const int nTrees = std::atoi(rest.c_str());
    for (int t = 0; t < nTrees; ++t) {
      if (!Expect(lines, k, "tree ", rest, error, file)) return false;
      std::istringstream ts(rest);
      int id = -1, nNodes = 0;
      ts >> id >> nNodes;
      if (id != t || nNodes <= 0) { error = file + ": tree records out of order"; return false; }
      std::vector<Node> nodes(nNodes);
      for (int j = 0; j < nNodes; ++j) {
        if (!Expect(lines, k, "node ", rest, error, file)) return false;
        std::istringstream ns(rest);
        std::string thr, leaf;
        int nid = -1;
        Node& n = nodes[j];
        ns >> nid >> n.input >> thr >> n.yes >> n.no >> n.missing >> leaf;
        if (ns.fail() || nid != j) { error = file + ": bad node record: " + rest; return false; }
        n.threshold = static_cast<float>(std::strtod(thr.c_str(), 0));
        n.leaf = static_cast<double>(static_cast<float>(std::strtod(leaf.c_str(), 0)));
        if (n.input >= nIn || (n.input >= 0 && (n.yes <= j || n.no <= j || n.missing <= j || n.yes >= nNodes || n.no >= nNodes ||
                                                n.missing >= nNodes))) {
          error = file + ": node out of range: " + rest;
          return false;
        }
      }
      mTrees.push_back(nodes);
    }
    if (k >= lines.size() || lines[k++] != "end" || k != lines.size()) { error = file + ": missing or misplaced \"end\""; return false; }
    if (mInputs.empty() || mTrees.empty()) { error = file + ": no inputs or no trees"; return false; }
    return true;
  }

  const std::vector<std::string>& Inputs() const { return mInputs; }
  const std::string& Learner() const { return mLearner; }
  size_t NumberOfTrees() const { return mTrees.size(); }

  // x: the input values (float, one per input, in the order of Inputs()).
  double Margin(const std::vector<Float_t>& x) const {
    double sum = 0.;
    for (size_t t = 0; t < mTrees.size(); ++t) {
      const std::vector<Node>& nodes = mTrees[t];
      int i = 0;
      while (nodes[i].input >= 0) {
        const float v = x[nodes[i].input];
        i = (v != v) ? nodes[i].missing : ((v < nodes[i].threshold) ? nodes[i].yes : nodes[i].no);
      }
      sum += nodes[i].leaf;
    }
    return sum + mBaseMargin;
  }

  double Score(const std::vector<Float_t>& x) const { return 2. / (1. + std::exp(-Margin(x))) - 1.; }

private:
  struct Node {
    int input;
    float threshold;
    int yes, no, missing;
    double leaf;
    Node() : input(-1), threshold(0.f), yes(-1), no(-1), missing(-1), leaf(0.) {}
  };
  static bool Expect(const std::vector<std::string>& lines, size_t& k, const std::string& prefix, std::string& rest,
                     std::string& error, const std::string& file) {
    if (k >= lines.size() || lines[k].compare(0, prefix.size(), prefix) != 0) {
      error = file + ": expected \"" + prefix + "...\"" + (k < lines.size() ? ", got \"" + lines[k] + "\"" : " at the end");
      return false;
    }
    rest = lines[k++].substr(prefix.size());
    return true;
  }

  std::string mLearner;
  double mBaseMargin;
  std::vector<std::string> mInputs;
  std::vector<std::vector<Node> > mTrees;
};

#endif
