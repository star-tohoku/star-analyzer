#ifndef CANDIDATE_BDT_SCORER_H
#define CANDIDATE_BDT_SCORER_H

#include "CandidateSource.h"
#include "CandidateTreeEnsemble.h"
#include "FlatYamlConfig.h"

#include "TMD5.h"
#include "TString.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// In-Maker BDT scoring for StCandidateFemtoMaker (rules: docs/ai/AGENT_RULES.md "Candidate BDT";
// procedure: docs/ai/skills/candidate-bdt.md).
//
// A model is one "<name>.model.yaml" file (written by common/macro/datadriven_bdt, action "model")
// plus the two tree-ensemble files (one per fold, CandidateTreeEnsemble.h) it points to. The Maker gets models through the YAML of the optional
// mainconf key "candidateBdt" (no key = no BDT, every score axis has one bin):
//   speciesKeys: k0s, xi_minus
//   outsideTrainingScore: -2               # score of candidates outside the training domain
//   <s>_model: bdt/<anaName>/<name>.model.yaml   # relative to config/, like every mainconf reference
//   <s>_scoreEdges: -3, -1, 0, 0.1, 0.2, 1 # score axis of the pair histograms; the first bin holds
//                                          # outsideTrainingScore only, the others cover [-1, 1]
//
// The training/application contract, enforced here:
// - Variables. Model expressions see exactly the branches of the candidate trees the Maker writes
//   (CandidateModelVariableNames: runId, eventId, cent9, mass, px, py, pz and the source variables),
//   filled from the same values (FillCandidateModelValues), so a model trained on those trees sees
//   identical numbers in the Maker. Derived quantities are expressions (sqrt(px*px+py*py)).
//   Names must be unique: a source variable may not be called like an event variable.
//   The trees also carry one flag per channel, sharedTrackRemoved_<channel> (the channel's shared-track rule
//   removed the candidate). The flags are not model variables: a model's selections define where the Maker
//   scores (the same score is used by every channel), so they must not depend on one channel's decision.
// - Inputs. The input expressions are read from the tree-ensemble files (records "input <i> <expression>"),
//   never repeated in a YAML file. Both folds must have the same inputs.
// - Folds. Only foldBy: runParity can be applied in the Maker: fold f was trained on runs with
//   runId % 2 == f, and a candidate of a run with parity p is scored with fold 1 - p.
// - Domain. A candidate that fails the model's eventSelection or preselection (the selections the
//   training sample was made with) gets outsideTrainingScore instead of a BDT score.
// Expressions: + - * / ; < <= > >= == != ; && || ! ; log10 log exp sqrt abs max min; numbers. Edge values
// as in ROOT 5.34 TTreeFormula: log/log10 of x <= 0 = 0, sqrt(|x|), exp(x) = exp(700) above 700 and 0 below -700, x/0 = 0.
// && and || may not be mixed without parentheses (TTreeFormula precedence differs between ROOT
// versions), and comparisons do not chain.

// Variable names a model may use, in the branch order of the candidate trees (without "score").
inline std::vector<std::string> CandidateModelVariableNames(const std::vector<std::string>& sourceVariables) {
  const char* event[] = {"runId", "eventId", "cent9", "mass", "px", "py", "pz"};
  std::vector<std::string> names(event, event + 7);
  names.insert(names.end(), sourceVariables.begin(), sourceVariables.end());
  return names;
}
// Index of the first source variable in the list above.
const size_t kCandidateModelFirstSourceVariable = 7;

// Values in the order of CandidateModelVariableNames. Floating values are the Float_t values the
// candidate trees store, so the Maker and a TTreeFormula on the trees evaluate the same numbers.
inline void FillCandidateModelValues(const FemtoSourceCandidate& c, Int_t runId, Int_t eventId, Int_t cent9,
                                     std::vector<Float_t>& floats, std::vector<double>& values) {
  floats.resize(kCandidateModelFirstSourceVariable + c.variables.size());
  floats[0] = static_cast<Float_t>(runId);
  floats[1] = static_cast<Float_t>(eventId);
  floats[2] = static_cast<Float_t>(cent9);
  floats[3] = c.massValue;
  floats[4] = c.px;
  floats[5] = c.py;
  floats[6] = c.pz;
  for (size_t i = 0; i < c.variables.size(); ++i) floats[kCandidateModelFirstSourceVariable + i] = c.variables[i];
  values.assign(floats.begin(), floats.end());
  // integers exactly (Float_t holds run and event numbers only up to 2^24)
  values[0] = runId;
  values[1] = eventId;
  values[2] = cent9;
}

// Tiny expression evaluator (compiled once to RPN).
class CandidateExpression {
public:
  bool Compile(const std::string& text, const std::vector<std::string>& names, std::string& error) {
    mText = text;
    mNames = names;
    mCode.clear();
    mPos = 0;
    mError.clear();
    ParseOr();
    SkipSpace();
    if (mError.empty() && mPos != mText.size()) mError = "unexpected '" + mText.substr(mPos) + "'";
    error = mError;
    return mError.empty();
  }

  // values[i] belongs to names[i] of Compile.
  double Evaluate(const std::vector<double>& values) const {
    std::vector<double> st;
    st.reserve(16);
    for (size_t i = 0; i < mCode.size(); ++i) {
      const Op& o = mCode[i];
      if (o.kind == kNumber) { st.push_back(o.value); continue; }
      if (o.kind == kVariable) { st.push_back(values[o.index]); continue; }
      if (o.kind == kUnary) {
        double& a = st.back();
        // log, sqrt, exp and division follow ROOT 5.34 TFormula/TTreeFormula:
        // log10(x <= 0) = log(x <= 0) = 0, sqrt(x) = sqrt(|x|), exp(x > 700) = exp(700), exp(x < -700) = 0, x / 0 = 0.
        if (o.index == 0) a = -a;
        else if (o.index == 1) a = (a > 0.) ? std::log10(a) : 0.;
        else if (o.index == 2) a = (a > 0.) ? std::log(a) : 0.;
        else if (o.index == 3) a = (a > 700.) ? std::exp(700.) : ((a < -700.) ? 0. : std::exp(a));
        else if (o.index == 4) a = std::sqrt(std::fabs(a));
        else if (o.index == 5) a = std::fabs(a);
        else a = (a == 0.) ? 1. : 0.;
        continue;
      }
      const double b = st.back();
      st.pop_back();
      double& a = st.back();
      switch (o.index) {
        case '+': a = a + b; break;
        case '-': a = a - b; break;
        case '*': a = a * b; break;
        case '/': a = (b == 0.) ? 0. : a / b; break;
        case 'M': a = a > b ? a : b; break;
        case 'm': a = a < b ? a : b; break;
        case '<': a = (a < b) ? 1. : 0.; break;
        case 'l': a = (a <= b) ? 1. : 0.; break;
        case '>': a = (a > b) ? 1. : 0.; break;
        case 'g': a = (a >= b) ? 1. : 0.; break;
        case '=': a = (a == b) ? 1. : 0.; break;
        case 'n': a = (a != b) ? 1. : 0.; break;
        case '&': a = (a != 0. && b != 0.) ? 1. : 0.; break;
        default: a = (a != 0. || b != 0.) ? 1. : 0.; break;
      }
    }
    return st.empty() ? 0. : st.back();
  }

  const std::string& Text() const { return mText; }

private:
  enum Kind { kNumber, kVariable, kUnary, kBinary };
  struct Op {
    Kind kind;
    double value;
    int index;
  };
  void Emit(Kind k, double v, int i) { Op o; o.kind = k; o.value = v; o.index = i; mCode.push_back(o); }
  void SkipSpace() { while (mPos < mText.size() && std::isspace(static_cast<unsigned char>(mText[mPos]))) ++mPos; }
  bool Accept(char c) {
    SkipSpace();
    if (mPos < mText.size() && mText[mPos] == c) { ++mPos; return true; }
    return false;
  }
  // Two-character token (&&, ||, <=, >=, ==, !=).
  bool Accept2(const char* t) {
    SkipSpace();
    if (mPos + 1 < mText.size() && mText[mPos] == t[0] && mText[mPos + 1] == t[1]) { mPos += 2; return true; }
    return false;
  }
  bool Peek(char c) {
    SkipSpace();
    return mPos < mText.size() && mText[mPos] == c;
  }
  // or := and ('||' and)* ; returns nothing. Mixing && and || at one level needs parentheses.
  void ParseOr() {
    const int andsFirst = ParseAnd();
    Int_t nOr = 0, nAnd = andsFirst;
    while (mError.empty() && Accept2("||")) {
      nAnd += ParseAnd();
      Emit(kBinary, 0., '|');
      ++nOr;
    }
    if (mError.empty() && nOr > 0 && nAnd > 0) mError = "mixed && and || without parentheses";
  }
  int ParseAnd() {
    ParseCompare();
    int n = 0;
    while (mError.empty() && Accept2("&&")) {
      ParseCompare();
      Emit(kBinary, 0., '&');
      ++n;
    }
    return n;
  }
  void ParseCompare() {
    ParseSum();
    if (!mError.empty()) return;
    char op = 0;
    if (Accept2("<=")) op = 'l';
    else if (Accept2(">=")) op = 'g';
    else if (Accept2("==")) op = '=';
    else if (Accept2("!=")) op = 'n';
    else if (Accept('<')) op = '<';
    else if (Accept('>')) op = '>';
    if (!op) return;
    ParseSum();
    Emit(kBinary, 0., op);
    if (mError.empty() && (Peek('<') || Peek('>') || Peek('=') || (Peek('!') && mPos + 1 < mText.size() && mText[mPos + 1] == '=')))
      mError = "chained comparison (use &&)";
  }
  void ParseSum() {
    ParseProduct();
    while (mError.empty()) {
      if (Accept('+')) { ParseProduct(); Emit(kBinary, 0., '+'); }
      else if (Accept('-')) { ParseProduct(); Emit(kBinary, 0., '-'); }
      else break;
    }
  }
  void ParseProduct() {
    ParseUnary();
    while (mError.empty()) {
      if (Accept('*')) { ParseUnary(); Emit(kBinary, 0., '*'); }
      else if (Accept('/')) { ParseUnary(); Emit(kBinary, 0., '/'); }
      else break;
    }
  }
  void ParseUnary() {
    if (Accept('-')) { ParseUnary(); Emit(kUnary, 0., 0); return; }
    if (Accept('+')) { ParseUnary(); return; }
    SkipSpace();
    if (mPos + 1 < mText.size() && mText[mPos] == '!' && mText[mPos + 1] != '=') {
      ++mPos;
      ParseUnary();
      Emit(kUnary, 0., 6);
      return;
    }
    ParsePrimary();
  }
  void ParsePrimary() {
    SkipSpace();
    if (mPos >= mText.size()) { mError = "unexpected end of expression"; return; }
    if (Accept('(')) {
      ParseOr();
      if (mError.empty() && !Accept(')')) mError = "missing ')'";
      return;
    }
    const char c = mText[mPos];
    if (std::isdigit(static_cast<unsigned char>(c)) || c == '.') {
      const char* begin = mText.c_str() + mPos;
      char* end = 0;
      const double v = std::strtod(begin, &end);
      mPos += end - begin;
      Emit(kNumber, v, 0);
      return;
    }
    if (!(std::isalpha(static_cast<unsigned char>(c)) || c == '_')) {
      mError = std::string("unexpected character '") + c + "'";
      return;
    }
    std::string id;
    while (mPos < mText.size() &&
           (std::isalnum(static_cast<unsigned char>(mText[mPos])) || mText[mPos] == '_'))
      id += mText[mPos++];
    if (Accept('(')) {
      const char* unary[] = {"", "log10", "log", "exp", "sqrt", "abs"};
      for (int u = 1; u < 6; ++u)
        if (id == unary[u]) {
          ParseOr();
          if (mError.empty() && !Accept(')')) mError = "missing ')' after " + id;
          Emit(kUnary, 0., u);
          return;
        }
      if (id == "max" || id == "min") {
        ParseOr();
        if (mError.empty() && !Accept(',')) { mError = id + " needs two arguments"; return; }
        ParseOr();
        if (mError.empty() && !Accept(')')) mError = "missing ')' after " + id;
        Emit(kBinary, 0., id == "max" ? 'M' : 'm');
        return;
      }
      mError = "unknown function '" + id + "'";
      return;
    }
    for (size_t i = 0; i < mNames.size(); ++i)
      if (mNames[i] == id) { Emit(kVariable, 0., static_cast<int>(i)); return; }
    mError = "unknown variable '" + id + "'";
  }

  std::string mText;
  std::vector<std::string> mNames;
  std::vector<Op> mCode;
  size_t mPos;
  std::string mError;
};

// One model file: "<name>.model.yaml" (flat YAML, every key required):
//   name, particle, foldBy (runParity), weights (pattern "<name>_f%d.trees.txt" with %d for the fold,
//   relative to the directory of the model file), eventSelection, preselection ("1" = none), trainingData,
//   trainingCandidates, recipe, gitRevision, learner, validation, created.
struct CandidateBdtModel {
  std::string path, name, particle, foldBy, weightsPattern, eventSelection, preselection;
  std::string trainingData, trainingCandidates, recipe, gitRevision, learner, validation, created;
  std::string weightFile[2], weightMd5[2];
  CandidateTreeEnsemble ensemble[2];
  std::vector<std::string> expressions;  // the input expressions (identical in both folds)

  bool Load(const std::string& modelPath) {
    path = modelPath;
    FlatYamlConfig c;
    if (!c.Load(modelPath)) return false;
    name = c.RequireString("name");
    particle = c.RequireString("particle");
    foldBy = c.RequireString("foldBy");
    const std::string weights = c.RequireString("weights");
    eventSelection = c.RequireString("eventSelection");
    preselection = c.RequireString("preselection");
    trainingData = c.RequireString("trainingData");
    trainingCandidates = c.RequireString("trainingCandidates");
    recipe = c.RequireString("recipe");
    gitRevision = c.RequireString("gitRevision");
    learner = c.RequireString("learner");
    validation = c.RequireString("validation");
    created = c.RequireString("created");
    c.RejectUnused();
    if (!c.Ok()) return false;
    if (foldBy != "runParity") {
      c.Fail("foldBy is '" + foldBy + "': only runParity models can be applied in the Maker (an event has a run "
             "number but no job-list index); retrain the model with fold_by = runParity");
      return false;
    }
    const std::string::size_type slash = modelPath.rfind('/');
    const std::string dir = (slash == std::string::npos) ? std::string(".") : modelPath.substr(0, slash);
    weightsPattern = (!weights.empty() && weights[0] == '/') ? weights : dir + "/" + weights;
    const std::string suffix = ".trees.txt";
    if (weightsPattern.find("%d") == std::string::npos || weightsPattern.size() < suffix.size() ||
        weightsPattern.compare(weightsPattern.size() - suffix.size(), suffix.size(), suffix) != 0) {
      c.Fail("weights must be a tree-ensemble pattern with %d for the fold, ending in .trees.txt: " + weights +
             " (models are written by common/macro/datadriven_bdt, action model)");
      return false;
    }
    for (Int_t fold = 0; fold < 2; ++fold) {
      weightFile[fold] = TString::Format(weightsPattern.c_str(), fold).Data();
      std::string error;
      if (!ensemble[fold].Load(weightFile[fold], error)) {
        std::cerr << "ERROR: " << path << ": " << error << std::endl;
        return false;
      }
      if (fold == 0) expressions = ensemble[0].Inputs();
      else if (ensemble[1].Inputs() != expressions) {
        std::cerr << "ERROR: " << path << ": the two folds have different inputs" << std::endl;
        return false;
      }
      TMD5* md5 = TMD5::FileChecksum(weightFile[fold].c_str());
      weightMd5[fold] = md5 ? md5->AsString() : "unreadable";
      delete md5;
    }
    return true;
  }

  // Provenance text written into the Maker output.
  std::string Describe() const {
    std::ostringstream o;
    o << "model " << path << "\n  name " << name << ", particle " << particle << ", foldBy " << foldBy
      << "\n  eventSelection " << eventSelection << "\n  preselection " << preselection
      << "\n  trainingData " << trainingData << "\n  trainingCandidates " << trainingCandidates
      << "\n  recipe " << recipe << ", gitRevision " << gitRevision << ", created " << created
      << "\n  learner " << learner << "\n  validation " << validation;
    for (Int_t f = 0; f < 2; ++f)
      o << "\n  fold " << f << " " << weightFile[f] << " md5 " << weightMd5[f] << ", " << ensemble[f].NumberOfTrees()
        << " trees (" << ensemble[f].Learner() << ")";
    for (size_t i = 0; i < expressions.size(); ++i) o << "\n  input " << i << " " << expressions[i];
    return o.str();
  }

};

class CandidateBdtScorer {
public:
  CandidateBdtScorer() : mOutside(-2.f), mOutsideCount(0) {}

  // modelVariableNames: CandidateModelVariableNames(source variables) of this species.
  bool Init(const std::string& speciesKey, const std::vector<std::string>& modelVariableNames,
            const std::string& modelPath, Float_t outsideTrainingScore) {
    mKey = speciesKey;
    mOutside = outsideTrainingScore;
    if (!mModel.Load(modelPath)) return false;
    std::string error;
    if (!mEventSelection.Compile(mModel.eventSelection, modelVariableNames, error) ||
        !mPreselection.Compile(mModel.preselection, modelVariableNames, error)) {
      std::cerr << "ERROR: " << modelPath << ": selection of " << speciesKey << ": " << error << std::endl;
      if ((mModel.eventSelection + mModel.preselection).find("sharedTrackRemoved_") != std::string::npos)
        std::cerr << "  The candidate-tree flags sharedTrackRemoved_<channel> may select a training sample but cannot be part "
                     "of a model's eventSelection/preselection: those define where the Maker scores, and a score that depends "
                     "on one channel's shared-track decision would couple the channels." << std::endl;
      return false;
    }
    mExpressions.assign(mModel.expressions.size(), CandidateExpression());
    mInputs.assign(mModel.expressions.size(), 0.f);
    for (size_t i = 0; i < mModel.expressions.size(); ++i)
      if (!mExpressions[i].Compile(mModel.expressions[i], modelVariableNames, error)) {
        std::cerr << "ERROR: " << modelPath << ": input '" << mModel.expressions[i] << "' of " << speciesKey << ": "
                  << error << std::endl;
        if (mModel.expressions[i].find("sharedTrackRemoved_") != std::string::npos)
          std::cerr << "  sharedTrackRemoved_<channel> is a candidate-tree flag, not a candidate property: exclude it from the "
                       "inputs (inputs.exclude of the recipe)." << std::endl;
        return false;
      }
    std::cout << "[CandidateBdtScorer] " << speciesKey << ": " << mModel.Describe() << std::endl;
    return true;
  }

  // values: FillCandidateModelValues. Candidates outside the training domain get the outside score;
  // the others are scored with the fold trained on the other run parity.
  Float_t Score(const std::vector<double>& values, Int_t runId) {
    if (mEventSelection.Evaluate(values) == 0. || mPreselection.Evaluate(values) == 0.) {
      ++mOutsideCount;
      return mOutside;
    }
    for (size_t i = 0; i < mExpressions.size(); ++i) mInputs[i] = static_cast<Float_t>(mExpressions[i].Evaluate(values));
    const Int_t fold = (runId % 2 == 0) ? 1 : 0;
    return static_cast<Float_t>(mModel.ensemble[fold].Score(mInputs));
  }

  const CandidateBdtModel& Model() const { return mModel; }
  Long64_t OutsideCount() const { return mOutsideCount; }

private:
  std::string mKey;
  CandidateBdtModel mModel;
  CandidateExpression mEventSelection, mPreselection;
  std::vector<CandidateExpression> mExpressions;
  std::vector<Float_t> mInputs;
  Float_t mOutside;
  Long64_t mOutsideCount;
};

#endif
