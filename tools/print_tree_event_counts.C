// Print "<events> <path>" for every reduced-tree ROOT file named in a list file. Used by
// script/submit_downstream_blocks.py to split a block into equal event slices. ROOT only, so it
// runs against the container's own ROOT and needs no STAR release.
#include <fstream>
#include <string>
#include "TFile.h"
#include "TTree.h"

void print_tree_event_counts(const char* listFile) {
  std::ifstream in(listFile);
  std::string path;
  while (std::getline(in, path)) {
    if (path.empty()) continue;
    TFile* f = TFile::Open(path.c_str());
    TTree* ev = f ? (TTree*)f->Get("FemtoEventTree") : 0;
    printf("COUNT %lld %s\n", ev ? ev->GetEntries() : -1LL, path.c_str());
    if (f) { f->Close(); delete f; }
  }
}
