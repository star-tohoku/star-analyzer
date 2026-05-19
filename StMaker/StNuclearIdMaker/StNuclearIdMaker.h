#ifndef StNuclearIdMaker_h
#define StNuclearIdMaker_h

#include "StMaker.h"
#include "TString.h"

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StPicoBTofPidTraits;
class HistManager;

class StNuclearIdMaker : public StMaker {
 public:
  StNuclearIdMaker(const char* name, StPicoDstMaker* picoMaker, const char* outName);
  virtual ~StNuclearIdMaker();

  virtual Int_t Init();
  virtual Int_t Make();
  virtual void Clear(Option_t* opt = "");
  virtual Int_t Finish();

  void WriteHistograms();

 private:
  StPicoDstMaker* mPicoDstMaker;
  StPicoDst* mPicoDst;
  TString mOutName;
  HistManager* m_histManager;
};

// Factory for CINT
StNuclearIdMaker* createStNuclearIdMaker(const char* name, StPicoDstMaker* picoMaker, const char* outName);
extern "C" void* createStNuclearIdMakerC(const char* name, void* picoMaker, const char* outName);

#endif
