#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
// Compiled replay binds every persisted STL branch through its RTTI dictionary.
#pragma link C++ class std::vector<float>+;
#pragma link C++ class std::vector<double>+;
#pragma link C++ class std::vector<int>+;
#pragma link C++ class std::vector<unsigned int>+;
#pragma link C++ class std::vector<TVector3>+;
#endif
