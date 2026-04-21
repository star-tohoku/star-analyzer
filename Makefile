# Makefile for StPhiMaker and StarAnaConfig (STAR analysis framework)
# Requires: starver SL24y (or SL24c), $STAR set
# Usage: ./script/setup.sh config/mainconf/main_<anaName>.yaml && make

ifeq ($(STAR),)
  $(error STAR environment variable not set. Run: ./script/setup.sh config/mainconf/main_<anaName>.yaml)
endif

# Directories - prefer x86_64 object tree when available
STAR_OBJ := $(STAR)/.sl74_x8664_gcc485
ifeq ($(wildcard $(STAR_OBJ)),)
  STAR_OBJ := $(STAR)/.sl73_x8664_gcc485
endif
ifeq ($(wildcard $(STAR_OBJ)),)
  STAR_OBJ := $(STAR)/.sl74_gcc485
endif
ifeq ($(wildcard $(STAR_OBJ)),)
  STAR_OBJ := $(STAR)/.sl73_gcc485
endif
STAR_INC_DIR := $(STAR_OBJ)/include
STAR_LIB_DIR := $(STAR_OBJ)/lib

STMAKER_DIR := StMaker/StPhiMaker
LIB_DIR := lib
YAML_CPP_DIR := src/third_party/yaml-cpp
YAML_CPP_BUILD := $(YAML_CPP_DIR)/build

# Compiler and flags
CXX := g++
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS := $(shell root-config --libs)

# STAR include and link - need both obj/include and StRoot for headers
STAR_INC := -I$(STAR_INC_DIR) -I$(STAR)/StRoot -I$(STAR)/StRoot/StPicoDstMaker -I$(STAR)/StRoot/StPicoEvent -I$(STAR)/StRoot/StMaker -I$(STAR)/StRoot/StarClassLibrary -I$(STAR)/StRoot/StBTofUtil
STAR_LDFLAGS := -L$(STAR_LIB_DIR) \
                -lStPicoDstMaker -lStPicoEvent \
                -lStarClassLibrary -lSt_base -lStChain -lStUtilities

# --- libStarAnaConfig (ConfigManager + YamlParser + cut configs) ---
STAR_ANA_CONFIG_SRCS := src/ConfigManager.cpp src/YamlParser.cpp src/HistManager.cpp src/kinematics.cpp \
  src/cuts/EventCutConfig.cpp src/cuts/TrackCutConfig.cpp src/cuts/PIDCutConfig.cpp \
  src/cuts/V0CutConfig.cpp src/cuts/PhiCutConfig.cpp src/cuts/LambdaCutConfig.cpp \
  src/cuts/Lambda1520CutConfig.cpp src/cuts/Sigma1385CutConfig.cpp src/cuts/MixingConfig.cpp
STAR_ANA_CONFIG_OBJS := $(addprefix $(LIB_DIR)/,$(notdir $(STAR_ANA_CONFIG_SRCS:.cpp=.o)))
CXXFLAGS_CONFIG := -O2 -Wall -fPIC -std=c++11 $(ROOTCFLAGS) -Iinclude -I$(YAML_CPP_DIR)/include
LDFLAGS_CONFIG := $(ROOTLDFLAGS) -shared -Wl,--whole-archive -L$(YAML_CPP_BUILD) -lyaml-cpp -Wl,--no-whole-archive

# --- libStPhiMaker (depends on libStarAnaConfig) ---
LIB_NAME := libStPhiMaker.so
CXXFLAGS_MAKER := -O2 -Wall -fPIC $(ROOTCFLAGS) -Iinclude $(STAR_INC)
LDFLAGS_MAKER := $(ROOTLDFLAGS) -shared -Wl,-rpath,$(STAR_LIB_DIR)
SRC := $(STMAKER_DIR)/StPhiMaker.cxx
OBJ := $(LIB_DIR)/StPhiMaker.o

# --- libStLambdaMaker (depends on libStarAnaConfig) ---
STLAMBDA_DIR := StMaker/StLambdaMaker
LIB_LAMBDA_NAME := libStLambdaMaker.so
SRC_LAMBDA := $(STLAMBDA_DIR)/StLambdaMaker.cxx
OBJ_LAMBDA := $(LIB_DIR)/StLambdaMaker.o

.PHONY: all clean

all: $(LIB_DIR)/libStarAnaConfig.so $(LIB_DIR)/$(LIB_NAME) $(LIB_DIR)/$(LIB_LAMBDA_NAME)

# Build yaml-cpp via CMake (static lib, must match STAR/ROOT 32-bit)
$(YAML_CPP_BUILD)/libyaml-cpp.a:
	@mkdir -p $(YAML_CPP_BUILD)
	@cd $(YAML_CPP_BUILD) && cmake .. -DYAML_CPP_BUILD_CONTRIB=OFF -DYAML_CPP_BUILD_TOOLS=OFF -DYAML_CPP_BUILD_TESTS=OFF -DYAML_BUILD_SHARED_LIBS=OFF -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DCMAKE_CXX_FLAGS="$(ROOTCFLAGS)" && $(MAKE) yaml-cpp

$(LIB_DIR):
	mkdir -p $(LIB_DIR)

# libStarAnaConfig.so (depends on yaml-cpp)
$(LIB_DIR)/libStarAnaConfig.so: $(LIB_DIR) $(YAML_CPP_BUILD)/libyaml-cpp.a $(STAR_ANA_CONFIG_OBJS)
	$(CXX) $(LDFLAGS_CONFIG) -o $@ $(STAR_ANA_CONFIG_OBJS) $(ROOTLIBS)

$(LIB_DIR)/ConfigManager.o: src/ConfigManager.cpp include/ConfigManager.h
	$(CXX) $(CXXFLAGS_CONFIG) -c src/ConfigManager.cpp -o $@
$(LIB_DIR)/YamlParser.o: src/YamlParser.cpp include/YamlParser.h
	$(CXX) $(CXXFLAGS_CONFIG) -c src/YamlParser.cpp -o $@
$(LIB_DIR)/EventCutConfig.o: src/cuts/EventCutConfig.cpp include/cuts/EventCutConfig.h
	$(CXX) $(CXXFLAGS_CONFIG) -c src/cuts/EventCutConfig.cpp -o $@
$(LIB_DIR)/TrackCutConfig.o: src/cuts/TrackCutConfig.cpp include/cuts/TrackCutConfig.h
	$(CXX) $(CXXFLAGS_CONFIG) -c src/cuts/TrackCutConfig.cpp -o $@
$(LIB_DIR)/PIDCutConfig.o: src/cuts/PIDCutConfig.cpp include/cuts/PIDCutConfig.h
	$(CXX) $(CXXFLAGS_CONFIG) -c src/cuts/PIDCutConfig.cpp -o $@
$(LIB_DIR)/V0CutConfig.o: src/cuts/V0CutConfig.cpp include/cuts/V0CutConfig.h
	$(CXX) $(CXXFLAGS_CONFIG) -c src/cuts/V0CutConfig.cpp -o $@
$(LIB_DIR)/PhiCutConfig.o: src/cuts/PhiCutConfig.cpp include/cuts/PhiCutConfig.h
	$(CXX) $(CXXFLAGS_CONFIG) -c src/cuts/PhiCutConfig.cpp -o $@
$(LIB_DIR)/LambdaCutConfig.o: src/cuts/LambdaCutConfig.cpp include/cuts/LambdaCutConfig.h
	$(CXX) $(CXXFLAGS_CONFIG) -c src/cuts/LambdaCutConfig.cpp -o $@
$(LIB_DIR)/Lambda1520CutConfig.o: src/cuts/Lambda1520CutConfig.cpp include/cuts/Lambda1520CutConfig.h
	$(CXX) $(CXXFLAGS_CONFIG) -c src/cuts/Lambda1520CutConfig.cpp -o $@
$(LIB_DIR)/Sigma1385CutConfig.o: src/cuts/Sigma1385CutConfig.cpp include/cuts/Sigma1385CutConfig.h
	$(CXX) $(CXXFLAGS_CONFIG) -c src/cuts/Sigma1385CutConfig.cpp -o $@
$(LIB_DIR)/MixingConfig.o: src/cuts/MixingConfig.cpp include/cuts/MixingConfig.h
	$(CXX) $(CXXFLAGS_CONFIG) -c src/cuts/MixingConfig.cpp -o $@
$(LIB_DIR)/HistManager.o: src/HistManager.cpp include/HistManager.h
	$(CXX) $(CXXFLAGS_CONFIG) -c src/HistManager.cpp -o $@
$(LIB_DIR)/kinematics.o: src/kinematics.cpp include/kinematics.h
	$(CXX) $(CXXFLAGS_CONFIG) -c src/kinematics.cpp -o $@

# libStPhiMaker.so (links against libStarAnaConfig)
$(LIB_DIR)/$(LIB_NAME): $(LIB_DIR)/libStarAnaConfig.so $(LIB_DIR) $(OBJ)
	$(CXX) $(LDFLAGS_MAKER) -o $@ $(OBJ) -L$(LIB_DIR) -lStarAnaConfig -Wl,-rpath,$(abspath $(LIB_DIR)) $(STAR_LDFLAGS) $(ROOTLIBS)

$(OBJ): $(SRC) $(STMAKER_DIR)/StPhiMaker.h include/HistManager.h
	$(CXX) $(CXXFLAGS_MAKER) -c $(SRC) -o $@

# libStLambdaMaker.so (links against libStarAnaConfig)
$(LIB_DIR)/$(LIB_LAMBDA_NAME): $(LIB_DIR)/libStarAnaConfig.so $(LIB_DIR) $(OBJ_LAMBDA)
	$(CXX) $(LDFLAGS_MAKER) -o $@ $(OBJ_LAMBDA) -L$(LIB_DIR) -lStarAnaConfig -Wl,-rpath,$(abspath $(LIB_DIR)) $(STAR_LDFLAGS) $(ROOTLIBS)

$(OBJ_LAMBDA): $(SRC_LAMBDA) $(STLAMBDA_DIR)/StLambdaMaker.h include/HistManager.h
	$(CXX) $(CXXFLAGS_MAKER) -c $(SRC_LAMBDA) -o $@

clean:
	rm -f $(LIB_DIR)/*.o $(LIB_DIR)/$(LIB_NAME) $(LIB_DIR)/$(LIB_LAMBDA_NAME) $(LIB_DIR)/libStarAnaConfig.so
	rm -rf $(YAML_CPP_BUILD)
