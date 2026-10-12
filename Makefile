# Builds build/orbitsi from src/ and the Escape orbit-counting library in extern/Escape.
#   make                      build
#   make BUILD_DIR=<dir>      build into another directory
#   make GCC_DIR=<dir>        use another GCC installation (C++17, OpenMP)
#   make test                 run the pytest suite for graphlet sizes 3, 4 and 5

# The compiler is given by path, so that no environment module has to be loaded. The -L makes
# the link use this GCC's libstdc++ and libgomp instead of the system's older ones, and the
# rpath lets the binary find them at run time.
GCC_DIR ?= /opt/gridware/depots/54e7fb3c/el8/pkg/compilers/gcc/14.1.0
CXX = $(GCC_DIR)/bin/g++
LDFLAGS = -L$(GCC_DIR)/lib64 -Wl,-rpath,$(GCC_DIR)/lib64

# -march=native: the binary is tuned to, and only runs on, the CPU family it is built on
#   (the search has an AVX2 code path).
# Never add -ffast-math or -Ofast: orbit counts are integers held in doubles and are compared
#   exactly.
# -flto=auto: link-time optimisation across src/ and the Escape sources.
# -gdwarf-4: older linkers (binutils 2.30) cannot read GCC 14's default DWARF 5.
BASE_CXXFLAGS = -std=c++17 -O3 -g -gdwarf-4 -fopenmp -march=native -flto=auto

APP_CXXFLAGS = $(BASE_CXXFLAGS) -Wall


EXECUTABLE = orbitsi

BUILD_DIR = build

SRC_DIR = src
INCLUDE_DIR = include
ESCAPE_DIR = extern/Escape

# -isystem: Escape's headers are included without their warnings.
INCLUDES = -I$(INCLUDE_DIR) -isystem extern

SRC_SOURCES    = $(wildcard $(SRC_DIR)/*.cpp)
ESCAPE_SOURCES = $(wildcard $(ESCAPE_DIR)/*.cpp)

SRC_OBJS    = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/src/%.o, $(SRC_SOURCES))
ESCAPE_OBJS = $(patsubst $(ESCAPE_DIR)/%.cpp, $(BUILD_DIR)/escape/%.o, $(ESCAPE_SOURCES))

ALL_OBJS = $(SRC_OBJS) $(ESCAPE_OBJS)

TARGET = $(BUILD_DIR)/$(EXECUTABLE)



all: $(TARGET)

$(TARGET): $(ALL_OBJS)
	@echo "Linking executable..."
	@mkdir -p $(@D)
	$(CXX) $(APP_CXXFLAGS) $(LDFLAGS) -o $@ $^
	@echo "Build successful! Executable is at: $(TARGET)"

$(BUILD_DIR)/src/%.o: $(SRC_DIR)/%.cpp $(wildcard $(INCLUDE_DIR)/*.h)
	@echo "Compiling [app]   : $<"
	@mkdir -p $(@D)
	$(CXX) $(APP_CXXFLAGS) $(INCLUDES) -c $< -o $@

# Escape is third-party code: compiled with the same optimisation flags, warnings off.
$(BUILD_DIR)/escape/%.o: $(ESCAPE_DIR)/%.cpp
	@echo "Compiling [escape]: $<"
	@mkdir -p $(@D)
	$(CXX) $(BASE_CXXFLAGS) -w $(INCLUDES) -c $< -o $@



test3: all
	@echo "\n--- Running Test Suite (pytest, graphlet-size=3) ---"
	GRAPHLET_SIZE=3 pytest test/

test4: all
	@echo "\n--- Running Test Suite (pytest, graphlet-size=4) ---"
	GRAPHLET_SIZE=4 pytest test/

test5: all
	@echo "\n--- Running Test Suite (pytest, graphlet-size=5) ---"
	GRAPHLET_SIZE=5 pytest test/

test: test3 test4 test5
	@echo "\n--- All Test Suites Passed ---"

# Removes the whole build directory, including any binaries saved in it.
clean:
	@echo "Cleaning up build artifacts..."
	rm -rf $(BUILD_DIR)

.PHONY: all clean test

