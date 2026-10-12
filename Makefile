

GCC_DIR ?= /opt/gridware/depots/54e7fb3c/el8/pkg/compilers/gcc/14.1.0
CXX = $(GCC_DIR)/bin/g++
LDFLAGS = -L$(GCC_DIR)/lib64 -Wl,-rpath,$(GCC_DIR)/lib64

BASE_CXXFLAGS = -std=c++17 -O3 -g -gdwarf-4 -fopenmp -march=native -flto=auto

APP_CXXFLAGS = $(BASE_CXXFLAGS) -Wall


EXECUTABLE = orbitsi

BUILD_DIR = build

SRC_DIR = src
INCLUDE_DIR = include
ESCAPE_DIR = extern/Escape

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

clean:
	@echo "Cleaning up build artifacts..."
	rm -rf $(BUILD_DIR)

.PHONY: all clean test

