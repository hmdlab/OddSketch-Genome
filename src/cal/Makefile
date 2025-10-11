# Makefile for compiling and running the C++ program

# コンパイラの設定
CXX = g++
CXXFLAGS = -std=c++17 -O3 -march=native -Wall -Wextra -fopenmp

BIN_CAL_JACCARD = cal_Jaccard
BIN_TRUE_JACCARD = true_jaccard

all: $(BIN_CAL_JACCARD) $(BIN_TRUE_JACCARD)

$(BIN_CAL_JACCARD): cal_Jaccard.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

$(BIN_TRUE_JACCARD): true_jaccard.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f $(BIN_CAL_JACCARD) $(BIN_TRUE_JACCARD)

.PHONY: all clean
