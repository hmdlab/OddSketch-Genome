CXX ?= clang++
CXXFLAGS := -std=c++17 -O3 -march=native -Wall -Wextra
CPPFLAGS := -I../include
LDFLAGS ?=
THREAD_LDFLAGS ?= -pthread

ODDSKETCH_BIN := oddsketch
TOOLS_BIN_DIR := ../experiments/tools/bin
TRUE_JACCARD_BIN := $(TOOLS_BIN_DIR)/true_jaccard
TRUE_INDEX_PAIRS_BIN := $(TOOLS_BIN_DIR)/true_index_pairs

ODDSKETCH_SRCS := apps/oddsketch_main.cpp core/oddsketch.cpp third_party/xxhash.cpp
TRUE_JACCARD_SRCS := ../experiments/tools/src/true_jaccard_main.cpp
TRUE_INDEX_PAIRS_SRCS := ../experiments/tools/src/true_index_pairs_main.cpp

.PHONY: all clean

all: $(ODDSKETCH_BIN) $(TRUE_JACCARD_BIN) $(TRUE_INDEX_PAIRS_BIN)

$(ODDSKETCH_BIN): $(ODDSKETCH_SRCS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(THREAD_LDFLAGS)

$(TRUE_JACCARD_BIN): $(TRUE_JACCARD_SRCS) | $(TOOLS_BIN_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(THREAD_LDFLAGS)

$(TRUE_INDEX_PAIRS_BIN): $(TRUE_INDEX_PAIRS_SRCS) | $(TOOLS_BIN_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(THREAD_LDFLAGS)

$(TOOLS_BIN_DIR):
	mkdir -p $(TOOLS_BIN_DIR)

clean:
	rm -f $(ODDSKETCH_BIN) $(TRUE_JACCARD_BIN) $(TRUE_INDEX_PAIRS_BIN)
