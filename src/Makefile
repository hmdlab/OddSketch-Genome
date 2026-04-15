CXX := clang++
CXXFLAGS := -std=c++17 -O3 -march=native -Wall -Wextra
CPPFLAGS := -I../include

ODDSKETCH_BIN := oddsketch
TRUE_JACCARD_BIN := cal/true_jaccard
TRUE_INDEX_PAIRS_BIN := cal/true_index_pairs

ODDSKETCH_SRCS := apps/oddsketch_main.cpp core/oddsketch.cpp third_party/xxhash.cpp
TRUE_JACCARD_SRCS := apps/true_jaccard_main.cpp
TRUE_INDEX_PAIRS_SRCS := apps/true_index_pairs_main.cpp

.PHONY: all clean

all: $(ODDSKETCH_BIN) $(TRUE_JACCARD_BIN) $(TRUE_INDEX_PAIRS_BIN)

$(ODDSKETCH_BIN): $(ODDSKETCH_SRCS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

$(TRUE_JACCARD_BIN): $(TRUE_JACCARD_SRCS) | cal
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

$(TRUE_INDEX_PAIRS_BIN): $(TRUE_INDEX_PAIRS_SRCS) | cal
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

cal:
	mkdir -p cal

clean:
	rm -f $(ODDSKETCH_BIN) $(TRUE_JACCARD_BIN) $(TRUE_INDEX_PAIRS_BIN)
