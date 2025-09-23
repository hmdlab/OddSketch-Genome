CXX := g++
CXXFLAGS := -std=c++17 -O3 -march=native -fopenmp -Wall -Wextra

TARGET := oddsketch
SRCS := oddsketch.cpp
OBJS := $(SRCS:.cpp=.o)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ -fopenmp

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
