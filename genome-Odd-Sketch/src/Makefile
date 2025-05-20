# ソースファイルのリスト
SOURCES = oddsketch.cpp xxhash.cpp

# オブジェクトファイルのリスト（.cpp を .o に変換）
OBJECTS = $(SOURCES:.cpp=.o)

# コンパイラとオプション
CXX = g++
CXXFLAGS = -Wall -std=c++11

# 実行ファイル名
TARGET = oddsketch

# コンパイルルール
$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $(TARGET)

# .cpp から .o への変換ルール
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# クリーンアップルール
clean:
	rm -f $(OBJECTS) $(TARGET)
