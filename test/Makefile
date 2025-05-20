# コンパイラとフラグ
CXX = g++
CXXFLAGS = -Wall -std=c++11

# 出力ファイル名
TARGET = test

# ソースファイル（test.cpp のみ）
SRC = test.cpp 

# インクルードディレクトリ（ヘッダーファイルが test ディレクトリ内にあるため）
INCLUDES = -I./

# コンパイルとリンク
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SRC) -o $(TARGET)

# クリーンアップ
clean:
	rm -f $(TARGET)
