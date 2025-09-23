# Makefile for compiling and running the C++ program

# コンパイラの設定
CXX = g++
CXXFLAGS = -std=c++17 -Wall

# ソースファイル名と出力ファイル名
SRC = cal_Jaccard.cpp
OUT = cal_Jaccard

# ターゲット: 実行ファイルのビルド
all: $(OUT)

# 実行ファイルをビルドするルール
$(OUT): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(OUT)

# 実行ファイルの削除（クリーン）
clean:
	rm -f $(OUT)

# 実行ファイルの実行
run: $(OUT)
	./$(OUT)

.PHONY: all clean run
