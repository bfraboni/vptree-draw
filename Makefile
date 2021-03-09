all: vptree

vptree: main.cpp
	g++ -o vptree main.cpp -O3 -std=c++17

clean:
	rm -rf vptree 