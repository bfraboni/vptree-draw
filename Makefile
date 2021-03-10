all: vptree

vptree: main.cpp
	g++ -o vptree main.cpp -O3 -std=c++17 -ICavalierContours/include/

clean:
	rm -rf vptree 