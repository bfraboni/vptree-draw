all: poster vptree compare

vptree: main.cpp
	clang++ -o vptree main.cpp -O3 -std=c++17 -ICavalierContours/include/ -Wall

poster: poster.cpp
	clang++ -o poster poster.cpp -O3 -std=c++17 -ICavalierContours/include/ -Wall

compare: compare.cpp
	clang++ -o compare compare.cpp -O3 -std=c++17 -ICavalierContours/include/ -Wall

clean:
	rm -f vptree compare poster