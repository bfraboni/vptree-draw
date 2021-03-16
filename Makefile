all: poster compare

poster: poster.cpp *.h
	clang++ -o poster poster.cpp -O3 -std=c++17 -ICavalierContours/include/ -Wall

compare: compare.cpp *.h
	clang++ -o compare compare.cpp -O3 -std=c++17 -ICavalierContours/include/ -Wall

clean:
	rm -f compare poster