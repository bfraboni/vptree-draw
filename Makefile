all: poster compare

poster: src/poster.cpp src/*.h
	clang++ -o poster src/poster.cpp -O3 -std=c++17 -ICavalierContours/include/ -Wall

compare: src/compare.cpp src/*.h
	clang++ -o compare src/compare.cpp -O3 -std=c++17 -ICavalierContours/include/ -Wall

test: test.cpp src/*.h
	clang++ -o test test.cpp -O3 -std=c++17 -ICavalierContours/include/ -Isrc/ -Wall

clean:
	rm -f compare poster test