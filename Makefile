all:
	g++ -std=c++11 -Wall -Wno-stringop-truncation src/Main.cpp -o Calculate.exe -I/usr/local/include -L/usr/local/lib -O3 -lgsl -lgslcblas