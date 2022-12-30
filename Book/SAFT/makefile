all: saft

saft: main.o img.o array.o bscan.o
	g++ -o saft main.o img.o array.o bscan.o
	mv saft ../

main.o: main.cpp  saft.h
	g++ -c main.cpp
img.o: img.cpp saft.h
	g++ -c img.cpp
array.o: array.cpp saft.h
	g++ -c array.cpp
bscan.o: bscan.cpp saft.h
	g++ -c bscan.cpp


