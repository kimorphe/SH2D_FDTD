all: wveq2d

wveq2d: domain.o field.o main.o cntrl.o wave1d.o source.o receiver.o wveq2d.h wave1d.h
	g++ -o wveq2d main.o domain.o field.o cntrl.o wave1d.o source.o  receiver.o

main.o: main.cpp wveq2d.h
	g++ -c main.cpp
field.o: field.cpp wveq2d.h
	g++ -c field.cpp
domain.o: domain.cpp wveq2d.h
	g++ -c domain.cpp
cntrl.o: cntrl.cpp wveq2d.h
	g++ -c cntrl.cpp

wave1d.o: wave1d.h wave1d.cpp
	g++ -c wave1d.cpp

source.o: source.cpp
	g++ -c source.cpp
receiver.o: receiver.cpp
	g++ -c receiver.cpp
