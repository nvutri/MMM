#!/bin/bash

function partA(){
	#Part A code
	echo "PART A"
	echo -e "TYPE\t N\tPAPI_FP_OPS\t\tTime Elapsed"
	for ((j=1; j<=2; j=j+1)); do
		./MATMULT $N $j A
	done
}
function partB(){
	#Part B code
	echo "PART B"
	echo -e "TYPE\t N\tPAPI_FP_OPS\t\tTime Elapsed"
	for ((j=1; j<=6; j=j+1)); do
		./MATMULT $N $j B
	done
}

function partC(){
	#Part C code
	echo "PART C"
	echo -e "NB\t N\tPAPI_FP_OPS\t\tTime Elapsed"
	for ((NB=16; NB<=80; NB=NB+4)); do
		./MATMULT $N $NB C
	done
}

function partD(){
	#Part D
	echo "PART D"
	echo -e "TYPE\t N\tPAPI_FP_OPS\t\tTime Elapsed"
	./MATMULT $N 64 D
}

function final(){
	echo "Final Fastest Version"
	echo -e "TYPE\t N\tPAPI_FP_OPS\t\tTime Elapsed"
	./MATMULT $N 64 F 
}

rm MATMULT
rm matMult.o

g++ -O2 -c matMult.cpp -I./papi-4.1.2.1/src
g++ -o MATMULT matMult.o -L./papi-4.1.2.1/src/ ./papi-4.1.2.1/src/libpapi.a

N=960
TYPE=6
MAX_NB=80

final