#!/bin/bash

rm MATMULT
rm matMult.o

g++ -c matMult.cpp -I./papi-4.1.2.1/src
g++ -o MATMULT matMult.o -L./papi-4.1.2.1/src/ ./papi-4.1.2.1/src/libpapi.a

N=1500
TYPE=4

echo -e "TYPE\tN\tPAPI_FP_OPS"
#for ((j=1; j<=TYPE; j=j+1)); do
./MATMULT $N 4
#./MATMULT $N 3
#done
