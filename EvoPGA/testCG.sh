#!/bin/bash

echo "Cleaning previous BactCG1.0 outputs..."
rm -rf BactCG1.0/result/
rm -rf BactCG1.0/seq/

echo "Copying input sequence folder to BactCG1.0/seq..."
cp -r ./input/seq ./BactCG1.0/

echo "Running CG alignment inside BactCG1.0..."
cd BactCG1.0
./CG seq nA_AG 0.7 0.7
cd ..

echo "Checking CG output..."
ls -lh BactCG1.0/result/out_mutbest_filt/