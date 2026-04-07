#!/bin/bash

#echo "Cleaning previous result files..."
#rm -rf results/*
#rm -rf BactCG1.0/result/
#rm -rf BactCG1.0/seq/

echo "Copying input files..."
cp -r ./input/seq ./BactCG1.0/
cp ./input/255_PG.txt ./results/
cp ./input/nA_AG.gbk ./results/

#echo "Running CG alignment program..."
#cd BactCG1.0
#./CG seq nA_AG 0.7 0.7
#cd ..

echo "Copying PGA program and alignment list..."
cp results/255_PG.txt BactCG1.0/result/out_mutbest_filt/
cp bin/PGA BactCG1.0/result/out_mutbest_filt/

echo "Running gbParse to generate tabular data..."

./bin/gbParse results/nA_AG.gbk > results/nA_AG_PGAG.tab.txt

cp results/nA_AG_PGAG.tab.txt BactCG1.0/result/out_mutbest_filt/

echo "Running PGA analysis..."
cd BactCG1.0/result/out_mutbest_filt/

./PGA 255_PG.txt nA_AG_PGAG.tab.txt PGAG nA_AG > nA.AG_PGAG_PGA.tab.txt

cd ../../../

echo "Copying final result to results/..."
cp BactCG1.0/result/out_mutbest_filt/nA.AG_PGAG_PGA.tab.txt results/

echo "Pipeline execution completed. Final result is in: results/nA.AG_PGAG_PGA.tab.txt"