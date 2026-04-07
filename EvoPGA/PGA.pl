system("rm -R *");
system("rm *");
system("rm -R  ../BactCG1.0/result/");
system("rm -R  ../BactCG1.0/seq/");
system("cp  ../test/seq  ../BactCG1.0/");
system("cp  ../test/26_PG.txt  ./");
system("cp  ../test/nA_AG.gbk  ./");
system("../BactCG1.0/CG  ../BactCG1.0/seq  nA_AG  0.7  0.7");
system("cp  ./26_PG.txt  ../BactCG1.0/result/out_mutbest_filt/");
system("cp  ../bin/PGA  ../BactCG1.0/result/out_mutbest_filt/");
system("../bin/gbParse  ./nA_AG.gbk  >nA_AG_PGAG.tab.txt");
system("cp  ./nA_AG_PGAG.tab.txt  ../BactCG1.0/result/out_mutbest_filt/");
system("../BactCG1.0/result/out_mutbest_filt/PGA  ../BactCG1.0/result/out_mutbest_filt/26_PG.txt  ../BactCG1.0/result/out_mutbest_filt/nA_AG_PGAG.tab.txt  PGAG  nA_AG   >../BactCG1.0/nA.AG_PGAG_PGA.tab.txt");
system("cp ../BactCG1.0/nA.AG_PGAG_PGA.tab.txt  ./");

