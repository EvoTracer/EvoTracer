#!/usr/bin/perl

my $dir_name = '';
my @dots = [];
my $fileName = '';
my $outfile = '';
my $filePrefix = '';
my $megafile = '';


$dir_name="output/CG_results/2.result/";


opendir(DIR, $dir_name."/") || die "Can't open directory $dir_name";
@dots = readdir(DIR);


foreach my $fileName(@dots){
 $fileName0 = $dir_name."/".$fileName;
 if($fileName =~ m/(\S+)\.fa/){
     $filePrefix = $1; 
     $outfile = "output/CG_results/3.gcg/".$filePrefix.".gcg";
     $megafile = "output/CG_results/3.mega/".$filePrefix.".meg";
     system("Tool/clustalw2 -infile=$fileName0 -type=protein -output=gcg -outfile=$outfile -pwmatrix=GONNET -pwgapopen=10 -pwgapext=0.1 -gapopen=10 -gapext=0.2 -gapdist=4 -align");
     system("perl Tool/gcg2meg.pl $outfile $megafile");
  }
}
closedir DIR;