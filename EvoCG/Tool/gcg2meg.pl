#gcg2meg.pl  <GCG_FILE>  <MEG_FILE>




#!/usr/bin/perl


my $flag = 0;
my $infile = $ARGV[0];
my $outfile = $ARGV[1];
my $seq = '';
my $line = '';
my %seq = ();
open(IN,$infile)||die"Cannot open $infile!";
open(OUT,">$outfile")||die"Cannot open $outfile!";
while(<IN>){
  if(/^\/\//) {
    $flag = 1;
  } elsif ($flag == 1) {
    $line = $_;
    if ($line =~ m/^(\S+)\s\s\s*(.*)/){
      $strain = $1;
      $seq = $2;
      $seq =~ s/\s//g;
      $seq{$strain} = $seq{$strain}.$seq;
    }
  }
}
close(IN);

print OUT "#mega\n";
print OUT "!Title CG;\n";
print OUT "!Format DataType=Protein indel=-;\n\n";

foreach my $k (keys %seq) {
  print OUT "#",$k,"\n";
  $seq{$k} =~ s/\./\-/g;
  print OUT $seq{$k},"\n\n";
}
close(OUT);
