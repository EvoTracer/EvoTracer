#!/usr/bin/perl
use strict;
use warnings;

my $inputFile = '';
my $seq = '';
my $genomeSeq = '';
my %genome = ();
my $ag = '';
my $sp = '';
my $st = 0;
my $en = 0;
my $seqLen = 0;

# Read the genome files
for (my $i = 0; $i <= 1; $i++) {
    $inputFile = $ARGV[$i] . ".fasta";
    open(IN, "<", $inputFile) or die "Cannot open $inputFile!\n";
    while (<IN>) {
        if (/^\>/) {
            next;
        } elsif (/(\S+)/) {
            $seq = $1;
            if (length($seq) > 0) {
                $genomeSeq .= $seq;
            }
        }
    }
    close(IN);
    $genome{$ARGV[$i]} = $genomeSeq;
    $genomeSeq = '';
    $seq = '';
}

# Read the positions file
open(IN1, "<", $ARGV[2]) or die "Cannot open $ARGV[2]!\n";
while (<IN1>) {
    if (/^(\S+)\t(\d+)\t(\d+)/) {
        $sp = $1;
        $st = $2;
        $en = $3;
        $seqLen = $en - $st + 1;
        $st = $st - 1;
        if (exists $genome{$sp}) {
            $ag .= substr($genome{$sp}, $st, $seqLen);
        } else {
            warn "Species $sp not found in genome data\n";
        }
    }
}
close(IN1);

# Print the result
print ">Salmonella orthologous_ancient genome\n";
for (my $j = 0; $j < length($ag); $j += 80) {
    print substr($ag, $j, 80), "\n";
}