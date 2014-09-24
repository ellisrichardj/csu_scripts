#!/usr/bin/perl -w

# The program fasta_stats_N50.pl calculates some basic stats (GC content, total length, number of gaps etc.) from a file with multiple nucleic acid sequences in fasta format (e.g. genomic contigs or scaffolds).
# Copyright (C) 2012 Minou Nowrousian

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details (http://www.gnu.org/licenses/).


use strict;
use warnings;

unless (@ARGV == 1) {
	die("program usage: fasta_stats_N50.pl inputfile.fasta\n");
}

unless (-e $ARGV[0]) {
	die("Can't open file $ARGV[0]: $!");
}

open INPUTFILE, "< $ARGV[0]" or die "Kann Datei nicht oeffnen: $!";

my $header_found = 0;  ## set variable to indicate that fasta header line (begins with >) was found
my $header;            ## contains the header line of the current fasta sequence
my $seq;
my $z = 0;             ## line counter

my $gc = 0;            ## total gc content of all reads
my $total_length = 0;  ## length of all reads
my $number_all_gaps = 0;   ## number of gaps of one or more undefined bases ("N" or "n") in all sequences
my $length_all_gaps = 0;   ## number of undefined bases ("N" or "n") in all sequences (add up length of gaps)
my @lengths;           ## array with the lengths of all fasta sequences

my @stats = "fasta header\tlength of sequence\tgc content in percent\tnumber of gaps\tlength of gaps\tpercent of gaps\n";

while (<INPUTFILE>) {      # parses through lines of inputfile

	my $line = $_;
	$z++;
	$line =~ s/\s+$//;     # remove all trailing whitespace, end of line, tabs etc.

	if ($line =~ /^>/) {
		if (defined($seq)) {
			my ($length, $gc_seq, $gaps, $length_gaps, $percent_gaps) = &stats($seq);
			push(@stats, "$header\t$length\t$gc_seq\t$gaps\t$length_gaps\t$percent_gaps\n");
		}
		$header = $line;
		$header =~ s/^>//;
		$header_found = 1;
		$seq = undef;
	} elsif ($line =~ /[^acgtnACGTN]/) {
		print "unexpected character (not a, c, g, t, n, A, C, G, T, N) in line $z, analysis will proceed anyway\n";
		$seq .= $line;
	} else {
		$seq .= $line;
	}

}

## process data for last sequence

if (defined($seq)) {
	my ($length, $gc_seq, $gaps, $length_gaps, $percent_gaps) = &stats($seq);
	push(@stats, "$header\t$length\t$gc_seq\t$gaps\t$length_gaps\t$percent_gaps\n");
}


## calculate N50

my @sorted_lengths = sort { $b <=> $a } @lengths; ## sort length of fasta sequences in reverse order (biggest first)

my $n = @sorted_lengths;
my $N50 = 0;
my $addup = 0;
my $halfsize = $total_length / 2;
until ( ($addup >= $halfsize) || ($n == 0) ) {
	$N50 = shift(@sorted_lengths);
	$addup += $N50;
	$n = @sorted_lengths;
}

## now put stats for all sequences in array for output

my $number_seqs = @lengths;
my @output = ("number of sequences: ", $number_seqs, "\n\n");
my $gc_content = sprintf "%.1f", ($gc/$total_length)*100;
push(@output, "GC content in percent\t" . $gc_content . "\n\n");
push(@output, "Total length of all sequences in bases\t" . $total_length . "\n\n");
push(@output, "N50 in bases\t" . $N50 . "\n\n");
push(@output, "number of gaps\t" . $number_all_gaps . "\n\n");
push(@output, "Total length of all gaps in bases\t" . $length_all_gaps . "\n\n");
my $percent_gaps = sprintf "%.1f", $length_all_gaps * 100 / $total_length;
push(@output, "Total gap length in percent\t" . $percent_gaps . "\n\n");

## now add stats for individual sequences to output and save output file

push(@output, @stats);

my $filename = $ARGV[0] . "_stats.txt";
open NEU, "> $filename" or die "Kann Datei nicht oeffnen: $!";
print NEU @output;
close NEU;


sub stats {

	my $seq = shift(@_);
	my $length = length($seq);
	$total_length += $length;
	push(@lengths, $length);
	
	## determine GC content
	
	my $gc_seq = 0;
	while ($seq =~ /[CG]/ig) {
		$gc_seq += 1;
	}	
	$gc += $gc_seq;
	$gc_seq = sprintf "%.1f", $gc_seq * 100 / $length;  ## GC content in percent

	## determine gaps
	
	my $gaps = 0;
	my $length_gaps = 0;
	my $percent_gaps = 0;
	
	if ($seq =~ /NN*/i) {
		my @partialseqs = split /NN*/i, $seq;
		$gaps = @partialseqs - 1;
		my $restseq = join "", @partialseqs;
		my $restlength = length($restseq);
		$length_gaps = $length - $restlength;
		$percent_gaps = sprintf "%.1f", ($length_gaps/$length)*100;
	} 
	
	$number_all_gaps += $gaps;
	$length_all_gaps += $length_gaps;
	
	return ($length, $gc_seq, $gaps, $length_gaps, $percent_gaps);

}



