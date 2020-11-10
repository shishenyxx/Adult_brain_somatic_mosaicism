#!/usr/bin/env perl
#######################################################################
## @version 0.34
##
## Change Log
##    Add swtich to print only the match lines
##
## @version 0.33
##
## change log
##    Auto use zcat to read files with extention .gz
##    Auto use bzcat to read files with extention .bz2
##  
## @version 0.32
##
## Change log 
##    Fix bugs in appendstring
##
## @version 0.31
##
## Change log 
##     Add the regexp support for the split separator
##
## @version 0.3
##
## change log
##     Add the connector for multi-keyword connection
##
## 
## @version 0.2
##
## change log
##     change parameter parser to Getopt::Long modular
##     Add the debug mode
##     Add one-to-more hit
##     Get ride of the bug for no hit
##                           
##                                  Sat May 14 22:45:29 CST 2005
#######################################################################

use warnings;
use strict;
use Getopt::Long;
use Getopt::Long qw(:config no_ignore_case bundling);

use vars qw($field_first $field_second $seperator_first $seperator_second $seperator_output $gap_num);
use vars qw( $file_first $file_second $seperator_connect $debug $top_num $append_string $onlymatch $help);

GetOptions(
	   "F=s"   => \$field_first,
	   "f=s"   => \$field_second,
       "S=s"   => \$seperator_first,
	   "s=s"   => \$seperator_second,
       "o=s"   => \$seperator_output,
	   "c=s"   => \$seperator_connect,
       "N=i"   => \$gap_num,
       "d"     => \$debug,
       "a=s"   => \$file_first,
	   "b=s"   => \$file_second,
	   "x=s"   => \$append_string,
	   "t=n"   => \$top_num,
	   "h"     => \$help,
	   "m"     => \$onlymatch,
           )   || usage();

usage() if defined($help);
usage() unless (defined($file_first) && defined($file_second));

my %index=();

$field_first="1" unless defined($field_first);
$field_second="1" unless defined($field_second);

$seperator_first="\t" unless defined($seperator_first);
$seperator_first=~s/^\"|\"$//g;

$seperator_second="\t" unless defined($seperator_second);
$seperator_second=~s/^\"|\"$//g;

$seperator_output="\t" unless defined($seperator_output);
$seperator_output=~s/^\"|\"$//g;

$seperator_connect="\t" unless defined($seperator_connect);
$seperator_connect=~s/^\"|\"$//g;

$gap_num=1 unless $gap_num;

#$append_string = "-" unless defined($append_string);
$append_string = "" unless defined($append_string);

if ( $file_first =~ /\.gz$/ ) {
	open(FILE1,"zcat '$file_first'|") or die "cannot open the file $file_first";
} elsif ( $file_first =~ /\.bz2$/ ) {
	open(FILE1,"bzcat '$file_first'|") or die "cannot open the file $file_first";
} else {
	open(FILE1,$file_first) or die "cannot open the file $file_first";
}

if ( $file_second =~ /\.gz$/ ) {
	open(FILE2,"zcat '$file_second'|") or die "cannot open the file $file_second";
} elsif ( $file_second =~ /\.bz2$/ ) {
	open(FILE2,"bzcat '$file_second'|") or die "cannot open the file $file_second";
} else {
	open(FILE2,$file_second) or die "cannot open the file $file_second";
}

my @keya=split(/\D+/,$field_first);
my @keyb=split(/\D+/,$field_second);

while(my $line=<FILE2>)
{
    chomp($line);
    my @items=split(/$seperator_second/,$line);

    my $key_file2="";

    for(my $i=0;$i<@keyb;$i++)
    {
	$key_file2.=$items[$keyb[$i]-1];
	$key_file2.="$seperator_connect" if $i<@keyb-1;
    }
	
    warn ">>> key2: [$key_file2] sep2: [$seperator_second]\n" if $debug;

    push(@{$index{$key_file2}},$line);
}

while(my $line=<FILE1>)
{
    chomp($line);
    my @items=split(/$seperator_first/,$line);

    my $key_file1="";

    for(my $i=0;$i<@keya;$i++)
    {
	$items[$keya[$i]-1]="" unless defined($items[$keya[$i]-1]);
	$key_file1.=$items[$keya[$i]-1];
	$key_file1.="$seperator_connect" if $i<@keya-1;
    }

    warn "<<< key1: [$key_file1] sep1: [$seperator_first]\n" if $debug;

    if ($key_file1 ne "" && exists($index{$key_file1})) {
	my $print_line = scalar(@{$index{$key_file1}});
	$print_line = $top_num if (defined($top_num) && $top_num < $print_line);
        for(my $i=0;$i<$print_line;$i++) {
	    print "$line$seperator_output";
	    print $index{$key_file1}->[$i]."\n";
        }
    }
    elsif ( !defined($onlymatch) ) {
    	print "$line$seperator_output";
	print "$append_string$seperator_output" x ($gap_num - 1) if $gap_num>1;
        print "$append_string\n";
    }

}
close(FILE1);
close(FILE2);

exit(0);
################################################
sub usage
{
    print "\nThis program join two file together used the common key\n\n";
    print "Usage:\n";
    print "$0 [-d] [-F 1,2] [-f 1,2] [-S \"\\t\"] [-s \"\\t\"] [-c \"\\t\"] [-N 1] [-m] [-x \" \"] [-t 1] [-o \"\\t\"] <-a file1> <-b file2>\n";
    print "  d: debug mode\n";
    print "  F: the field group of key in file1, separate by ',' [default : 1]\n";
    print "  f: the field group of key in file2, separate by ',' [default : 1]\n";
    print "  c: the connector of keylist, used when connect with many keywords [default :\\t]\n";
    print "  S: the seperator of file1, regexp support warning \\, quote with ' [default : \\t]\n";
    print "  s: the seperator of file2, regexp support warning \\, quote with ' [default : \\t]\n";
    print "  o: the connector for two file [default : \\t]\n";
    print "  N: the number of \\t if lost data in file2 [default : 1]\n";
    print "  m: output the match lines, support the unique lines\n";
    print "  x: append string for lost data [default : SPACE]\n";
    print "  t: the top number left if there is multi-item exists in file2 for one item in file1\n";
    print "  file1, file2:  use zcat to open the file with extention name .gz and\n";
    print "                 use bzcat to open the file with extention name .bz2\n\n";
    exit(0);
}
