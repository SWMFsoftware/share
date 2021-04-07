#!/usr/bin/perl -s
# Usage: uptodate.pl [-v] FILE DIR
# returns empty string if file is newer than any file in directory DIR
# returns "UPDATE" otherwise
# -v is for verbose output.
my $verbose = $v;

use strict;

my $file = $ARGV[0] or die "ERROR in uptodate.pl: Missing file name\n";
my $dir  = ($ARGV[1] or ".");
print "file=$file dir=$dir\n" if $verbose;

if(-f $file){
    my $filetime = -M $file;
    my $latest = `ls -tR $dir | head -1`; chop $latest;
    my $latesttime = -M $latest;
    print "$file modtime $filetime, $latest modtime=$latesttime: ",
	$filetime < $latesttime,"\n" if $verbose;
    if( $filetime <= $latesttime){exit};
}
print "UPDATE\n";

