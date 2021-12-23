#!/usr/bin/perl
use strict;

my $pattern = $ARGV[0];
my $qstat   = '/PBS/bin/qstat -u $USER';
my $qdel    = '/PBS/bin/qdel';

if(not $pattern){
    print "
Usage: watch.pfe.pl PATTERN >& watch.log &

Watch jobs matching pattern. If any of them starts to run, delete the others.
The following jobs are available now:

",`$qstat`,"

Use PATTERN to select a subset of these, e.g.

watch.pfe.pl SWMF >& watch.log &
";
    exit;
}

my @results;
my $running;
my $user = $ENV{USER};
LOOP:{
    @results = `$qstat | grep $pattern`;
    my $ids;
    foreach (@results){
	/^(\d+).*([A-Z]) +\S+/;
	print "id=$1 status=$2: $_";
	$ids .= " $1";
	$running = $1 if $2 eq "R";
    }
    print "-------------------------\n";
    if($running){
	$ids =~ s/ $running//; # remove running ID
      QDEL:{
	  # Delete all other jobs
	  print "--- qdel $ids ---\n";
	  `$qdel $ids`;
	  # Check if the delete succeeded
	  sleep 2;
	  @results = `$qstat | grep $pattern`;
	  $ids = "";
	  foreach (@results){
	      /^(\d+).*([A-Z]) +(\S+)/;
	      print "id=$1: $_";
	      $ids .= " $1" unless $1 eq $running;
	  }
	  # Delete jobs again if there are any jobs left
	  redo QDEL if $ids;
	}
	last LOOP;
    }
    sleep 5;
    redo;
}

print `$qstat`;

print "Finished watch, job $running is running\n";
