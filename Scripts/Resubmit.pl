#!/usr/bin/perl -s

my $jobfile = $ARGV[0];
if(not $jobfile or $h or $help){
    print "
Usage: resubmit.frontera.pl jobfile >& resubmit.log & 
";
    exit 0;
}
use strict;

my $Success1 = "SWMF.SUCCESS";
my $Success2 = "BATSRUS.SUCCESS";
my $Done1    = "SWMF.DONE";
my $Done2    = "BATSRUS.DONE";

unlink $Done1, $Done2;                    # remove stop indicators
while(not -f $Done1 and not -f $Done2){   # loop until done
    unlink $Success1, $Success2;          # remove success indicators
    `sbatch $jobfile`;                    # submit job
    print "$jobfile submitted on ", `date`;
    sleep 10 while not (-f $Success1 or -f $Success2); 
    print "job finished successfully on ", `date`;
    sleep 60;                             # wait for the job to quit
    'Restart.pl';                         # process restart files
    if(-f "PARAM.in.restart" and not -f "PARAM.in.start"){
	`mv PARAM.in PARAM.in.start`;
	`cp PARAM.in.restart PARAM.in`;
    }
}
print "job reached final time on ", `date`;
exit 0;
