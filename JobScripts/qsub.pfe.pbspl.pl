#!/usr/bin/perl -s

my $NoWatch=$n or $nowatch;

use strict;

my $script = shift(@ARGV);
my $name   = shift(@ARGV);
my @machine = @ARGV;


# Default for machine types
if(not @machine){
    if($script =~ /nvidia/){
	@machine = ('sky_gpu', 'cas_gpu');
    }else{
	@machine = ('Ivy', 'San', 'Has', 'Bro', 'Sky_ele', 'Bro_ele', 'Cas_ait');
    }
}

my $qsub  = "qsub.pfe.pbspl.pl";

if(not $script or $script =~ /\-+h/i){
    print "
Usage: $qsub [-n] SCRIPT [NAME [MACHINE1 [MACHINE2 [MACHINE3]]]] ...

Submit generic job script to multiple machine types.
Use a unique NAME argument to identify the jobs.
Only the first four characters of the NAME are used. Default NAME is the 
last 4 characters of the directory name where the job is submitted from.

If no machine is specified, then typically 7 jobs will be submitted for
7 machine types: 'Ivy, San, Has, Bro, Bro_ele, Sky_ele, Cas_Ait'. 
If the name of the script contains 'nvidia', then two jobs are submitted
for 'sky_gpu' and 'cas_gpu'. The abbreviations stand for Ivybridge,
Sandybridge, Haswell, Broadwell, Skylake, Cascadelake CPU names and
Aitken and Electra machine names.
Otherwise, the job will be submitted for the listed machines that
can be chosen from the 9 names above.

Unless the -n (or -nowatch) flag is used and more than one job scripts
were submitted, the code starts watching qstat for all jobs matching
the NAME argument to make sure that when any of the jobs start to run,
the others get deleted with qdel. The output is piped into
watch.NAME.log.

Note you can add or delete jobs with matching NAME while $qsub is running.

Example:

$qsub ensemble.job Event03
";
    exit;
}

# Default for job ID
($name) = (`pwd` =~ /(....)$/) if not $name;

# Read original script into $text
print "$qsub reading $script\n";
my $text;
open(SCRIPT, $script) or die "Could not open $script\n";
$text = join("", <SCRIPT>);
close SCRIPT;

# Copy original script
my $machine;
$script =~ s/.*\///; # remove path (save into local directory) 
$script .= ".$name"; # add name

foreach $machine (@machine){    
    my $fileout = "$script.$machine";
    print "creating $fileout\n";

    # Change name of the job to show machine name
    $text =~ s/^(#PBS -N).*/$1 $name$machine/m;

    # Comment out all machine selections
    $text =~ s/^#+ *(PBS -l.*model=.*)$/### $1/igm;

    # Uncomment the line for model=$machine
    $text =~ s/^### (PBS -l.*model=$machine)$/#$1/im;
    
    # Change the name of the resubmit script
    $text =~ s/^qsub .*$/qsub $fileout/m;

    # Change the event name
    $text =~ s/^(setenv EVENT) Event../$1 $name/m if $name =~ /^Event\d\d$/;

    open(SCRIPT, ">$fileout") or die "Could not open $fileout\n";
    print SCRIPT $text;
    close SCRIPT;
}

# submit jobs;
foreach $machine (@machine){
    print "qsub $script.$machine\n";
    `/PBS/bin/qsub $script.$machine`;
}

exit 0 if $NoWatch or @machine < 2;

my $watch = "watch.$name.log";
print "Start watching jobs. See $watch\n";

# Continue running in the background
if(fork() > 0){exit 0;} 

use POSIX "setsid"; 
POSIX::setsid();
$SIG{HUP}="IGNORE";

open STDOUT, ">$watch" or die "Could not open $watch\n";
open STDERR, ">&STDOUT";

my $pattern = $name;
my $qstat   = '/PBS/bin/qstat -u $USER';
my $qdel    = '/PBS/bin/qdel';

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

close(STDOUT);
close(STDERR);

exit 0;
