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
my $watch = "watch.pfe.pbspl.pl";

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

Unless the -n (or -nowatch) flag is used, the code starts $watch with
the NAME argument to make sure that when any of the jobs start to run, 
the others get deleted with qdel. The output is piped into watch.log. 

Note you can add or delete jobs with matching NAME while $watch is running.

Example:

$qsub job.long Mars
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
my @script;

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

    open(SCRIPT, ">$fileout") or die "Could not open $fileout\n";
    print SCRIPT $text;
    close SCRIPT;
}

# submit jobs;
foreach $machine (@machine){
    print "qsub $script.$machine\n";
    `/PBS/bin/qsub $script.$machine`;
}

# start $watch in the background
unless($NoWatch){
    print "$watch $name >& watch.log\n";
    exec("./$watch $name > watch.log 2>&1") unless fork();
}

exit 0;
