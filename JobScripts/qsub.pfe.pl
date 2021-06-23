#!/usr/bin/perl -s

my $NoWatch=$n or $nowatch;

use strict;

my $script = shift(@ARGV);
my $name   = shift(@ARGV);
my @machine = @ARGV;


# Default for machine types
@machine = ('Ivy', 'San', 'Has', 'Bro', 'Sky_ele', 'Bro_ele', 'Cas_ait') if not @machine;


if(not $script or $script =~ /\-+h/i){
    print "
Usage: qsub.pfe.pl [-n] SCRIPT [NAME [MACHINE1 [MACHINE2 [MACHINE3]]]] ...

Submit generic job script to multiple machine types.
Use a unique NAME argument to identify the jobs.
Only the first four characters of the NAME are used. Default NAME is the 
last 4 characters of the directory name where the job is submitted from.
If no machine is specified, 7 jobs will be submitted for the 7 machine types:
IvyBridge, SandyBridge, Haswell, Broadwell, Broadwell_Electra, Skylake_Electra,
and Cascadelake_Aitken. 
Otherwise,the job will be submitted for the listed machines.

The machine name(s) must be provided strictly as:
@machine

Unless the -n (or -nowatch) flag is used, the code starts watch.pfe.pl with
the NAME argument to make sure that when any of the jobs start to run, 
the others get deleted with qdel. The output is piped into watch.log. 

Note you can add or delete jobs with matching NAME while watch.pfe.pl is running.

Example:

qsub.pfe.pl job.long Mars
";
    exit;
}

# Default for job ID
($name) = (`pwd` =~ /(....)$/) if not $name;

# Read original script into $text
print "qsub.pfe.pl reading $script\n";
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
    $text =~ s/^(if.*qsub) .*$/$1 $fileout/m;

    open(SCRIPT, ">$fileout") or die "Could not open $fileout\n";
    print SCRIPT $text;
    close SCRIPT;
}

# submit jobs;
foreach $machine (@machine){
    print "qsub $script.$machine\n";
    `/PBS/bin/qsub $script.$machine`;
}

# start watch.pfe.pl in the background
unless($NoWatch){
    print "watch.pfe.pl $name >& watch.log\n";
    exec("./watch.pfe.pl $name > watch.log 2>&1") unless fork();
}

exit 0;
