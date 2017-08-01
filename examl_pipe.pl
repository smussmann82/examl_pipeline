#! /usr/bin/perl

use warnings;
use strict;
use Getopt::Std;
use Data::Dumper;
use Cwd;
use Cwd 'abs_path';
use List::Util 'shuffle';
use File::Path;
use File::Copy;
use File::Spec;
#use File::Copy qw( move );
#use File::Copy::Recursive qw( dirmove );

# kill program and print help if no command line arguments were given
if( scalar( @ARGV ) == 0 ){
  &help;
  die "Exiting program because no command line options were used.\n\n";
}

# take command line arguments
my %opts;
getopts( 'b:B:hp:t:', \%opts );

# if -h flag is used, or if no command line arguments were specified, kill program and print help
if( $opts{h} ){
  &help;
  die "Exiting program because help flag was used.\n\n";
}

# declare variables
my $dir = getcwd(); # get current working directory
chomp( $dir );

# parse the command line
my( $phy, $bin, $ntrees, $nboot ) = &parsecom( \%opts );
my $pfile = "raxml_commands"; #file to hold commands that will be sent to raxml
my $bfile = "boot_commands";
my $name = "StartingTree"; #basename for starting tree files
my $streedir = "starting_trees";
my $bootdir = "boot_phylip";
my $bootstartdir = "boot_starting_trees";
my $bootbindir = "boot_binary";
my $pars_con_trees = "starting.trees";

# create input / output directories
mkpath( [$streedir], 1, 0755 ) or die "Couldn't make directory to hold raxml StartingTree files in $dir: $!";
mkpath( [$bootdir], 1, 0755 ) or die "Couldn't make directory to hold raxml StartingTree files in $dir: $!";
mkpath( [$bootstartdir], 1, 0755 ) or die "Couldn't make directory to hold raxml StartingTree files in $dir: $!";
mkpath( [$bootbindir], 1, 0755 ) or die "Couldn't make directory to hold raxml StartingTree files in $dir: $!";

=pod
mkpath( [$phylip], 1, 0755 ) or die "Couldn't make directory to hold phylip files in $dir: $!";
mkpath( [$grouped_loci], 1, 0755 ) or die "Couldn't make directory to hold grouped locus files in $dir: $!";
mkpath( [$astral], 1, 0755 ) or die "Couldn't make directory to hold astral input and output in $dir: $!";
=cut

&binaryConvert( $phy, $bin ); #make initial binary conversion
&startParallel( $ntrees, $phy, $name, $pfile, "$dir/$streedir" ); #generate starting tree commands
&gnuParallel( $pfile ); #generate the actual starting trees

system("cat $dir/$streedir/RAxML_parsimonyTree.$name.* > $pars_con_trees"); # concatenate starting trees into a single file

&runExaml($pars_con_trees, $bin, "T1"); # run EXaML


&makeBootTrees( $nboot, $phy, "$dir/$bootdir" ); #generate RAxML commands for  pseudoreplicate phylip files for bootstrapping

#generate starting tree commands for all bootstrap replicates
&makeBootTreeCommands( "$dir/$bootdir", $ntrees, $bfile, "$dir/$bootbindir", "$dir/$bootstartdir" );
&gnuParallel( $bfile );

opendir( WD, $dir ) or die "Can't open $dir, d-bag: $!\n\n";

my @contents = readdir( WD );
foreach my $file( @contents ){
	if( $file =~ /^boot_\d+/ ){
		my $dest = "$dir/$bootbindir/$file";
		my $current = "$dir/$file";
		move( $current, $dest );
	}
}

closedir WD;


for( my $i=0; $i<$nboot; $i++ ){
	my $output = "bootfile.$i";
	open( OUT, '>', $output ) or die "Can't open $output, d-bag: $!\n\n";

	for( my $j=0; $j<$ntrees; $j++ ){
		my $input = "$dir/$bootstartdir/RAxML_parsimonyTree.bootstart.$i.$j";
		open( IN, "$input" ) or die "Can't open $input, d-bag: $!\n\n";
		while( my $line = <IN> ){
			chomp $line;
			print OUT $line, "\n";
		}
		close IN;
	}

	close OUT;
}

for( my $i=0; $i<$nboot; $i++ ){
	&runExaml( "bootfile.$i", "$dir/$bootbindir/boot_$i", "BOOT$i"  );
}

system("cat ExaML_result.BOOT* > all_boot_trees.txt");

system("raxmlSSE3-serial -f b -t ExaML_result.T1 -z all_boot_trees.txt -m GTRCAT -s $phy -n final");

exit;

#####################################################################################################
############################################ Subroutines ############################################
#####################################################################################################

# subroutine to print help
sub help{
  
  print "\nexaml_pipe.pl is a perl script developed by Steven Michael Mussmann\n";
  print "This script executes the examl pipeline using the phylip file from pyRAD.\n\n";
  print "To report bugs send an email to mussmann\@email.uark.edu\n";
  print "When submitting bugs please include all input files, options used for the program, and all error messages that were printed to the screen\n\n";
  print "Program Options:\n";
  print "\t\t[ -b | -B | -h | -p | -t ]\n\n";
  
  print "\t-b:\tSpecify the number of bootstrap replicates to be conducted\n";
  print "\t\tOptional.  Default=100.  Set to 0 to produce a tree without bootstrapping\n\n";

  print "\t-B:\tUse this flag to change the basename of the binary input files.\n";
  print "\t\tOptional.  Default=spd.\n\n";

  print "\t-h:\tUse this flag to display this help message.\n";
  print "\t\tProgram will exit after displaying the help message.\n\n";

  print "\t-p:\tUse this flag to display the input phylip file name.\n";
  print "\t\tProgram will exit if name is not provided.\n\n";
  
  print "\t-t:\tUse this flag to specify a number of starting trees.\n";
  print "\t\tOptional.  Default=32.  Must be positive integer\n\n";

}

#####################################################################################################
# subroutine to parse the command line options

sub parsecom{ 
  
	my( $params ) =  @_;
	my %opts = %$params;

	# set default values for command line arguments
	my $phy = $opts{p} or die "\nNo phylip file specified.\n\n";	
	my $bin = $opts{B} || "spd";
	my $ntrees = $opts{t} || "32";
	my $nboot = $opts{b} || "100";

	return( $phy, $bin, $ntrees, $nboot );

}

#####################################################################################################
# subroutine to convert phylip file to binary

sub binaryConvert{
	my( $file, $bin ) = @_;

	system("parse-examl -s $file -m DNA -n $bin");

}

#####################################################################################################
# subroutine to write commands for generating starting trees to a file that will be passed to gnu parallel

sub startParallel{
	my( $ntrees, $file, $name, $pfile, $outdir ) = @_; #$ntrees = number of starting trees, $file = phylip input file, $name = basename for output files, $pfile is text file that will contain raxml commands

	open( OUT, '>', $pfile ) or die "Can't open $pfile, d-bag:$!\n\n";

	for( my $i=0; $i<$ntrees; $i++ ){
		my $rand = int(rand(2147483647));
		print OUT "raxmlSSE3-serial -y -m GTRCAT -p $rand -s $file -n $name.$i -w $outdir\n";
	}

	close OUT;

}

#####################################################################################################
# subroutine to run a file through gnu parallel

sub gnuParallel{

	my( $file ) = @_;
#	my $unique = "unique-nodelist.txt";
#	my $nodefile = $ENV{'PBS_NODEFILE'};
#	system( "sort -u $nodefile > $unique" );
#	system("cat $file | parallel --sshloginfile $unique");
	system("cat $file | parallel");

}

#####################################################################################################
# subroutine to run a file through gnu parallel

sub runExaml{
	my( $parsconcat, $bin, $exaName ) = @_;
	
	my $np = 0;

	my $nodefile = $ENV{'PBS_NODEFILE'};
	open( NODE, $nodefile ) or die "Can't open $nodefile, d-bag: $!\n\n";
	while( my $line = <NODE> ){
		$np++;
	}
	close NODE;
	
	print $np, "\n";

	system("mpirun -np $np -machinefile $nodefile examlAVX -t $parsconcat -m PSR -s $bin.binary -n $exaName");

}
#####################################################################################################
#

sub makeBootTrees{
	my( $nboot, $file, $outdir ) = @_;
	
	my $rand = int(rand(2147483647));

	system( "raxmlSSE3-serial -# $nboot -b $rand -f j -m GTRCAT -s $file -n REPS -w $outdir" );

}

#####################################################################################################
# subroutine to generate bootstrap tree commands

sub makeBootTreeCommands{
	my( $dir, $ntrees, $commandfile, $outdir_bin, $outdir_trees  ) = @_;


	# get list of bootstrap files
	opendir( WD, $dir ) or die "Can't open $dir, d-bag: $!\n\n";

	open( OUT, '>', $commandfile ) or die "Can't open $commandfile, d-bag: $!\n\n";

	my @contents = readdir( WD );
	foreach my $file( @contents ){
		if( $file =~ /BS(\d+)$/ ){
			print $file, "\n";
			print OUT "parse-examl -s $dir/$file -m DNA -n boot_$1\n";
		}
	}

	foreach my $file( @contents ){
		if( $file =~ /BS(\d+)$/ ){
			for( my $i=0; $i<$ntrees; $i++ ){
				my $rand = int(rand(2147483647));
				print OUT "raxmlSSE3-serial -y -m GTRCAT -p $rand -s $dir/$file -n bootstart.$1.$i -w $outdir_trees\n";
			}
		}
	}

	close OUT;

	closedir WD;

}

#####################################################################################################
