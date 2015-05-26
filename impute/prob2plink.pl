#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use List::Util qw[min max];


###############################################################
##################### DOCUMENTATION ###########################

## Author: Yun Li
## Date: 2010-11-29
## Goal: convert MaCH output prob+info file into PLINK input dosage dat file
my $location = $0;
my $title = "prob2plink version 0.0.1";

print "$title\n";
print "\t\@: $location\n";
print "\tby: Yun Li(ylwtx\@umich.edu)\n\n";

## Input Parameters
##

## Input Files
## (1) prob
## (2) info (only need 1. SNP name; 2. AL1; and 3. AL2)

## Output Files
## .plink_dat
## .fam

################# END OF DOCUMENTATION ########################
###############################################################
my $starttime = time();

my $MIN_NUM_OPTS = 2;
my $MIN_INFO_FIELDS = 3;

if($#ARGV<($MIN_NUM_OPTS*2 - 1))
{
	&usage();
	die "Failed to parse parameters\n\n";
}

my %opts = ();

	# Default Optional Options
$opts{o} = "recoded.plink";

Getopt::Long::GetOptions(\%opts,qw(
	prob=s
	info=s
	o=s
)) || die "Failed to parse options\n\n";

&printPar();

my @SNPs = ();
my @alleleOnes = ();
my @alleleTwos = ();
my $numSNPs = 0;

my @probs = ();
my @famIDs = ();
my @subIDs = ();
my $numPeople = 0;
&main();

my $endtime = time();
my $lapsetime = $endtime - $starttime;
print "\nAnalysis took $lapsetime seconds\n\n";

# begin sub
sub usage
{
	print "\n";
	print "Usage:\n\t";
	print "-prob\t MaCH output prob file (Required) \n\t";
	print "-info\t MaCH output info file (Required) \n\t";
	print "-o\t Output prefix (Optional, default=\"recoded.plink\") \n\t";
	print "\n";
}

sub printPar
{
	print "\n";
	print "Parameters in Effect:\n";
	print "\t MaCH output prob file \n\t\t-prob '$opts{prob}'\n";
	print "\t MaCH output info file \n\t\t-info '$opts{info}'\n";
	print "\t Output prefix \n\t\t-o '$opts{o}'\n";
	print "\n";
}

sub rm
{
	if(-e $_[0])
	{
		print "WARNING: $_[0] existed and removed!!\n\n";
		system "/bin/rm -rf $_[0]";
	}
}

#==================================================================
# Subroutine:
#   Open($filename, $handle)
#
# Usage:
#   Open('myfile.txt', \*IN);
#   Open($pedfile, \*PED);
#
#   Then refernce the handles <IN> and <PED> as usual.
#
# Description:
#   Open the file 'filename'. The file can be a normal file
#   or one compressed with gzip or bzip2.
#   Function does not return if the open fails.
#
# Arguments:
#   $filename - path to file to open
#   $handle - GLOB name of handle to use. Include '\*'
#
#==================================================================
sub Open {
    my ($filename, $handle) = @_;

    # Handle gzipped data
    if ($filename =~ /\.gz$/) {
        my $cmd = "gunzip -c $filename";
        open($handle, "$cmd |") ||
            die "Unable to invoke command '$cmd': $!\n";
        if (eof($handle)) { die "Cmd failed: '$cmd'\n"; }
        return 1;
    }
    # Handle bzip2 data
    if ($filename =~ /\.bz2$/) {
        my $cmd = "bunzip2 -c $filename";
        open($handle, "$cmd |") ||
            die "Unable to invoke command '$cmd': $!\n";
        if (eof($handle)) { die "Cmd failed: '$cmd'\n"; }
        return 2;
    }
    # Handle normal file
    open($handle, $filename) ||
        die "Unable to open file '$filename': $!\n";
    return 3;
}

sub main
{
  &loadInfo();
  &loadProb(); #may not scale well with huge dataset, since loading all into RAM first;
  &output();
}

sub loadInfo
{
  open(IN,$opts{info}) || die "cannot open file $opts{info}\n\n";
  # sanity check
  my $headerline = <IN>;
  chomp($headerline);
  my @headerarr = (split(/\s+/, $headerline));
  my $numFields = scalar(@headerarr);
  if ($numFields < $MIN_INFO_FIELDS)
  {
    print "ERROR: info file must have at least three fields: SNP, Al1 and Al2\n\n";
    exit;
  }
  
  if ($headerarr[0] ne "SNP" || $headerarr[1] ne "Al1" || $headerarr[2] ne "Al2")
  {
    print "ERROR: the first three fields of the info file must be: SNP, Al1 and Al2 in that order\n\n";
    exit;
  }
  
  # read SNPs
  while (my $line = <IN>)
  {
    chomp($line);
    my @linearr = (split(/\s+/, $line));
    if (scalar(@linearr) == 0) {next;}
    if (scalar(@linearr) != $numFields)
    {
      print "ERROR: the following line in info file contains different number of fields than specified by ther header line\n";
      print "       $line\n\n";
      exit;
    }
    
    push @SNPs, $linearr[0];
    push @alleleOnes, $linearr[1];
    push @alleleTwos, $linearr[2];
  }
  
  $numSNPs = scalar(@SNPs);
  print "Read $numSNPs SNPs from info file\n";
  
  close(IN);
}

sub loadProb
{
  Open($opts{prob},\*IN) || die "cannot open file $opts{prob}\n\n";
  my $expFieldsPerLine = $numSNPs * 2 + 2;
  my $linenum = 0;
  while(my $line = <IN>)
  {
    $linenum ++;
    chomp($line);
    my @linearr = (split(/\s+/,$line));
    if (scalar(@linearr) == 0) {next;}
    
    if (scalar(@linearr) != $expFieldsPerLine)
    {
      print "ERROR: line $linenum in prob file contains different number of fields than expected from info file\n\n";
      exit;
    }
    
    my @IDs = (split(/->/,$linearr[0]));
    
    ## could print out the IDs directly to output .fam file without storing in RAM ##
    push @famIDs, $IDs[0];
    push @subIDs, $IDs[1];
    
    ## store probs in RAM ##
    for (my $j = 2; $j <= $#linearr; $j++)
    {
      push @probs, $linearr[$j];
    }
    
  }
  close(IN);
  $numPeople = scalar(@famIDs);
  print "Read $numPeople people from prob file\n";
}

sub output
{
  ######################
  ## output .fam file ##
  ######################
  my $ofam = $opts{o}.".fam";
  &rm($ofam);
  print "Generating .fam output ...\n";
  open(OUT, ">>$ofam") || die "cannot create file $ofam\n\n";
  for (my $i = 0; $i < $numPeople; $i ++)
  {
    print OUT "$famIDs[$i] $subIDs[$i] 0 0 1 -9\n";
  }
  close(OUT);
  
  ############################
  ## output .plink_dat file ##
  ############################
  my $odat = $opts{o}.".plink_dat";
  &rm($odat);
  print "Generating .plink_dat output ...\n";
  open(OUT, ">>$odat") || die "cannot create file $odat\n\n";
  
  ## print header line
  print OUT "SNP\tA1\tA2\t";
  foreach (my $i = 0; $i < $numPeople; $i ++)
  {
    print OUT "$famIDs[$i]\t$subIDs[$i]\t";
  } 
  print OUT "\n";
  # end of printing header line
  
  # print one SNP per line
  for (my $j = 0; $j < $numSNPs; $j ++)
  {
    print OUT "$SNPs[$j]\t$alleleOnes[$j]\t$alleleTwos[$j]\t";
    for (my $i = 0; $i < $numPeople; $i ++)
    {
      my $indexone = $i * $numSNPs * 2 + $j * 2;
      my $indextwo = $indexone + 1;
      print OUT "$probs[$indexone]\t$probs[$indextwo]\t";
    } # end of person $i
    print OUT "\n";
  } # end of SNP j
  # end of printing one SNP per line
  close(OUT);
}





