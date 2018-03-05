#!/usr/bin/perl
#############################################
#                                           #
# Copyright (C) 2015                        #
#                                           #
# AptaIT GmbH                               #
# Zenettistr. 11                            #
# 80337 Munich                              #
# Germany                                   #
#                                           #
# Author                                    #
# Dr. Carsten Groeber                       #
# Implerstr. 48a                            #
# 81371 Munich                              #
# Germany                                   #
#                                           #
# Proprietary and confidential information. #
# Redistribution, modification and/or usage #
# is prohibited. All Rights Reserved.       #
#                                           #
#############################################

use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use IO::Uncompress::AnyUncompress;
use File::Spec;
use File::Temp;
use DBI;
use Tk;

local $| = 1;

######################
# Global data tables #
######################

# List of possible constant amino acids at end of V-region
my @globalVRegWordArray = ("YFCA", "YLCA", "YICA", "YRCA", "YYCA", "YLCS", "YICS", "YLCT", "YICV", "YLCV", "YHCA", "YQCA", "YLCG", "YFCG", "YLRA", "YLYS", "YYCL", "YFCV", "YYCI");
my $globalVRegWordNumber = scalar(@globalVRegWordArray);
for (my $i = 0; $i < $globalVRegWordNumber; $i++) {
  $globalVRegWordArray[$i] =~ s/X/./g;
}

# List of possible constant amino acids at front of J-region
my @globalJRegWordArray = ("FAXG", "FGXG", "WGXG", "VGXG", "LGXG");
my $globalJRegWordNumber = scalar(@globalJRegWordArray);
for (my $i = 0; $i < $globalJRegWordNumber; $i++) {
  $globalJRegWordArray[$i] =~ s/X/./g;
}

# Table of phred scores needed to decode fastq sequence quality
my @globalPhredScoresArray;
my $globalNumberPhredScores = 93;
my $globalBaseCharPhredScores = 33;
for (my $i = 0; $i < $globalNumberPhredScores; $i++) {
  $globalPhredScoresArray[$i] = 10.0**($i/-10.0);
}

# Hashes with key region-bases and value array of numbers of all sequences (starting from 0)
# in which region-bases occur, number of entries in arrays are total counts for each region. 
my %globalOptBarcodeSequence;
my %globalVRegionSequence;
my %globalCDR3RegionSequence;
my %globalJRegionSequence;

# Corresponding hashes if second chain is present
my %globalOptScndChainOptBarcodeSequence;
my %globalOptScndChainVRegionSequence;
my %globalOptScndChainCDR3RegionSequence;
my %globalOptScndChainJRegionSequence;

# Arrays of unique region-bases, index of array is the number of the sequence
# (starting from 0) in which the corresponding region-bases were found first.
my @globalUniqueOptBarcode;
my @globalUniqueVRegion;
my @globalUniqueCDR3Region;
my @globalUniqueJRegion;

# Array names from name data base, index of array is the number of the sequence
# (starting from 0) in which the corresponding region-bases were found first.
my @globalVRegionName;
my @globalJRegionName;

# Arrays of number of sequence (starting from 0) for which the region-bases of
# the corresponding sequence were found first, i.e. to query region-bases for
# sequences use $globalUniqueRegion[$globalSequenceRegion[$numberOfSequence]]
# and for names use $globalRegionName[$globalSequenceRegion[$numberOfSequence]].
my @globalSequenceOptBarcode;
my @globalSequenceVRegion;
my @globalSequenceCDR3Region;
my @globalSequenceJRegion;

# Optional hash with unique combinations of regions
my %globalRegionConnectivity;

# Additional information on data tables
my $globalOptBarcodeMaxLength = -1;
my $globalVRegionMaxLength = -1;
my $globalCDR3RegionMaxLength = -1;
my $globalJRegionMaxLength = -1;
my $globalOptBarcodeMinLength = -1;
my $globalVRegionMinLength = -1;
my $globalCDR3RegionMinLength = -1;
my $globalJRegionMinLength = -1;

# Number objects pushed to data base at once
my $globalMaxTermsInsertTable = 1;
# Seems to may be larger 1 for linux only
if ($^O =~ m/linux/i) {
  $globalMaxTermsInsertTable = 100;
}

# Including constant flanking amino acids
my $globalVRegionMinLen = 5;
my $globalJRegionMinLen = 5;

# Flags whether barcode marked by constant region exists
my $globalOptBarcodePresent = 0;

# If value > 0 number of sequence optional second chain
# from second file starts with (i.e. M first and N second)
# or if value == 0 flag that sequences 0 and 1, 2 and 3, ...
# were coupled chains in data file (i.e. 1:1 connection)
my $globalOptScndChainFirstSequence = -1;

# Total number of unique combinations of regions
my $globalTotalNumberConnect = 0;

# Total number sequences parsed
my $globalTotalNumberSequence = 0;

###################################
# Function to re-init global data #
###################################

sub initGlobalDataTables {
  if (scalar(%globalOptBarcodeSequence)) {
    %globalOptBarcodeSequence = ();
  }
  if (scalar(%globalVRegionSequence)) {
    %globalVRegionSequence = ();
  }
  if (scalar(%globalCDR3RegionSequence)) {
    %globalCDR3RegionSequence = ();
  }
  if (scalar(%globalJRegionSequence)) {
    %globalJRegionSequence = ();
  }
  if (scalar(%globalOptScndChainOptBarcodeSequence)) {
    %globalOptScndChainOptBarcodeSequence = ();
  }
  if (scalar(%globalOptScndChainVRegionSequence)) {
    %globalOptScndChainVRegionSequence = ();
  }
  if (scalar(%globalOptScndChainCDR3RegionSequence)) {
    %globalOptScndChainCDR3RegionSequence = ();
  }
  if (scalar(%globalOptScndChainJRegionSequence)) {
    %globalOptScndChainJRegionSequence = ();
  }
  if (scalar(@globalUniqueOptBarcode) > 0) {
    @globalUniqueOptBarcode = ();
  }
  if (scalar(@globalUniqueVRegion) > 0) {
    @globalUniqueVRegion = ();
  }
  if (scalar(@globalUniqueCDR3Region) > 0) {
    @globalUniqueCDR3Region = ();
  }
  if (scalar(@globalUniqueJRegion) > 0) {
    @globalUniqueJRegion = ();
  }
  if (scalar(@globalVRegionName) > 0) {
    @globalVRegionName = ();
  }
  if (scalar(@globalJRegionName) > 0) {
    @globalJRegionName = ();
  }
  if (scalar(@globalSequenceOptBarcode) > 0) {
    @globalSequenceOptBarcode = ();
  }
  if (scalar(@globalSequenceVRegion) > 0) {
    @globalSequenceVRegion = ();
  }
  if (scalar(@globalSequenceCDR3Region) > 0) {
    @globalSequenceCDR3Region = ();
  }
  if (scalar(@globalSequenceJRegion) > 0) {
    @globalSequenceJRegion = ();
  }
  if (scalar(%globalRegionConnectivity)) {
    %globalRegionConnectivity = ();
  }
  $globalOptBarcodeMaxLength = -1;
  $globalVRegionMaxLength = -1;
  $globalCDR3RegionMaxLength = -1;
  $globalJRegionMaxLength = -1;
  $globalOptBarcodeMinLength = -1;
  $globalVRegionMinLength = -1;
  $globalCDR3RegionMinLength = -1;
  $globalJRegionMinLength = -1;
  $globalOptBarcodePresent = 0;
  $globalOptScndChainFirstSequence = -1;
  $globalTotalNumberConnect = 0;
  $globalTotalNumberSequence = 0;
  return;
}

##################################
# Functions to query data tables #
##################################

sub barcodeBases {
  my $numberOfSequence = shift;
  if (($globalOptBarcodePresent != 0) && ($numberOfSequence >= 0) && ($numberOfSequence < $globalTotalNumberSequence)) {
    return $globalUniqueOptBarcode[$globalSequenceOptBarcode[$numberOfSequence]];
  }
  return "";
}

sub vRegionBases {
  my $numberOfSequence = shift;
  if (($numberOfSequence >= 0) && ($numberOfSequence < $globalTotalNumberSequence)) {
    return $globalUniqueVRegion[$globalSequenceVRegion[$numberOfSequence]];
  }
  return "";
}

sub cdr3RegionBases {
  my $numberOfSequence = shift;
  if (($numberOfSequence >= 0) && ($numberOfSequence < $globalTotalNumberSequence)) {
    return $globalUniqueCDR3Region[$globalSequenceCDR3Region[$numberOfSequence]];
  }
  return "";
}

sub jRegionBases {
  my $numberOfSequence = shift;
  if (($numberOfSequence >= 0) && ($numberOfSequence < $globalTotalNumberSequence)) {
    return $globalUniqueJRegion[$globalSequenceJRegion[$numberOfSequence]];
  }
  return "";
}

sub vRegionName {
  my $sequence = shift;
  my $optChainNumber = shift;
  if (defined($optChainNumber) && ($optChainNumber != 0)) {
    if (exists($globalOptScndChainVRegionSequence{$sequence})) {
      return $globalVRegionName[$globalSequenceVRegion[@{$globalOptScndChainVRegionSequence{$sequence}}[0]]];
    }
  } else {
    if (exists($globalVRegionSequence{$sequence})) {
      return $globalVRegionName[$globalSequenceVRegion[@{$globalVRegionSequence{$sequence}}[0]]];
    }
  }
  return "";
}

sub jRegionName {
  my $sequence = shift;
  my $optChainNumber = shift;
  if (defined($optChainNumber) && ($optChainNumber != 0)) {
    if (exists($globalOptScndChainJRegionSequence{$sequence})) {
      return $globalJRegionName[$globalSequenceJRegion[@{$globalOptScndChainJRegionSequence{$sequence}}[0]]];
    }
  } else {
    if (exists($globalJRegionSequence{$sequence})) {
      return $globalJRegionName[$globalSequenceJRegion[@{$globalJRegionSequence{$sequence}}[0]]];
    }
  }
  return "";
}

sub barcodeNumber {
  my $sequence = shift;
  my $optChainNumber = shift;
  if ($globalOptBarcodePresent != 0) {
    if (defined($optChainNumber) && ($optChainNumber != 0)) {
      if (exists($globalOptScndChainOptBarcodeSequence{$sequence})) {
	return scalar(@{$globalOptScndChainOptBarcodeSequence{$sequence}});
      }
    } else {
      if (exists($globalOptBarcodeSequence{$sequence})) {
	return scalar(@{$globalOptBarcodeSequence{$sequence}});
      }
    }
  }
  return -1;
}

sub vRegionNumber {
  my $sequence = shift;
  my $optChainNumber = shift;
  if (defined($optChainNumber) && ($optChainNumber != 0)) {
    if (exists($globalOptScndChainVRegionSequence{$sequence})) {
      return scalar(@{$globalOptScndChainVRegionSequence{$sequence}});
    }
  } else {
    if (exists($globalVRegionSequence{$sequence})) {
      return scalar(@{$globalVRegionSequence{$sequence}});
    }
  }
  return -1;
}

sub cdr3RegionNumber {
  my $sequence = shift;
  my $optChainNumber = shift;
  if (defined($optChainNumber) && ($optChainNumber != 0)) {
    if (exists($globalOptScndChainCDR3RegionSequence{$sequence})) {
      return scalar(@{$globalOptScndChainCDR3RegionSequence{$sequence}});
    }
  } else {
    if (exists($globalCDR3RegionSequence{$sequence})) {
      return scalar(@{$globalCDR3RegionSequence{$sequence}});
    }
  }
  return -1;
}

sub jRegionNumber {
  my $sequence = shift;
  my $optChainNumber = shift;
  if (defined($optChainNumber) && ($optChainNumber != 0)) {
    if (exists($globalOptScndChainJRegionSequence{$sequence})) {
      return scalar(@{$globalOptScndChainJRegionSequence{$sequence}});
    }
  } else {
    if (exists($globalJRegionSequence{$sequence})) {
      return scalar(@{$globalJRegionSequence{$sequence}});
    }
  }
  return -1;
}

sub barcodeSequences {
  my $sequence = shift;
  my $optChainNumber = shift;
  if ($globalOptBarcodePresent != 0) {
    if (defined($optChainNumber) && ($optChainNumber != 0)) {
      if (exists($globalOptScndChainOptBarcodeSequence{$sequence})) {
	return \@{$globalOptScndChainOptBarcodeSequence{$sequence}};
      }
    } else {
      if (exists($globalOptBarcodeSequence{$sequence})) {
	return \@{$globalOptBarcodeSequence{$sequence}};
      }
    }
  }
  return [];
}

sub vRegionSequences {
  my $sequence = shift;
  my $optChainNumber = shift;
  if (defined($optChainNumber) && ($optChainNumber != 0)) {
    if (exists($globalOptScndChainVRegionSequence{$sequence})) {
      return \@{$globalOptScndChainVRegionSequence{$sequence}};
    }
  } else {
    if (exists($globalVRegionSequence{$sequence})) {
      return \@{$globalVRegionSequence{$sequence}};
    }
  }
  return [];
}

sub cdr3RegionSequences {
  my $sequence = shift;
  my $optChainNumber = shift;
  if (defined($optChainNumber) && ($optChainNumber != 0)) {
    if (exists($globalOptScndChainCDR3RegionSequence{$sequence})) {
      return \@{$globalOptScndChainCDR3RegionSequence{$sequence}};
    }
  } else {
    if (exists($globalCDR3RegionSequence{$sequence})) {
      return \@{$globalCDR3RegionSequence{$sequence}};
    }
  }
  return [];
}

sub jRegionSequences {
  my $sequence = shift;
  my $optChainNumber = shift;
  if (defined($optChainNumber) && ($optChainNumber != 0)) {
    if (exists($globalOptScndChainJRegionSequence{$sequence})) {
      return \@{$globalOptScndChainJRegionSequence{$sequence}};
    }
  } else {
    if (exists($globalJRegionSequence{$sequence})) {
      return \@{$globalJRegionSequence{$sequence}};
    }
  }
  return [];
}

sub numberBarcode {
  if ($globalOptBarcodePresent != 0) {
    return scalar(@globalUniqueOptBarcode);
  }
  return 0;
}

sub numberVRegion {
  return scalar(@globalUniqueVRegion);
}

sub numberCdr3Region {
  return scalar(@globalUniqueCDR3Region);
}

sub numberJRegion {
  return scalar(@globalUniqueJRegion);
}

##################################
# Function to build connectivity #
##################################

sub buildConnectivity {
  for (my $i = 0; ($i < $globalTotalNumberSequence) && (($globalOptScndChainFirstSequence != 0) || (($i + 1) < $globalTotalNumberSequence)); $i++) {
    my $key = "";
    if ($globalOptScndChainFirstSequence > 0) {
      if ($i < $globalOptScndChainFirstSequence) {
	$key .= "ALPHA ";
      } else {
	$key .= "BETA ";
      }
      if ($globalOptBarcodePresent != 0) {
	$key .= barcodeBases($i) . " ";
      }
    } elsif ($globalOptScndChainFirstSequence == 0) {
      if ($globalOptBarcodePresent != 0) {
	$key .= barcodeBases($i) . " ";
      }
      $key .= "ALPHA ";
    } elsif ($globalOptBarcodePresent != 0) {
      $key .= barcodeBases($i) . " ";
    }
    $key .= vRegionBases($i) . " " . cdr3RegionBases($i) . " " . jRegionBases($i);
    if ($globalOptScndChainFirstSequence == 0) {
      $i++;
      $key .= " BETA " . vRegionBases($i) . " " . cdr3RegionBases($i) . " " . jRegionBases($i);
    }
    if (!exists($globalRegionConnectivity{$key})) {
      $globalRegionConnectivity{$key} = 0;
    }
    $globalRegionConnectivity{$key}++;
  }
  $globalTotalNumberConnect = scalar(keys(%globalRegionConnectivity));
  return;
}

##################################
# Functions to print data tables #
##################################

sub vRegionOfCdr3 {
  my $sequence = shift;
  my $optChainNumber = shift;
  my $sequenceArray;
  if (defined($optChainNumber) && ($optChainNumber != 0)) {
    $sequenceArray = cdr3RegionSequences($sequence, $optChainNumber);
  } else {
    $sequenceArray = cdr3RegionSequences($sequence);
  }
  if (scalar(@$sequenceArray) > 0) {
    my %regionHash;
    for (my $i = 0; $i < scalar(@$sequenceArray); $i++) {
      my $bases = vRegionBases($sequenceArray->[$i]);
      if (exists($regionHash{$bases})) {
	$regionHash{$bases} += 1;
      } else {
	$regionHash{$bases} = 1;
      }
    }
    my $result = "";
    foreach my $key (sort{$regionHash{$b} <=> $regionHash{$a}} keys(%regionHash)) {
      if ($result ne "") {
	$result .= "\n";
      }
      $result .= $key;
      my $name;
      if (defined($optChainNumber) && ($optChainNumber != 0)) {
	$name = vRegionName($key, $optChainNumber);
      } else {
	$name = vRegionName($key);
      }
      if ($name ne "") {
	$result .= " (" . $name . ")";
      }
      $result .= " " . $regionHash{$key};
    }
    return $result;
  }
  return "";
}

sub jRegionOfCdr3 {
  my $sequence = shift;
  my $optChainNumber = shift;
  my $sequenceArray;
  if (defined($optChainNumber) && ($optChainNumber != 0)) {
    $sequenceArray = cdr3RegionSequences($sequence, $optChainNumber);
  } else {
    $sequenceArray = cdr3RegionSequences($sequence);
  }
  if (scalar(@$sequenceArray) > 0) {
    my %regionHash;
    for (my $i = 0; $i < scalar(@$sequenceArray); $i++) {
      my $bases = jRegionBases($sequenceArray->[$i]);
      if (exists($regionHash{$bases})) {
	$regionHash{$bases} += 1;
      } else {
	$regionHash{$bases} = 1;
      }
    }
    my $result = "";
    foreach my $key (sort{$regionHash{$b} <=> $regionHash{$a}} keys(%regionHash)) {
      if ($result ne "") {
	$result .= "\n";
      }
      $result .= $key;
      my $name;
      if (defined($optChainNumber) && ($optChainNumber != 0)) {
	$name = jRegionName($key, $optChainNumber);
      } else {
	$name = jRegionName($key);
      }
      if ($name ne "") {
	$result .= " (" . $name . ")";
      }
      $result .= " " . $regionHash{$key};
    }
    return $result;
  }
  return "";
}

sub printBarcodeSequences {
  if ($globalOptBarcodePresent != 0) {
    print "*MESSAGE* Barcode bases vs. original valid sequence number:\n";
    foreach my $key (keys(%globalOptBarcodeSequence)) {
      if ($globalOptScndChainFirstSequence >= 0) {
	print "ALPHA ";
      }
      print "$key";
      my $sequences = barcodeSequences($key);
      for (my $i = 0; $i < scalar(@$sequences); $i++) {
	print " ", $sequences->[$i] + 1;
      }
      print "\n";
    }
    if ($globalOptScndChainFirstSequence >= 0) {
      foreach my $key (keys(%globalOptScndChainOptBarcodeSequence)) {
	print "BETA $key";
	my $sequences = barcodeSequences($key, 1);
	for (my $i = 0; $i < scalar(@$sequences); $i++) {
	  print " ", $sequences->[$i] + 1;
	}
	print "\n";
      }
    }
    print "\n";
  }
  return;
}

sub printVRegionSequences {
  print "*MESSAGE* V-region bases vs. original valid sequence number:\n";
  foreach my $key (keys(%globalVRegionSequence)) {
    if ($globalOptScndChainFirstSequence >= 0) {
      print "ALPHA ";
    }
    print "$key";
    my $sequences = vRegionSequences($key);
    my $name = vRegionName($key);
    if ($name ne "") {
      print " (", $name, ")";
    }
    for (my $i = 0; $i < scalar(@$sequences); $i++) {
      print " ", $sequences->[$i] + 1;
    }
    print "\n";
  }
  if ($globalOptScndChainFirstSequence >= 0) {
    foreach my $key (keys(%globalOptScndChainVRegionSequence)) {
      print "BETA $key";
      my $sequences = vRegionSequences($key, 1);
      my $name = vRegionName($key, 1);
      if ($name ne "") {
	print " (", $name, ")";
      }
      for (my $i = 0; $i < scalar(@$sequences); $i++) {
	print " ", $sequences->[$i] + 1;
      }
      print "\n";
    }
  }
  print "\n";
  return;
}

sub printCDR3RegionSequences {
  print "*MESSAGE* CDR3-region bases vs. original valid sequence number:\n";
  foreach my $key (keys(%globalCDR3RegionSequence)) {
    if ($globalOptScndChainFirstSequence >= 0) {
      print "ALPHA ";
    }
    print "$key";
    my $sequences = cdr3RegionSequences($key);
    for (my $i = 0; $i < scalar(@$sequences); $i++) {
      print " ", $sequences->[$i] + 1;
    }
    print "\n";
  }
  if ($globalOptScndChainFirstSequence >= 0) {
    foreach my $key (keys(%globalOptScndChainCDR3RegionSequence)) {
      print "BETA $key";
      my $sequences = cdr3RegionSequences($key, 1);
      for (my $i = 0; $i < scalar(@$sequences); $i++) {
	print " ", $sequences->[$i] + 1;
      }
      print "\n";
    }
  }
  print "\n";
  return;
}

sub printJRegionSequences {
  print "*MESSAGE* J-region bases vs. original valid sequence number:\n";
  foreach my $key (keys(%globalJRegionSequence)) {
    if ($globalOptScndChainFirstSequence >= 0) {
      print "ALPHA ";
    }
    print "$key";
    my $sequences = jRegionSequences($key);
    my $name = jRegionName($key);
    if ($name ne "") {
      print " (", $name, ")";
    }
    for (my $i = 0; $i < scalar(@$sequences); $i++) {
      print " ", $sequences->[$i] + 1;
    }
    print "\n";
  }
  if ($globalOptScndChainFirstSequence >= 0) {
    foreach my $key (keys(%globalOptScndChainJRegionSequence)) {
      print "BETA $key";
      my $sequences = jRegionSequences($key, 1);
      my $name = jRegionName($key, 1);
      if ($name ne "") {
	print " (", $name, ")";
      }
      for (my $i = 0; $i < scalar(@$sequences); $i++) {
	print " ", $sequences->[$i] + 1;
      }
      print "\n";
    }
  }
  print "\n";
  return;
}

sub printConnectivity {
  if (!scalar(%globalRegionConnectivity) || ($globalTotalNumberConnect == 0)) {
    buildConnectivity();
  }
  print "*MESSAGE* Unique combinations of regions with total counts:\n";
  my $number = 1;
  foreach my $combination (keys(%globalRegionConnectivity)) {
    my $numberSequence = $globalRegionConnectivity{$combination};
    print "Combination ", $number++;
    my $optScndChain = "";
    if ($globalOptScndChainFirstSequence == 0) {
      $optScndChain = $combination;
      $optScndChain =~ s/^.+ BETA //;
      $combination =~ s/ BETA.+$//;
    }
    my $optChainIdentifier = "";
    if ($globalOptScndChainFirstSequence > 0) {
      $combination =~ s/^([^ ]*) //;
      $optChainIdentifier = $1;
      print " ", $optChainIdentifier;
    }
    if ($globalOptBarcodePresent != 0) {
      $combination =~ s/^([^ ]*) //;
      print " ", $1;
    }
    if ($globalOptScndChainFirstSequence == 0) {
      $combination =~ s/^([^ ]*) //;
      $optChainIdentifier = $1;
      print " ", $optChainIdentifier;
    }
    my ($vBases, $cdr3Bases, $jBases) = split(' ', $combination);
    print " ", $vBases;
    my $vName = "";
    if ($optChainIdentifier eq "BETA") {
      $vName = vRegionName($vBases, 1);
    } else {
      $vName = vRegionName($vBases);
    }
    if ($vName ne "") {
      print " (", $vName, ")";
    }
    print " ", $cdr3Bases, " ", $jBases;
    my $jName = "";
    if ($optChainIdentifier eq "BETA") {
      $jName = jRegionName($jBases, 1);
    } else {
      $jName = jRegionName($jBases);
    }
    if ($jName ne "") {
      print " (", $jName, ")";
    }
    if ($optScndChain ne "") {
      my ($optScndChainVBases, $optScndChainCdr3Bases, $optScndChainJBases) = split(' ', $optScndChain);
      print " BETA ", $optScndChainVBases;
      my $optScndChainVName = vRegionName($optScndChainVBases, 1);
      if ($optScndChainVName ne "") {
	print " (", $optScndChainVName, ")";
      }
      print " ", $optScndChainCdr3Bases, " ", $optScndChainJBases;
      my $optScndChainJName = jRegionName($optScndChainJBases, 1);
      if ($optScndChainJName ne "") {
	print " (", $optScndChainJName, ")";
      }
    }
    print " ", $numberSequence . "\n";
  }
  print "\n";
  return;
}

sub printDataTable {
  print "*MESSAGE* Sequence V-, CDR3- and J-regions with total counts:\n";
  for (my $i = 0; $i < $globalTotalNumberSequence; $i++) {
    print "Sequence ", $i + 1;
    if ($globalOptScndChainFirstSequence > 0) {
      if ($i < $globalOptScndChainFirstSequence) {
	print " ALPHA";
      } else {
	print " BETA";
      }
    } elsif ($globalOptScndChainFirstSequence == 0) {
      if (($i % 2) == 0) {
	print " ALPHA";
      } else {
	print " BETA";
      }
    }
    if ($globalOptBarcodePresent != 0) {
      my $totalNumberBarcode = 0;
      if ((($globalOptScndChainFirstSequence > 0) && ($i >= $globalOptScndChainFirstSequence)) ||
	  (($globalOptScndChainFirstSequence == 0) && (($i % 2) != 0))) {
	$totalNumberBarcode = barcodeNumber(barcodeBases($i), 1);
      } else {
	$totalNumberBarcode = barcodeNumber(barcodeBases($i));
      }
      print " ", barcodeBases($i), " (", $totalNumberBarcode, ")";
    }
    my $totalNumberVRegion = 0;
    my $totalNumberCDR3Region = 0;
    my $totalNumberJRegion = 0;
    if ((($globalOptScndChainFirstSequence > 0) && ($i >= $globalOptScndChainFirstSequence)) ||
	(($globalOptScndChainFirstSequence == 0) && (($i % 2) != 0))) {
      $totalNumberVRegion = vRegionNumber(vRegionBases($i), 1);
      $totalNumberCDR3Region = cdr3RegionNumber(cdr3RegionBases($i), 1);
      $totalNumberJRegion = jRegionNumber(jRegionBases($i), 1);
    } else {
      $totalNumberVRegion = vRegionNumber(vRegionBases($i));
      $totalNumberCDR3Region = cdr3RegionNumber(cdr3RegionBases($i));
      $totalNumberJRegion = jRegionNumber(jRegionBases($i));
    }
    print " ", vRegionBases($i), " (", $totalNumberVRegion, ") ", 
      cdr3RegionBases($i), " (", $totalNumberCDR3Region, ") ", 
	jRegionBases($i), " (", $totalNumberJRegion, ")\n";
  }
  print "\n";
  return;
}

####################
# Helper functions #
####################

sub fuzzyDiff {
  my $string1 = shift;
  my $string2 = shift;
  my $maxNumWrong = shift;

  if (length($string1) != length($string2)) {
    die "*ERROR* Cannot fuzzy compare strings with different length\n";
  }

  my $numberWrong = 0;
  for (my $i = 0; ($numberWrong <= $maxNumWrong) && ($i < length($string1)); $i++) {
    if (substr($string1, $i, 1) ne substr($string2, $i, 1)) {
      $numberWrong++;
    }
  }

  return $numberWrong;
}

sub fuzzyMatchOne {
  my $testedString = shift;
  my $fuzzyString = shift;

  if (length($testedString) != length($fuzzyString)) {
    die "*ERROR* Cannot fuzzy compare strings with different length\n";
  }

  if (length($testedString) < 5) {
    die "*ERROR* Cannot fuzzy compare strings of length smaller 5 with 1 error\n";
  }

  my $halfLen = length($fuzzyString) / 2;
  my $half1 = substr($fuzzyString, 0, $halfLen);
  my $half2 = $fuzzyString;
  $half2 =~ s/^$half1//;
  if ($testedString =~ m/^$half1/) {
    my $minusLen = length($half2) * -1;
    if (fuzzyDiff(substr($testedString, $minusLen), $half2, 1) <= 1) {
      return 0;
    }
  } elsif ($testedString =~ m/$half2$/) {
    if (fuzzyDiff(substr($testedString, 0, $halfLen), $half1, 1) <= 1) {
      return 0;
    }
  }

  return 1;
}

sub fuzzyMatchTwo {
  my $testedString = shift;
  my $fuzzyString = shift;

  if (length($testedString) != length($fuzzyString)) {
    die "*ERROR* Cannot fuzzy compare strings with different length\n";
  }

  if (length($testedString) < 8) {
    die "*ERROR* Cannot fuzzy compare strings of length smaller 8 with 2 error\n";
  }

  my $third2Len = length($fuzzyString) / 3;
  my $third1Len = $third2Len + ((length($fuzzyString) % 3) / 2);
  my $third1 = substr($fuzzyString, 0, $third1Len);
  my $third2 = substr($fuzzyString, $third1Len, $third2Len);
  my $third3 = $fuzzyString;
  $third3 =~ s/^$third1//;
  $third3 =~ s/^$third2//;
  if ($testedString =~ m/^$third1/) {
    my $numberWrong = fuzzyDiff(substr($testedString, $third1Len, $third2Len), $third2, 2);
    if ($numberWrong <= 2) {
      my $minusLen = length($third3) * -1;
      if (($numberWrong + fuzzyDiff(substr($testedString, $minusLen), $third3, 2)) <= 2) {
	return 0;
      }
    }
  } elsif ($testedString =~ m/$third3$/) {
    my $numberWrong = fuzzyDiff(substr($testedString, 0, $third1Len), $third1, 2);
    if ($numberWrong <= 2) {
      if (($numberWrong + fuzzyDiff(substr($testedString, $third1Len, $third2Len), $third2, 2)) <= 2) {
	return 0;
      }
    }
  } elsif (substr($testedString, $third1Len, $third2Len) =~ m/^$third2$/) {
    my $numberWrong = fuzzyDiff(substr($testedString, 0, $third1Len), $third1, 2);
    if ($numberWrong <= 2) {
      my $minusLen = length($third3) * -1;
      if (($numberWrong + fuzzyDiff(substr($testedString, $minusLen), $third3, 2)) <= 2) {
	return 0;
      }
    }
  }

  return 1;
}

sub printProgress {
  my $currentNumber = shift;
  my $totalNumber = shift;

  if (($totalNumber >= 10) && ($currentNumber % ($totalNumber / 10)) == 0) {
    printf "%.0f%%", ($currentNumber / $totalNumber) * 100;
  } elsif (($totalNumber >= 100) && ($currentNumber % ($totalNumber / 100)) == 0) {
    print ".";
  }

  return;
}

####################################################################
# Functions to translate and build reverse complement of sequences #
####################################################################

sub reverseComplements {
  my $bases = shift;
  my $result = "";
  my $numBases = length($bases);
  for (my $i = 0; $i < $numBases; $i++) {
    my $base = uc(substr($bases, $numBases - $i - 1, 1));
    if ($base eq "A") {
      $result .= "T";
    } elsif ($base eq "C") {
      $result .= "G";
    } elsif ($base eq "G") {
      $result .= "C";
    } elsif ($base eq "T") {
      $result .= "A";
    } else {
      $result .= "N";
    }
  }
  return $result;
}

sub translateBasesAminoForward {
  my $base = shift;
  $base = uc($base);
  if (($base eq "GCT") || ($base eq "GCC") || ($base eq "GCA") || ($base eq "GCG")) {
    return "A";
  } elsif (($base eq "TGT") || ($base eq "TGC")) {
    return "C";
  } elsif (($base eq "GAT") || ($base eq "GAC")) {
    return "D";
  } elsif (($base eq "GAA") || ($base eq "GAG")) {
    return "E";
  } elsif (($base eq "TTT") || ($base eq "TTC")) {
    return "F";
  } elsif (($base eq "GGT") || ($base eq "GGC") || ($base eq "GGA") || ($base eq "GGG")) {
    return "G";
  } elsif (($base eq "CAT") || ($base eq "CAC")) {
    return "H";
  } elsif (($base eq "ATT") || ($base eq "ATC") || ($base eq "ATA")) {
    return "I";
  } elsif (($base eq "AAA") || ($base eq "AAG")) {
    return "K";
  } elsif (($base eq "TTA") || ($base eq "TTG") || ($base eq "CTT") || ($base eq "CTC") || ($base eq "CTA") || ($base eq "CTG")) {
    return "L";
  } elsif ($base eq "ATG") {
    return "M";
  } elsif (($base eq "AAT") || ($base eq "AAC")) {
    return "N";
  } elsif (($base eq "CCT") || ($base eq "CCC") || ($base eq "CCA") || ($base eq "CCG")) {
    return "P";
  } elsif (($base eq "CAA") || ($base eq "CAG")) {
    return "Q";
  } elsif (($base eq "CGT") || ($base eq "CGC") || ($base eq "CGA") || ($base eq "CGG") || ($base eq "AGA") || ($base eq "AGG")) {
    return "R";
  } elsif (($base eq "TCT") || ($base eq "TCC") || ($base eq "TCA") || ($base eq "TCG") || ($base eq "AGT") || ($base eq "AGC")) {
    return "S";
  } elsif (($base eq "ACT") || ($base eq "ACA") || ($base eq "ACC") || ($base eq "ACG")) {
    return "T";
  } elsif (($base eq "GTT") || ($base eq "GTC") || ($base eq "GTA") || ($base eq "GTG")) {
    return "V";
  } elsif ($base eq "TGG") {
    return "W";
  } elsif (($base eq "TAC") || ($base eq "TAT")) {
    return "Y";
  }
  return "X";
}

##############################################
# Functions to parse human TCR sequence data #
##############################################

sub localFileWithPath {
  my $pathScript = shift;
  my $fileName = shift;
  my $systemVolume = "";
  my $directories = "";
  my $optFileTmp = "";
  ($systemVolume, $directories, $optFileTmp) = File::Spec->splitpath($pathScript);
  return File::Spec->catfile(File::Spec->splitdir($directories), $fileName);
}

sub bareFileName {
  my $inputFile = shift;
  my $systemVolume = "";
  my $directories = "";
  my $bareFileName = "";
  ($systemVolume, $directories, $bareFileName) = File::Spec->splitpath($inputFile);
  return $bareFileName;
}

sub mergePairedEnd {
  my $pathScript = shift;
  my $inputFile1 = shift;
  my $inputFile2 = shift;
  my $mergeTool = shift;
  my $numThread = shift;

  my $mergeToolExec = $mergeTool;
  if ($mergeToolExec eq "") {
    $mergeToolExec = localFileWithPath($pathScript, "pear");
  }

  if (! -x $mergeToolExec) {
    die "*ERROR* Cannot find paired end merge tool $mergeToolExec\n";
  }

  # Do not care for compressed files, pear handels those itself
  my $bareFile1 = bareFileName($inputFile1);
  my $bareFile2 = bareFileName($inputFile2);

  my $assembledBase = $bareFile1 . "_" . $bareFile2;
  my $execStrg = $mergeToolExec . " -f " . $inputFile1 . " -r " . $inputFile2 . " -o " . $assembledBase;
  if ($numThread > 0) {
    $execStrg .= " -j " . $numThread;
  }
  system($execStrg);

  my $assembledFile = $assembledBase . ".assembled.fastq";
  if (! -s $assembledFile) {
    die "*ERROR* Paired end output $assembledFile does not exist or is not readable\n";
  }

  return $assembledFile;
}

sub discardBadQuality {
  my $phredScores = shift;

  # Discard sequence if probability for error summed over any 10 characters exceeds 1.0
  my @lastPhredScores;
  my $maximumError = 1.0;
  my $numberChecked = 0;
  my $checkLength = 10;
  my $checkIndex = 0;
  my $sumPErr = 0.0;

  for (my $i = 0; $i < length($phredScores); $i++) {
    my $index = ord(substr($phredScores, $i, 1)) - $globalBaseCharPhredScores;
    if ($index < 0) {
      return 1;
    } elsif ($index > $globalNumberPhredScores) {
      $index = $globalNumberPhredScores - 1;
    }
    $lastPhredScores[$checkIndex] = $index;
    $sumPErr += $globalPhredScoresArray[$index];
    $numberChecked++;
    $checkIndex++;
    if ($numberChecked >= $checkLength) {
      if ($sumPErr >= $maximumError) {
	return 1;
      } else {
	if ($checkIndex == $checkLength) {
	  $sumPErr -= $globalPhredScoresArray[$lastPhredScores[0]];
	} else {
	  $sumPErr -= $globalPhredScoresArray[$lastPhredScores[$checkIndex]];
	}
      }
    }
    if ($checkIndex == $checkLength) {
      $checkIndex = 0;
    }
  }
  if (($numberChecked < $checkLength) &&
      ($sumPErr >= ($maximumError * $numberChecked / $checkLength))) {
    return 1;
  }

  return 0;
}

sub parseLineFrames {
  my $firstSequenceBarcode = shift;
  my $firstSequenceVRegion = shift;
  my $firstSequenceCDR3Region = shift;
  my $firstSequenceJRegion = shift;
  my $chainDelimiterConstReg = shift;
  my $wellBarcodeConstReg = shift;
  my $inputDataLine = shift;
  my $optChainNumber = shift;

  my $foundData = 0;

  for (my $offSet = 0; ($foundData == 0) && ($offSet < 3); $offSet++) {
    my $translatedLine = "";
    for (my $pos = $offSet; $pos < (length($inputDataLine) - 2); $pos += 3) {
      my $bases = substr($inputDataLine, $pos, 3);
      $translatedLine .= translateBasesAminoForward($bases);
    }
    if ($translatedLine !~ m/X/) {
      my $barcode = "";
      if ($globalOptBarcodePresent != 0) {
	if ($translatedLine =~ m/^(.+)$wellBarcodeConstReg/) {
	  $barcode = $1;
	  $translatedLine =~ s/^.+$wellBarcodeConstReg//;
	}
      }
      if (($globalOptBarcodePresent == 0) || ($barcode ne "")) {
	my $foundVLen = -1;
	for (my $i = 0; ($foundVLen < 0) && ($i < $globalVRegWordNumber); $i++) {
	  if ($translatedLine =~ m/$globalVRegWordArray[$i]/) {
	    $translatedLine =~ s/($globalVRegWordArray[$i])/$1 /;
	    $foundVLen = $+[1];
	  }
	}
	if ($foundVLen >= $globalVRegionMinLen) {
	  my $foundJLen = -1;
	  my $foundJPos = -1;
	  for (my $i = 0; ($foundJPos < 0) && ($i < $globalJRegWordNumber); $i++) {
	    if ($translatedLine =~ m/$globalJRegWordArray[$i]/) {
	      $translatedLine =~ s/($globalJRegWordArray[$i])/ $1/;
	      $foundJPos = $-[1];
	      $foundJLen = length($translatedLine) - $foundJPos - 1;
	    }
	  }
	  if ($foundJPos > $foundVLen) {
	    my $scndChainFoundVLen = -1;
	    my $scndChainFoundJLen = -1;
	    if (($chainDelimiterConstReg ne "") && ($translatedLine =~ m/^[^ ]+ [^ ]+ [^ ]+$chainDelimiterConstReg[^ ]+$/)) {
	      $translatedLine =~ s/^([^ ]+ [^ ]+ [^ ]+)$chainDelimiterConstReg([^ ]+)$/$1 $2/;
	      $foundJLen = $+[1] - $foundJPos + 1;
	      for (my $i = 0; ($scndChainFoundVLen < 0) && ($i < $globalVRegWordNumber); $i++) {
		if ($translatedLine =~ m/^[^ ]+ [^ ]+ [^ ]+ [^ ]+$globalVRegWordArray[$i]/) {
		  $translatedLine =~ s/^([^ ]+ [^ ]+ [^ ]+ )([^ ]+$globalVRegWordArray[$i])/$1$2 /;
		  $scndChainFoundVLen = $+[2] - $-[2] + 1;
		}
	      }
	      for (my $i = 0; ($scndChainFoundJLen < 0) && ($i < $globalJRegWordNumber); $i++) {
		if ($translatedLine =~ m/^[^ ]+ [^ ]+ [^ ]+ [^ ]+ [^ ]+$globalJRegWordArray[$i]/) {
		  $translatedLine =~ s/^([^ ]+ [^ ]+ [^ ]+ [^ ]+ [^ ]+)($globalJRegWordArray[$i])/$1 $2/;
		  $scndChainFoundJLen = length($translatedLine) - $-[2] - 1;
		}
	      }
	    }
	    if ((($chainDelimiterConstReg eq "") || (($scndChainFoundVLen >= $globalVRegionMinLen) && ($scndChainFoundJLen >= $globalJRegionMinLen))) &&
		($foundJLen >= $globalJRegionMinLen)) {
	      my @region = split(' ', $translatedLine);
	      # In order to prevent empty CDR3s
	      if ((scalar(@region) == 3) || (scalar(@region) == 6)) {
		if (defined($optChainNumber) && ($optChainNumber > 0)) {
		  if (($globalOptBarcodePresent != 0) && ($barcode ne "")) {
		    if (exists($globalOptScndChainOptBarcodeSequence{$barcode})) {
		      push(@{$globalOptScndChainOptBarcodeSequence{$barcode}}, $globalTotalNumberSequence);
		    } else {
		      $globalOptScndChainOptBarcodeSequence{$barcode} = [$globalTotalNumberSequence];
		      push(@globalUniqueOptBarcode, $barcode);
		      $firstSequenceBarcode->{$barcode} = scalar(@globalUniqueOptBarcode) - 1;
		      my $regLen = length($barcode);
		      if (($globalOptBarcodeMaxLength < 0) || ($regLen > $globalOptBarcodeMaxLength)) {
			$globalOptBarcodeMaxLength = $regLen;
		      }
		      if (($globalOptBarcodeMinLength < 0) || ($regLen < $globalOptBarcodeMinLength)) {
			$globalOptBarcodeMinLength = $regLen;
		      }
		    }
		    push(@globalSequenceOptBarcode, $firstSequenceBarcode->{$barcode});
		  }
		  if (exists($globalOptScndChainVRegionSequence{$region[0]})) {
		    push(@{$globalOptScndChainVRegionSequence{$region[0]}}, $globalTotalNumberSequence);
		  } else {
		    $globalOptScndChainVRegionSequence{$region[0]} = [$globalTotalNumberSequence];
		    push(@globalUniqueVRegion, $region[0]);
		    $firstSequenceVRegion->{$region[0]} = scalar(@globalUniqueVRegion) - 1;
		    push(@globalVRegionName, "");
		    my $regLen = length($region[0]);
		    if (($globalVRegionMaxLength < 0) || ($regLen > $globalVRegionMaxLength)) {
		      $globalVRegionMaxLength = $regLen;
		    }
		    if (($globalVRegionMinLength < 0) || ($regLen < $globalVRegionMinLength)) {
		      $globalVRegionMinLength = $regLen;
		    }
		  }
		  push(@globalSequenceVRegion, $firstSequenceVRegion->{$region[0]});
		  if (exists($globalOptScndChainCDR3RegionSequence{$region[1]})) {
		    push(@{$globalOptScndChainCDR3RegionSequence{$region[1]}}, $globalTotalNumberSequence);
		  } else {
		    $globalOptScndChainCDR3RegionSequence{$region[1]} = [$globalTotalNumberSequence];
		    push(@globalUniqueCDR3Region, $region[1]);
		    $firstSequenceCDR3Region->{$region[1]} = scalar(@globalUniqueCDR3Region) - 1;
		    my $regLen = length($region[1]);
		    if (($globalCDR3RegionMaxLength < 0) || ($regLen > $globalCDR3RegionMaxLength)) {
		      $globalCDR3RegionMaxLength = $regLen;
		    }
		    if (($globalCDR3RegionMinLength < 0) || ($regLen < $globalCDR3RegionMinLength)) {
		      $globalCDR3RegionMinLength = $regLen;
		    }
		  }
		  push(@globalSequenceCDR3Region, $firstSequenceCDR3Region->{$region[1]});
		  if (exists($globalOptScndChainJRegionSequence{$region[2]})) {
		    push(@{$globalOptScndChainJRegionSequence{$region[2]}}, $globalTotalNumberSequence);
		  } else {
		    $globalOptScndChainJRegionSequence{$region[2]} = [$globalTotalNumberSequence];
		    push(@globalUniqueJRegion, $region[2]);
		    $firstSequenceJRegion->{$region[2]} = scalar(@globalUniqueJRegion) - 1;
		    push(@globalJRegionName, "");
		    my $regLen = length($region[2]);
		    if (($globalJRegionMaxLength < 0) || ($regLen > $globalJRegionMaxLength)) {
		      $globalJRegionMaxLength = $regLen;
		    }
		    if (($globalJRegionMinLength < 0) || ($regLen < $globalJRegionMinLength)) {
		      $globalJRegionMinLength = $regLen;
		    }
		  }
		  push(@globalSequenceJRegion, $firstSequenceJRegion->{$region[2]});
		} else {
		  if (($globalOptBarcodePresent != 0) && ($barcode ne "")) {
		    if (exists($globalOptBarcodeSequence{$barcode})) {
		      push(@{$globalOptBarcodeSequence{$barcode}}, $globalTotalNumberSequence);
		    } else {
		      $globalOptBarcodeSequence{$barcode} = [$globalTotalNumberSequence];
		      push(@globalUniqueOptBarcode, $barcode);
		      $firstSequenceBarcode->{$barcode} = scalar(@globalUniqueOptBarcode) - 1;
		      my $regLen = length($barcode);
		      if (($globalOptBarcodeMaxLength < 0) || ($regLen > $globalOptBarcodeMaxLength)) {
			$globalOptBarcodeMaxLength = $regLen;
		      }
		      if (($globalOptBarcodeMinLength < 0) || ($regLen < $globalOptBarcodeMinLength)) {
			$globalOptBarcodeMinLength = $regLen;
		      }
		    }
		    push(@globalSequenceOptBarcode, $firstSequenceBarcode->{$barcode});
		  }
		  if (exists($globalVRegionSequence{$region[0]})) {
		    push(@{$globalVRegionSequence{$region[0]}}, $globalTotalNumberSequence);
		  } else {
		    $globalVRegionSequence{$region[0]} = [$globalTotalNumberSequence];
		    push(@globalUniqueVRegion, $region[0]);
		    $firstSequenceVRegion->{$region[0]} = scalar(@globalUniqueVRegion) - 1;
		    push(@globalVRegionName, "");
		    my $regLen = length($region[0]);
		    if (($globalVRegionMaxLength < 0) || ($regLen > $globalVRegionMaxLength)) {
		      $globalVRegionMaxLength = $regLen;
		    }
		    if (($globalVRegionMinLength < 0) || ($regLen < $globalVRegionMinLength)) {
		      $globalVRegionMinLength = $regLen;
		    }
		  }
		  push(@globalSequenceVRegion, $firstSequenceVRegion->{$region[0]});
		  if (exists($globalCDR3RegionSequence{$region[1]})) {
		    push(@{$globalCDR3RegionSequence{$region[1]}}, $globalTotalNumberSequence);
		  } else {
		    $globalCDR3RegionSequence{$region[1]} = [$globalTotalNumberSequence];
		    push(@globalUniqueCDR3Region, $region[1]);
		    $firstSequenceCDR3Region->{$region[1]} = scalar(@globalUniqueCDR3Region) - 1;
		    my $regLen = length($region[1]);
		    if (($globalCDR3RegionMaxLength < 0) || ($regLen > $globalCDR3RegionMaxLength)) {
		      $globalCDR3RegionMaxLength = $regLen;
		    }
		    if (($globalCDR3RegionMinLength < 0) || ($regLen < $globalCDR3RegionMinLength)) {
		      $globalCDR3RegionMinLength = $regLen;
		    }
		  }
		  push(@globalSequenceCDR3Region, $firstSequenceCDR3Region->{$region[1]});
		  if (exists($globalJRegionSequence{$region[2]})) {
		    push(@{$globalJRegionSequence{$region[2]}}, $globalTotalNumberSequence);
		  } else {
		    $globalJRegionSequence{$region[2]} = [$globalTotalNumberSequence];
		    push(@globalUniqueJRegion, $region[2]);
		    $firstSequenceJRegion->{$region[2]} = scalar(@globalUniqueJRegion) - 1;
		    push(@globalJRegionName, "");
		    my $regLen = length($region[2]);
		    if (($globalJRegionMaxLength < 0) || ($regLen > $globalJRegionMaxLength)) {
		      $globalJRegionMaxLength = $regLen;
		    }
		    if (($globalJRegionMinLength < 0) || ($regLen < $globalJRegionMinLength)) {
		      $globalJRegionMinLength = $regLen;
		    }
		  }
		  push(@globalSequenceJRegion, $firstSequenceJRegion->{$region[2]});
		}
		$globalTotalNumberSequence++;
		if (scalar(@region) == 6) {
		  if (($globalOptBarcodePresent != 0) && ($barcode ne "")) {
		    if (exists($globalOptScndChainOptBarcodeSequence{$barcode})) {
		      push(@{$globalOptScndChainOptBarcodeSequence{$barcode}}, $globalTotalNumberSequence);
		    } else {
		      $globalOptScndChainOptBarcodeSequence{$barcode} = [$globalTotalNumberSequence];
		      push(@globalUniqueOptBarcode, $barcode);
		      $firstSequenceBarcode->{$barcode} = scalar(@globalUniqueOptBarcode) - 1;
		      my $regLen = length($barcode);
		      if (($globalOptBarcodeMaxLength < 0) || ($regLen > $globalOptBarcodeMaxLength)) {
			$globalOptBarcodeMaxLength = $regLen;
		      }
		      if (($globalOptBarcodeMinLength < 0) || ($regLen < $globalOptBarcodeMinLength)) {
			$globalOptBarcodeMinLength = $regLen;
		      }
		    }
		    push(@globalSequenceOptBarcode, $firstSequenceBarcode->{$barcode});
		  }
		  if (exists($globalOptScndChainVRegionSequence{$region[3]})) {
		    push(@{$globalOptScndChainVRegionSequence{$region[3]}}, $globalTotalNumberSequence);
		  } else {
		    $globalOptScndChainVRegionSequence{$region[3]} = [$globalTotalNumberSequence];
		    push(@globalUniqueVRegion, $region[3]);
		    $firstSequenceVRegion->{$region[3]} = scalar(@globalUniqueVRegion) - 1;
		    push(@globalVRegionName, "");
		    my $regLen = length($region[3]);
		    if (($globalVRegionMaxLength < 0) || ($regLen > $globalVRegionMaxLength)) {
		      $globalVRegionMaxLength = $regLen;
		    }
		    if (($globalVRegionMinLength < 0) || ($regLen < $globalVRegionMinLength)) {
		      $globalVRegionMinLength = $regLen;
		    }
		  }
		  push(@globalSequenceVRegion, $firstSequenceVRegion->{$region[3]});
		  if (exists($globalOptScndChainCDR3RegionSequence{$region[4]})) {
		    push(@{$globalOptScndChainCDR3RegionSequence{$region[4]}}, $globalTotalNumberSequence);
		  } else {
		    $globalOptScndChainCDR3RegionSequence{$region[4]} = [$globalTotalNumberSequence];
		    push(@globalUniqueCDR3Region, $region[4]);
		    $firstSequenceCDR3Region->{$region[4]} = scalar(@globalUniqueCDR3Region) - 1;
		    my $regLen = length($region[4]);
		    if (($globalCDR3RegionMaxLength < 0) || ($regLen > $globalCDR3RegionMaxLength)) {
		      $globalCDR3RegionMaxLength = $regLen;
		    }
		    if (($globalCDR3RegionMinLength < 0) || ($regLen < $globalCDR3RegionMinLength)) {
		      $globalCDR3RegionMinLength = $regLen;
		    }
		  }
		  push(@globalSequenceCDR3Region, $firstSequenceCDR3Region->{$region[4]});
		  if (exists($globalOptScndChainJRegionSequence{$region[5]})) {
		    push(@{$globalOptScndChainJRegionSequence{$region[5]}}, $globalTotalNumberSequence);
		  } else {
		    $globalOptScndChainJRegionSequence{$region[5]} = [$globalTotalNumberSequence];
		    push(@globalUniqueJRegion, $region[5]);
		    $firstSequenceJRegion->{$region[5]} = scalar(@globalUniqueJRegion) - 1;
		    push(@globalJRegionName, "");
		    my $regLen = length($region[5]);
		    if (($globalJRegionMaxLength < 0) || ($regLen > $globalJRegionMaxLength)) {
		      $globalJRegionMaxLength = $regLen;
		    }
		    if (($globalJRegionMinLength < 0) || ($regLen < $globalJRegionMinLength)) {
		      $globalJRegionMinLength = $regLen;
		    }
		  }
		  push(@globalSequenceJRegion, $firstSequenceJRegion->{$region[5]});
		  $globalTotalNumberSequence++;
		}
		$foundData = 1;
	      }
	    }
	  }
	}
      }
    }
  }

  return $foundData;
}

sub parseOneLine {
  my $firstSequenceBarcode = shift;
  my $firstSequenceVRegion = shift;
  my $firstSequenceCDR3Region = shift;
  my $firstSequenceJRegion = shift;
  my $chainDelimiterConstReg = shift;
  my $wellBarcodeConstReg = shift;
  my $barcodeDataSet = shift;
  my $inputDataLine = shift;
  my $optChainNumber = shift;

  my $foundData = 0;

  my $orgInputDataLine = "";
  my $wrongBarcode1 = 0;
  if ($barcodeDataSet ne "") {
    if ($inputDataLine !~ /^$barcodeDataSet/) {
      $wrongBarcode1 = 1;
    } else {
      $orgInputDataLine = $inputDataLine;
      $inputDataLine =~ s/^$barcodeDataSet//;
    }
  }
  if ($wrongBarcode1 == 0) {
    $foundData = parseLineFrames($firstSequenceBarcode, $firstSequenceVRegion, $firstSequenceCDR3Region, $firstSequenceJRegion, 
				 $chainDelimiterConstReg, $wellBarcodeConstReg, $inputDataLine, $optChainNumber);
  }

  my $wrongBarcode2 = 0;
  if ($foundData == 0) {
    my $revCmpLine = "";
    if ($orgInputDataLine ne "") {
      $revCmpLine = reverseComplements($orgInputDataLine);
    } else {
      $revCmpLine = reverseComplements($inputDataLine);
    }
    if ($barcodeDataSet ne "") {
      if ($revCmpLine !~ /^$barcodeDataSet/) {
	$wrongBarcode2 = 1;
      } else {
	$revCmpLine =~ s/^$barcodeDataSet//;
      }
    }
    if ($wrongBarcode2 == 0) {
      $foundData = parseLineFrames($firstSequenceBarcode, $firstSequenceVRegion, $firstSequenceCDR3Region, $firstSequenceJRegion, 
				   $chainDelimiterConstReg, $wellBarcodeConstReg, $revCmpLine, $optChainNumber);
    }
  }

  if (($wrongBarcode1 != 0) && ($wrongBarcode2 != 0)) {
    return 1;
  }

  return 0;
}

sub parseDataFile {
  my $inputFile = shift;
  my $optScndFile = shift;
  my $chainDelimiterConstReg = shift;
  my $wellBarcodeConstReg = shift;
  my $barcodeDataSet = shift;

  my %firstSequenceBarcode;
  my %firstSequenceVRegion;
  my %firstSequenceCDR3Region;
  my %firstSequenceJRegion;
  my $wrongBarcode = 0;
  my $scoreNext = 0;
  my $discarded = 0;
  my $newLine = "";
  my $comment = "";
  my $number = 0;
  my $lines = "";

  my $numberFiles = 1;
  if ($optScndFile ne "") {
    $numberFiles = 2;
  } elsif ($chainDelimiterConstReg ne "") {
    $globalOptScndChainFirstSequence = 0;
  }

  for (my $numberOfFile = 0; $numberOfFile < $numberFiles; $numberOfFile++) {
    my $nameOfFile = $inputFile;
    if ($numberOfFile == 1) {
      $globalOptScndChainFirstSequence = $globalTotalNumberSequence;
      $nameOfFile = $optScndFile;
      print "\n";
    }

    print "*MESSAGE* Parsing sequence data from file $nameOfFile\n";

    my $cnt = new IO::Uncompress::AnyUncompress($nameOfFile) || die "*ERROR* Cannot open file $nameOfFile\n";
    my $totalNumberLines = 0;
    while(<$cnt>) {
      $totalNumberLines++;
    }
    $cnt->close();

    my $inp = new IO::Uncompress::AnyUncompress($nameOfFile) || die "*ERROR* Cannot open file $nameOfFile\n";
    my $lineNumber = 0;
    while (($newLine = <$inp>)) {
      chomp($newLine);
      $newLine = uc($newLine);
      if ($scoreNext != 0) {
	if (discardBadQuality($newLine) != 0) {
	  $discarded++;
	  $lines = "";
	}
	$scoreNext = 0;
      } elsif (($newLine =~ m/^>/) || ($newLine =~ m/^@/)) {
	if ($lines ne "") {
	  $wrongBarcode += parseOneLine(\%firstSequenceBarcode, \%firstSequenceVRegion, \%firstSequenceCDR3Region, \%firstSequenceJRegion, 
					$chainDelimiterConstReg, $wellBarcodeConstReg, $barcodeDataSet, $lines, $numberOfFile);
	}
	$scoreNext = 0;
	$lines = "";
	$number++;
      } elsif ($newLine =~ m/^\+/) {
	$scoreNext = 1;
      } else {
	$lines .= $newLine;
      }
      $lineNumber++;
      printProgress($lineNumber, $totalNumberLines);
    }
    if ($lines ne "") {
      $wrongBarcode += parseOneLine(\%firstSequenceBarcode, \%firstSequenceVRegion, \%firstSequenceCDR3Region, \%firstSequenceJRegion, 
				    $chainDelimiterConstReg, $wellBarcodeConstReg, $barcodeDataSet, $lines, $numberOfFile);
    }
    $inp->close();
  }

  print "\n*MESSAGE* Parsed $number sequences, ";
  if ($discarded > 0) {
    print "discarded ", $discarded, " due to bad quality (", $number - $discarded, " remain), ";
  }
  if ($wrongBarcode > 0) {
    print "discarded ", $wrongBarcode, " due to wrong barcode (", $number - $discarded - $wrongBarcode, " remain), ";
  }
  print "found $globalTotalNumberSequence valid ones\n";
  if ($globalOptBarcodePresent != 0) {
    print "*MESSAGE* Data base contains ", numberBarcode(), " barcodes (length ", $globalOptBarcodeMinLength, " to ", $globalOptBarcodeMaxLength, ")\n";
  }
  print "*MESSAGE* Data base contains ", numberVRegion(), " V-region (length ", $globalVRegionMinLength, " to ", $globalVRegionMaxLength, ")\n";
  print "*MESSAGE* Data base contains ", numberCdr3Region(), " CDR3-region (length ", $globalCDR3RegionMinLength, " to ", $globalCDR3RegionMaxLength, ")\n";
  print "*MESSAGE* Data base contains ", numberJRegion(), " J-region (length ", $globalJRegionMinLength, " to ", $globalJRegionMaxLength, ")\n\n";

  return;
}

###########################################################
# Functions to annotate V/J-regions using name data bases #
###########################################################

sub parseNameDB {
  my $namesArray = shift;
  my $sequenceArray = shift;
  my $inputFile = shift;

  my $inp = new IO::Uncompress::AnyUncompress($inputFile) || die "*ERROR* Cannot open file $inputFile\n";

  my %alreadyStored;
  my $sequence = "";
  my $newLine = "";
  my $number = 0;
  my $name = "";

  while (($newLine = <$inp>)) {
    chomp($newLine);
    $newLine = uc($newLine);
    if ($newLine =~ m/^>[ ]*([^ ]+)/) {
      my $tmp = $1;
      if (($name ne "") && ($sequence ne "") && !exists($alreadyStored{$name})) {
	$alreadyStored{$name} = 1;
	$namesArray->[$number] = $name;
	$sequenceArray->[$number] = $sequence;
	$number++;
      }
      $sequence = "";
      $name = $tmp;
    } else {
      $sequence .= $newLine;
    }
  }

  if (($name ne "") && ($sequence ne "") && !exists($alreadyStored{$name})) {
    $namesArray->[$number] = $name;
    $sequenceArray->[$number] = $sequence;
    $number++;
  }

  $inp->close();

  return $number;
}

sub annotateVRegion {
  my $inputFile = shift;

  print "*MESSAGE* Parsing V-region name data base $inputFile\n";

  my @namesArray;
  my @sequenceArray;
  my $totalNumber = parseNameDB(\@namesArray, \@sequenceArray, $inputFile);
  
  print "*MESSAGE* Parsed $totalNumber name data base entries\n";

  if ($totalNumber > 0) {
    my $numberRegion = numberVRegion();
    print "*MESSAGE* Annotating ", $numberRegion, " unique V-regions\n";
    my $annotated = 0;
    for (my $i = 0; $i < $numberRegion; $i++) {
      my $bases = $globalUniqueVRegion[$i];
      my $lengthBases = length($bases);
      my $minusLen = $lengthBases * -1;
      my $match = 0;
      for (my $j = 0; ($match == 0) && ($j < $totalNumber); $j++) {
	my $lengthArray = length($sequenceArray[$j]); 
	if ((($lengthArray >= $lengthBases) && ($bases eq substr($sequenceArray[$j], $minusLen))) ||
	    (($lengthArray < $lengthBases) && (substr($bases, $lengthArray * -1) eq $sequenceArray[$j]))) {
	  $globalVRegionName[$i] = $namesArray[$j];
	  $match = 1;
	}
      }
      # First character is highly likely wrong...
      if (($match == 0) && ($lengthBases > 5)) {
	$bases =~ s/^.//;
	$lengthBases = length($bases);
	$minusLen = $lengthBases * -1;
	for (my $j = 0; ($match == 0) && ($j < $totalNumber); $j++) {
	  my $lengthArray = length($sequenceArray[$j]); 
	  if ((($lengthArray >= $lengthBases) && ($bases eq substr($sequenceArray[$j], $minusLen))) ||
	      (($lengthArray < $lengthBases) && (substr($bases, $lengthArray * -1) eq $sequenceArray[$j]))) {
	    $globalVRegionName[$i] = $namesArray[$j];
	    $match = 1;
	  }
	}
      }
      if (($match == 0) && ($lengthBases > 5)) {
	for (my $j = 0; ($match == 0) && ($j < $totalNumber); $j++) {
	  my $lengthArray = length($sequenceArray[$j]); 
	  if ((($lengthArray >= $lengthBases) && (fuzzyMatchOne($bases, substr($sequenceArray[$j], $minusLen)) == 0)) ||
	      (($lengthArray < $lengthBases) && ($lengthArray > 5) && (fuzzyMatchOne(substr($bases, $lengthArray * -1), $sequenceArray[$j]) == 0))) {
	    $globalVRegionName[$i] = $namesArray[$j];
	    $match = 1;
	  }
	}
      }
      if (($match == 0) && ($lengthBases > 10)) {
	for (my $j = 0; ($match == 0) && ($j < $totalNumber); $j++) {
	  my $lengthArray = length($sequenceArray[$j]); 
	  if ((($lengthArray >= $lengthBases) && (fuzzyMatchTwo($bases, substr($sequenceArray[$j], $minusLen)) == 0)) ||
	      (($lengthArray < $lengthBases) && ($lengthArray > 10) && (fuzzyMatchTwo(substr($bases, $lengthArray * -1), $sequenceArray[$j]) == 0))) {
	    $globalVRegionName[$i] = $namesArray[$j];
	    $match = 1;
	  }
	}
      }
      if ($match != 0) {
	$annotated++;
      }
      printProgress($i + 1, $numberRegion);
    }
    print "\n*MESSAGE* Annotated ", $annotated, " unique V-regions\n";
  }
  print "\n";

  return;
}

sub annotateJRegion {
  my $inputFile = shift;

  print "*MESSAGE* Parsing J-region name data base $inputFile\n";

  my @namesArray;
  my @sequenceArray;
  my $totalNumber = parseNameDB(\@namesArray, \@sequenceArray, $inputFile);
  
  print "*MESSAGE* Parsed $totalNumber name data base entries\n";

  if ($totalNumber > 0) {
    my $numberRegion = numberJRegion();
    print "*MESSAGE* Annotating ", $numberRegion, " unique J-regions\n";
    my $annotated = 0;
    for (my $i = 0; $i < $numberRegion; $i++) {
      my $bases = $globalUniqueJRegion[$i];
      my $lengthBases = length($bases);
      my $match = 0;
      for (my $j = 0; ($match == 0) && ($j < $totalNumber); $j++) {
	my $lengthArray = length($sequenceArray[$j]); 
	if ((($lengthArray >= $lengthBases) && ($bases eq substr($sequenceArray[$j], 0, $lengthBases))) ||
	    (($lengthArray < $lengthBases) && (substr($bases, 0, $lengthArray) eq $sequenceArray[$j]))) {
	  $globalJRegionName[$i] = $namesArray[$j];
	  $match = 1;
	}
      }
      # Last character is highly likely wrong...
      if (($match == 0) && ($lengthBases > 5)) {
	$bases =~ s/.$//;
	$lengthBases = length($bases);
	for (my $j = 0; ($match == 0) && ($j < $totalNumber); $j++) {
	  my $lengthArray = length($sequenceArray[$j]); 
	  if ((($lengthArray >= $lengthBases) && ($bases eq substr($sequenceArray[$j], 0, $lengthBases))) ||
	      (($lengthArray < $lengthBases) && (substr($bases, 0, $lengthArray) eq $sequenceArray[$j]))) {
	    $globalJRegionName[$i] = $namesArray[$j];
	    $match = 1;
	  }
	}
      }
      if (($match == 0) && ($lengthBases > 5)) {
	for (my $j = 0; ($match == 0) && ($j < $totalNumber); $j++) {
	  my $lengthArray = length($sequenceArray[$j]); 
	  if ((($lengthArray >= $lengthBases) && (fuzzyMatchOne($bases, substr($sequenceArray[$j], 0, $lengthBases)) == 0)) ||
	      (($lengthArray < $lengthBases) && ($lengthArray > 5) && (fuzzyMatchOne(substr($bases, 0, $lengthArray), $sequenceArray[$j]) == 0))) {
	    $globalJRegionName[$i] = $namesArray[$j];
	    $match = 1;
	  }
	}
      }
      if (($match == 0) && ($lengthBases > 10)) {
	for (my $j = 0; ($match == 0) && ($j < $totalNumber); $j++) {
	  my $lengthArray = length($sequenceArray[$j]); 
	  if ((($lengthArray >= $lengthBases) && (fuzzyMatchTwo($bases, substr($sequenceArray[$j], 0, $lengthBases)) == 0)) ||
	      (($lengthArray < $lengthBases) && ($lengthArray > 10) && (fuzzyMatchTwo(substr($bases, 0, $lengthArray), $sequenceArray[$j]) == 0))) {
	    $globalJRegionName[$i] = $namesArray[$j];
	    $match = 1;
	  }
	}
      }
      if ($match != 0) {
	$annotated++;
      }
      printProgress($i + 1, $numberRegion);
    }
    print "\n*MESSAGE* Annotated ", $annotated, " unique J-regions\n";
  }
  print "\n";

  return;
}

###################################
# Functions to export data tables #
###################################

sub exportBarcode {
  my $dbh = shift;

  if ($globalOptBarcodePresent != 0) {
    my $sth = $dbh->prepare("DROP TABLE IF EXISTS barcode");
    $sth->execute();

    my $query = "CREATE TABLE barcode (id INTEGER PRIMARY KEY";
    if ($globalOptScndChainFirstSequence >= 0) {
      $query .= ", chain VARCHAR";
    }
    $query .= ", bases VARCHAR, number INTEGER)";
    $sth = $dbh->prepare($query);
    $sth->execute();

    $query = "";
    my $number = 0;
    my $numberRegion = numberBarcode();
    print "*MESSAGE* Exporting ", $numberRegion, " unique barcodes\n";
    foreach my $region (keys(%globalOptBarcodeSequence)) {
      if ($query eq "") {
	$query = "INSERT INTO barcode (";
	if ($globalOptScndChainFirstSequence >= 0) {
	  $query .= "chain, ";
	}
	$query .= "bases, number) VALUES ";
      } else {
	$query .= ",";
      }
      $query .= "(";
      if ($globalOptScndChainFirstSequence >= 0) {
	$query .= "\"ALPHA\",";
      }
      $query .= "\"" . $region . "\"," . barcodeNumber($region) . ")";
      if (($number % $globalMaxTermsInsertTable) == 0) {
	$sth = $dbh->prepare($query);
	$sth->execute();
	$query = "";
      }
      printProgress(++$number, $numberRegion);
    }
    if ($query ne "") {
      $sth = $dbh->prepare($query);
      $sth->execute();
      $query = "";
    }
    if ($globalOptScndChainFirstSequence >= 0) {
      foreach my $region (keys(%globalOptScndChainOptBarcodeSequence)) {
	if ($query eq "") {
	  $query = "INSERT INTO barcode (chain, bases, number) VALUES ";
	} else {
	  $query .= ",";
	}
	$query .= "(\"BETA\",\"" . $region . "\"," . barcodeNumber($region, 1) . ")";
	if (($number % $globalMaxTermsInsertTable) == 0) {
	  $sth = $dbh->prepare($query);
	  $sth->execute();
	  $query = "";
	}
	printProgress(++$number, $numberRegion);
      }
      if ($query ne "") {
	$sth = $dbh->prepare($query);
	$sth->execute();
      }
    }
    print "\n\n";
  }

  return;
}

sub exportVRegion {
  my $dbh = shift;

  my $sth = $dbh->prepare("DROP TABLE IF EXISTS vregion");
  $sth->execute();

  my $query = "CREATE TABLE vregion (id INTEGER PRIMARY KEY";
  if ($globalOptScndChainFirstSequence >= 0) {
    $query .= ", chain VARCHAR";
  }
  $query .= ", bases VARCHAR, name VARCHAR, number INTEGER)";
  $sth = $dbh->prepare($query);
  $sth->execute();

  $query = "";
  my $number = 0;
  my $numberRegion = numberVRegion();
  print "*MESSAGE* Exporting ", $numberRegion, " unique V-regions\n";
  foreach my $region (keys(%globalVRegionSequence)) {
    my $name = vRegionName($region);
    if ($name eq "") {
      $name = "undef";
    }
    if ($query eq "") {
      $query = "INSERT INTO vregion (";
      if ($globalOptScndChainFirstSequence >= 0) {
	$query .= "chain, ";
      }
      $query .= "bases, name, number) VALUES ";
    } else {
      $query .= ",";
    }
    $query .= "(";
    if ($globalOptScndChainFirstSequence >= 0) {
      $query .= "\"ALPHA\",";
    }
    $query .= "\"" . $region . "\",\"" . $name . "\"," . vRegionNumber($region) . ")";
    if (($number % $globalMaxTermsInsertTable) == 0) {
      $sth = $dbh->prepare($query);
      $sth->execute();
      $query = "";
    }
    printProgress(++$number, $numberRegion);
  }
  if ($query ne "") {
    $sth = $dbh->prepare($query);
    $sth->execute();
    $query = "";
  }
  if ($globalOptScndChainFirstSequence >= 0) {
    foreach my $region (keys(%globalOptScndChainVRegionSequence)) {
      my $name = vRegionName($region, 1);
      if ($name eq "") {
	$name = "undef";
      }
      if ($query eq "") {
	$query = "INSERT INTO vregion (chain, bases, name, number) VALUES ";
      } else {
	$query .= ",";
      }
      $query .= "(\"BETA\",\"" . $region . "\",\"" . $name . "\"," . vRegionNumber($region, 1) . ")";
      if (($number % $globalMaxTermsInsertTable) == 0) {
	$sth = $dbh->prepare($query);
	$sth->execute();
	$query = "";
      }
      printProgress(++$number, $numberRegion);
    }
    if ($query ne "") {
      $sth = $dbh->prepare($query);
      $sth->execute();
    }
  }
  print "\n\n";

  return;
}

sub exportCDR3Region {
  my $dbh = shift;

  my $sth = $dbh->prepare("DROP TABLE IF EXISTS cdr3region");
  $sth->execute();

  my $query = "";
  $query = "CREATE TABLE cdr3region (id INTEGER PRIMARY KEY";
  if ($globalOptScndChainFirstSequence >= 0) {
    $query .= ", chain VARCHAR";
  }
  $query .= ", bases VARCHAR, number INTEGER, vregion VARCHAR, jregion VARCHAR)";
  $sth = $dbh->prepare($query);
  $sth->execute();

  $query = "";
  my $number = 0;
  my $numberRegion = numberCdr3Region();
  print "*MESSAGE* Exporting ", $numberRegion, " unique CDR3-regions\n";
  foreach my $region (keys(%globalCDR3RegionSequence)) {
    if ($query eq "") {
      $query = "INSERT INTO cdr3region (";
      if ($globalOptScndChainFirstSequence >= 0) {
	$query .= "chain, ";
      }
      $query .= "bases, number, vregion, jregion) VALUES ";
    } else {
      $query .= ",";
    }
    $query .= "(";
    if ($globalOptScndChainFirstSequence >= 0) {
      $query .= "\"ALPHA\",";
    }
    $query .= "\"" . $region . "\"," . cdr3RegionNumber($region) . ",\"" . vRegionOfCdr3($region) . "\",\"" . jRegionOfCdr3($region) . "\")";
    if (($number % $globalMaxTermsInsertTable) == 0) {
      $sth = $dbh->prepare($query);
      $sth->execute();
      $query = "";
    }
    printProgress(++$number, $numberRegion);
  }
  if ($query ne "") {
    $sth = $dbh->prepare($query);
    $sth->execute();
    $query = "";
  }
  if ($globalOptScndChainFirstSequence >= 0) {
    foreach my $region (keys(%globalOptScndChainCDR3RegionSequence)) {
      if ($query eq "") {
	$query = "INSERT INTO cdr3region (chain, bases, number, vregion, jregion) VALUES ";
      } else {
	$query .= ",";
      }
      $query .= "(\"BETA\",\"" . $region . "\"," . cdr3RegionNumber($region, 1) . ",\"" . vRegionOfCdr3($region, 1) . "\",\"" . jRegionOfCdr3($region, 1) . "\")";
      if (($number % $globalMaxTermsInsertTable) == 0) {
	$sth = $dbh->prepare($query);
	$sth->execute();
	$query = "";
      }
      printProgress(++$number, $numberRegion);
    }
    if ($query ne "") {
      $sth = $dbh->prepare($query);
      $sth->execute();
    }
  }
  print "\n\n";

  return;
}

sub exportJRegion {
  my $dbh = shift;

  my $sth = $dbh->prepare("DROP TABLE IF EXISTS jregion");
  $sth->execute();

  my $query = "";
  $query = "CREATE TABLE jregion (id INTEGER PRIMARY KEY";
  if ($globalOptScndChainFirstSequence >= 0) {
    $query .= ", chain VARCHAR";
  }
  $query .= ", bases VARCHAR, name VARCHAR, number INTEGER)";
  $sth = $dbh->prepare($query);
  $sth->execute();

  $query = "";
  my $number = 0;
  my $numberRegion = numberJRegion();
  print "*MESSAGE* Exporting ", $numberRegion, " unique J-regions\n";
  foreach my $region (keys(%globalJRegionSequence)) {
    my $name = jRegionName($region);
    if ($name eq "") {
      $name = "undef";
    }
    if ($query eq "") {
      $query = "INSERT INTO jregion (";
      if ($globalOptScndChainFirstSequence >= 0) {
	$query .= "chain, ";
      }
      $query .= "bases, name, number) VALUES ";
      } else {
      $query .= ",";
    }
    $query .= "(";
    if ($globalOptScndChainFirstSequence >= 0) {
      $query .= "\"ALPHA\",";
    }
    $query .= "\"" . $region . "\",\"" . $name . "\"," . jRegionNumber($region) . ")";
    if (($number % $globalMaxTermsInsertTable) == 0) {
      $sth = $dbh->prepare($query);
      $sth->execute();
      $query = "";
    }
    printProgress(++$number, $numberRegion);
  }
  if ($query ne "") {
    $sth = $dbh->prepare($query);
    $sth->execute();
    $query = "";
  }
  if ($globalOptScndChainFirstSequence >= 0) {
    foreach my $region (keys(%globalOptScndChainJRegionSequence)) {
      my $name = jRegionName($region, 1);
      if ($name eq "") {
	$name = "undef";
      }
      if ($query eq "") {
	$query = "INSERT INTO jregion (chain, bases, name, number) VALUES ";
      } else {
	$query .= ",";
      }
      $query .= "(\"BETA\",\"" . $region . "\",\"" . $name . "\"," . jRegionNumber($region, 1) . ")";
      if (($number % $globalMaxTermsInsertTable) == 0) {
	$sth = $dbh->prepare($query);
	$sth->execute();
	$query = "";
      }
      printProgress(++$number, $numberRegion);
    }
    if ($query ne "") {
      $sth = $dbh->prepare($query);
      $sth->execute();
    }
  }
  print "\n\n";

  return;
}

sub exportConnectivity {
  my $dbh = shift;

  my $sth = $dbh->prepare("DROP TABLE IF EXISTS connectivity");
  $sth->execute();

  my $query = "CREATE TABLE connectivity (id INTEGER PRIMARY KEY";
  if ($globalOptScndChainFirstSequence  > 0) {
    $query .= ", chain VARCHAR";
  }
  if ($globalOptBarcodePresent != 0) {
    $query .= ", barcode VARCHAR";
  }
  if ($globalOptScndChainFirstSequence  == 0) {
    $query .= ", \"alpha vregion bases\" VARCHAR, \"alpha vregion name\" VARCHAR, \"alpha cdr3region\" VARCHAR, \"alpha jregion bases\" VARCHAR, \"alpha jregion name\" VARCHAR, ";
    $query .= "\"beta vregion bases\" VARCHAR, \"beta vregion name\" VARCHAR, \"beta cdr3region\" VARCHAR, \"beta jregion bases\" VARCHAR, \"beta jregion name\" VARCHAR, number INTEGER)";
  } else {
    $query .= ", \"vregion bases\" VARCHAR, \"vregion name\" VARCHAR, \"cdr3region\" VARCHAR, \"jregion bases\" VARCHAR, \"jregion name\" VARCHAR, number INTEGER)";
  }
  $sth = $dbh->prepare($query);
  $sth->execute();

  if (!scalar(%globalRegionConnectivity) || ($globalTotalNumberConnect == 0)) {
    buildConnectivity();
  }

  $query = "";
  my $number = 0;
  my $numberConnect = scalar(keys(%globalRegionConnectivity));
  print "*MESSAGE* Exporting ", $numberConnect, " unique combinations\n";
  foreach my $combination (keys(%globalRegionConnectivity)) {
    my $fullCombination = $combination;
    if ($query eq "") {
      $query = "INSERT INTO connectivity (";
      if ($globalOptScndChainFirstSequence > 0) {
	$query .= "chain, ";
      }
      if ($globalOptBarcodePresent != 0) {
	$query .= "barcode, ";
      }
      if ($globalOptScndChainFirstSequence == 0) {
	$query .= "\"alpha vregion bases\", \"alpha vregion name\", \"alpha cdr3region\", \"alpha jregion bases\", \"alpha jregion name\", ";
	$query .= "\"beta vregion bases\", \"beta vregion name\", \"beta cdr3region\", \"beta jregion bases\", \"beta jregion name\", number) VALUES ";
      } else {
	$query .= "\"vregion bases\", \"vregion name\", \"cdr3region\", \"jregion bases\", \"jregion name\", number) VALUES ";
      }
    } else {
      $query .= ",";
    }
    $query .= "(";
    my $optScndChain = "";
    if ($globalOptScndChainFirstSequence == 0) {
      $optScndChain = $combination;
      $optScndChain =~ s/^.+ BETA //;
      $combination =~ s/ BETA.+$//;
    }
    my $optChainIdentifier = "";
    if ($globalOptScndChainFirstSequence > 0) {
      $combination =~ s/^([^ ]*) //;
      $optChainIdentifier = $1;
      $query .= "\"" . $optChainIdentifier . "\",";
    }
    if ($globalOptBarcodePresent != 0) {
      $combination =~ s/^([^ ]*) //;
      $query .= "\"" . $1 . "\",";
    }
    if ($globalOptScndChainFirstSequence == 0) {
      $combination =~ s/^([^ ]*) //;
      $optChainIdentifier = $1;
    }
    my ($vBases, $cdr3Bases, $jBases) = split(' ', $combination);
    $query .= "\"" . $vBases . "\",";
    my $vName = "";
    if ($optChainIdentifier eq "BETA") {
      $vName = vRegionName($vBases, 1);
    } else {
      $vName = vRegionName($vBases);
    }
    if ($vName eq "") {
      $vName = "undef";
    }
    $query .= "\"" . $vName . "\",\"" . $cdr3Bases . "\",\"" . $jBases . "\",";
    my $jName = "";
    if ($optChainIdentifier eq "BETA") {
      $jName = jRegionName($jBases, 1);
    } else {
      $jName = jRegionName($jBases);
    }
    if ($jName eq "") {
      $jName = "undef";
    }
    $query .= "\"" . $jName . "\",";
    if ($optScndChain ne "") {
      my ($optScndChainVBases, $optScndChainCdr3Bases, $optScndChainJBases) = split(' ', $optScndChain);
      $query .= "\"" . $optScndChainVBases . "\",";
      my $optScndChainVName = vRegionName($optScndChainVBases, 1);
      if ($optScndChainVName eq "") {
	$optScndChainVName = "undef";
      }
      $query .= "\"" . $optScndChainVName . "\",\"" . $optScndChainCdr3Bases . "\",\"" . $optScndChainJBases . "\",";
      my $optScndChainJName = jRegionName($optScndChainJBases, 1);
      if ($optScndChainJName eq "") {
	$optScndChainJName = "undef";
      }
      $query .= "\"" . $optScndChainJName . "\",";
    }
    $query .= $globalRegionConnectivity{$fullCombination} . ")";
    if (($number % $globalMaxTermsInsertTable) == 0) {
      $sth = $dbh->prepare($query);
      $sth->execute();
      $query = "";
    }
    printProgress(++$number, $numberConnect);
  }
  if ($query ne "") {
    $sth = $dbh->prepare($query);
    $sth->execute();
  }
  print "\n\n";

  return;
}

sub exportDataTable {
  my $dbfile = shift;
  my $connectitvity = shift;

  # Type of data base is defined here
  my $dsn = "dbi:SQLite:dbname=$dbfile";
  my $user = "";
  my $password = "";
  my $dbh = DBI->connect($dsn, $user, $password, {
    PrintError       => 0,
    RaiseError       => 1,
    AutoCommit       => 1,
    FetchHashKeyName => 'NAME_lc',
  });

  if ($globalOptBarcodePresent != 0) {
    exportBarcode($dbh);
  }

  exportVRegion($dbh);
  exportCDR3Region($dbh);
  exportJRegion($dbh);

  if ($connectitvity != 0) {
    exportConnectivity($dbh);
  }

  $dbh->disconnect();

  return;
}

################################
# Function to process data set #
################################

sub processDataSet {
  my $pathScript = shift;
  my $inputFile = shift;
  my $optPairedInputFile = shift;
  my $optScndChainFile = shift;
  my $optPairedScndChainFile = shift;
  my $chainDelimiterConstReg = shift;
  my $wellBarcodeConstReg = shift;
  my $barcodeDataSet = shift;
  my $noConnectivity = shift;
  my $mergeTool = shift;
  my $numThread = shift;
  my $vNameDB = shift;
  my $jNameDB = shift;
  my $nameDB = shift;
  my $verbose = shift;

  # Reinit global data for successive runs
  initGlobalDataTables();

  # Check whether constant region before barcode is given and barcodes are present
  if ($wellBarcodeConstReg ne "") {
    $globalOptBarcodePresent = 1;
  }

  # Check conflicting two files and directly linked paired chains approaches
  if (($chainDelimiterConstReg ne "") && ($optScndChainFile ne "")) {
    die "*ERROR* Paired chains cannot be in two files and directly linked at the same time\n";
  }

  # Check existence of name data bases
  if ((($vNameDB ne "") || ($jNameDB ne "")) && ($nameDB ne "")) {
    die "*ERROR* Please specfiy separate V- and J-region name data bases or one name data base for both\n";
  }
  if (($vNameDB ne "") && (! -s $vNameDB)) {
    die "*ERROR* V-region name data base $vNameDB does not exist or is not readable\n";
  }
  if (($jNameDB ne "") && (! -s $jNameDB)) {
    die "*ERROR* J-region name data base $jNameDB does not exist or is not readable\n";
  }
  if (($nameDB ne "") && (! -s $nameDB)) {
    die "*ERROR* Name data base $nameDB does not exist or is not readable\n";
  }

  # Check existence of sequence data file
  if ($inputFile eq "") {
    die "*ERROR* Please specify input file\n";
  }
  if (! -s $inputFile) {
    die "*ERROR* Input file $inputFile does not exist or is not readable\n";
  }

  # Check whether more sequence data files are present and merge if needed
  if ($optPairedInputFile ne "") {
    if (! -s $optPairedInputFile) {
      die "*ERROR* Input file $optPairedInputFile does not exist or is not readable\n";
    }
    $inputFile = mergePairedEnd($pathScript, $inputFile, $optPairedInputFile, $mergeTool, $numThread);
  }
  if ($optScndChainFile ne "") {
    if (! -s $optScndChainFile) {
      die "*ERROR* Input file $optScndChainFile does not exist or is not readable\n";
    }
    if ($optPairedScndChainFile ne "") {
      if (! -s $optPairedScndChainFile) {
	die "*ERROR* Input file $optPairedScndChainFile does not exist or is not readable\n";
      }
      $optScndChainFile = mergePairedEnd($pathScript, $optScndChainFile, $optPairedScndChainFile, $mergeTool, $numThread);
    }
  }

  # Parse sequence data file
  parseDataFile($inputFile, $optScndChainFile, $chainDelimiterConstReg, $wellBarcodeConstReg, $barcodeDataSet);

  # Check if default name database exists
  if ($nameDB eq "") {
    $nameDB = localFileWithPath($pathScript, "nameDataBaseVJC.txt");
    if (! -s $nameDB) {
      $nameDB = "";
    }
  }

  # Annotate V-region sequences
  if ($vNameDB ne "") {
    annotateVRegion($vNameDB);
  } elsif ($nameDB ne "") {
    annotateVRegion($nameDB);
  }

  # Annotate J-region sequences
  if ($jNameDB ne "") {
    annotateJRegion($jNameDB);
  } elsif ($nameDB ne "") {
    annotateJRegion($nameDB);
  }

  # Export data tables
  my $dbfile = $inputFile . ".sqlite";
  if ($optScndChainFile ne "") {
    $dbfile = $inputFile . "_" . bareFileName($optScndChainFile) . ".sqlite";
  }
  my $connectivity = 1;
  if ($noConnectivity != 0) { 
    $connectivity = 0;
  } else {
    buildConnectivity();
  }
  exportDataTable($dbfile, $connectivity);

  # Print out data parsed if needed
  if ($verbose != 0) {
    if ($globalOptBarcodePresent != 0) {
      printBarcodeSequences();
    }
    printVRegionSequences();
    printCDR3RegionSequences();
    printJRegionSequences();
    if ($connectivity != 0) {
      printConnectivity();
    }
    printDataTable();
  }

  return;
}

#################
# Main function #
#################

# Local variables
my $inputFile = "";
my $optPairedInputFile = "";
my $optScndChainFile = "";
my $optPairedScndChainFile = "";
my $chainDelimiterConstReg = "";
my $wellBarcodeConstReg = "";
my $barcodeDataSet = "";
my $noConnectivity = 0;
my $mergeTool = "";
my $numThread = 0;
my $vNameDB = "";
my $jNameDB = "";
my $nameDB = "";
my $verbose = 0;
my $help = 0;
my $option;

# Parse options
Getopt::Long::Configure("pass_through");
$option = GetOptions('-help|h' => \$help);
$option = GetOptions('-verbose|v' => \$verbose);
$option = GetOptions('-chainDelimiterConstReg|cdcr=s' => \$chainDelimiterConstReg);
$option = GetOptions('-wellBarcodeConstReg|wbcr=s' => \$wellBarcodeConstReg);
$option = GetOptions('-barcodeDataSet|bcds=s' => \$barcodeDataSet);
$option = GetOptions('-noConnectivity|nc' => \$noConnectivity);
$option = GetOptions('-mergeTool|mt=s' => \$mergeTool);
$option = GetOptions('-numberThread|nt=i' => \$numThread);
$option = GetOptions('-vRegNameDB|vrndb=s' => \$vNameDB);
$option = GetOptions('-jRegNameDB|jrndb=s' => \$jNameDB);
$option = GetOptions('-nameDB|ndb=s' => \$nameDB);

# Print help if needed
if ($help != 0) {
  print "Usage: aperimTCRKit.pl [options] [inputFile(s)]\n\n";
  print "Files: - One file with TCR-data of one chain (with\n";
  print "         or without constant region of barcodes given)\n\n";
  print "       - Two files without constant region of barcodes\n";
  print "         given are treated as paired end files and merged\n\n";
  print "       - Two files with constant region of barcodes given\n";
  print "         are treated as alpha- and beta-chain data files\n\n";
  print "       - Four files with constant region of barcodes given\n";
  print "         are treated as paired end files and merged and are\n";
  print "         used as alpha- and beta-chain data files then\n\n";
  print "Options: -mergeTool|mt <mergeToolPearExecutable>\n";
  print "         (default is \"pear\" in installation directory)\n";
  print "         -numberThread|nt <numberThreadsForMerging>\n\n";
  print "         -vRegNameDB|vrndb <vRegionNameDataBaseFile>\n";
  print "         -jRegNameDB|jrndb <jRegionNameDataBaseFile>\n";
  print "         or -nameDB|ndb <nameDataBase> (default is file\n";
  print "         \"nameDataBaseVJC.txt\" in installation directory)\n\n";
  print "         -chainDelimiterConstReg|bccr <constantAminoRegion>\n";
  print "         (constant region that delimits optional connected\n";
  print "         alpha/beta-chains: V-CDR3-J-constReg-V-CDR3-J)\n\n";
  print "         -wellBarcodeConstReg|bccr <constantAminoRegion>\n";
  print "         (constant region that delimits optional barcode\n";
  print "         of well and TCR-data: barcode-constReg-V-CDR3-J)\n\n";
  print "         -barcodeDataSet|bcds <nucleicAcidBarcode>\n";
  print "         (optional constant region at beginning of lines\n";
  print "         to identify data set for multiplexed sequencing)\n\n";
  print "         -noConnectivity|nc (skip export of all individual\n";
  print "         combinations of sequences of different regions)\n\n";
  print "         -verbose|v (print TCR-data to screen too)\n\n";
  print "         -help|h (this message)\n\n";
  exit(0);
}

# Map inputs to variables or raise GUI
if (scalar(@ARGV) > 0) {
  $inputFile = $ARGV[0];
  if (($wellBarcodeConstReg ne "") && (scalar(@ARGV) > 3)) {
    $optPairedInputFile = $ARGV[1];
    $optScndChainFile = $ARGV[2];
    $optPairedScndChainFile = $ARGV[2];
  } elsif (($wellBarcodeConstReg ne "") && (scalar(@ARGV) > 1)) {
    $optScndChainFile = $ARGV[1];
  } elsif (scalar(@ARGV) > 1) {
    $optPairedInputFile = $ARGV[1];
  }
  processDataSet(abs_path($0), $inputFile, $optPairedInputFile, $optScndChainFile, $optPairedScndChainFile,
		 $chainDelimiterConstReg, $wellBarcodeConstReg, $barcodeDataSet, $noConnectivity, 
		 $mergeTool, $numThread, $vNameDB, $jNameDB, $nameDB, $verbose);
} else {
  my $mainWindow = MainWindow->new();
  $mainWindow->title("AptaIT APERIM TCR-Kit");
  $mainWindow->protocol("WM_DELETE_WINDOW" => sub { $mainWindow->destroy(); });
  my $inputFileLabel = $mainWindow->Label(-text => "Alpha- or beta-chain\nTCR-data input file");
  my $inputFileEntry = $mainWindow->Entry(-width => 40, -relief => "sunken", -bd => 2, -textvariable => \$inputFile);
  my $inputFileButton = $mainWindow->Button(-text => "Browse", -command => sub {
					      $inputFile = $mainWindow->getOpenFile();
					    });
  my $optPairedInputFileLabel = $mainWindow->Label(-text => "Corresponding\npaired end file");
  my $optPairedInputFileEntry = $mainWindow->Entry(-width => 40, -relief => "sunken", -bd => 2, -textvariable => \$optPairedInputFile);
  my $optPairedInputFileButton = $mainWindow->Button(-text => "Browse", -command => sub {
						       $optPairedInputFile = $mainWindow->getOpenFile();
						     });
  my $optScndChainFileLabel = $mainWindow->Label(-text => "Beta-chain\nTCR-data input file");
  my $optScndChainFileEntry = $mainWindow->Entry(-width => 40, -relief => "sunken", -bd => 2, -textvariable => \$optScndChainFile);
  my $optScndChainFileButton = $mainWindow->Button(-text => "Browse", -command => sub {
						     $optScndChainFile = $mainWindow->getOpenFile();
						   });
  my $optPairedScndChainFileLabel = $mainWindow->Label(-text => "Corresponding\npaired end file");
  my $optPairedScndChainFileEntry = $mainWindow->Entry(-width => 40, -relief => "sunken", -bd => 2, -textvariable => \$optPairedScndChainFile);
  my $optPairedScndChainFileButton = $mainWindow->Button(-text => "Browse", -command => sub {
							   $optPairedScndChainFile = $mainWindow->getOpenFile();
							 });
  my $wellBarcodeConstRegLabel = $mainWindow->Label(-text => "Constant region between\nwell barcode and TCR-data");
  my $wellBarcodeConstRegEntry = $mainWindow->Entry(-width => 40, -relief => "sunken", -bd => 2, -textvariable => \$wellBarcodeConstReg);
  my $inlineConstRegLabel = $mainWindow->Label(-text => "Constant region\nbetween linked chains");
  my $inlineConstRegEntry = $mainWindow->Entry(-width => 40, -relief => "sunken", -bd => 2, -textvariable => \$chainDelimiterConstReg);
  my $logFileOut = $mainWindow->Scrolled("Text", -wrap => "none");
  my $scndChainInlineLastState = 0;
  my $scndChainInline = 0;
  my $scndChainFile = 0;
  my $pairedEnd = 0;
  sub adaptMainWindow {
    if (($scndChainFile != 0) && ($scndChainInline != 0)) {
      if ($scndChainInlineLastState != 0) {
	$scndChainInline = 0;
      } else {
	$scndChainFile = 0;
      }
    }
    if (($pairedEnd != 0) && ($^O !~ m/linux/i)) {
      $logFileOut->delete("0.0", "end");
      $logFileOut->insert("end", "Merge tool \"pear\" is present for linux only,\nso please merge input files manually...\n");
      $pairedEnd = 0;
    }
    if (($pairedEnd != 0) && ($scndChainFile != 0) && ($scndChainInline == 0)) {
      $inputFileLabel->configure(-text => "Alpha-chain\nTCR-data input file");
      $optPairedInputFileLabel->grid(-row => 1, -column => 0);
      $optPairedInputFileEntry->grid(-row => 1, -column => 1, -columnspan => 2);
      $optPairedInputFileButton->grid(-row => 1, -column => 3);
      $optScndChainFileLabel->grid(-row => 2, -column => 0);
      $optScndChainFileEntry->grid(-row => 2, -column => 1, -columnspan => 2);
      $optScndChainFileButton->grid(-row => 2, -column => 3);
      $optPairedScndChainFileLabel->grid(-row => 3, -column => 0);
      $optPairedScndChainFileEntry->grid(-row => 3, -column => 1, -columnspan => 2);
      $optPairedScndChainFileButton->grid(-row => 3, -column => 3);
      $wellBarcodeConstRegLabel->grid(-row => 4, -column => 0, -columnspan => 1);
      $wellBarcodeConstRegEntry->grid(-row => 4, -column => 1, -columnspan => 2);
      $inlineConstRegLabel->gridForget();
      $inlineConstRegEntry->gridForget();
      $scndChainInlineLastState = 0;
      $chainDelimiterConstReg = "";
    } elsif (($pairedEnd != 0) && ($scndChainFile == 0) && ($scndChainInline != 0)) {
      $inputFileLabel->configure(-text => "Linked-chains\nTCR-data input file");
      $optPairedInputFileLabel->grid(-row => 1, -column => 0);
      $optPairedInputFileEntry->grid(-row => 1, -column => 1, -columnspan => 2);
      $optPairedInputFileButton->grid(-row => 1, -column => 3);
      $optScndChainFileLabel->gridForget();
      $optScndChainFileEntry->gridForget();
      $optScndChainFileButton->gridForget();
      $optPairedScndChainFileLabel->gridForget();
      $optPairedScndChainFileEntry->gridForget();
      $optPairedScndChainFileButton->gridForget();
      $wellBarcodeConstRegLabel->gridForget();
      $wellBarcodeConstRegEntry->gridForget();
      $inlineConstRegLabel->grid(-row => 4, -column => 0, -columnspan => 1);
      $inlineConstRegEntry->grid(-row => 4, -column => 1, -columnspan => 2);
      $scndChainInlineLastState = 1;
      $optScndChainFile = "";
      $optPairedScndChainFile = "";
      $wellBarcodeConstReg = "";
    } elsif (($pairedEnd != 0) && ($scndChainFile == 0) && ($scndChainInline == 0)) {
      $inputFileLabel->configure(-text => "Alpha- or beta-chain\nTCR-data input file");
      $optPairedInputFileLabel->grid(-row => 1, -column => 0);
      $optPairedInputFileEntry->grid(-row => 1, -column => 1, -columnspan => 2);
      $optPairedInputFileButton->grid(-row => 1, -column => 3);
      $optScndChainFileLabel->gridForget();
      $optScndChainFileEntry->gridForget();
      $optScndChainFileButton->gridForget();
      $optPairedScndChainFileLabel->gridForget();
      $optPairedScndChainFileEntry->gridForget();
      $optPairedScndChainFileButton->gridForget();
      $wellBarcodeConstRegLabel->gridForget();
      $wellBarcodeConstRegEntry->gridForget();
      $inlineConstRegLabel->gridForget();
      $inlineConstRegEntry->gridForget();
      $scndChainInlineLastState = 0;
      $optScndChainFile = "";
      $optPairedScndChainFile = "";
      $wellBarcodeConstReg = "";
      $chainDelimiterConstReg = "";
    } elsif (($pairedEnd == 0) && ($scndChainFile != 0) && ($scndChainInline == 0)) {
      $inputFileLabel->configure(-text => "Alpha-chain\nTCR-data input file");
      $optPairedInputFileLabel->gridForget();
      $optPairedInputFileEntry->gridForget();
      $optPairedInputFileButton->gridForget();
      $optScndChainFileLabel->grid(-row => 1, -column => 0);
      $optScndChainFileEntry->grid(-row => 1, -column => 1, -columnspan => 2);
      $optScndChainFileButton->grid(-row => 1, -column => 3);
      $optPairedScndChainFileLabel->gridForget();
      $optPairedScndChainFileEntry->gridForget();
      $optPairedScndChainFileButton->gridForget();
      $wellBarcodeConstRegLabel->grid(-row => 2, -column => 0, -columnspan => 1);
      $wellBarcodeConstRegEntry->grid(-row => 2, -column => 1, -columnspan => 2);
      $inlineConstRegLabel->gridForget();
      $inlineConstRegEntry->gridForget();
      $scndChainInlineLastState = 0;
      $optPairedInputFile = "";
      $optPairedScndChainFile = "";
      $chainDelimiterConstReg = "";
    } elsif (($pairedEnd == 0) && ($scndChainFile == 0) && ($scndChainInline != 0)) {
      $inputFileLabel->configure(-text => "Linked-chains\nTCR-data input file");
      $optPairedInputFileLabel->gridForget();
      $optPairedInputFileEntry->gridForget();
      $optPairedInputFileButton->gridForget();
      $optScndChainFileLabel->gridForget();
      $optScndChainFileEntry->gridForget();
      $optScndChainFileButton->gridForget();
      $optPairedScndChainFileLabel->gridForget();
      $optPairedScndChainFileEntry->gridForget();
      $optPairedScndChainFileButton->gridForget();
      $wellBarcodeConstRegLabel->gridForget();
      $wellBarcodeConstRegEntry->gridForget();
      $inlineConstRegLabel->grid(-row => 3, -column => 0, -columnspan => 1);
      $inlineConstRegEntry->grid(-row => 3, -column => 1, -columnspan => 2);
      $scndChainInlineLastState = 1;
      $optPairedInputFile = "";
      $optScndChainFile = "";
      $optPairedScndChainFile = "";
      $wellBarcodeConstReg = "";
    } else {
      $inputFileLabel->configure(-text => "Alpha- or beta-chain\nTCR-data input file");
      $optPairedInputFileLabel->gridForget();
      $optPairedInputFileEntry->gridForget();
      $optPairedInputFileButton->gridForget();
      $optScndChainFileLabel->gridForget();
      $optScndChainFileEntry->gridForget();
      $optScndChainFileButton->gridForget();
      $optPairedScndChainFileLabel->gridForget();
      $optPairedScndChainFileEntry->gridForget();
      $optPairedScndChainFileButton->gridForget();
      $wellBarcodeConstRegLabel->gridForget();
      $wellBarcodeConstRegEntry->gridForget();
      $inlineConstRegLabel->gridForget();
      $inlineConstRegEntry->gridForget();
      $scndChainInlineLastState = 0;
      $optPairedInputFile = "";
      $optScndChainFile = "";
      $optPairedScndChainFile = "";
      $wellBarcodeConstReg = "";
      $chainDelimiterConstReg = "";
      $scndChainFile = 0;
      $scndChainInline = 0;
      $pairedEnd = 0;
    }
    return;
  }
  my $scndChainFileButton = $mainWindow->Checkbutton(-text => "Alpha/beta-chains\nin two files (m:n)", -variable => \$scndChainFile, -onvalue => 1, -offvalue => 0, -command => [\&adaptMainWindow]);
  my $scndChainInlineButton = $mainWindow->Checkbutton(-text => "Alpha/beta-chains\nlinked inline (1:1)", -variable => \$scndChainInline, -onvalue => 1, -offvalue => 0, -command => [\&adaptMainWindow]);
  my $pairedEndButton = $mainWindow->Checkbutton(-text => "Paired\nend run", -variable => \$pairedEnd, -onvalue => 1, -offvalue => 0, -command => [\&adaptMainWindow]);
  my $startButton = $mainWindow->Button(-text => "  Start   ", -command => sub {  
					  $logFileOut->delete("0.0", "end");
					  if (($inputFile ne "") && (-s $inputFile) &&
					      (($pairedEnd == 0) || (($optPairedInputFile ne "") && (-s $optPairedInputFile))) &&
					      (($scndChainFile == 0) || (($optScndChainFile ne "") && (-s $optScndChainFile))) &&
					      ((($pairedEnd == 0) || ($scndChainFile == 0)) || (($optPairedScndChainFile ne "") && (-s $optPairedScndChainFile))) &&
					      (($scndChainFile == 0) || ($wellBarcodeConstReg ne "")) &&
					      (($scndChainInline == 0) || ($chainDelimiterConstReg ne ""))) {
					    my $tempFileHandle = File::Temp->new();
					    local *STDOUT = $tempFileHandle;
					    processDataSet(abs_path($0), $inputFile, $optPairedInputFile, $optScndChainFile, $optPairedScndChainFile,
							   $chainDelimiterConstReg, $wellBarcodeConstReg, $barcodeDataSet, $noConnectivity, 
							   $mergeTool, $numThread, $vNameDB, $jNameDB, $nameDB, $verbose);
					    seek($tempFileHandle, 0, 0);
					    while (defined(my $line = <$tempFileHandle>)) {
					      $logFileOut->insert("end", $line);
					      $logFileOut->yview("end");
					    }
					    close($tempFileHandle);
					  } else {
					    if (($inputFile eq "") || (! -s $inputFile) ||
						(($pairedEnd != 0) && (($optPairedInputFile eq "") || (! -s $optPairedInputFile))) ||
						(($scndChainFile != 0) && (($optScndChainFile eq "") || (! -s $optScndChainFile))) ||
						(($pairedEnd != 0) && ($scndChainFile != 0) && (($optPairedScndChainFile eq "") || (! -s $optPairedScndChainFile)))) {
					      $logFileOut->insert("end", "Please specify existing input file(s)\n");
					    }
					    if (($scndChainFile != 0) && ($wellBarcodeConstReg eq "")) {
					      $logFileOut->insert("end", "Please specify constant region between well barcode and TCR-data\n");
					    }
					    if (($scndChainInline != 0) && ($chainDelimiterConstReg eq "")) {
					      $logFileOut->insert("end", "Please specify constant region between linked alpha- and beta-chains\n");
					    }
					  }
					});
  # my $quitButton = $mainWindow->Button(-text => "Quit", -command => sub { 
  #                                        $mainWindow->destroy();
  #                                       });
  $inputFileLabel->grid(-row => 0, -column => 0);
  $inputFileEntry->grid(-row => 0, -column => 1, -columnspan => 2);
  $inputFileButton->grid(-row => 0, -column => 3);
  $scndChainFileButton->grid(-row => 5, -column => 0);
  $scndChainInlineButton->grid(-row => 5, -column => 1);
  $pairedEndButton->grid(-row => 5, -column => 2);
  $startButton->grid(-row => 5, -column => 3);
  # $startButton->grid(-row => 6, -column => 0);
  # $quitButton->grid(-row => 6, -column => 3);
  $logFileOut->grid(-row => 6, -column => 0, -columnspan => 4);
  $mainWindow->MainLoop();
}

exit(0);

