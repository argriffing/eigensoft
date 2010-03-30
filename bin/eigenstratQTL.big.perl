#!/usr/bin/perl

use Getopt::Std ;
use File::Basename ;

$MAXSIZE = 2000000000;

### process flags
my @flaglist = ("i","j","p","l","o");
$x = @ARGV;
for($n=0; $n<$x; $n++)
{
  foreach $flag (@flaglist) 
  {
    if($ARGV[$n] eq "-$flag") { $specified{$flag} = 1; }
  }
}
foreach $flag ("i","j","p","o")
{
  unless($specified{$flag}) { die("OOPS -$flag flag not specified"); }
}
getopts('i:j:p:l:o:',\%opts);
$i = $opts{"i"}; 
$j = $opts{"j"}; 
$p = $opts{"p"}; 
$l = 10; if($specified{"l"}) { $l = $opts{"l"}; }
$o = $opts{"o"}; 

open(OUT,">$o") || die("COF");
print OUT ("eigenstratQTL.big.perl program run using parameters\n");
print OUT (" -i $i\n");
print OUT (" -j $j\n");
print OUT (" -p $p\n");
print OUT (" -l $l\n"); 
print OUT (" -o $i\n");
print OUT ("\n");
print OUT ("Chisq EIGENSTRAT\n");

$NSAMPLES = `wc -l $p`; $NSAMPLES=$NSAMPLES+0;
$maxnSNP = int($MAXSIZE / $NSAMPLES);
open(IN,$i) || die("COF");
$R = 0;
while(1) # process at most $maxnSNP at a time
{
  $thisi = "$i.MAX$maxnSNP"; $thiso = "$o.MAX$maxnSNP";
  open(THISI,">$thisi") || die("COF");
  $m = 0;
  while($line = <IN>)
  {
    print THISI ("$line");
    $m++;
    if($m == $maxnSNP) { last; }
  }
  close(THISI);
  if($m == 0) { system("rm $thisi"); system("rm $thiso"); last; }

  $command = "eigenstratQTL";
  $command .= " -i $thisi";
  $command .= " -j $j";
  $command .= " -p $p";
  $command .= " -l $l";
  $command .= " -o $thiso";
  print("$command\n");
  system("$command");

  open(THISO,$thiso) || die("COF");
  for($n=0; $n<8; $n++) { $line = <THISO>; }
  while($line = <THISO>) { print OUT ("$line"); }
  close(THISO);
 
  $R++;
}

