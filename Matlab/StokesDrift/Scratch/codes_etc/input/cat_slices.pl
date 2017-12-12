#! /usr/bin/perl
use strict;
use warnings;

{         # main

 my $numArgs = $#ARGV + 1;
 if($numArgs != 2 )
 {
  print "usage: ./cat_slices.pl islice numprocs \n";
  last;
 }
 
 my $islice=$ARGV[0];
 my $numprocs=$ARGV[1];
 my $cmd;
 
 
#------------------------------------------
#  loop through each of the time slices
#------------------------------------------
 for( my $tid=$islice-1; $tid <= $islice-1; $tid++ ) {
  
  $cmd = "rm -f slice*.nc";
  system($cmd);
  
  $cmd = "./concat_across_z.pl $tid $numprocs";
  print "$cmd \n";
  system($cmd);
  
  
 }
 
 
 
}  # end of main
