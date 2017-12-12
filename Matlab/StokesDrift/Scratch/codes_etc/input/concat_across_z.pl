#! /usr/bin/perl
use strict;
use warnings;

{# main
  my $numArgs = $#ARGV + 1;
   if($numArgs != 2 )
   {
    print "usage: ./concat_across_z.pl time_index np \n";
    last;
   }

  #my $NCKS='/opt/local/bin/ncks';
  #my $NCPDQ='/opt/local/bin/ncpdq';
  #my $NCRCAT='/opt/local/bin/ncrcat';
  my $NCKS='/usr/bin/ncks';
  my $NCPDQ='/usr/bin/ncpdq';
  my $NCRCAT='/usr/bin/ncrcat';
  my $np=$ARGV[1];
  my $datadir='../../output/2D/';
  my $fnroot='YZplane_0-';
  my $istep=sprintf("%06d",$ARGV[0]);
  my ($pid,$ncfile,$tmpfile,$newfile,$cid,$cmd);

  #=====================================================
  #  Loop over pid from 0 to np-1 
  #=====================================================
  for($pid=0;$pid<$np;$pid++){
  
    $cid = sprintf("%03d", $pid);
    $ncfile=$datadir . $fnroot . $cid . ".nc";
    $tmpfile= "tmp" . $cid . ".nc";
    
    $cmd="$NCKS -O -d timedimension,$istep $ncfile $tmpfile";
    system($cmd);
    #print "$cmd \n";
 
    if($pid==0){
     $cmd="$NCPDQ -O -a kdimension,timedimension $tmpfile $tmpfile >> /dev/null"; 
    }
    else{
     $cmd="$NCPDQ -O -a kdimension,timedimension $tmpfile $tmpfile >> /dev/null";
    }
   
    system($cmd);       
  }
  
  $newfile="slice_$istep" . ".nc";
  $cmd="$NCRCAT -O " . 'tmp*.nc ' . "$newfile";
  system($cmd);
  
  $cmd='rm -f tmp*.nc';
  system($cmd);
  
  
} # end main
