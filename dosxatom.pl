#;-*- Perl -*-
#--------ALBERTO ROLDAN-----------2-2011---
#
#usage:file with 2 colums

eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;


$numargs=$#ARGV+1;
$file=$ARGV[0];
$natoms=$ARGV[1];

if ($numargs lt 1 ) { print "\n ---  Two-columns file and number of atoms \n \n" ; exit 0 ;};
if ($numargs lt 2 ) { print "\n ---  Number of atoms \n \n" ; exit 0 ;};

open IN, $file; 
open OUT,">>$ARGV[0].xatom";
   while (<IN>) {
       @cf=split(/\s+/,$_);
     if ($cf[1] =~ /[0-9]+/) {
        $y=@cf[2]/$natoms;
        @point=("@cf[1]  $y");
       print OUT "@point\n";
     }else{ print OUT "\n"; };
   };
close OUT;
close IN ;
