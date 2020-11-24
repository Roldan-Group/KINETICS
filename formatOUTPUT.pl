#;-*- Perl -*-
#
########   Alberto Roldan
#
#
#
eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
use Scalar::Util qw(looks_like_number);
#-------------------------------------------------------------------------------------------------------------------- prelimimar considerations
    $numargs=$#ARGV+1;
  if ($numargs lt 1 ) { $file="OUTPUT.kinetics.in.dat";
     if (!-e $file) { while (!-e $file) { print " ---  Introduce an existing file --- \n"; $file=<STDIN>; }; };
  }elsif( $numargs eq 1 ) { $file=$ARGV[0];
     if (!-e $file) { while (!-e $file) { print " ---  Introduce a valid file --- \n"; $file=<STDIN>; }; }; };

open IN, $file;
  while (<IN>) { $in.=$_;};
close IN;
    @in=split(/\n/,$in);
    $i=0;
open OUT, ">$file.formated";
  while (@in[$i]) {
      @line=split(/\s+/,@in[$i]); @Nline=();
    foreach $l (@line){ if (($l) or ($l =~ 0.0)) { push(@Nline,$l); }; }; @line=@Nline; @Nline=();
    if (looks_like_number(@line[0])){
        printf OUT "%.4f\t%.3f\t%.2f\t%.5f",@line[0..3];
       for($j=4; $j<=$#line; $j++) { 
             @number=split(//,@line[$j]);
             $exp='no';
          foreach $n (@number) { if ($n eq "e") {$exp='yes'; };};
          if ($exp eq 'no') {
#              $k=1;
#            while ((@line[$j] < 0) or ($k < $i)) {
#               @line2=split(/\s+/,@in[$i-$k]); @Nline=();
#              foreach $l (@line2){ if (($l) or ($l =~ 0)) { push(@Nline,$l); }; }; @line2=@Nline; @Nline=();
#              if (@line2[$j] >= 0){ @line[$j]=@line2[$j]; $k=$i; }else{ $k++; };
#            }; # while
            printf OUT "\t%.15f",@line[$j];
          }else{ @number=split(/e/,@line[$j]);
            if (@number[0] < 0) { @line[$j]="\t1e@number[1]"; };
#----------------------------------------------------------------------------------------------- proper by time demanding
#              $k=1;
#            while ((@number[0] < 0) or ($k < $i)) {
#               @line2=split(/\s+/,@in[$i-$k]); @Nline=();
#              foreach $l (@line2){ if (($l) or ($l =~ 0)) { push(@Nline,$l); }; }; @line2=@Nline; @Nline=();
#                @number2=split(/e/,@line2[$j]);
#              if (@number2[0] > 0.0){ @number[0]=@number2[0]; @line[$j]=@line2[$j]; $k=$i;}else{ $k++; };
#              }; # while
#---------------------------------------------------------------------------------------------
          print OUT "\t@line[$j]"; };
        
       }; #for
      print OUT "\n";
     }else{ 
        foreach $l (@line){ if (($l) and ($l ne "\"")) { print OUT "\t$l"; }; }; @line=@Nline; @Nline=();
        print OUT "\n"; }; 
   $i++; }; #while
close OUT;

