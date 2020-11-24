#;-*- Perl -*-
#--------ALBERTO ROLDAN-----------07/12---
#
# this program sums the second column of differente files 
#  and divides by the numbre of files=columns
#
# INPUT = list of files
# OUTPUT= sum_averae.dat file with 2 columns


eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 @argv:q' if 0;

#------------------------------------------------------------------------------------ inputs
#$numargs=$#ARGV+1;
# if ($#ARGV lt 2 ) { print "\n ---  Please, provide at least two files  --- \n "; exit 0 ; };
 foreach $file (@ARGV) {
    if (-z $file) { print "\n ---  The $file file is empty ! --- \n \n" ; exit 0 ; 
     }else{ print "\n  ---  Loaded file is $file  ---\n\n"; }; };
#------------------------------------------------------------------------------------ read
     @Xcolumn=0; @Ycolumn=0;
 for ($file=0; $file<=$#ARGV; $file++) { 
     $input=(); 
     $i=0; $j=0;
   open IN, $ARGV[$file]; 
     while (<IN>) {
       $input.=$_;
       @lines=split(/\n/,$input);
       @comment=();
       @line=split(/\s+/,@lines[$i]);
     foreach $l (@line) { if ($l) { push(@Nline,$l); }; };
       @line=@Nline; @Nline=();
       @comment=split(//,@line[0]);
     if ($comment[1] =~ /[0-9]+/) {
          @Xcolumn[$j]=@line[0];
          @Ycolumn[$j]=@Ycolumn[$j]+@line[1];
           $j++;
      }else{
          @Xcolumn[$j]=" ";
          @Ycolumn[$j]=" ";
            $j++;
     }; # if number
       $i++;
     }; # while lines
   close(IN);
 }; # for files
#------------------------------------------------------------------------------------ ponderation
#------------------------------------------------------------------------------------ print
  $k=0;
 open OUT, ">>SUM_AVERAGE.dat";
 for ($k=0; $k<=$j-1; $k++) {
     if (@Ycolumn[$k] =~ /[0-9]+/) { @Ycolumn[$k]=@Ycolumn[$k]/$#ARGV; 
         printf OUT " %+.5f %+.5E\n",@Xcolumn[$k],$Ycolumn[$k];
      }else{ print OUT "\n"; };
  };
 close;
