#;-*- Perl -*-
#
#   INPUT: CONTCAR_files 
# 	adds the different CONTCAR files in a square-like
#       the program squentially attaches images along the X axis until $nx
#
#   OUTPUT: file.xy.pov
#
########   Alberto 01-2014

eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#-----------------------------------------------------------------------------------------------------
  if ($#ARGV lt 1 ) { print " ---  Introduce at least 2 files  ---\n"; exit 0; };
#----------------------------------------------------------------------------------------------------- images
    $nx=''; $ny='';
 if ( $ARGV[0] =~ /^[0-9]+$/) { $nx=shift(@ARGV); $ny=int(($#ARGV+1)/$nx); } 
  elsif ( $ARGV[$#ARGV] =~ /^[0-9]+$/) { $nx=pop(@ARGV); $ny=(($#ARGV+1)/$nx); };
  if (! $nx) {  $n=int(sqrt($#ARGV+1));
   if ($n*$n == $#ARGV+1) { $nx=$n; $ny=$n; }else{ $ny=$n; $nx=$n+1; }; };
 for ($i=0; $i<=$#ARGV; $i++) {
  if (!-e $ARGV[$i]) { print " ---  The file $ARGV[$i] does not exist ---\n"; exit 0; };};
    printf "---  The total number of images is %d ---\n",$#ARGV+1;
    print "---  The output is a $nx x $ny images ---\n";
#----------------------------------------------------------------------------------------------------- Nvectors
    @NX=''; @NY=''; @NZ=''; $filename=''; @file=''; $system='';
      $x=1; $y=1;
#----------------------------------------------------------------------------------------------------- READ
 foreach $arg (@ARGV) { $file=(); 
  if ($arg) {
    $filename=$arg;
  open (IN,$filename); while (<IN>) { $file.= $_; }; close (IN);
      @file=split(/\n/,$file);
#----------------------------------------------------------------------------------------------------- selective
     @selec=split(//,@file[7]);
     foreach $p (@selec) {if (($p) or ($p eq 0)) { push(@Nline,$p); }; }; @selec=@Nline; @Nline=();
   if ((@selec[0] eq 'S') or (@selec[0] eq 's')) { $s=1; }else{ $s=0; };
#----------------------------------------------------------------------------------------------------- direct
   if ($s = 1){ $d=8; }else{ $d=7; };
     @direct=split(//,@file[$d]);
     foreach $p (@direct) {if (($p) or ($p eq 0)) { push(@Nline,$p); }; }; @direct=@Nline; @Nline=();
   if (@direct[0] eq 'D') { system ("\$HOME/software/cartesian.pl $filename");
        $system="$filename.cart"; 
   }else{ $system=$filename; };
#----------------------------------------------------------------------------------------------------- where
# alberto
#   print "i=$i....x=$x-nx=$nx--y=$y-ny=$ny\n";
   if ($i < $nx*$y+1) { ($HEAD,$ATOMS)=&attach($system,$x,$y); $x++; $y=1;
    }else{ ($HEAD,$ATOMS)=&attach($system,$x,$y); $y++; };
     @HEAD=@$HEAD; @ATOMS=@$ATOMS;
   if (-e "$filename.cart") { system ("rm $filename.cart"); };
 }; }; # for files
#----------------------------------------------------------------------------------------------------- WRITE
  open OUT, ">attached";
    print OUT "@HEAD\n  1.0000\n";
      printf OUT "     %.15f  %.15f  %.15f\n",@NX;
      printf OUT "     %.15f  %.15f  %.15f\n",@NY;
      printf OUT "     %.15f  %.15f  %.15f\n",@NZ;
    print OUT "  @HEAD\n  @ATOMS\n";
      printf OUT "Selective Dynamics\nCartesian\n";
  close OUT;
   foreach $H (@HEAD) { system ("cat $H.positions.temp >> attached; rm $H.positions.temp"); };
#----------------------------------------------------------------------------------------------------- to .pov
#  system ("\~/software/vasp2pov.pl attached 0" );

###############################################################################################################################
#----------------------------------------------------------------------------------------------------------- to attach
 sub attach {
         ($system,$x,$y)=@_;
       $file=''; @file='';
#----------------------------------------------------------------------------------------------------- READ
  open (IN,$system); while (<IN>) { $file.= $_; }; close (IN);
     @file=split(/\n/,$file);
#----------------------------------------------------------------------------------------------------- param
    foreach $l (@file[1]) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; $param=@Nline[0]; @Nline=();
#----------------------------------------------------------------------------------------------------- vectors
     @line2=split(/\s+/,@file[2]);
     @line3=split(/\s+/,@file[3]);
     @line4=split(/\s+/,@file[4]);
   foreach $l (@line2) { if ($l) { $L=$l*$param; push(@Nline,$L); }; }; @Xvec=@Nline; @Nline=();
   foreach $l (@line3) { if ($l) { $L=$l*$param; push(@Nline,$L); }; }; @Yvec=@Nline; @Nline=();
   foreach $l (@line4) { if ($l) { $L=$l*$param; push(@Nline,$L); }; }; @Zvec=@Nline; @Nline=();
#----------------------------------------------------------------------------------------------------- Nvectors
   for ($j=0; $j<=2; $j++) { @NX[$j]=@Xvec[$j]*$x; @NY[$j]=@Yvec[$j]*$y; };
       @NZ=@Zvec;
#   for ($j=0; $j<=2; $j++) { $Nvec[$j]=@NX[$j]+@NY[$j]+@NZ[$j]; };
#----------------------------------------------------------------------------------------------------- atoms
     @line5=split(/\s+/,@file[5]);
   foreach $l (@line5) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; @head=@Nline; @Nline=();
     @line6=split(/\s+/,@file[6]);
   foreach $l (@line6) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; @iatoms=@Nline; @Nline=();
     $natoms=0; $nspecies=0; %atoms=();
   foreach $a (@head) { $atoms{$a}=$iatoms[$nspecies]; $natoms+=$iatoms[$nspecies]; $nspecies++; };
     
#----------------------------------------------------------------------------------------------------- selective
     @selec=split(//,@file[7]);
   foreach $l (@selec) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; @selec=@Nline; @Nline=();
   if ((@selec[0] eq 'S') or (@selec[0] eq 's')) { $s=1; }else{ $s=0; };
#----------------------------------------------------------------------------------------------------- positions
     $k=1; $d=7+$s;
 foreach $h (@head) {
  open OUT, ">>$h.positions.temp";
    for ($p=$k; $p<=$atoms{$h}+$k-1; $p++) { 
       @position=split(/\s+/,@file[$d+$p]);
     foreach $pst (@position) { if (($pst) or ($pst eq 0)) { push(@Nline,$pst); }; }; @pos=@Nline; @Nline=();
	     $Npos[0]=$pos[0]*$param+(@Xvec[0]*($x)+@Yvec[0]*($y)+@Zvec[0]);
             $Npos[1]=$pos[1]*$param+(@Xvec[1]*($x)+@Yvec[1]*($y)+@Zvec[1]);
      if ($pos[4]) { printf OUT "  %.15f  %.15f  %.15f  @pos[3..5]\n",$Npos[0],$Npos[1],$pos[2];
      }else{         printf OUT "  %.15f  %.15f  %.15f\n",$Npos[0],$Npos[1],$pos[2]; };

#      if ($h eq "Au"){ print "$system-->$h=$atoms{$h}-$d+$p=@pos\n"; };

        @pos='';
    }; # for p=natoms
  close OUT;
      $k=$k+$atoms{$h};  
      $addATOMS{$h}=$addATOMS{$h}+$atoms{$h};
 }; # foreach
      $k=0;
   foreach $h (@head) {  $exist=0; foreach $addH (@addHEAD) { if ($h eq $addH) { $exist=1; }; }; 
    if ($exist == 0) { push(@addHEAD,$h); };
      $addATOMS[$k]=$addATOMS[$k]+$atoms{$h}; $k++; }; # foreach head
       $file=''; @file='';
   return (\@addHEAD,\@addATOMS);
 }; #sub attach
######################################################################################################

