#;-*- Perl -*-
#
#   Give the angle respect to <0,0,1> plane of a three bobies molecules (H2O)
#
#
#
#   INPUT: CONTCAR_files + who?
#   OUTPUT: orientation.txt
#
########   Alberto 01-2014

eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
use Math::Trig;
#-----------------------------------------------------------------------------------------------------
  if (!$ARGV[0]) { $ARGV[0]="CONTCAR"; };
  if ( $ARGV[0] =~ /^[0-9]+$/) { print " ---  Introduce a CONTCAR-format files  ---\n"; exit 0; }
    elsif (!-e $ARGV[0]) { print " ---  The file $ARGV[0] does not exist ---\n"; exit 0; }
    else{ $filename=$ARGV[0]; print " ---  Input file is $filename   ---\n"; };
#----------------------------------------------------------------------------------------------------- who?
#                                                                              position in the CONTCAR file -- Default 2+1
  if (($ARGV[1]) and ( $ARGV[1] =~ /^[0-9]+$/)) { $who=$ARGV[1]; }else{ $who=2; };
#----------------------------------------------------------------------------------------------------- READ
  open (IN,$filename); while (<IN>) { $file.= $_; }; close (IN);
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
#----------------------------------------------------------------------------------------------------- atoms
     @line5=split(/\s+/,@file[5]);
   foreach $l (@line5) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; @head=@Nline; @Nline=();
     @line6=split(/\s+/,@file[6]);
   foreach $l (@line6) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; @atoms=@Nline; @Nline=();
     my $natoms=0; ($natoms+=$_) for @atoms;
     my $nspecies=-1; ($nspecies++) for @atoms;
#----------------------------------------------------------------------------------------------------- who?
     $whoatom=@head[$who];
    print " ---  The three bodies system is centred in $whoatom ---\n";
#----------------------------------------------------------------------------------------------------- selective
     @selec=split(//,@file[7]);
   if ((@selec[0] eq 'S') or (@selec[0] eq 's')) { $s=1; }else{ $s=0; };
#----------------------------------------------------------------------------------------------------- selective
     @selec=split(//,@file[7]);
   if ((@selec[0] eq 'S') or (@selec[0] eq 's')) { $s=1; }else{ $s=0; };
#----------------------------------------------------------------------------------------------------- direct
   if ($s = 1){ $d=8; }else{ $d=7; };
     @direct=split(//,@file[$d]);
   if (@direct[0] eq 'D') { system ("\$HOME/software/cartesian.pl $filename");
        $system="$filename.cart"; 
   }else{ $system=$filename; };
#----------------------------------------------------------------------------------------------------- READ
   @file='';
  open (IN2,$system); while (<IN2>) { $file2.= $_; }; close (IN2);
      @file2=split(/\n/,$file2);
#----------------------------------------------------------------------------------------------------- 3-bodies
    $many=0;
   while ($n < $who) { $many=$many+@atoms[$n]; $n++; }; $many=$many+8;
    $i=$many+1; $k=0; 
   while ($i <= $many+@atoms[$who] ) {
        @line=split(/\s+/,@file2[$i]);
      foreach $l (@line) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; @pos=@Nline; @Nline=();
        @Aposition=@pos[0..2];
# -----------------------------------------------------------------------------------------------------3-bodies elements
     for ($j=$many+@atoms[$who]+$k+1; $j<=$many+@atoms[$who]+$k+2; $j++) {
        @line=split(/\s+/,@file2[$j]);
      foreach $l (@line) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; @pos2=@Nline; @Nline=();
          @Bposition=@pos2[0..2];
         @V=(@Bposition[0]-@Aposition[0],@Bposition[1]-@Aposition[1],@Bposition[2]-@Aposition[2]);
          $Vmodul=sqrt(@V[0]**2+@V[1]**2+@V[2]**2); 
          $Wmodul=@V[2];
        if ( $Vmodul != 0 ) { if ($Wmodul == 0){ $sin=0; }else{ $sin=$Wmodul/$Vmodul; }; };
         $angle=asin($sin)*360/(2*pi);
        push(@angles,$angle);
#print "A=$i;B=$j;zA=@Aposition[2];zB=@Bposition[2];Vm=$Vmodul;Wm=$Wmodul;sin=$sin\n";
      }; #for $j
        $k+=2;
    $i++; 
   }; # while
      $sum=0; $k=0;
     foreach $a (@angles) { $sum=$sum+$a; $k++; }
      $ave=$sum/$k;
#----------------------------------------------------------------------------------------------------- PRINT
  open OUT, ">>orientation.txt";
     print OUT "\n\n# Angle of 3-bodies centered on @head[$who] respect to <0,0,1> plane\n\n";
     print OUT " number of 3-bodies = @atoms[$who]\n";
     print OUT " list of angles ($k):\n";
           foreach $a (@angles) { printf OUT " %.2f",$a; }; print OUT "\n";
     printf OUT " average of angles = %.2f\n",$ave;
     print OUT "\n\n";
  close;
#-----------------------------------------------------------------------------------------------------
system("rm $system");
