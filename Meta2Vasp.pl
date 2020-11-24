#;-*- Perl -*-
#
#
########   Alberto 09-12-2010

eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
use Math::Trig;
#---------------------------------------------------------------------------------------------------------------------- FILE
  $filename=$ARGV[0];
  $numargs=$#ARGV+1;
if ($numargs lt 1 ) { print "\n ---  Provide a proper file  --- \n\n "; exit 0; };

open IN, $filename;
   while (<IN>) { $file.= $_; };
close (IN);
     @file=split(/\n/,$file);
#---------------------------------------------------------------------------------------------------------------------- HEAD
     $header=$file[0];
     $param=$file[1];
#---------------------------------------------------------------------------------------------------------------------- VECTORS
     @Xv=split(/\s+/,@file[2]);
     @Yv=split(/\s+/,@file[3]);
     @Zv=split(/\s+/,@file[4]);
      $j=0;
   foreach $xv (@Xv) { if (!$xv) { }else{ @Avec[$j]=$xv; $j++; }; };
      $j=0;
   foreach $yv (@Yv) { if (!$yv) { }else{ @Bvec[$j]=$yv; $j++; }; };
      $j=0;
   foreach $zv (@Zv) { if (!$zv) { }else{ @Cvec[$j]=$zv; $j++; }; };
#---------------------------------------------------------------------------------------------------------------------- ATOMS
     @at=split(/\s+/, @file[5]); my $natoms=0; ($natoms+=$_) for @at;
     $j=0;
   foreach $a (@at) { if (!$a) { }else{ @atoms[$j]=$a; $j++; }; };
#---------------------------------------------------------------------------------------------------------------------- CARTESIAN
     $cartesian=$file[6];
#====================================================================================================================== PRINT
open OUT, ">>$filename.tmp";
   print OUT "$header\n$param\n";
   printf OUT "  %3.10f  %3.10f  %3.10f\n",@Avec[1],@Avec[2],@Avec[0];
   printf OUT "  %3.10f  %3.10f  %3.10f\n",@Bvec[1],@Bvec[2],@Bvec[0];
   printf OUT "  %3.10f  %3.10f  %3.10f\n",@Cvec[1],@Cvec[2],@Cvec[0];
   print OUT "   @atoms\n$cartesian\n";
close OUT;
#---------------------------------------------------------------------------------------------------------------------- POSITIONS
     $i=1;
   while (@file[$i+6]) {
       @P=split(/\s+/,@file[$i+6]);
        $j=0;
     foreach $pos (@P) { if (!$pos) { }else{ @position[$j]=$pos; $j++; }; };
#---------------------------------------------------------------------------------------------------------------- correct POSITION
      $tmp=@position[0];
      @position[0]=@position[1];
      @position[1]=@position[2];
      @position[2]=$tmp;
#---------------------------------------------------------------------------------------------------------------- upside down
#---------------------------------------------------------------------------------------------------------- xy reflexion
   $Mx[0]=1; $Mx[1]=0; $Mx[2]=0;
   $My[0]=0; $My[1]=1; $My[2]=0;
   $Mz[0]=0; $Mz[1]=0; $Mz[2]=-1;
     $newPos[0]=(@position[0]*$Mx[0]+@position[1]*$My[0]+@position[2]*$Mz[0])-abs(@Cvec[1]);
     $newPos[1]=(@position[0]*$Mx[1]+@position[1]*$My[1]+@position[2]*$Mz[1])-abs(@Cvec[2]);
     $newPos[2]=(@position[0]*$Mx[2]+@position[1]*$My[2]+@position[2]*$Mz[2])+@Cvec[0];
#====================================================================================================================== PRINT
open OUT, ">>$filename.tmp";
   printf OUT "  %3.10f  %3.10f  %3.10f\n",@newPos[0..2];
close OUT;
      $i++;
    }; # while


