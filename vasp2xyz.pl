#;-*- Perl -*-
#
#
#   It vonverts a POSCAR/CONTCAR file in a XYZ file
#
#####   Alberto 12-04-2011

eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
use Math::Trig;

#______________________________________________________________ VASP VERSION
#   $vasp=5;     # for version 4.*
   $vasp=6;     # for version 5.*
#______________________________________________________________
  $filename=$ARGV[0];
  $numargs=$#ARGV+1;
if ($numargs lt 1 ) {
     $filename="POSCAR";
   if (-s $filename) { print "\n ---  The loaded file is ".$filename." --- \n\n ";
    }else{ print "\n ---  Provide a proper file  --- \n\n "; exit 0; };
};

#---------------------------------------------------------------------------------------------------- inbox.pl

##  system("inbox.pl $filename");

#----------------------------------------------------------------------------------------------------- read
  if (-e "$filename.box") { $filename2="$filename.box"; } else { $filename2=$filename; };
 
open (IN,$filename2);
   while (<IN>) { $file.= $_; };
close (IN);
     @file=split(/\n/,$file);
     $header=$file[0];
     $param=$file[1];
     @Xvec=split(/\s+/,@file[2]);
     @Yvec=split(/\s+/,@file[3]);
     @Zvec=split(/\s+/,@file[4]);
     $j=0;
   foreach $c (@Xvec) { if ($c) { @Xvec[$j]=$c; $j++; }; };
     $j=0;
   foreach $c (@Yvec) { if ($c) { @Yvec[$j]=$c; $j++; }; };
     $j=0;
   foreach $c (@Zvec) { if ($c) { @Zvec[$j]=$c; $j++; }; };
   for ($i=0; $i<=2; $i++) {
      @Xvec[$i]=@Xvec[$i]*$param;
      @Yvec[$i]=@Yvec[$i]*$param;
      @Zvec[$i]=@Zvec[$i]*$param;
   };
#---------------------------------------------------------------------------------------------------------------------- ATOMS
   if ($vasp == 6) { $header=@file[5]; };
     @atoms=split(/\s+/, @file[$vasp]); my $natoms=0;    ($natoms+=$_) for @atoms;
     $j=0;
   foreach $a (@atoms) { if (!$a) { }else{ @atoms[$j]=$a; $j++; }; };
     @species=split(/  /,@file[$vasp]); my $nspecies=-1; ($species++ ) for @species;
     @types=split(/\s+/,$header);
     $j=0;
   foreach $t (@types) { if (!$t) { }else{ @types[$j]=$t; $j++; }; };
#----------------------------------------------------------------------------------------------------------------------
open OUT, ">>$filename.xyz";
   print OUT "$natoms\nAlberto vasp2xyz.pl file for VASP 5.*\n";
     @selec=split(//,@file[$vasp+1]);
   if ((@selec[0] eq 'S') or (@selec[0] eq 's')) { $n=1; }else{ $n=0; };
     @direc=split(//,@file[$vasp+1+$n]);
   if ((@direc[0] eq 'D') or (@direc[0] eq 'd')) {
#---------------------------------------------------------------------------------------------------------------------- VECTORS
#														     TRIGONOMETRY
     $XxY=@Xvec[0]*@Yvec[0]+@Xvec[1]*@Yvec[1]+@Xvec[2]*@Yvec[2]; # --- omega x to y scalar product
     $XxZ=@Xvec[0]*@Zvec[0]+@Xvec[1]*@Zvec[1]+@Xvec[2]*@Zvec[2]; # --- beta= x to z
     $YxZ=@Yvec[0]*@Zvec[0]+@Yvec[1]*@Zvec[1]+@Yvec[2]*@Zvec[2]; # --- alfa= y to z
      $Mx=sqrt(@Xvec[0]**2+@Xvec[1]**2+@Xvec[2]**2);  # --- M=modulus
      $My=sqrt(@Yvec[0]**2+@Yvec[1]**2+@Yvec[2]**2);
      $Mz=sqrt(@Zvec[0]**2+@Zvec[1]**2+@Zvec[2]**2);
   $cos1=$YxZ/($Mz*$My);   # --- cosine (alfa)
   $cos2=$$XxZ/($Mx*$Mz);  # --- cosine (Beta)
   $cos3=$XxY/($Mx*$My);
     $omega=acos($cos3);   # --- omega angle
   $sin3=sin($omega);   # --- sine (omega)
#-------------
     $C=($cos1-$cos2*$cos3)/$sin3;
     $Vcte=sqrt(1-$cos1**2-$cos2**2-$cos3**2+2*$cos1*$cos2*$cos3)/$sin3;
#     $Vcte=sqrt(1-$cos1**2-$cos2**2-$cos3**2+2*$Mx*$My*$cos1*$cos2*$cos3)/$sin3; ###### probably this one if Z has a or b component
#-------------
       $sum=@atoms[0];
       $k=0;
      for ($i=1; $i<=$natoms; $i++) {
         @pos=split(/\s+/,$file[$i+$vasp+1+$n]);
         $j=0;
       foreach $p (@pos) { if ($p) { @pos[$j]=$p; $j++; }; };
#---------------------------------------------------------------------------------------------------
           @newpos[0]=@pos[0]*$Mx+@pos[1]*$cos3*$My+@pos[2]*$cos2*$Mz;
           @newpos[1]=@pos[1]*$sin3*$My+@pos[2]*$C*$Mz;
           @newpos[2]=@pos[2]*$Vcte*$Mz;
#---------------------------------------------------------------------------------------------------
         if ($i > $sum) { $k++; $sum=$sum+@atoms[$k]; };
         printf OUT "@types[$k]\t% .7f\t% .7f\t% .7f\n",@newpos[0..2]; 
      };
#.......................................................................................................................................
   }elsif ((@direc[0] eq 'C') or (@direc[0] eq 'c')) {
       $j=0; $sum=@atoms[0];
      for ($i=1; $i<=$natoms; $i++) {
           @pos=split(/\s+/,$file[$i+$vasp+1+$n]);
   if (!@pos[0]) { @pos[0]=@pos[1]; @pos[1]=@pos[2]; @pos[2]=@pos[3]; };
           @pos[0]=@pos[0]*$param; @pos[1]=@pos[1]*$param; @pos[2]=@pos[2]*$param;
         if ($i > $sum) { $j++; $sum=$sum+@atoms[$j]; };
         printf OUT "@types[$j]\t% .6f\t% .6f\t% .6f\n",@pos[0..2];
      };
   }else{ print "\n --- This file doesn't contains direct or cartesian line after atoms position --- \n\n"; exit 0; };
close OUT;

if (-e "$filename.box") { system("rm $filename.box"); };
