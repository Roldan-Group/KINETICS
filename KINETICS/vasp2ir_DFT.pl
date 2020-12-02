#! /usr/bin/perl

#-------------
# VASPIR code
#-------------
#
# Beta    0.1 (Jun 2008)
# Version 0.2 (Feb 2009)
#
# By Javier Carrasco
# By Alberto Roldan 4-2011 
#
# Comments:
# -- Tested for VASP 4.6.28
# -- Version 0.2 movie of each eigenvector added 
# -- for VASP 5.2.*
# 
# How to use:
#-------------
# perl vasp2ir.pl OUTCAR.file
#
# Outputs:
# -- IRCAR: full output
# -- IRSPECTRA: spectra (XY format)  
# -- IR.agr: spectrum for xmgrace
# desactivated::-- XXXX-mode.xyz: movie of the XXXX cm-1 eigenvector 
#
# Requirements:
#---------------
# -- Previous frequency calculation with VASP (OUTCAR file needed)
# -- In the VASP input (INCAR), DIPOL correction required (LDIPOL=TRUE, IDIPOL=3)
# -- In the VASP input (INCAR), DIPOL correction required (LDIPOL=FALSE, IDIPOL=4) ------------ to CHECK
#
# Uses:
#-------
# (1) Calculate IR intensities
#
# Implementation details
#------------------------
#
# Based on: A. Varcarcel et al., J. Phys. Chem. B 108, 18297 (2004)
#           &
#           A. Varcarcel PhD Thesis (Chapter 2)
#
# The VASP code already computes the dipole moment components at each
# nuclear configuration used for the construction of the Hessian
# matrix. The VASPIR code computes a numerical estimate of the
# dipole moment derivatives, Ddipol(z)/Dz, on the basis of the atomic
# Cartesian displacements. Note that, because only those vibrational
# modes that give rise an oscillating dipolar moment
# perpendicular to the surface are active, we have considered only
# the z component of the dipole moment. The calculated Ddipol(z)/Dz
# values are then projected onto the basis of normal modes to
# obtain the dynamic dipole moments of the vibrational normal
# modes, Ddipol(z)/DQ_k. The square of the latter is directly related to
# the IR intensities of fundamental bands. It is important to point
# out that vibrational frequencies for the normal modes are mass
# sensitive and, consequently, the intensities of IR bands may also
# change when substituted isotopic species spectra are examined.
#
#        /         \ 2   /                                        \ 2 
#       | Ddipol(z) |   |                    P_ki                  |
# I_k = |-----------| = | Sum_{i=1}^{3N} ----------- Ddipol(z)/Dz_i| 
#       |    DQ_k   |   |                 sqrt(m_i)                |
#        \         /     \                                        /
#
# P_ki / sqrt(m_i) : mass weighted coordinate matrix of the normal
# mode k 
#
# In the current implementation only derivatives of the z-component of
# the dipolmoment with respect to Z coordinate (X,Y neglected) are 
# considered, estimating such value from  the slope of a linear 
# regression employing the displacements around the considered 
# atom along Z     

#
# PARAMETERS DEFINED BY THE USER:
#
$frmin=0.00;       # Spectrum scale MIN
$frmax=4000;      # Spectrum scale MAX
$specint=1.0;     # Spacing between calculated spectrum (in cm-1) 
$smear=0.01;      # Smearing for gaussian smoothing (1/f --> 100 cm-1 --> 0.01)

# -------------------------------------------
#   sqr()  -  Return the square of a number.
# -------------------------------------------

sub sqr {$_[0] * $_[0];}

open (OUT, ">IRCAR");
 print OUT "\n";
 print OUT "********************************************************************\n";
 print OUT "*                              VASPIR                              *\n";
 print OUT "********************************************************************\n";
 print OUT "\n";
 print OUT "   Input file:     $ARGV[0]\n";
 print OUT "   Initial date:  ";
close (OUT);
   system ("date >> IRCAR");
# print "\n";
#---------------------------------------------------------------------------------------------------------------------------
# print " (1/3) Extracting data from $ARGV[0]\n";
#---------------------------------------------------------------------------------------------------------------------------
# Reading data from OUTCAR files:
# (0) considered displacement ($step) and number of them ($nstep)
# (1) dipolmoment ($dm[$k]) for each configuration: (#mnv) * (nstep) + 1 
#      -- $dm[1]: reference (opt geometry)
#      -- $dm[2]: for atom 1, x + ($step) 
#      -- $dm[3]: for atom 1, x - ($step) 
#      -- $dm[4]: for atom 1, y + ($step) 
#      -- ...
#      -- $dm[7]: for atom 1, z - ($step) 
#      -- $dm[8]: for atom 2, x + ($step) 
#      -- ...
# (2) coordinates of each displacement 
# (3) mass weighted coordinate matrix elements of each mnv 
#
#
#---------------------------------------------------------------------------------------------------------------------------
   system("grep LDIPOL $ARGV[0] > tmpldipol");
   open (IN1, tmpldipol);
      while (<IN1>) { if ($.>0) { @target=split(/\s+/,$_);  $ldipol="@target[3]"; $kk++; }; };
   close (IN1);
   system("rm -f tmpldipol");
#---------------------------------------------------------------------------------------------------------------------------
   system("grep IDIPOL $ARGV[0] > tmpidipol");
   open (IN1, tmpidipol);
      while (<IN1>) { if ($.>0) { @target=split(/\s+/,$_);  $idipol=@target[3]; $kk++; }; };
   close (IN1);
   system("rm -f tmpidipol"); 
#---------------------------------------------------------------------------------------------------------------------------
   system("grep NFREE $ARGV[0] > tmpcore1"); 
   open (IN, tmpcore1);
     while (<IN>) { if ($.>0) { @target=split (' '); $nstep=$target[2];  $k++; }; };
   close(IN);
   if (($ibrion eq 7) or ($ibrion eq 8)) {  if ($nstep eq 0) { $nstep=2; }; };
   system("rm -f tmpcore1");
#print "...nstep=$nstep...\n";
#---------------------------------------------------------------------------------------------------------------------------
   system("grep 'POTIM  =' $ARGV[0] > tmpcore1");
   open (IN, tmpcore1);
     while (<IN>) { if ($.>0) { @target=split (' '); $step=$target[2]; $k++; }; };
   close(IN);
   system("rm -f tmpcore1"); 
#print "...potim=$step...\n";
#---------------------------------------------------------------------------------------------------------------------------
   if ($ldipol eq "T") {
       system("grep -B 35 'potential at core' $ARGV[0] > tmpcore1; cp tmpcore1 caca1");
      system("grep 'dipolmoment' tmpcore1 > tmpcore2; cp tmpcore2 caca"); 
      system("rm -f tmpcore1");
   }elsif ($ldipol eq "F") {
      system("grep 'dipolmoment' $ARGV[0] > tmpcore2;"); };
     $k=0;
   open (IN, tmpcore2);
       while (<IN>) { if ($.>0) { @target=split (' '); 
                   $xdm[$k]=$target[1];
                   $ydm[$k]=$target[2];
                   $zdm[$k]=$target[3]; $k++; }; };
#print "...k=$k...\n";
   close(IN);
        $configs=$k;   # TOTAL NUMBER OF CONFIGURATIONS INCLUDING (x0, y0, z0)
     $actat=$configs/(3*$nstep); # ACTIVE ATOMS ALLOWED TO RELAX
   system("rm -f tmpcore2 caca caca1");
#print "..configs=$configs..\n...actat=$actat...\n";
#---------------------------------------------------------------------------------------------------------------------------
   system("grep NIONS $ARGV[0] > tmpcore1");
   open (IN, tmpcore1);
     while (<IN>) {if ($.>0) { @target=split (' ');  $nions=$target[11];  }; };
   close(IN);
   system("rm -f tmpcore1");
#print "...nions=$nions...\n";
#------------------------------------------------------------------------------------------------------- General output info
   open (OUT, ">>IRCAR");
     print OUT "\n";
     print OUT "**************************** Parameters ****************************\n";
     print OUT "\n";
     printf OUT "   Total number of atoms:            %6.0f (NIONS)\n",$nions;
     printf OUT "   Active atoms displaced:           %6.0f \n",$actat;
     printf OUT "   Number of displacements / corrd:  %6.0f (NFREE)\n",$nstep;
     printf OUT "   Total number of configurations:   %6.0f (NFREE*NIONS*3+1)\n",$configs; 
     printf OUT "   Magnitude of the displacements:   %6.4f (POTIM)\n",$step;
   close (OUT);
#---------------------------------------------------------------------------------------------------------------------------
     $matches=$nions+1;
   system ("grep -A $matches POSITION $ARGV[0] | grep -v '\\-\\-' | grep -v POSITION >> tmpcore1");
     $i=1; $j=0; # $j=1 configuration (x0, y0, z0). x+ (j=2), x- (j=3), y+(j=4)... 
   open (IN, tmpcore1);
     while (<IN>) {
       if ($.>0) { @target=split (' ');  $xdisp[$j][$i]=$target[0];
                                         $ydisp[$j][$i]=$target[1];
                                         $zdisp[$j][$i]=$target[2];
#print "...$i...x=$xdisp[$j][$i] y=$ydisp[$j][$i] z=$zdisp[$j][$i]\n";
          if ($i>=$nions){ $i=0; $j++; }; 
       }; $i++;
     };
   close(IN);
   system("rm -f tmpcore1");
   open (OUT, ">>IRCAR");
     print OUT "\n";
     print OUT "********************** Minimum configuration ***********************\n";
     print OUT "\n";
       $w=1;
     for ($w=1;$w<=$nions;$w++){
         printf OUT "   %11.7f %11.7f %11.7f\n",$xdisp[0][$w],$ydisp[0][$w],$zdisp[0][$w]; # freq/atom
     };
   close (OUT);
#---------------------------------------------------------------------------------------------------------------------------
#print " (2/3) Computing Intensities\n";
#---------------------------------------------------------------------------------------------------------------------------
     $matches=$nions+1;
   system ("grep -A $matches THz $ARGV[0] | grep -v '\\-\\-' | grep -v 'dx'> tmpcore1");
     $line=($nions+1)*($actat*3); $i=1; $l=0;
open IN, "tmpcore1";
   while (<IN>) { $file.= $_; };
close (IN);
     @file=split(/\n/,$file);
   while ($l <= $line) {
        @target=split(/\s+/,@file[$l]);
     foreach $t (@target) { if (($t) or ($t == 0)) { push(@Ntarget,$t); }; }; @target=@Ntarget; @Ntarget=();
       if ($target[2] eq 'f'){ $ftag[$i]=@target[1];
                             $fcm[$i]=@target[8];
                             $fmev[$i]=@target[10];
                             $i++;  $j=1;
     }elsif ($target[2] eq 'f/i=') { break;
     }elsif ($target[2] ne 'f') {
                $dx[$j][$i-1]=$target[4];
                $dy[$j][$i-1]=$target[5];
                $dz[$j][$i-1]=$target[6];
               $j++;
      }; $l++;
    };
system("rm -f tmpcore1"); 
#______________________________________________________________________________________________________________________________________________
open (OUT1, ">tmpdipol");
     $j=1; $i=1; $kx=1;$kxx=1; $ky=1;$kyy=1; $kz=1;$kzz=1; $tag1=1; $atag[1]=f; # cebador
    for ($j=1;$j<=$configs-1; $j++){
      for ($i=1;$i<=$nions; $i++){
       if (($tag1/$nstep) > 3) {  $kxx=1;  $kyy=1;  $kzz=1;  $tag1=1; };
        $tag1++;
       if (($xdisp[$j][$i])==($xdisp[0][$i])&&($ydisp[$j][$i])==($ydisp[0][$i])&&($zdisp[$j][$i])==($zdisp[0][$i])){
           $atag[$i]=f; 
        }else{ if ($i ne @activeatoms[$#activeatoms]) { push(@activeatoms,$i); };
        if ($xdisp[$j][$i] > $xdisp[0][$i]+0.0001){
           printf OUT1 "%3.0f 1 1 : (x0) :  %2.5f \t %2.5f %2.5f %2.5f\n",$i, $xdisp[0][$i],$xdm[0],$ydm[0], $zdm[0];
           printf OUT1 "%3.0f 2 1 : (x0+) :  %2.5f \t %2.5f %2.5f %2.5f\n",$i, $xdisp[$j][$i],$xdm[$j],$ydm[$j], $zdm[$j];
        }elsif ($xdisp[$j][$i] < $xdisp[0][$i]-0.0001){
           printf OUT1 "%3.0f 3 1 : (x0-) :  %2.5f \t %2.5f %2.5f %2.5f\n",$i, $xdisp[$j][$i],$xdm[$j],$ydm[$j], $zdm[$j];
        }elsif ($ydisp[$j][$i] > $ydisp[0][$i]+0.0001){
           printf OUT1 "%3.0f 1 2 : (y0) :  %2.5f \t %2.5f %2.5f %2.5f\n",$i, $ydisp[0][$i],$xdm[0],$ydm[0], $zdm[0];
           printf OUT1 "%3.0f 2 2 : (y0+) :  %2.5f \t %2.5f %2.5f %2.5f\n",$i, $ydisp[$j][$i],$xdm[$j],$ydm[$j], $zdm[$j];
        }elsif ($ydisp[$j][$i] < $ydisp[0][$i]-0.0001){
           printf OUT1 "%3.0f 3 2 : (y0-) :  %2.5f \t %2.5f %2.5f %2.5f\n",$i, $ydisp[$j][$i],$xdm[$j],$ydm[$j], $zdm[$j];
        }elsif ($zdisp[$j][$i] > $zdisp[0][$i]+0.0001){
           printf OUT1 "%3.0f 1 3 : (z0) :  %2.5f \t %2.5f %2.5f %2.5f\n",$i, $zdisp[0][$i],$xdm[0],$ydm[0], $zdm[0];
           printf OUT1 "%3.0f 2 3 : (z0+) :  %2.5f \t %2.5f %2.5f %2.5f\n",$i, $zdisp[$j][$i],$xdm[$j],$ydm[$j], $zdm[$j];
        }elsif ($zdisp[$j][$i] < $zdisp[0][$i]-0.0001){
           printf OUT1 "%3.0f 3 3 : (z0-) :  %2.5f \t %2.5f %2.5f %2.5f\n",$i, $zdisp[$j][$i],$xdm[$j],$ydm[$j], $zdm[$j];
       };

 }; }; };
close(OUT1);
#------------------------------------------------------------------------------------------------------- General output info
open (OUT, ">>IRCAR");
   print OUT "\n";
   print OUT "*********************** Displacements summary **********************\n";
   print OUT "\n";
   print OUT "< atom | displac. tag | 1-X, 2-Y, 3-Z | . |  coord | x y z-dipolmoment >\n";
   print OUT "\n";
close (OUT);
  system ("cat tmpdipol >> IRCAR");
#-------------------------------------------------------------------------------------------------------
open (IN, tmpdipol);
   while (<IN>) { if ($. > 0) { @target=split (' ');
        $a=$target[0]; $b=$target[1]; $c=$target[2]; $r[$a][$b][$c]=$target[6];
        $xdipol[$a][$b][$c]=$target[7]; $ydipol[$a][$b][$c]=$target[8]; $zdipol[$a][$b][$c]=$target[9]; 
};};

close(IN);
  system ("rm -f tmpdipol");
open (OUT2, ">>IRCAR");
   print OUT2 "\n";
   print OUT2 "******************* Ddipol/Dz (linear regression) ******************\n";
    $i=1; $j=1; $k=1; $n=$nstep+1;
#   for ($j=1;$j<=$nions;$j++){
 foreach $j (@activeatoms){
      print OUT2 "\n";
      printf OUT2 " Atom %4.0f\n",$j; 
      print OUT2 "--------------------------------------------------------------------\n";
     for ($i=1;$i<=3;$i++) {
        if ($i == 1) { print OUT2 " Displacement along X:\n"; };
        if ($i == 2) { print OUT2 " Displacement along Y:\n"; };
        if ($i == 3) { print OUT2 " Displacement along Z:\n"; };
          $xsumx=0; $xsumx2=0; $xsumxy=0; $xsumy=0; $xsumy2=0;
          $ysumx=0; $ysumx2=0; $ysumxy=0; $ysumy=0; $ysumy2=0;
          $zsumx=0; $zsumx2=0; $zsumxy=0; $zsumy=0; $zsumy2=0;
        for ($k=1;$k<=$n;$k++){
           printf OUT2 " %10.6f \t %5.6f %5.6f %5.6f\n",$r[$j][$k][$i],$xdipol[$j][$k][$i],$ydipol[$j][$k][$i],$zdipol[$j][$k][$i];
#------------------------------------------------------------------------------------------------------- Z
          if ($xdipol[$j][$k][$i] != 0) {
              $xsumx = $xsumx + $r[$j][$k][$i];
              $xsumx2 = $xsumx2 + $r[$j][$k][$i] * $r[$j][$k][$i];
              $xsumxy = $xsumxy + $r[$j][$k][$i] * $xdipol[$j][$k][$i];
              $xsumy  = $xsumy + $xdipol[$j][$k][$i];
              $xsumy2 = $xsumy2 + $xdipol[$j][$k][$i] * $xdipol[$j][$k][$i];
          };
          if ($ydipol[$j][$k][$i] != 0) {
              $ysumx = $ysumx + $r[$j][$k][$i];
              $ysumx2 = $ysumx2 + $r[$j][$k][$i] * $r[$j][$k][$i];
              $ysumxy = $ysumxy + $r[$j][$k][$i] * $ydipol[$j][$k][$i];
              $ysumy  = $ysumy + $ydipol[$j][$k][$i];
              $ysumy2 = $ysumy2 + $ydipol[$j][$k][$i] * $ydipol[$j][$k][$i];
          };
          if ($zdipol[$j][$k][$i] != 0) {
              $zsumx = $zsumx + $r[$j][$k][$i];             
              $zsumx2 = $zsumx2 + $r[$j][$k][$i] * $r[$j][$k][$i];      
              $zsumxy = $zsumxy + $r[$j][$k][$i] * $zdipol[$j][$k][$i];    
              $zsumy  = $zsumy + $zdipol[$j][$k][$i];       
              $zsumy2 = $zsumy2 + $zdipol[$j][$k][$i] * $zdipol[$j][$k][$i]; 
          };
       };
       if ($xsumy != 0) {
           $xm[$j][$i]=($n*$xsumxy-$xsumx*$xsumy)/($n*$xsumx2-sqr($xsumx));
           $xcoef[$j][$i]=($xsumy*$xsumx2-$xsumx*$xsumxy)/($n*$xsumx2-sqr($xsumx));
           $xrfac[$j][$i]=($xsumxy-$xsumx*$xsumy/$n)/sqrt(($xsumx2-sqr($xsumx)/$n)*($xsumy2-sqr($xsumy)/$n));
           $xrfac2[$j][$i]=$xrfac[$j][$i]*$xrfac[$j][$i];
       }else{ $xm[$j][$i]=0; $xcoef[$j][$i]=0; $zrfac2[$j][$i]=0; };
       if ($ysumy != 0) {
           $ym[$j][$i]=($n*$ysumxy-$ysumx*$ysumy)/($n*$ysumx2-sqr($ysumx));
           $ycoef[$j][$i]=($ysumy*$ysumx2-$ysumx*$ysumxy)/($n*$ysumx2-sqr($ysumx));
           $yrfac[$j][$i]=($ysumxy-$ysumx*$ysumy/$n)/sqrt(($ysumx2-sqr($ysumx)/$n)*($ysumy2-sqr($ysumy)/$n));
           $yrfac2[$j][$i]=$yrfac[$j][$i]*$yrfac[$j][$i];
       }else{ $ym[$j][$i]=0; $ycoef[$j][$i]=0; $yrfac2[$j][$i]=0; };         
       if ($zsumy != 0) {
           $zm[$j][$i]=($n*$zsumxy-$zsumx*$zsumy)/($n*$zsumx2-sqr($zsumx));
           $zcoef[$j][$i]=($zsumy*$zsumx2-$zsumx*$zsumxy)/($n*$zsumx2-sqr($zsumx));
           $zrfac[$j][$i]=($zsumxy-$zsumx*$zsumy/$n)/sqrt(($zsumx2-sqr($zsumx)/$n)*($zsumy2-sqr($zsumy)/$n));             
           $zrfac2[$j][$i]=$zrfac[$j][$i]*$zrfac[$j][$i];
       }else{ $zm[$j][$i]=0; $zcoef[$j][$i]=0; $zrfac2[$j][$i]=0; };
          printf OUT2 " Fit: Ddipol/Dx = %9.6f  (B = %9.6f; R^2 = %8.6f)\n",$xm[$j][$i],$xcoef[$j][$i],$xrfac2[$j][$i];
          printf OUT2 " Fit: Ddipol/Dy = %9.6f  (B = %9.6f; R^2 = %8.6f)\n",$ym[$j][$i],$ycoef[$j][$i],$yrfac2[$j][$i];
          printf OUT2 " Fit: Ddipol/Dz = %9.6f  (B = %9.6f; R^2 = %8.6f)\n",$zm[$j][$i],$zcoef[$j][$i],$zrfac2[$j][$i];  
          print OUT2 "--------------------------------------------------------------------\n";
   }; };

close(OUT2); 
   for ($i=1;$i<=($actat*3);$i++){ $int[$i]=0; $Xint=0;$Yint=0;$Zint=0;
     foreach $j (@activeatoms){
       if ($idipol == 1) {
          $int[$i]=$int[$i]+$dx[$j][$i]*$xm[$j][1]+$dy[$j][$i]*$xm[$j][2]+$dz[$j][$i]*$xm[$j][3];  #------------- INTENs....
       }elsif ($idipol == 2) {
          $int[$i]=$int[$i]+$dx[$j][$i]*$ym[$j][1]+$dy[$j][$i]*$ym[$j][2]+$dz[$j][$i]*$ym[$j][3];  
       }elsif ($idipol == 3) {
          $int[$i]=$int[$i] + ($dx[$j][$i]*$zm[$j][1]) + ($dy[$j][$i]*$zm[$j][2]) + ($dz[$j][$i]*$zm[$j][3]);  
       }elsif ($idipol == 4) {
            $Xint=$dx[$j][$i]*$xm[$j][1]+$dy[$j][$i]*$xm[$j][2]+$dz[$j][$i]*$xm[$j][3];
            $Yint=$dx[$j][$i]*$ym[$j][1]+$dy[$j][$i]*$ym[$j][2]+$dz[$j][$i]*$ym[$j][3];  
            $Zint=$dx[$j][$i]*$zm[$j][1]+$dy[$j][$i]*$zm[$j][2]+$dz[$j][$i]*$zm[$j][3];  
          $int[$i]=$int[$i]+$Xint+$Yint+$Zint;
       };
#       for ($k=1;$k<=5;$k++){
#         $xmov[$i][$j][$k]=$xmode[$i][$j]+$dxmov[$i][$j]*($k/5);
#         $ymov[$i][$j][$k]=$ymode[$i][$j]+$dymov[$i][$j]*($k/5);
#         $zmov[$i][$j][$k]=$zmode[$i][$j]+$dzmov[$i][$j]*($k/5);
#     };
      };
       $int2[$i] = sqr($int[$i]);
#----------------------------------------------------------------------------------------------------- frequencies files to videos
#   open (OUT3, ">f$ftag[$i]_$fcm[$i]-mode.xyz");
#      for ($k=1;$k<=5;$k++){ print OUT3 "$nions\n";
#         print OUT3 "frame $k\n";
#        for ($kk=1;$kk<=$nions;$kk++){
#           print OUT3 "H $xmov[$i][$kk][$k] $ymov[$i][$kk][$k] $zmov[$i][$kk][$k]\n";
#      }; };
#   close(OUT3);
    };

#  $int2max=$int2[1];  
  foreach $intens (@int2) { if (!$int2max) { $int2max=$intens; }elsif ($intens > $int2max){ $int2max=$intens; }; };
############################################################################################################################
##------------------------------------------------------------------------------------------------------------------------- Added on Alberto Version 
  for ($i=1;$i<=($actat*3);$i++){ $int2s[$i]=($int2[$i]/$int2max); };
#------------------------------------------------------------------------------------------------------- General output info
open (OUT, ">>IRCAR");
    print OUT "\n";
    print OUT "***************** IR intensities for each frequency ****************\n";
    print OUT "    (Normalized with respect to the most intense signal and *1)\n"; 
    print OUT "\n";
    print OUT " mnv     f (cm-1)      f (meV)   Int (a.u.)\n"; 
    print OUT "--------------------------------------------------------------------\n";
   for ($i=1;$i<=($actat*3);$i++){ printf OUT "%4.0f  %11.4f  %11.4f    %9.4f\n",$i,$fcm[$i],$fmev[$i],$int2s[$i]; };
close (OUT);
##------------------------------------------------------------------------------------------------------------------------- Smearing
#print " (3/3) Generating spectrum\n";

open (OUT, ">>IRSPECTRA");
   for ($fr=$frmin; $fr<=$frmax; ($fr=$fr+$specint)){
   $spec=0;
     for ($i=1;$i<=($actat*3);$i++){                     
        $spec=$spec+$int2s[$i]*exp(-$smear*($fcm[$i]-$fr)**2);  
#$k1=$fcm[$i]*$fr; $k2=($fcm[$i]*$fr)**2; $k3=-$smear*($fcm[$i]*$fr)**2; $k4 =exp(-$smear*($fcm[$i]*$fr)**2); print OUT "$k1 $k2 $k3 $k4\n";
     };
     if ($spec>100){ $spec2=$spec/2; }else{ $spec2=$spec; };                #----------------------------------------------------------- Alberto Version
      printf OUT "%.2f %.8f\n",$fr,$spec2;
   };
close (OUT);
#------------------------------------------------------------------------------------------------------- General output info
open (OUT, ">>IRCAR");
   print OUT "\n";
   print OUT "********************************************************************\n";
   print OUT "\n";
   print OUT " Final date: ";
   system ("date >> IRCAR");
   print OUT "\n";
#   print OUT "Beendet.\n";
   print OUT "\n";
close (OUT);
#   print "\n";
#------------------------------------------------------------------------------------------------ agr.file

# system("mkdir ./intensities");
#print " --- Generating IR.agr file ---\n";
#system ("cp  /home/alberto/software/VASP2IR_SOURCE/IR.agr tmp0.agr");
#   open (OUT, ">>tmp1.agr");
#      for ($i=0;$i<=($actat*3);$i++){
#       if (($fcm[$i] > 500) and ($fcm[$i] < 3500)) {
#         print OUT "@ with string\n@    string on\n@    string loctype world\n@    string g0\n";
#         printf OUT "@    string  %5.2f, %1.2f\n",$fcm[$i],$int2s[$i]+0.01;
#         print OUT "@    string color 1\n@    string rot 90\n@    string font 0\n@    string just 0\n";
#         printf OUT "@    string char size 1.000000\n@    string def \"%d\"\n",$fcm[$i];
#        }; };
#      print OUT "@    s0 comment \"IRSPECTRA\"\n@ target G0.S0\n@ type xy\n# ------------------------------- DATA\n";
#    close (OUT);

#system ("cat tmp0.agr tmp1.agr >> tmp.agr; cat tmp.agr IRSPECTRA >> ./intensities/IR.agr");
#-------------------------------------------------------------------------------------------------------

#system ("rm tmp*.agr ; mv IRCAR IRSPECTRA ./intensities");
#print "  Done.\n";

