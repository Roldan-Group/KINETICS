#;-*- Perl -*-
#
#
########   Alberto 11-2016
#
## input --> OUTCAR CONTCAR N (slab atoms)
#   OUTPUT: Distances.txt 
#
#
#
eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
use Math::Trig;
#---------------------------------------------------------------------------
  $numargs=$#ARGV+1; $slab=$ARGV[0];  
     if ("CONTCAR") { $filename="CONTCAR"; };
  print "\n ---  The loaded file is $filename ---\n"; 
     if ($ARGV[0]) {$slab=$ARGV[0]; } else { $slab=45; };
  print " ---  The atoms in the slab is $slab ---\n\n";
#---------------------------------------------------------------------------
  system ("vasp2xyz.pl $filename");
    $filename="$filename.xyz";
open (IN,$filename);
   while (<IN>) { $file.= $_; };
close (IN);
     @file=split(/\n/,$file);
#---------------------------------------------------------------------------
        $z='', $preZ=''; $origin=0.0; $scale=2.0;
#------------------------------------------------------------------------------------------------------ Distance
   for ($i=2; $i<=$slab+2; $i++) { 
     @line=split(/\s+/,@file[$i]);
     foreach $l (@line) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; @line=@Nline; @Nline=();
       if (@line[3] ge $origin-$scale) { push(@preZ,@line[3]); $origin=@line[3]; }; };  #for
   foreach $z (@preZ) { if ($z ge $origin-$scale) { push(@surfaceAtoms,$z); }; };
#---------------------------------------------------------------------------
   for ($i=$slab+2; $i<=$slab+7; $i++) { 
     @line=split(/\s+/,@file[$i]);
     foreach $l (@line) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; @line=@Nline; @Nline=();
     push(@ringAtoms,$line[3]); };  #for
#---------------------------------------------------------------------------
     $SumSurf=0;
   for ($n=0; $n<=$#surfaceAtoms; $n++) { $SumSurf+=@surfaceAtoms[$n]; }; $AverageSurf=$SumSurf/($#surfaceAtoms+1);
   for ($n=0; $n<=$#ringAtoms; $n++) { $SumRing+=@ringAtoms[$n]; }; $AverageRing=$SumRing/($#ringAtoms+1);
     $distance=$AverageRing-$AverageSurf;
#------------------------------------------------------------------------------------------------------ Angle
      $i=$slab+8; $nH=0; $numberH=5;
   while ($nH < $numberH) {
      @line=split(/\s+/,@file[$i]);
     foreach $l (@line) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; @line=@Nline; @Nline=();
     if (@line[0] ne "H") { $i++; }else{  $zH=$line[3]-$AverageRing;
      system("sed -n \"/neighbor/,/LATTYP/p\" OUTCAR | grep \"  $i  \" > temp");
open (IN2,"temp");
   while (<IN2>) { $temp.= $_; };
close (IN2);
     @temp=split(/\n/,$temp);
      @dist=split(/\s+/,@temp[0]); 
      foreach $d (@dist) { if (($d) or ($d eq 0)) { push(@Nline,$d); }; }; @dist=@Nline; @Nline=();
        $distC=@dist[5];
        $angle=asin($zH/$distC)*360/(2*pi);
        push(@angles,$angle); $nH++; $SumAngles+=$angle; $i++; }; # if H
   }; # while
      $AverageAngle=$SumAngles/$numberH;
#---------------------------------------------------------------------------
  open OUT, ">>Distances.dat";
    printf OUT "Zatoms - Zsurface = distance | angle\n";
    printf OUT "%.5f - %.5f = %.5f %.2f",$AverageRing,$AverageSurf,$distance,$AverageAngle;
  close OUT;

#---------------------------------------------------------------------------
    system("rm CONTCAR.xyz temp");
#---------------------------------------------------------------------------
    system("more Distances.dat");
