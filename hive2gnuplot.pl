#;-*- Perl -*-
#
#
##################  Alberto Roldan 2-2012  ###############
#
# 1. the program reads hive 3D object output & POSCAR-type file

eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
use Math::Trig;
#------------------------------------------------------------------------------------ input
   $numargs=$#ARGV+1;
  if (-s $ARGV[0]) { $filename=$ARGV[0];
     }else{ print "\n ---  Provide a proper file  --- \n\n "; exit 0; };
  if (-s $ARGV[1]) { $filename2=$ARGV[1];
   }elsif (-s "CONTCAR") { $filename2="CONTCAR"; 
   }elsif (-s "POSCAR") { $filename2="POSCAR";
   }elsif (-s "PARCHG") { $filename2="PARCHG";
   }elsif (-s "HPARCHG") { $filename2="HPARCHG";
  }else{ print "\n ---  Provide a second file to check cell geometry --- \n\n "; exit 0; };
  if ($filename) { print "\n  ---  Loaded file is $filename  --- \n "; };
  if ($filename2) { print " ---  Loaded file is $filename2  --- \n "; }; 
#------------------------------------------------------------------------------------ read file 2
open IN2, $filename2;
   while (<IN2>) { $file2.= $_; };
close (IN2);
     @file2=split(/\n/,$file2);
#------------------------------------------------------------------------------------ head
     $header=$file2[0];
   foreach $l (@file2[1]) { if ($l) { push(@Nline,$l); }; }; $param=@Nline[0]; @Nline=();
     @line2=split(/\s+/,@file2[2]);
     @line3=split(/\s+/,@file2[3]);
     @line4=split(/\s+/,@file2[4]);
   foreach $l (@line2) { if ($l) { push(@Nline,$l); }; }; @Xvec=@Nline; @Nline=();
   foreach $l (@line3) { if ($l) { push(@Nline,$l); }; }; @Yvec=@Nline; @Nline=();
   foreach $l (@line4) { if ($l) { push(@Nline,$l); }; }; @Zvec=@Nline; @Nline=();
       $Yveccos=@Yvec[0]/(sqrt(@Yvec[0]**2+@Yvec[1]**2+@Yvec[2]**2));
       $Yvecsen=@Yvec[1]/(sqrt(@Yvec[0]**2+@Yvec[1]**2+@Yvec[2]**2));
   if (($Yveccos != 0) and (@Yvec[1] != 0)) { $square="no";
       print " ---  Unit cell is not XY-square  --- \n ";
    }else{ print " ---  Unit cell is XY-square  --- \n "; };
#------------------------------------------------------------------------------------ read file 1
open IN, $filename;
   while (<IN>) { $file.= $_; };
close (IN);
     @file=split(/\n/,$file);
#------------------------------------------------------------------------------------ head
       $pathfile=$file[0];
       @line1=split(/\s+/,$file[1]);
    foreach $l (@line1) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; @line1=@Nline; @Nline=();
       $points=@line1[0];
#------------------------------------------------------------------------------------ points
       $j=0;
    for ($i=4; $i<=$points+3; $i++) {
        @point=split(/\s+/,$file[$i]);
       foreach $p (@point) { if (($p) or ($p eq 0)) { push(@Npoint,$p); }; };
       if ($square eq "no") {
            $Npoint[0]=sprintf("%.6f",($Npoint[0]/$Yveccos+$Npoint[1]*$Yveccos)); 
            $Npoint[1]=sprintf("%.6f",$Npoint[1]*$Yvecsen);
       };
         $Npoint[1]=$Npoint[1]*(-1);
         @points[$j]="$Npoint[0] $Npoint[1] $Npoint[2]"; @Npoint=(); $j++;
    }; # for
#------------------------------------------------------------------------------------ sort by value
#    @sorted=(sort { $a cmp $b } @points);
@sorted=@points; ###########################################
#------------------------------------------------------------------------------------ print .dat 
open OUT, ">>$filename.dat";
    $k=0; $coord=1;
   while ($sorted[$k]) {
       @line=split(/\s+/,@sorted[$k]);
    foreach $l (@line) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; @line=@Nline; @Nline=();
         $value=$line[$coord];
      $kk=$k;
     while ($line[$coord] eq $value) {
       printf OUT "  %.6f    %.6f   %.6f\n", @line[0],@line[1],@line[2];
          $kk++;
        @line=split(/\s+/,$sorted[$kk]);
       foreach $l (@line) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; @line=@Nline; @Nline=();
     };
     print OUT "\n";
     $k=$kk;
   };
close OUT;
#------------------------------------------------------------------------------------ print .plt
open OUT, ">>$filename.plt";
        @MIN=split(/\s+/,@sorted[$0]);
        @MAX=split(/\s+/,@sorted[$n-1]);
       $Xmin=@MIN[0]; $Xmax=@MAX[0];
       $Ymin=@MIN[1]; $Ymax=@MAX[1];
   print OUT "#!/usr/bin/gnuplot -persist\n";
   print OUT " set terminal postscript eps enhanced color solid \"Helvetica\" 40\n";
   print OUT "set view map\nset surface\n";
   print OUT "set encoding iso_8859_1\nset size ratio -1 5,5\nset ticslevel 0\n";
   print OUT "set noxtics ##border in scale 1,0.5 mirror norotate  offset character 0, 0, 0\n";
##   print OUT "set xtics 0.5  norangelimit\n";
   print OUT "set noytics ##border in scale 1,0.5 mirror norotate  offset character 0, 0, 0\n";
##   print OUT "set ytics 0.5  norangelimit\n";
   print OUT "set noztics\nset nox2tics\nset noy2tics\nset nocbtics\n";
   print OUT "set xlabel \"A\"\n";      #\"A /\\305\"\n"; # to Angstroms
   print OUT "set xlabel  offset character 0, 0, 0 font \"\" textcolor lt -1 norotate\n";
   print OUT "set xrange [ $Xmin : $Xmax ] noreverse nowriteback\n";
   print OUT "set ylabel \"B\"\n";     #\"B /\\305\"\n"; # to Angstroms
   print OUT "set ylabel  offset character 0, 0, 0 font \"\" textcolor lt -1 rotate by -270\n";
   print OUT "set yrange [ $Ymin : $Ymax] noreverse nowriteback\n";
   print OUT "set pm3d explicit at s\nset pm3d scansautomatic\n";
   print OUT "set pm3d interpolate 1,1 flush begin noftriangles nohidden3d corners2color mean\n";
   print OUT "set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB\n";
   print OUT "set palette rgbformulae 7, 5, 15\nset colorbox default\n";
   print OUT "set colorbox vertical origin screen 0.9, 0.2, 0 size screen 0.05, 0.6, 0 front bdefault\n";
   print OUT "set loadpath\nset fontpath\nset fit noerrorvariables\n";
#GNUTERM = \"wxt\"\n";
   print OUT "splot \"$filename.dat\" w pm3d\n\n";
   print OUT "	set output \"./$filename.eps\"; rep; q\n";
close OUT;

system("gnuplot < $filename.plt > gnu.txt");


print "\n#######################################\n";
print "                 gimp:\n";
print "  Filters/Blur/Selective_Gaussian_Blur\n";
print "           Blur radius = 5.0\n";
print "           Max. Delta  = 100\n";
print "#######################################\n\n";

