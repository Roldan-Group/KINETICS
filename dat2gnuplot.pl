#;-*- Perl -*-
#
#
##################  Alberto Roldan 2-2012  ###############

eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
use Scalar::Util qw(looks_like_number);
#------------------------------------------------------------------------------------------------

   $filename=$ARGV[0];
   $numargs=$#ARGV+1;
  if ($numargs lt 1 ) { print "\n ---  Provide a proper file  --- \n\n "; exit 0; 
  }else{ print "\n ---  Loaded file is $filename  --- \n\n "; };

open IN, $filename;
   while (<IN>) { $file.= $_; };
close (IN);
open OUT, ">>$filename.gnu";
     @file=split(/\n/,$file); 
   foreach $f (@file) {
       @line=split(/\s+/,$f);
     foreach $l (@line){ if (($l) or ($l =~ /^-?(?:\d+(?:\.\d*)?&\.\d+)$/) or ($l =~ 0)) { push(@Nline,$l); }; }; @line=@Nline; @Nline=();
    if (looks_like_number($line[0])) {
     if (!$value) { print OUT "\n"; $value="$line[0]"; };
     if (($line[0] == $value) or ($line[0] eq "0")) {
#        print OUT " @line\n";
#             if (@line[2] < 0.0) { @line[2]=0.0; }; ############################################################################################################ MODIFY DATA TO 0
#        printf OUT " %.5f %.5f %.5f\n",@line[0],@line[1],@line[2];
       if (!$Xmin){$Xmin=$line[0]; }elsif($Xmin<$line[0]){$Xmin=$line[0];};
       if (!$Xmax){$Xmax=$line[0]; }elsif($Xmax>$line[0]){$Xmax=$line[0];};
       if (!$Ymin){$Ymin=$line[1]; }elsif($Ymin<$line[1]){$Ymin=$line[1];};
       if (!$Ymax){$Ymax=$line[1]; }elsif($Ymax>$line[1]){$Ymax=$line[1];};
        printf OUT "@line[0]  @line[1]  @line[2]\n";
     }else{ undef($value);
     if (!$value) { print OUT "\n"; $value="$line[0]"; };
     if (($line[0] == $value) or ($line[0] eq "0")) {
#        print OUT " @line\n";
#             if (@line[2] < 0.0) { @line[2]=0.0; }; ############################################################################################################ MODIFY DATA TO 0
#        printf OUT " %.5f %.5f %.5f\n",@line[0],@line[1],@line[2];
       if (!$Xmin){$Xmin=$line[0]; }elsif($Xmin<$line[0]){$Xmin=$line[0];};
       if (!$Xmax){$Xmax=$line[0]; }elsif($Xmax>$line[0]){$Xmax=$line[0];};
       if (!$Ymin){$Ymin=$line[1]; }elsif($Ymin<$line[1]){$Ymin=$line[1];};
       if (!$Ymax){$Ymax=$line[1]; }elsif($Ymax>$line[1]){$Ymax=$line[1];};
        printf OUT "@line[0]  @line[1]  @line[2]\n";
 }; }; }; }; 

close OUT;
#------------------------------------------------------------------------------------ print .plt
open OUT, ">>$filename.plt";
   print OUT "#!/usr/bin/gnuplot -persist\n";
   print OUT " set terminal postscript eps enhanced color solid \"Arial\" 80\n";
   print OUT "#set view 60, 350, 1, 1\nset view 62, 198, 1, 1\n#set view 62, 10, 1, 1\nset surface\n";
   print OUT "set encoding iso_8859_1\nset size ratio -1 5,5\nset ticslevel 0\n";
   print OUT "set xtics border in scale 5,5 mirror norotate offset 0,-0.3\n";
   print OUT "set xtics   norangelimit\n";
   print OUT "set ytics border in scale 5,5 mirror norotate offset character 0, 0, 0\n";
   print OUT "set ytics   norangelimit\n";
   print OUT "set ztics border in scale 5,5 mirror norotate offset character 0, 0, 0\n";
   print OUT "set ztics   norangelimit\n";
   print OUT "set nox2tics\nset noy2tics\nset nocbtics\n";
   print OUT "set xlabel \"pH\"\n";      #\"A /\\305\"\n"; # to Angstroms
   print OUT "set xlabel  offset 0,-1 font \"\" textcolor lt -1 rotate parallel\n";
   print OUT "set xrange [ $Xmin : $Xmax ] noreverse nowriteback\n";
   print OUT "set ylabel \"Temperature /K\"\n";     #\"B /\\305\"\n"; # to Angstroms
   print OUT "set ylabel  offset 0,-1 font \"\" textcolor lt -1 rotate parallel\n";
   print OUT "set yrange [ $Ymin : $Ymax] noreverse nowriteback\n";
   print OUT "set zlabel \"coverage /ML\"\n";
   print OUT "set zlabel  offset 0,-1 font \"\" textcolor lt -1 rotate parallel\n";
   print OUT "set zrange [ : ] noreverse nowriteback\n";
   print OUT "set pm3d explicit at s\nset pm3d scansautomatic\n";
   print OUT "set pm3d interpolate 1,1 flush begin noftriangles nohidden3d corners2color mean\n";
   print OUT "set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB\n";
   print OUT "set palette rgbformulae 7, 5, 15\nset colorbox default\n";
   print OUT "set colorbox vertical origin screen 0.9, 0.2, 0 size screen 0.05, 0.6, 0 front bdefault\n";
   print OUT "set loadpath\nset fontpath\nset fit noerrorvariables\n";
#GNUTERM = \"wxt\"\n";
   print OUT "splot \"$filename.gnu\" w pm3d\n\n";
   print OUT "	set output \"./$filename.eps\"; rep; q\n";
close OUT;

system("gnuplot < $filename.plt > gnu.txt");
