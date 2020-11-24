#;-*- Perl -*-
#
#
##################  Alberto Roldan 2-2012  ###############
#
# 1. the program reads Gdis povray outputs

eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
#------------------------------------------------------------------------------------ input
$numargs=$#ARGV+1;
 if ($numargs lt 1 ) { print "\n ---  Please, provide the input file  --- \n "; exit 0 ; }; 
 if (-z $ARGV[0]) { print "\n ---  The ".$ARGV[0]." file is empty ! --- \n \n" ; exit 0 ; };
   $filename=$ARGV[0];
  print "\n  ---  Loaded file is $filename  ---\n\n";
#------------------------------------------------------------------------------------ output
 $outfile="Morphology.pov";
 open OUT, ">>$outfile";
    print OUT ' //__Alberto_Roldan_martinez_____xyz2pov'."\n".' // tags: +ua( sin fondo)';
    print OUT ' +FN (formato .PNG) +H413 (heigh) +W550 (width)'."\n";
    print OUT ' // orthographic'."\n\n";
   system("head -18 $filename >> $outfile");
    print OUT "\n #declare cylin_rad = 0.005;\n";
    print OUT "\n#declare morphology = object{ union{\n";
 close OUT;
#------------------------------------------------------------------------------------ read
 system("sed -e\"1,21d\" -e \"/#declare/d\" -e \"s/<//g\" -e \"s/,//g\" -e \"s/>//g\" $filename >> $outfile.tmp");
open IN, "$outfile.tmp";
    while (<IN>) { $file.=$_;};
close IN;
      @lines=split(/\n/,$file);
#------------------------------------------------------------------------------------ polygon0
     $i=0; @line=split(/\s+/,@lines[0]);
      foreach $l (@line) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; 
        @line=@Nline; @Nline=();
   while (@line[0] ne "polygon") {$i++; };
#------------------------------------------------------------------------------------ polygon
 open OUT, ">>$outfile";
   while (@lines[$i]) {
       @line=split(/\s+/,@lines[$i]);
      foreach $l (@line) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; };
          @line=@Nline; @Nline=();
          $nvertice=@line[2];
#------------------------------------------------------------------------------------ points
      for ($j=0; $j<=$nvertice-1; $j++) {
       @line=split(/\s+/,@lines[$j+$i+1]);
         foreach $l (@line) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; 
           @line=@Nline; @Nline=();
           @vertex[$j]=("< @line[0], @line[1], @line[2] >");
      };
     print OUT " polygon { $nvertice\n";
      for ($k=0; $k<=$nvertice-2; $k++) { print OUT "   @vertex[$k],\n"; };
     print OUT "   @vertex[$nvertice-1]\n";
     print OUT "     texture { ";
     print OUT "pigment { color rgbt <0.70, 0.70, 0.70, 0.10> }\n";
     print OUT "           finish { Phong_Shiny } } no_shadow }\n";
      for ($k=0; $k<=$nvertice-2; $k++) {
     print OUT "   cylinder { @vertex[$k], @vertex[$k+1] cylin_rad\n";
     print OUT "     texture { ";
     print OUT "pigment { color rgbt <0.0, 0.0, 0.0, 0.0> }\n";
      print OUT "           finish { Phong_Shiny } } no_shadow }\n";
      };
     print OUT "   cylinder { @vertex[$nvertice-1], @vertex[0] cylin_rad\n";
     print OUT "     texture { ";
     print OUT "pigment { color rgbt <0.0, 0.0, 0.0, 0.0> }\n";
     print OUT "           finish { Phong_Shiny } } no_shadow }\n";
      $i=$i+$nvertice+3+1;
   };
     print OUT " matrix <  1.0, 0.0, 0.0,\n";
     print OUT "          0.0, 1.0, 0.0,\n";
     print OUT "          0.0, 0.0, 1.0,\n";
     print OUT "          0, 0, 0 >\n";
     print OUT "  } translate < 0.0, 0.0, 0.0> rotate < 0.0, 0.0, 0.0> }\n\n";
    print OUT " #declare axis_size = 0.080;\n";
    print OUT " #declare axis_arrows = object { union {\n";
    print OUT " cylinder {< 0.0, 0.0, 0.0> , < 0.0, 0.0, 10.0> 0.3 texture { pigment {Blue} }  no_shadow  }\n";
    print OUT "     cone {<0.0, 0.0, 10.0>0.8< 0.0, 0.0, 11.0> 0.0 texture { pigment {Blue} }  no_shadow  }\n";
    print OUT " cylinder {< 0.0, 0.0, 0.0> , < 0.0, -10.0, 0.0> 0.3 texture { pigment {Red} }  no_shadow  }\n";
    print OUT "     cone {<0.0, -10.0, 0.0>0.8< 0.0, -11.0, 0.0> 0.0 texture { pigment {Red} }  no_shadow  }\n";
    print OUT " cylinder {< 0.0, 0.0, 0.0> , < 10.0, 0.0, 0.0> 0.3 texture { pigment {Green } }  no_shadow  }\n";
    print OUT "     cone {<10.0, 0.0, 0.0>0.8< 11.0, 0.0, 0.0> 0.0 texture { pigment {Green} }  no_shadow  }\n";
    print OUT " } matrix <  axis_size, 0.0, 0.0,\n";
    print OUT "          0.0, axis_size, 0.0,\n";
    print OUT "          0.0, 0.0, axis_size,\n";
    print OUT "          0, 0, 0 >\n";
    print OUT " }\n";
    print OUT " #declare axis_labels = object { union {\n";
    print OUT ' text { ttf "timrom.ttf" "Z" 0.15, 0 texture { pigment {Blue} }  no_shadow rotate <90,0,0> translate < 0, 0.0, 11.5> }'."\n";
    print OUT ' text { ttf "timrom.ttf" "Y" 0.15, 0 texture { pigment {Red} }  no_shadow rotate <90,0,0> translate < 0, -11.5, 0.0> }'."\n";
    print OUT ' text { ttf "timrom.ttf" "X" 0.15, 0 texture { pigment {Green} }  no_shadow rotate <90,0,0> translate < 11.5, 0, 0.0> }'."\n";
    print OUT " } matrix <  axis_size, 0.0, 0.0,\n";
    print OUT "          0.0, axis_size, 0.0,\n";
    print OUT "          0.0, 0.0, axis_size,\n";
    print OUT "          0, 0, 0 >\n";
    print OUT " }\n";
    print OUT "#declare axis = object{ union { object { axis_arrows } object { axis_labels } rotate < -90.0, 0.0, 0.0> }};\n\n";
    print OUT " object { axis translate < 0.0, 0.0, 0.0 >  rotate < 0.0, 0.0, 0.0>}\n";
    print OUT " object { morphology translate < 0.0, 0.0, 0.0 >  rotate < 0.0, 0.0, 0.0>\n}\n";
 close OUT;
#------------------------------------------------------------------------------------

system(" rm Morphology.pov.tmp");


