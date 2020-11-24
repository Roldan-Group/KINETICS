#;-*- Perl -*-
#
#
##################  Alberto Roldan 10-2013  ###############

eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;

#----------------------------------------------------------------------------------------------------
     $filename=$ARGV[0];
     $numargs=$#ARGV+1;
   if (!-e $filename) { print "\n ---  Provide a valid file  ---\n "; exit 0; }; 
   print "\n ---  Loaded file is $filename  ---\n ";
#---------------------------------------------------------------------------------------------------- print
 open OUT, ">>$filename.pov0";
     print OUT "//__Alberto_Roldan_martinez_____vmd2pov.pl\n";
     print OUT "// tags: +ua (sin fondo) +FN (formato .PNG) +H (heigh) +W (width)\n";
     print OUT "// tags: +KFF (final frame) +KFI (initial frame) +KF (final clock) +KI (initial clock) +KC (cyclic on)\n";
     print OUT "// orthographic\n\n";
     print OUT " #include \"colors.inc\"\n #include \"finish.inc\"\n #include \"metals.inc\"\n\n";
     print OUT "//camera {orthographic location <-1.4,0.75,-1.15> right -2/2*x angle 45.000 look_at <-1.4,0.75,0> rotate <0,0,0>} // side rotate 90\n";
     print OUT "camera {orthographic location <-1.2,0.3,-2> right -2/2*x angle 45.000 look_at <-1.2,0.3,0> rotate <0,0,0>} // top\n";
     print OUT "light_source { <0.0, 0.0, -10.0> rgb<1.0,1.0,1.0> rotate < 0.0, 0.0, 0.0 >}\n";
     print OUT "background { color rgb < 1.000, 1.000, 1.000 > }\n\n";
 close OUT;
#----------------------------------------------------------------------------------------------------
     $mesh=0; $i=0;
 open IN, $filename;
   while (<IN>) { $file.=$_; };
 close (IN);
     @line=split(/\n/,$file); $file=();
     $mesh=0; $atom=0; $i=0;
#---------------------------------------------------------------------------------------------------- mesh
     open OUT, ">>$filename.mesh"; print OUT " #declare charge = object{ union{\n"; close OUT;
#----------------------------------------------------------------------------------------------------
   while (@line[$i]) {
        @checkmesh=split(/\s+/,@line[$i]);
     if (@checkmesh[0] eq "mesh2") { $mesh=1; $atom=0; };
        @checkline=split(/</,@line[$i]);
     if ((@checkline[0] eq "VMD_sphere(") or (@checkline[0] eq "VMD_cylinder(")) { $mesh=0; $atom=1; };
#---------------------------------------------------------------------------------------------------- mesh
     if (($mesh == 1) and (@checkline[0] eq "  VMDC(")) {
         open OUT, ">>$filename.mesh"; print OUT "   texture { pigment { rgbt <0.000,0.000,1.000,0.500> }}\n"; close OUT;
      }elsif (($mesh == 1) and (@checkline[0] ne "  VMDC(")) {
         open OUT, ">>$filename.mesh"; print OUT "@line[$i]\n"; close OUT; };
#---------------------------------------------------------------------------------------------------- atoms
     if ($atom == 1) { open OUT, ">>atoms"; print OUT "@line[$i]\n"; close OUT; };
    @line[$i]=();
    $i++;
   }; #while
#---------------------------------------------------------------------------------------------------- mesh
     open OUT, ">>$filename.mesh"; 
       print OUT "}\nmatrix <1.0, 0.0, 0.0,\n        0.0, 1.0, 0.0,\n        0.0, 0.0, 1.0,\n";
       print OUT "        0.0, 0.0, 0.0>\ntranslate <0.0, 0.0, 0.0> rotate <0.0, 0.0, 0.0> }\n\n";
     close OUT;

   system("\$HOME/software/VMD2POV_SOURCE/vmd2pov.sh atoms; rm atoms; mv atoms.pov1 $filename.pov1");
   print " --- declare Ncharge and/or Pcharge ---\n\n";

 system("~/software/VMD2POV_SOURCE/vmd2mesh.sh $filename");


