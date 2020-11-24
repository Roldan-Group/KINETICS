#;-*- Perl -*-
#
#   INPUT: CONTCAR_file slab_atom_numbers
# 	replicates the catalyst surface by 3 (Default) using multiplecell
#	adds the molecules in one cell to the replicated surface
#
#   OUTPUT: file.x3.pov
#
########   Alberto 10-2012

eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
use Math::Trig;

#-----------------------------------------------------------------------------------------------------
    $numargs=$#ARGV+1; $filename=$ARGV[0]; $slabatoms=$ARGV[1]; 
  if ($numargs lt 1 ) {
           print " ---  Introduce the file name --- \n ";
              $filename=<STDIN>; chomp($filename);
     if (!-e $filename) {
       while (!-e $filename) {
           print " ---  Introduce an existing file --- \n ";
              $filename=<STDIN>; chomp($filename);
       }; };
           print " ---  Introduce the number of atoms in the slab --- \n ";
              $slabatoms=<STDIN>; chomp($slabatoms);
      if ($slabatoms !~ /^[0-9]+$/) {
        while ($slabatoms !~ /^[0-9]+$/) {
           print " ---  Introduce the number of atoms in the slab --- \n ";
              $slabatoms=<STDIN>; chomp($slabatoms);
    }; }; 
  }elsif ($numargs eq 1 ) {
    if ($ARGV[0] =~ /^[0-9]+$/) { $slabatoms=$ARGV[0]; 
      if (-e "CONTCAR") { $filename="CONTCAR";
       }else{ print " ---  Introduce the file name --- \n ";
              $filename=<STDIN>; chomp($filename);
        if (!-e $filename) {
          while (!-e $filename) {
           print " ---  Introduce an existing file --- \n ";
              $filename=<STDIN>; chomp($filename);
       }; }; };
     }elsif($ARGV[0] !~ /^[0-9]+$/) { 
      if (-e $ARGV[0]) { $filename=$ARGV[0]; };     
           print " ---  Introduce the number of atoms in the slab --- \n ";
              $slabatoms=<STDIN>; chomp($slabatoms);
      if ($slabatoms !~ /^[0-9]+$/) {
       while ($slabatoms !~ /^[0-9]+$/) {
           print " ---  Introduce the number of atoms in the slab --- \n ";
              $slabatoms=<STDIN>; chomp($slabatoms);
    }; }; }; 
  }elsif ($numargs gt 1 ) {
    if ($ARGV[0] !~ /^[0-9]+$/) {
      if (-e $ARGV[0]) { $filename=$ARGV[0]; $slabatoms=$ARGV[1];
       }elsif (!-e $filename) {
        while (!-e $filename) {
           print " ---  Introduce an existing file --- \n ";
              $filename=<STDIN>; chomp($filename);
       }; };
    }elsif($ARGV[0] =~ /^[0-9]+$/) { $slabatoms=$ARGV[0]; $filename=$ARGV[1]; };
  };
           print " ---  The loaded file is $filename --- \n ";
           print "---  The slab contains $slabatoms atoms --- \n ";
#----------------------------------------------------------------------------------------------------- multicell
  if ($ARGV[2]) { $multiple=$ARGV[2];
    print "---  The output is a $multiple x $multiple times the original cell --- \n";
   }else{
    print "---  How many replications? --- \n ";
       $multiple=<STDIN>; chomp($multiple);
     if ($multiple == 0) { $multiple=1; };      
    print "---  The output is a $multiple x $multiple times the original cell --- \n"; };
#----------------------------------------------------------------------------------------------------- read
open (IN,$filename);
   while (<IN>) { $file.= $_; };
close (IN);
     @file=split(/\n/,$file);
#----------------------------------------------------------------------------------------------------- head & param
     $header=$file[0];
   foreach $l (@file[1]) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; $param=@Nline[0]; @Nline=();
#----------------------------------------------------------------------------------------------------- vectors
     @line2=split(/\s+/,@file[2]);
     @line3=split(/\s+/,@file[3]);
     @line4=split(/\s+/,@file[4]);
   foreach $l (@line2) { if ($l) { push(@Nline,$l); }; }; @Xvec=@Nline; @Nline=();
   foreach $l (@line3) { if ($l) { push(@Nline,$l); }; }; @Yvec=@Nline; @Nline=();
   foreach $l (@line4) { if ($l) { push(@Nline,$l); }; }; @Zvec=@Nline; @Nline=();
#----------------------------------------------------------------------------------------------------- atoms
     @line5=split(/\s+/,@file[5]);
   foreach $l (@line5) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; @header=@Nline; @Nline=();
     @line6=split(/\s+/,@file[6]);
   foreach $l (@line6) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; @atoms=@Nline; @Nline=();
     my $natoms=0; ($natoms+=$_) for @atoms;
     my $nspecies=-1; ($nspecies++) for @atoms;
      $slabspecies=0; $tmpatoms=0;
   while ($tmpatoms < $slabatoms) { $tmpatoms=$tmpatoms+@atoms[$slabspecies]; $slabspecies++; };
    $slabspecies=$slabspecies-1;
#----------------------------------------------------------------------------------------------------- selective
     @selec=split(//,@file[7]);
   if ((@selec[0] eq 'S') or (@selec[0] eq 's')) { $n=1; }else{ $n=0; };
#----------------------------------------------------------------------------------------------------- surf
open OUT, ">$filename.surf";
      print OUT "@header\n$param\n";
      printf OUT "  %.15f  %.15f  %.15f\n",@Xvec[0],@Xvec[1],@Xvec[2];
      printf OUT "  %.15f  %.15f  %.15f\n",@Yvec[0],@Yvec[1],@Yvec[2];
      printf OUT "  %.15f  %.15f  %.15f\n",@Zvec[0],@Zvec[1],@Zvec[2];
      print OUT "  @header[0..$slabspecies]\n  @atoms[0..$slabspecies]\n";
   if ($n == 1) { print OUT "@file[7]\n"; }; 
      print OUT "@file[7+$n]\n";
   for ($i=1; $i<=$slabatoms; $i++) { print OUT "@file[7+$n+$i]\n"; };
close OUT;
#----------------------------------------------------------------------------------------------------- surf.x3
system ("\$HOME/software/multicell.pl $filename.surf $multiple >>kk-temp.$filename");
#----------------------------------------------------------------------------------------------------- read2
open (IN2,"$filename.surf.x$multiple");
   while (<IN2>) { $file2.= $_; };
close (IN2);
     @file2=split(/\n/,$file2);
     @line6=split(/\s+/,@file2[6]);
   foreach $l (@line6) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; @Tatoms=@Nline; @Nline=();
     my $nTatoms=0; ($nTatoms+=$_) for @Tatoms;
   for ($s=$slabspecies+1; $s<=$nspecies; $s++) { push(@Tatoms,$atoms[$s]); };
#----------------------------------------------------------------------------------------------------- print multiple
open OUT, ">$filename.x$multiple";
      print OUT "@header\n";
    for ($i=1; $i<=4; $i++) { print OUT "@file2[$i]\n"; };
      print OUT "  @header\n  @Tatoms\n";
   if ($n == 1) { print OUT "@file2[7]\n"; };
      print OUT "@file2[7+$n]\n";
   for ($i=1; $i<=$nTatoms; $i++) { print OUT "@file2[7+$n+$i]\n"; };
close OUT;
#----------------------------------------------------------------------------------------------------- cartesian
system ("\$HOME/software/cartesian.pl $filename >>kk-temp.$filename");
    $molecatoms=$natoms-$slabatoms;
system ("tail -$molecatoms $filename.cart >> $filename.x$multiple");
system ("\$HOME/software/vasp2pov.sh $filename.x$multiple 0; mv $filename.x$multiple.x1.pov tmppov" );
system ("sed \"s/declare ncells=1/declare ncells=$multiple/g\" tmppov >> $filename.x$multiple.x1.pov ;");
system ("rm $filename.cart $filename.surf.x$multiple $filename.surf kk-temp.$filename tmppov");

