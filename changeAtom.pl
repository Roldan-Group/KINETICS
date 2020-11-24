#;-*- Perl -*-
#
#
##################  Alberto Roldan 9-2013  ###############

eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;

#---------------------------------------------------------------------------------------------------- default
   $changebyNatom=8;
#----------------------------------------------------------------------------------------------------
   $filename=$ARGV[0];
   $element=$ARGV[1];
   $swapPosition=$ARGV[2];
   $numargs=$#ARGV+1;
  if ($numargs lt 1 ) {
   if (-e "CONTCAR") { $filename="CONTCAR"; 
   }else{
    print "\n ---  Provide a valid file and Element symbol ---\n "; exit 0; };
   if (-e $filename) { print "\n ---  Provide a valid file and Element symbol ---\n "; exit 0; }; };
   if (!$element) { $element="X"; };
   if (!$swapPosition) { $swapPosition=3; };

  print "\n ---  Loaded file is $filename  ---\n ";
  print "---  Swaping Element is $element  ---\n ";
  print "---  The $element position is $swapPosition in $filename.changed file ---\n ";

#----------------------------------------------------------------------------------------------------
 open IN, $filename;
   while (<IN>) { $file.= $_; };
 close (IN);
    @file=split(/\n/,$file);
#---------------------------------------------------------------------------------------------------- atoms
    @line=split(/\s+/,@file[5]);
    $n=1;
      foreach $l (@line) { if ($l) {
        if (($n < $swapPosition) or ($n > $swapPosition)) {
           push(@Nline,$l);
         }else{ push(@Nline,$element); push(@Nline,$l);  }; $n++; }; }; @atoms=@Nline; @Nline=();
#---------------------------------------------------------------------------------------------------- Natoms
    @line=split(/\s+/,@file[6]);
    $n=0; $changed="no"; $Tnatoms=0;
      foreach $l (@line) { if ($l) { $Tnatoms=$Tnatoms+$l;
        if (($changebyNatom <= $Tnatoms) and ($changed eq "no")) { $l=$l-1; $changed="yes"; };
        if (($n < $swapPosition-1) or ($n > $swapPosition-1)) { $Natoms{$atoms[$n]}=$l;
         }elsif ($n = $swapPosition-1) { $Natoms{$atoms[$n]}=1; $Natoms{$atoms[$n+1]}=$l; $n++; };
            $n++;
       }; }; 
#---------------------------------------------------------------------------------------------------- Selective
    @line=split(/\s+/,@file[7]);
      foreach $l (@line) { if ($l) { push(@Nline,$l); }; }; @line=@Nline; @Nline=();
        @selective=split(//,@line[0]);
      if ((@selective[0] eq "S") or (@selective[0] eq "s")) { $s=1; }else{ $s=0; };
#---------------------------------------------------------------------------------------------------- Newatom
    @NewatomPosition=@file[$changebyNatom+7+$s];
    delete @file[$changebyNatom+7+$s];
#---------------------------------------------------------------------------------------------------- print
 open OUT, ">>$filename.changed";
     print OUT "@atoms\n";
    for ($i=1; $i<=4; $i++) { print OUT "@file[$i]\n"; };
     print OUT "   @atoms\n   ";
    foreach $n (@atoms) { print OUT "$Natoms{$n} "; }; print OUT "\n";
     $m0=7+$s+$changebyNatom-1;
    for ($i=7; $i<=$m0; $i++) { print OUT "@file[$i]\n"; };
      $sumNatoms=0;
    for ($n=0; $n<=$swapPosition-2; $n++) { $sumNatoms=$sumNatoms+$Natoms{$atoms[$n]}; };
     $m1=7+$s+$sumNatoms+1;
    for ($i=$m0+2; $i<=$m1; $i++) { print OUT "@file[$i]\n"; };
     print OUT "@NewatomPosition\n";
      $j=$m1+1;
    while (@file[$j]) { print OUT "@file[$j]\n"; $j++; };
 close OUT;
