#;-*- Perl -*-
#
#
#   multicell CONTCAR Times_of_replicated_cell
#
########   Alberto 07-2012

eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
use Math::Trig;

  $numargs=$#ARGV+1; $filename=$ARGV[0]; $ntimes=$ARGV[1]; 
     if ($numargs lt 1 ) { $filename="CONTCAR"; };
  print "\n ---  The loaded file is $filename ---\n"; 
     if ($numargs le 1 ) { $ntimes=3; };
  print " ---  The output is a $ntimes x $ntimes times the original cell ---\n"; 

#--------------------------------------------------------------------------------------------------------------------------------- odd or even
 if ($ntimes%2 == 1) { $odd=1; }else{ $odd=0; };
#--------------------------------------------------------------------------------------------------------------------------------- times
  if ($odd == 1) { $Ptimes=($ntimes-1)/2; $Ntimes=$Ptimes;
  }elsif ($odd == 0) { if ($ntimes > 2) { $Ptimes=$ntimes/2; $Ntimes=$Ptimes-1; 
                       }elsif ($ntimes = 2) { $Ptimes=$ntimes/2; $Ntimes=0; };
  };
#----------------------------------------------------------------------------------------------------- cartesian
open (IN,$filename);
   while (<IN>) { $file0.= $_; };
close (IN);
     @file0=split(/\n/,$file0);
     @selec=split(//,$file0[7]);
    if (($selec[0] eq "S") or (@selec eq "s")) { $n=0; }else{ $n=1; };
     @cartesian=split(//,$file0[8-$n]);
    if (($cartesian[0] eq "D") or (@cartesian eq "d")) { system("cartesian.pl $filename");
      $filename2="$filename.cart"; 
     }elsif (($cartesian[0] eq "C") or (@cartesian eq "c")) { $filename2=$filename; };
#----------------------------------------------------------------------------------------------------- read
open (IN,$filename2);
   while (<IN>) { $file.= $_; };
close (IN);
     @file=split(/\n/,$file);
     $header=$file[0];
   foreach $l (@file[1]) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; $param=@Nline[0]; @Nline=();
     @line2=split(/\s+/,@file[2]);
     @line3=split(/\s+/,@file[3]);
     @line4=split(/\s+/,@file[4]);
   foreach $l (@line2) { if ($l) { push(@Nline,$l); }; }; @Xvec=@Nline; @Nline=();
   foreach $l (@line3) { if ($l) { push(@Nline,$l); }; }; @Yvec=@Nline; @Nline=();
   foreach $l (@line4) { if ($l) { push(@Nline,$l); }; }; @Zvec=@Nline; @Nline=();
#---------------------------------------------------------------------------------------------------------------------- ATOMS
     @line5=split(/\s+/,@file[5]);
   foreach $l (@line5) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; @header=@Nline; @Nline=();
     @line6=split(/\s+/,@file[6]);
   foreach $l (@line6) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; @atoms=@Nline; @Nline=();
     my $natoms=0; ($natoms+=$_) for @atoms;
     my $nspecies=-1; ($species++ ) for @atoms;
#---------------------------------------------------------------------------------------------------------------------- Selective
     @selec=split(//,@file[7]);
   if ((@selec[0] eq 'S') or (@selec[0] eq 's')) { $n=1; }else{ $n=0; };
#---------------------------------------------------------------------------------------------------------------------- PRINT
  open OUT, ">>$filename.x$ntimes";
    print OUT "@header\n$param\n";
    printf OUT "  %.15f  %.15f  %.15f\n",@Xvec[0]*$ntimes,@Xvec[1]*$ntimes,@Xvec[2]*$ntimes;
    printf OUT "  %.15f  %.15f  %.15f\n",@Yvec[0]*$ntimes,@Yvec[1]*$ntimes,@Yvec[2]*$ntimes;
    printf OUT "  %.15f  %.15f  %.15f\n",@Zvec[0]*$ntimes,@Zvec[1]*$ntimes,@Zvec[2];
    print OUT "@header\n";
  foreach $a (@atoms) { printf OUT "  %d",$a*($ntimes**2); }; print OUT "\n";
  if ($n eq 1) { print OUT "Selective Dynamics\nCartesian\n";
   }else{ print OUT "Cartesian\n"; };
#---------------------------------------------------------------------------------------------------------------------- CARTESIAN
#------------------------------------------------------------------------------------------------------------------------------- cell replication
    for ($i=1; $i<=$natoms; $i++) {
       @lines=(); @lines=split(/\s+/,$file[$i+7+$n]);
     foreach $l (@lines) { if ($l) { push(@Nline,$l); }; }; @pos=@Nline; @Nline=();
      for ($xn=0; $xn<=$Ptimes; $xn++) {
        for ($yn=0; $yn<=$Ptimes; $yn++) {
            $newpos='';
             @newpos[0]=$pos[0]+@Xvec[0]*$xn+@Yvec[0]*$yn;
             @newpos[1]=$pos[1]+@Xvec[1]*$xn+@Yvec[1]*$yn;
         if ($n eq 1) {
            printf OUT "%2.15f  %2.15f  %2.15f @pos[3] @pos[4] @pos[5]\n",@newpos[0],@newpos[1],@pos[2];
         }else{ 
            printf OUT "%2.15f  %2.15f  %2.15f\n",@newpos[0],@newpos[1],@pos[2]; };
       }; };
      for ($xn=1; $xn<=$Ntimes; $xn++) {
         for ($yn=0; $yn<=$Ptimes; $yn++) {
            $newpos='';
             @newpos[0]=$pos[0]-@Xvec[0]*$xn+@Yvec[0]*$yn;
             @newpos[1]=$pos[1]-@Xvec[1]*$xn+@Yvec[1]*$yn;
         if ($n eq 1) {
            printf OUT "%2.15f  %2.15f  %2.15f @pos[3] @pos[4] @pos[5]\n",@newpos[0],@newpos[1],@pos[2];
         }else{
            printf OUT "%2.15f  %2.15f  %2.15f\n",@newpos[0],@newpos[1],@pos[2]; };
       }; };
      for ($xn=0; $xn<=$Ptimes; $xn++) {
         for ($yn=1; $yn<=$Ntimes; $yn++) {
            @newpos='';
             @newpos[0]=$pos[0]+@Xvec[0]*$xn-@Yvec[0]*$yn;
             @newpos[1]=$pos[1]+@Xvec[1]*$xn-@Yvec[1]*$yn;
         if ($n eq 1) {
            printf OUT "%2.15f  %2.15f  %2.15f @pos[3] @pos[4] @pos[5]\n",@newpos[0],@newpos[1],@pos[2];
         }else{
            printf OUT "%2.15f  %2.15f  %2.15f\n",@newpos[0],@newpos[1],@pos[2]; };
       }; };
      for ($xn=1; $xn<=$Ntimes; $xn++) {
         for ($yn=1; $yn<=$Ntimes; $yn++) {
            @newpos='';
             @newpos[0]=$pos[0]-@Xvec[0]*$xn-@Yvec[0]*$yn;
             @newpos[1]=$pos[1]-@Xvec[1]*$xn-@Yvec[1]*$yn;
         if ($n eq 1) {
            printf OUT "%2.15f  %2.15f  %2.15f @pos[3] @pos[4] @pos[5]\n",@newpos[0],@newpos[1],@pos[2];
         }else{
            printf OUT "%2.15f  %2.15f  %2.15f\n",@newpos[0],@newpos[1],@pos[2]; };
       }; }; }; #--> for natoms
close OUT;

if ((@cartesian[0] eq 'D') or ($cartesian[0] eq 'd')) { system("rm $filename2");  };


