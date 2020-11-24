#;-*- Perl -*-
#
########   Alberto 04-2011

eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;

$numargs=$#ARGV+1;
$file=$ARGV[0];
if ($numargs lt 1 ){ if (-s "OUTCAR"){$file="OUTCAR"; }else{print "\n ---  Please, provide a valid OUTCAR file after the executable --- \n\n"; exit 0; };};


#============================================================================================================================== MOLCAR FILE
open OUT, ">>MOLCAR";
   print OUT " [Molden Format]\n [GEOCONV]\n energy\n";
close OUT;
#-----------------------------------------------------------------------------------------------------------energies
   system("grep -A 4 '  FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)' $file >> tmp");
   system('grep "energy(sigma->0" tmp >> tmp2');
 open IN, "tmp2";
    while (<IN>) {$ein.= $_ ;};
 close IN ;
      @Ein=split(/\n/,$ein);
      $i=0;
open OUT, ">>MOLCAR";
   while (@Ein[$i]) { @E=split(/\s+/,@Ein[$i]); print OUT " @E[7]\n"; $i++; };
#-----------------------------------------------------------------------------------------------------------geometries
   print OUT " [GEOMETRIES] XYZ\n";
close OUT;
   system("grep 'POSCAR =' $file >> tmp3");
 open IN, "tmp3";
    while (<IN>) {$pos.= $_ ;};
 close IN ;
      @POS=split(/\=/,$pos);
      @atoms=split(/\s+/,$POS[1]);
   system ("grep 'ions per type' $file >> tmp4");
 open IN, "tmp4";
    while (<IN>) {$ions.= $_ ;};
 close IN ;
      @IONS=split(/\=/,$ions);
      @nions=split(/\s+/,$IONS[1]);
      $Tnions=0; ($Tnions+=$_) for @nions;
      $tmp=$Tnions+1;
   system("grep -A$tmp 'POSITION' $file | cut -b1-45 >>tmp5");
 open IN, "tmp5";
    while (<IN>) { $in.=$_ ;}
 close IN;
       @input=split(/\n/,$in);
       $k=0;
   while (@input[$k]) { @ini=split(//,@input[$k]);
open OUT, ">>MOLCAR";
      if (($ini[1] ne 'P') and ($ini[1] ne '-')) {
 print OUT " $Tnions\n\n";              #______________________________________________________________________  total number of atoms
               $j=1; $l=1;
           for ($i=1; $i<=$Tnions;) {
                 if ($l <= @nions[$j]) {
 print OUT " @atoms[$j] @input[$k]\n";              #____________________________________________________  atom type and position
                      $l++; $k++; $i++;
                 } elsif ($l > @nions[$j]) { $j++; $l=1; };
          }; #--> for
      }else{ $k++; };
close OUT;
   }; #--> while
#============================================================================================================================== MOLCAR FILE END
   system("~/software/MOLDEN/molden5.2/gmolden MOLCAR");
   system("rm tmp tmp2 tmp3 tmp4 tmp5 MOLCAR");
##############################################################################################################
