#;-*- Perl -*-
#
#
########   Alberto 11-2016
#		   05-2017
#		   11-2018
#
## input --> POSCAR/CONTCAR N Element  ... (N= slab atoms)
#
# >>> CAUTION <<< it will fail to calculate neighbours if the entire cluster is not in the cell
#
eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;
 use Math::Trig;
 use Term::ANSIColor qw( colored );
#---------------------------------------------------------------------------
  $numargs=$#ARGV+1;   
     if ($ARGV[0]) { $filename=$ARGV[0]; };
  print "\n ---  The loaded file is $filename ---\n"; 
     if ($ARGV[1]) {$slab0=$ARGV[1]; } else { $slab0=128; };
  print " ---  The atoms in the slab is $slab0 ---\n";
    if ($ARGV[2]) {$element=$ARGV[2]; } else { $element="Au"; };
  print " ---  The searched atoms are $element ---\n\n";


###############################################################################################################

 $threshold=&diameter($element)*2;   #------------> re-defines the margin for interface atoms
 $scale=0.5; 	    #------------> defines the range for surface top atoms

###############################################################################################################

 #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> FIX PATHWAY FOR asp2surf.pl <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   system ("NP_vasp2surf.pl $filename $slab0 3; vasp2xyz.pl $filename.x3");
    $filename="$filename.x3.xyz"; $extraSatoms=$slab0*8; $slab=$slab0*9;
open (IN,$filename);
   while (<IN>) { $file.= $_; };
close (IN);
     @file=split(/\n/,$file); @line=split(/\s+/,@file[0]);
     foreach $l (@line) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; $Tatoms=@Nline[0]; @Nline=();  #----- Total atoms
#----------------------------------------------------------------------------------------------------------------
       $z='', $preZ=''; $origin=0.0; 
   for ($i=1; $i<=$slab; $i++) { @line=split(/\s+/,@file[$i+1]);
     foreach $l (@line) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; @line=@Nline; @Nline=();
     if ($line[3] > $origin) { $origin=$line[3];  }; };  #for --------------------------------------------------->  looking for the highest Z of the surface
	$counter=1;
   for ($i=1; $i<=$slab; $i++) { @line=split(/\s+/,@file[$i+1]); 
     foreach $l (@line) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; @line=@Nline; @Nline=();
     if ($line[3] >= $origin-$scale) { push(@atomsINsurface,$i);
	     push(@surfaceAtoms,$line[3]); push(@printSurface,"$line[0]$counter z=$line[3]"); }; $counter++; };
#----------------------------------------------------------------------------------------------------------------
     $SumSurf=0;
       for ($n=0; $n<=$#surfaceAtoms; $n++) { $SumSurf=$SumSurf+$surfaceAtoms[$n]; }; $AverageSurf=$SumSurf/($#surfaceAtoms+1); # ---> average Z of the surface

##################################################################################
# surfaceAtoms = height of the atoms from the surface at the uppermost layer
# AverageSurf = average height of the uppermost surface layer
# atomsINsurface = line number in filename of surface atoms
#
# clusterAtoms = line number in filename of cluster atoms
# ZclusterAtoms = Z coordinate of the atoms from the cluster
# ZAveragecluster = average height of the cluster
# ZAveragedistance = main distance from the cluster to the surface
# Zdistance = distance of the cluster atoms to the surface
################################################################i###################

#----------------------------------------------------------------------------------------------------------------- cluster atoms
   for ($i=$slab+1; $i<=$Tatoms; $i++) { @line=split(/\s+/,@file[$i+1]);
        foreach $l (@line) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; @line=@Nline; @Nline=();
             push(@clusterAtoms,$i); push(@ZclusterAtoms,$line[3]); };  #for
   for ($n=0; $n<=$#clusterAtoms; $n++) { $SumRing+=$ZclusterAtoms[$n]; $Zdistance[$n]=$ZclusterAtoms[$n]-$AverageSurf;};
           $ZAveragecluster=$SumRing/($#ZclusterAtoms); $ZAveragedistance=$ZAveragecluster-$AverageSurf;
#----------------------------------------------------------------------------------------------------------- Neighbors
           @Asurfatomdist=(); @Aclusteratomdist=(); $dsum=0; @interfaceAtoms=();
       for ($n=0; $n<=$#clusterAtoms; $n++) { $a=$clusterAtoms[$n]-$extraSatoms; @line=split(/\s+/,@file[$clusterAtoms[$n]+1]);
            $dist=(); @surfatomdist=(); @clusteratomdist=();
          if (@line[3] <= $AverageSurf+$threshold) { #---------------------------------------------------------->  the interface
		  $dsum=$dsum+($line[3]-$AverageSurf); 
              for ($m=0; $m<=$#atomsINsurface; $m++) { $b=@atomsINsurface[$m];
		      @line2=split(/\s+/,@file[@atomsINsurface[$m]+1]);  
                   $dist=sqrt((@line2[1]-@line[1])**2+(@line2[2]-@line[2])**2+(@line2[3]-@line[3])**2);
 		   $bond=&diameter($line[0])+&diameter($line2[0])*2;   #--------------------------------------------------> to define neighbors
   		   if ($dist <= $bond) { push(@surfatomdist,"$dist,$a,$b"); push(@interfaceAtoms,$a); };};};
   	@tmp=(); @tmp2=(); @tmp=sort {$a <=> $b} @surfatomdist ; @surfatomdist=@tmp; $Adist=0;
   	foreach $sa (@surfatomdist) { @tmp=split(/,/,$sa); 
   		if (!@tmp2) { push(@tmp2,"$tmp[1],$tmp[2],$tmp[0]"); $saa=$tmp[0]+$scale; $Adist=$tmp[0];
		}else{ if ($tmp[0] <= $saa) {  push(@tmp2,"$tmp[1],$tmp[2],$tmp[0]"); $Adist=($Adist+$tmp[0])/2; };};};
		$interfaceClusterAtomDist[$n]="@tmp2"; 

           for ($m=0; $m<=$#clusterAtoms;  $m++) { if ($n != $m) { $b=$clusterAtoms[$m]; $dist=''; #------------> the cluster
		  @line2=split(/\s+/,@file[@clusterAtoms[$m]+1]);
		  foreach $l (@line2) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; @line2=@Nline; @Nline=();
                  $dist=sqrt((@line2[1]-@line[1])**2+(@line2[2]-@line[2])**2+(@line2[3]-@line[3])**2);
		  $bond=&diameter($line[0])+&diameter($line2[0]);
		  $b1=&diameter($line[0]);
             if ($dist <= $bond) { push(@clusteratomdist,"$a,$b,$dist"); };};}; 
             $ClusterAtomDist[$n]="@clusteratomdist";
      }; # for n
	foreach $l (@interfaceClusterAtomDist) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; @interfaceClusterAtomDist=@Nline; @Nline=();
        foreach $l (@ClusterAtomDist) { if (($l) or ($l eq 0)) { push(@Nline,$l); }; }; @ClusterAtomDist=@Nline; @Nline=();
	foreach $a (@interfaceAtoms) { $new=1; foreach $b (@Nline) { if ($a eq $b) { $new=0; };}; 
		if ($new eq 1) { push(@Nline,$a); };}; @interfaceAtoms=@Nline; @Nline=();

#---------------------------------------------------------------------------
  open OUT, ">>clusterANDsurfaceHeight.dat";
    printf OUT "\nThe atoms defining the surface with average height of %.3f Angstroms are:\n",$AverageSurf;
         foreach $pS (@printSurface) { print OUT "   $pS\n"; };
  close OUT;	 
#---------------------------------------------------------------------------
  open OUT, ">>Distances.dat";
    printf OUT "\nThe atoms defining the surface with average height of %.3f Angstroms are in clusterANDsurfaceHeight.dat\n",$AverageSurf;
    printf OUT "Number of $element atoms above the surface: %d\n",$#clusterAtoms+1;
    printf OUT "   The heights of the $element cluster atoms are in clusterANDsurfaceHeight.dat\n";
  close OUT;
  open OUT, ">>clusterANDsurfaceHeight.dat"; print OUT "\n";
       for ($n=0; $n<=$#clusterAtoms; $n++) {
    printf OUT "   heights of $element %d(%d) from the surface average: %.3f Angstroms\n",$clusterAtoms[$n]-$extraSatoms,$n+1,@Zdistance[$n]; }; #for n
  close OUT;
  open OUT, ">>Distances.dat";
    printf OUT "Total number of $element atoms at the interface is %d :\n\t",$#interfaceClusterAtomDist+1;
       foreach $i (@interfaceClusterAtomDist) { @Asum=split(/,/,$i); 
		if (@Asum) { printf OUT " %d(%d)",$Asum[0],$Asum[0]-$slab0; };}; print OUT "\n";
    print OUT "   Coordination details on ClusterCoordination.dat\n";
   close OUT;
   open OUT, ">> ClusterCoordination.dat";
    printf OUT "Total number of $element atoms at the interface is %d :",$#interfaceClusterAtomDist+1;
       foreach $i (@interfaceClusterAtomDist) { @Asum=split(/,/,$i);
                if (@Asum) { printf OUT " %d(%d)",$Asum[0],$Asum[0]-$slab0; };}; print OUT "\n";
    printf OUT "   The average distance between the surface and interface cluster atoms is %.3f Angstroms $k\n",($dsum)/($#interfaceClusterAtomDist+1);

	$neighInterface=0; $neighInterfaceCluster=0; $neighCluster=0;
       for ($n=0; $n<=$#clusterAtoms; $n++) { $a=$clusterAtoms[$n]-$extraSatoms; @A=();

	foreach $i (@interfaceClusterAtomDist) { @IFatoms=split(/\s+/,$i);
		if ($#IFatoms < 1) { @tmp=split(/,/,$i); 
			if ($a eq $tmp[0]) { push(@A,"$tmp[0],$tmp[1],$tmp[2]"); $neighInterface++; }; 
		}else{ foreach $IFa (@IFatoms) { @tmp=split(/,/,$IFa);
			if ($a eq $tmp[0]) { push(@A,"$tmp[0],$tmp[1],$tmp[2]"); $neighInterface++; };};};}; 
	foreach $i (@ClusterAtomDist) { @CLatoms=split(/\s+/,$i);
		if ($#CLatoms < 1) {  @tmp=split(/,/,$i); 
			if ($a eq $tmp[0]) { push(@A,"$tmp[0],$tmp[1],$tmp[2]"); $neighCluster++;
				foreach $j (@interfaceAtoms) { $b=$j+$extraSatoms; if ($b eq $tmp[1]) { $neighInterfaceCluster++; };};};		       
		}else{ foreach $CLa (@CLatoms) { @tmp=split(/,/,$CLa); 
			if ($a eq $tmp[0]) { push(@A,"$tmp[0],$tmp[1],$tmp[2]"); $neighCluster++;
				foreach $j (@interfaceAtoms) { $b=$j+$extraSatoms; if ($b eq $tmp[1]) { $neighInterfaceCluster++; };};};};};};

    printf OUT "The $element $a (%d) has %d neighbors\n",$n+1,$#A+1; #------------> surface + INcluster neigh
       foreach $i (@A) { @Asum=split(/,/,$i); 
			if ($Asum[1] <= $slab) {
    printf OUT "   surface neighbor %4s is at a distance of %.3f Angstroms\n",$Asum[1],$Asum[2]; 
			}else{
    printf OUT "   cluster neighbor %4s (%d) is at a distance of %.3f Angstroms\n",$Asum[1]-$extraSatoms,$Asum[1]-$slab,$Asum[2]; };};
       }; # for n
  close OUT;
  open OUT, ">>Distances.dat";
    printf OUT "The average distance between $element atoms at the interface (%d) and binding sites (%d) on the surface is %.3f Angstroms\n",$#interfaceClusterAtomDist+1,$neighInterface,$Adist;
    printf OUT "The average neighbors ($neighInterface) at the interface (only $element) with the support is %.3f (considering %d $element)\n", $neighInterface/($#interfaceClusterAtomDist+1), $#interfaceClusterAtomDist+1;
    printf OUT "The average neighbors ($neighInterfaceCluster) at the interface (only $element) within the cluster is %.3f (considering %d $element)\n",   $neighInterfaceCluster/($#interfaceClusterAtomDist+1), $#interfaceClusterAtomDist+1;
    printf OUT "The average neighbors ($neighCluster) of the entire cluster is %.3f (considering %d $element)\n",   $neighCluster/($#clusterAtoms+1), $#clusterAtoms+1;
	print OUT "\n";

  close OUT;

#---------------------------------------------------------------------------
    system("rm $filename; more Distances.dat");

#=================================================================================================================
   sub diameter {
        $El=''; ($El)=@_;

# radii depends on the periodic table row and atomic number + a arbitrari paramenter

@symbolradius{"Va"}=$row/10+0.2*$row-3E-03*(26-18);
@symbolradius{"X"}=$row/10+0.2*$row-3E-03*(26-18);

#---------------- ROW 1 --------------
@symbolradius{"H"}=0.3-3E-03; $longBond{"H"}=0.30;
#---------------- ROW 2 --------------
$row=2;
@symbolradius{"Li"}=$row/10+0.2*$row-3E-03*(3-2);
@symbolradius{"Be"}=$row/10+0.2*$row-3E-03*(4-2);
@symbolradius{"C"}=$row/10+0.2*$row-3E-03*(6-2); $longBond{"C"}=0.18;
@symbolradius{"N"}=$row/10+0.2*$row-3E-03*(7-2); $longBond{"N"}=0.15;
@symbolradius{"O"}=$row/10+0.2*$row-3E-03*(8-2)+0.65;
@symbolradius{"O2"}=$row/10+0.2*$row-3E-03*(8-2)+0.65;
@symbolradius{"Ow"}=$row/10+0.2*$row-3E-03*(8-2)+0.65;   # Oxygen from H2O
#---------------- ROW 3 --------------
$row=3;
@symbolradius{"Na"}=$row/10+0.2*$row-3E-03*(11-10); $longBond{"Na"}=0.10;
@symbolradius{"Mg"}=$row/10+0.2*$row-3E-03*(12-10)-0.8;
@symbolradius{"Si"}=$row/10+0.2*$row-3E-03*(14-10);
@symbolradius{"S"}=$row/10+0.2*$row-3E-03*(16-10); $longBond{"S"}=0.10;
@symbolradius{"Cl"}=$row/10+0.2*$row-3E-03*(17-10);
#---------------- ROW 4 --------------
$row=4;
@symbolradius{"K"}=$row/10+0.2*$row-3E-03*(19-18); $longBond{"K"}=0.25;
@symbolradius{"Ca"}=$row/10+0.2*$row-3E-03*(20-18); $longBond{"Ca"}=0.25;
@symbolradius{"Ti"}=$row/10+0.2*$row-3E-03*(22-18); $longBond{"Ti"}=0.25;
@symbolradius{"V"}=$row/10+0.2*$row-3E-03*(23-18); $longBond{"V"}=0.25;
@symbolradius{"Cr"}=$row/10+0.2*$row-3E-03*(24-18); $longBond{"Cr"}=0.10;
@symbolradius{"Fe"}=$row/10+0.2*$row-3E-03*(26-18); $longBond{"Fe"}=0.05;
@symbolradius{"Ni"}=$row/10+0.2*$row-3E-03*(28-18)+0.35;
@symbolradius{"Cu"}=$row/10+0.2*$row-3E-03*(29-18);
@symbolradius{"Zn"}=$row/10+0.2*$row-3E-03*(30-18);
@symbolradius{"Se"}=$row/10+0.2*$row-3E-03*(34-18);
#---------------- ROW 5 --------------
$row=5;
@symbolradius{"Rb"}=$row/10+0.2*$row-3E-03*(37-36);
@symbolradius{"Sr"}=$row/10+0.2*$row-3E-03*(38-36);
@symbolradius{"Y"}=$row/10+0.2*$row-3E-03*(39-36)-0.05;
@symbolradius{"Zr"}=$row/10+0.2*$row-3E-03*(40-36)-0.20;
@symbolradius{"Mo"}=$row/10+0.2*$row-3E-03*(42-36);
@symbolradius{"Ru"}=$row/10+0.2*$row-3E-03*(44-36);
@symbolradius{"Pd"}=$row/10+0.2*$row-3E-03*(46-36);
@symbolradius{"Ag"}=$row/10+0.2*$row-3E-03*(47-36);
@symbolradius{"Cd"}=$row/10+0.2*$row-3E-03*(48-36);
@symbolradius{"Sn"}=$row/10+0.2*$row-3E-03*(50-36)+0.30;
@symbolradius{"Sb"}=$row/10+0.2*$row-3E-03*(51-36)+0.30;
#---------------- ROW 6 --------------
$row=6;
@symbolradius{"Cs"}=$row/10+0.2*$row-3E-03*(55-54);
@symbolradius{"Ba"}=$row/10+0.2*$row-3E-03*(56-54);
@symbolradius{"Ce"}=$row/10+0.2*$row-3E-03*(58-54);
@symbolradius{"W"}=$row/10+0.2*$row-3E-03*(74-54)-0.30;
@symbolradius{"Pt"}=$row/10+0.2*$row-3E-03*(78-54)-0.30;
@symbolradius{"Au"}=$row/10+0.2*$row-3E-03*(79-54)-0.10;

#------------------------------------------
     $diam=@symbolradius{$El};
     return ($diam);
   }; #--> sub diameter
#=================================================================================================================




