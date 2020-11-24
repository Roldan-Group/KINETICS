#;-*- Perl -*-

#--------ALBERTO ROLDAN-----------8-11-2010---
#

## system ('vi '.$ARGV[0]);
eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;

#---------------------------------------------------------
#$vasp=5;   #version 4.*
$vasp=6;   #version 5.*
#---------------------------------------------------------


open IN, $ARGV[0];
  while (<IN>) { $raw.=$_;};
close IN;
    @raw=split(/\n/ , $raw);
    @array=split(/\s+/,$raw[$vasp]); my $natoms = 0;  ($natoms+=$_) for @array ;
open OUT, ">$ARGV[0].tmp" ;
    for ($i=0;$i<=$vasp;$i++){ print OUT "@raw[$i]\n";  };
@selective=split(//,$raw[$vasp+1]);
    if ($selective[0] eq 'S' or $selective[0] eq 's'){ $k=$vasp+2; } else { $k=$vasp+1; };
for ($j=$k;$j<=$natoms+$k;$j++){ print OUT " @raw[$j] \n"; };
 close OUT;

   @elements=('H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe');
   @nelem=split(/\s+/,@elements);  my $nelements = 0;  ($nelements+=$_) for @nelem ;
#   @line=split(/\s+/,$raw[$vasp-1]);  #version 4.*
   @line=split(/\s+/,$raw[$vasp-1]); #version 5.*
     $j=0;
   foreach $el (@line) { if ($el) { @lineatoms[$j]=$el; $j++; }; };
   $ka=1;
foreach $katom (@lineatoms) { for ($e=0; $e<=$nelements; $e++){ if (@elements[$e] eq $katom) { @ne[$ka]=$e+1; @xsf[$ka]=("-$ka @ne[$ka]"); $ka++; }; }; };

print "@xsf\n";

 system ("/home/alberto/.xsfconvert/bin/v2xsf $ARGV[0].tmp @xsf"); 
 system ("xcrysden --xsf $ARGV[0].tmp.xsf.gz") ;
 system ("rm $ARGV[0].tmp $ARGV[0].tmp.*") ;

