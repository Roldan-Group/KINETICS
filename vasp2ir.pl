#;-*- Perl -*-
#
########   Alberto 21-10-2011
eval '(exit $?0)' && eval 'exec perl -S $0 ${1+"$@"}' && eval 'exec perl -S $0 $argv:q' if 0;

#-------------------------------------------------------------------------------------------- file
   if ($ARGV[0]) { $file=$ARGV[0]; }else{ $file=OUTCAR; };
   print "\n ---  The loaded file is $file  ---\n "; 

#-------------------------------------------------------------------------------------------- ibrion

   system("grep IBRION $file > tmpibrion");
   open (IN0, tmpibrion);
      while (<IN0>) { $tmp.= $_ ; };
   close (IN0);
        @target=split(/\s+/,$tmp);  $ibrion=@target[3];
   system("rm -f tmpibrion");

print "\t---  IBRION is $ibrion  ---\n";
#-------------------------------------------------------------------------------------------- idipol
   system("grep IDIPOL $file > tmpidipol");
   open (IN1, tmpidipol);
      while (<IN1>) { if ($.>0) { @target=split(/\s+/,$_);  $idipol=@target[3]; $kk++; }; };
   close (IN1);
   system("rm -f tmpidipol");

#print "...idipol=$idipol...\n";

#-------------------------------------------------------------------------------------------- vasp2ir
   if (($ibrion eq 5) or ($ibrion eq 6)) {
print " ---  Calling vasp2ir_DFT.pl  ---\n\n";
       system("/home/alberto/Software/VASP2IR_SOURCE/vasp2ir_DFT.pl $file");
	system("python3 /home/alberto/Software/VASP2IR_SOURCE/vasp_raman.py");
   }elsif (($ibrion eq 7) or ($ibrion eq 8)) {
print " ---  Calling vasp2ir_Perturbation.sh  ---\n\n";
       system("/home/alberto/software/VASP2IR_SOURCE/vasp2ir_Perturbation.sh $file");
      if ($idipol eq 4) {
print "\n ---  Calling vasp2ir_Perturbation_3D.sh  ---\n\n";
       system("/home/alberto/Software/VASP2IR_SOURCE/vasp2ir_Perturbation_3D.sh $file");
print "\n ---  Applaying Gaussian Smearing  ---\n";
        system("/home/alberto/Software/VASP2IR_SOURCE/Smearing.pl ./intensities/results/exact.reform.res.txt");
        system("mv ./IRSPECTRA ./intensities/");
      }else{
print "\n ---  Applaying Gaussian Smearing  ---\n";
        system("/home/alberto/Software/VASP2IR_SOURCE/Smearing.pl ./intensities/results/exact.res.txt");
        system("mv ./IRSPECTRA ./intensities/");
      };
   }else{ print " ---  It is not a frequencies $file  --- \n"; };


