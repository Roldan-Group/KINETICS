sigma ()
{
	echo -ne "$1\t"
	grep -i sigma $1 | tail -1 | cut -b68-78

}


#-------------------------------------------------------------------------------- running jobs
    if [ $@ ]; then
        a=$@
    else
        a=$(find . -name "*OUTCA*" -print)
   fi 

for i in $a; do
	sigma $i
done
