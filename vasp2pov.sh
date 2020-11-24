#!/bin/bash



#--------------------------------------------------------------------------------------------------------------------------------------------------------------
        MgO=4.195157

	lattice=$MgO
#--------------------------------------------------------------------------------------------------------------------------------------------------------------

#    number=^[+-]?[0-9]+([.][0-9]+)?$

    if [ $1 ]; then
            if [[ "$1" != *([+-])*([0-9])*(.)*([0-9]) ]]; then
                	file1=$1
		else
			n=$1
	    fi
    fi

    if [ $2 ]; then
            if [[ "$2" != *([+-])*([0-9])*(.)*([0-9]) ]] && [ $file1 ]; then
			file2=$2
                elif [[ "$2" = *([+-])*([0-9])*(.)*([0-9]) ]] && [ ! $n ]; then
                	n=$2
            fi
    fi

    if [ $3 ]; then
            if [[ "$3" != *([+-])*([0-9])*(.)*([0-9]) ]] && [ ! $file2 ]; then
                	file2=$3
                elif [[ "$3" = *([+-])*([0-9])*(.)*([0-9]) ]] && [ ! $n ]; then
                	n=$3
            fi	
    fi    

        if [ $4 ]; then
            if [[ "$4" != *([+-])*([0-9])*(.)*([0-9]) ]] && [ ! $file2 ]; then
                	file2=$4
		elif [[ "$4" = *([+-])*([0-9])*(.)*([0-9]) ]] && [ ! $n2 ]; then
                	n2=$4
            fi
    fi


    if [ $file1 ] && [ $file2 ]; then
            if [ $n ] && [ $n2 ]; then
    		    ~/software/VASP2POV_SOURCE/interpolate.pl $file1 $file2 $n $n2
		    cd $file1$file2
		    list=$(ls)
#		    ~/software/VASP2POV_SOURCE/executer.sh ${list[@]}
            elif [ $n ] && [ ! $n2 ]; then
    		    ~/software/VASP2POV_SOURCE/interpolate.pl $file1 $file2 $n $lattice
		    cd $file1$file2
                    list=$(ls)
                    ~/software/VASP2POV_SOURCE/executer.sh ${list[@]}
	    elif [ $n ]; then
		   ~/software/VASP2POV_SOURCE/interpolate.pl $file1 $file2 10 $lattice
		   cd $file1$file2
                   list=$(ls)
                   ~/software/VASP2POV_SOURCE/executer.sh ${list[@]}
            fi
    elif [ $file1 ] && [ ! $file2 ]; then
	    if [ $n ]; then
		    ~/software/VASP2POV_SOURCE/vasp2pov.pl $file1 $n
	    else
		    ~/software/VASP2POV_SOURCE/vasp2pov.pl $file1
	    fi
    fi

