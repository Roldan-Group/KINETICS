#/bin/sh


    if [ $1 ]; then
	    file=$1
    else
	    file=Data.dat
    fi

	python3 /home/alberto/software/OTHER/NeuralNetwork/Training.py $file

