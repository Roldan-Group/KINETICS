#/bin/sh



  dir=$(pwd)

	file=$(find . -name "x*" -type d)
	for i in $file;
       	do 
		cd $i
		python3 /home/alberto/Software/OTHER/NeuralNetwork/GetData.py $dir
		if [ -f "data.dat" ]; then
			cat data.dat >> ../Data.tmp
			rm data.dat
			mv labels.txt ../.
        python3 ~/Software/OTHER/NeuralNetwork/CatStructure.py

        cat Eadh.dat >> ../Eadh.dat; rm Eadh.dat
#        cat cc_Ecoh.dat >> ../cc_Ecoh.dat; rm cc_Ecoh.dat
        cat Eb.dat >> ../Eb.dat; rm Eb.dat
        cat Etotal.dat >> ../Etotal.dat; rm Etotal.dat
#        cat Esurf.dat >> ../Esurf.dat; rm Esurf.dat

		fi
		cd ..
       	done
	cat labels.txt Data.tmp >> Data.dat
	rm labels.txt Data.tmp
  python3 /home/alberto/Software/OTHER/NeuralNetwork/Pre_Data_Collection/Param.py Data.dat
  cat Param.dat >> ../Param.dat
  python3 /home/alberto/Software/OTHER/NeuralNetwork/Pre_Data_Collection/get_minima.py Data.dat
  cat min_Param.dat >> ../min_Param.dat
  python3 /home/alberto/Software/OTHER/NeuralNetwork/Pre_Data_Collection/Interpolate_MP_xy.py Data.dat
  cat xy_fit.dat >> ../xy_fit.dat
  python3 /home/alberto/Software/OTHER/NeuralNetwork/Pre_Data_Collection/Interpolate_MP_xy_2.py Data.dat
  cat 2xy_fit.dat >> ../2xy_fit.dat
	python3 /home/alberto/Software/OTHER/NeuralNetwork/Plot.py Data.dat
#	python3 /home/alberto/plot_e_adh.py Eadh.dat

# to normalise -- implemented before Ebinding
#	python3 /home/alberto/software/OTHER/NeuralNetwork/ENormalisation.py Data.dat
#	python3 /home/alberto/software/OTHER/NeuralNetwork/Plot.py NData.dat


