for a in 1 2 3 4 5; do mkdir $a; done
cd Scripts
python datapreprocessor.py > ../metrics.txt
for a in 1 2 3 4 5; 
	do echo >> ../metrics.txt;
	cp training.rb validation.rb ../$a;
	cd ../$a;
	pwd >> ../metrics.txt;
	mv train.csv train.csv.bak;
	cut -d ',' -f 2- train.csv.bak > train.csv;
	mv val.csv val.csv.bak; 
	cut -d ',' -f 2- val.csv.bak > val.csv; 
	ruby training.rb >> ../metrics.txt; 
	ruby validation.rb >> ../metrics.txt;
	echo >> ../metrics.txt;
done
cd ..
