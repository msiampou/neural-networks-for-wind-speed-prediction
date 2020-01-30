cd ../app
make
cd build/
./cluster -i ../../datasets/new_representations.csv -c ../../config/cluster.conf -o ../../results/new_representations_12_112.csv --complete --init random --assign lloyd --update mean
cd ..
make clean
