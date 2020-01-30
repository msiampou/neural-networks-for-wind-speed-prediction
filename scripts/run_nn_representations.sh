cd ../app
make
cd build/
./cluster -i ../../datasets/nn_representations.csv -c ../../config/cluster.conf -o ../../results/nn_representations_4_112.csv --complete --init random --assign lloyd --update mean
cd ..
make clean
