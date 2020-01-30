import argparse
import util
import time
import os

import numpy as np
import pandas as pd

from keras.models import load_model
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error

MODEL_PATH = '../models/WindDenseNN1.h5'
OUTPUT_PATH = '../results/new_representations.csv'

def make_args_parser():
    # create an ArgumentParser object
    parser = argparse.ArgumentParser(usage='%(prog)s [-h] [-i <input file>]')
    # fill parser with information about program arguments
    parser.add_argument('-i', '--input', action='store',
                        default=None,  metavar='',
                        help='relative path to the input file')
    # return an ArgumentParser object
    return parser.parse_args()

def print_args(args):
    print("Running with the following configuration")
    # get the __dict__ attribute of args using vars() function
    args_map = vars(args)
    for key in args_map:
        print('\t', key, '-->', args_map[key])
    # print an empty line
    print()

def main():
    # parse and print arguments
    args = make_args_parser()
    print_args(args)
    # load input file as well as the timestamps
    X, ts = util.load_file(args.input)
    # load and compile model
    model = util.load_and_compile_model(MODEL_PATH)
    # print model summary
    print("------------------------- Model Summary -------------------------")
    model.summary()
    # Get new model having as output the first hidden of our pretrained model
    new_model = util.get_intermediate_layer_model(model, 'dense_1')
    # print model summary
    print("--------------- Intermediate Layer Model Summary ----------------")
    new_model.summary()
    # predict new model output
    y_pred = new_model.predict(X, batch_size=32)
    # Open file an write from scratch
    with open(OUTPUT_PATH, 'w') as file:
        df = pd.concat([ts, pd.DataFrame(y_pred)], axis=1)
        df.to_csv(file, index=False, header=False, sep='\t', encoding='utf-8')

if __name__ == '__main__':
    main()
