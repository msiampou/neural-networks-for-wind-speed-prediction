import argparse
import util
import os

import numpy as np
import pandas as pd

MODEL_PATH = '../models/WindDenseNN1.h5'
ACTUAL_PATH = '../datasets/actual.csv'
PREDICTED_PATH = '../results/predicted.csv'

def make_args_parser():
    # Create an ArgumentParser object
    parser = argparse.ArgumentParser(usage='%(prog)s [-h] [-i <input file>]')
    # Fill parser with information about program arguments
    parser.add_argument('-i', '--input', action='store',
                        default=None,  metavar='',
                        help='relative path to the input file')
    # Return an ArgumentParser object
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
    # Parse and print arguments
    args = make_args_parser()
    print_args(args)
    # Load input file as well as the timestamps
    X, ts = util.load_file(args.input)
    # Load and compile model
    model = util.load_and_compile_model(MODEL_PATH)
    # Print model summary
    print("------------------------- Model Summary -------------------------")
    model.summary()
    # Predict model output
    y_pred = model.predict(X, batch_size=32)
    # Load actual values
    y_true, _ = util.load_file(ACTUAL_PATH)
    # Compute MAE, MAPE and MSE
    mae = util.compute_mae(X, y_true, model)
    mape = util.compute_mape(y_true, y_pred, model)
    mse = util.compute_mse(X, y_true, model)
    # Check if file already exists. If so delete it
    if os.path.exists(PREDICTED_PATH):
        os.remove(PREDICTED_PATH)
    # Open file an write from scratch
    with open(PREDICTED_PATH, 'a') as file:
        file.write('MAE: {:.4f}\tMAPE: {:.4f}\tMSE: {:.4f}\n'.format(mae,mape,mse))
        df = pd.concat([ts, pd.DataFrame(y_pred)], axis=1)
        df.to_csv(file, index=False, header=False, sep='\t', encoding='utf-8')

if __name__ == '__main__':
    main()
