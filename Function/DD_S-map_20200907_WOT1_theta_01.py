# Import modules
import numpy as np
import os
from sklearn.linear_model import ElasticNetCV
import pandas as pd
import time
from joblib import Parallel, delayed
import multiprocessing
from optparse import OptionParser
from local_functions_20200813 import make_weights, weight_data, select_frequent_dominate_genera, Elastic_net_fitting


# Define functions
def make_Smap_input(abund, target_otu, uncontinuous, uncontinous_loc):
    print('Process data for otu No. %s' % str(target_otu + 1))
    # Make input for the elastic_net
    abund = np.matrix(abund)
    block = np.append(abund[1:, target_otu], abund[0:-1, ], axis=1)
    # Delete uncontinuous data
    if uncontinuous == True:
        block = np.delete(block, uncontinous_loc, axis=0)
    # Scaling the input
    ## Each time series is normalized to have a mean of 0 and standard deviation of 1 before analysis with S-maps
    block = (block - np.average(block, axis=0)) / np.std(block, axis=0)

    return block


def predict_weight_possible_shift(abund, empirical_states, target_otu, target_abund, direction, theta,
                                  weighted_posshifts, target_states, DD_analyzed_states):
    print('Predict weighted possible shift at state No. %s' % str(target_abund + 1))
    # Make input for the elastic_net using both time series and empirical data
    ## Predict the possible shift
    ### Select the wanted genera according to the time series data
    wanted_empirical_state = pd.DataFrame(empirical_states, columns=abund.columns)
    wanted_empirical_state.fillna(0, inplace=True)
    ### Combine the time series and empirical data except state at time t and t+1
    if target_abund == abund.shape[0] - 1:
        abund == abund
    else:
        abund.drop([target_abund+1])
    wanted_abund = abund.append(wanted_empirical_state)

    ### Find possible shift
    target_state = abund.iloc[target_abund, :]
    if direction == 'increase':
        possible_state = np.where(wanted_abund.iloc[:, target_otu] > target_state.iloc[target_otu])[0]
    if direction == 'decrease':
        possible_state = np.where(wanted_abund.iloc[:, target_otu] < target_state.iloc[target_otu])[0]
    if direction == 'constant':
        possible_state = np.where(wanted_abund.iloc[:, target_otu] == target_state.iloc[target_otu])[0]
        possible_state = np.delete(possible_state, np.where(possible_state == target_abund))
    if len(possible_state > 0):
        possible_shifts = wanted_abund.iloc[possible_state, :]
        ### Calculate the weights
        E_dist = np.array(np.sqrt(np.sum((possible_shifts - target_state) ** 2, axis=1)))
        w = np.array(make_weights(E_dist, theta))
        ### Predict and collect weighted possible shift
        weighted_posshift = np.dot(w, possible_shifts) / sum(w)
        weighted_posshifts.append(weighted_posshift)
        target_states.append(target_state)
        DD_analyzed_states.append(target_abund)
    else:
        print('No possible shift in direction %s' % direction)

    return weighted_posshifts, target_states, DD_analyzed_states

def main(interest_otu):
    # Density dependent regularized S-map
    ## Make block data from the selected abundance data
    print('Process data for otu No. %s' % str(interest_otu + 1))
    for direction in ['increase', 'decrease']:
        # Start making input
        weighted_posshifts = []
        target_states = []
        DD_analyzed_states = []
        for target_abund in range(abund.shape[0]):
            weighted_posshifts, target_states, DD_analyzed_states = predict_weight_possible_shift(abund,
                                                                                                  empirical_states,
                                                                                                  interest_otu,
                                                                                                  target_abund,
                                                                                                  direction, theta,
                                                                                                  weighted_posshifts,
                                                                                                  target_states,
                                                                                                  DD_analyzed_states)
        weighted_posshifts = np.matrix(weighted_posshifts)
        target_states = np.matrix(target_states)
        for target_otu in range(weighted_posshifts.shape[1]):
            block = np.append(weighted_posshifts[:, target_otu], target_states, axis=1)

            ##Output analyzed state numbers
            print('/'.join([output_dir_DD, direction, '_'.join([str(interest_otu), str(target_otu), 'analyzed_states.csv'])]))
            pd.DataFrame(data=DD_analyzed_states).to_csv(
                '/'.join([output_dir_DD, direction, '_'.join([str(interest_otu), str(target_otu), 'analyzed_states.csv'])]),
                encoding='utf-8')

            ## Scaling the input
            ## Each time series is normalized to have a mean of 0 and standard deviation of 1 before analysis with S-maps
            block = (block - np.average(block, axis=0)) / np.std(block, axis=0)

            ## Use elastic_net to infer Jacobian matrices from the block
            Elastic_net_fitting(block, target_otu, interest_otu, theta, train_len,
                                cv, iteration, l_grid, '/'.join([output_dir_DD, direction]))


if __name__ == '__main__':
    # Imnput data and setting
    parse = OptionParser()
    parse.add_option('-I', '--input', dest='input', default='../inputs/month_sample.csv')
    parse.add_option('-O', '--output', dest='output', default='../outputs')
    (options, args) = parse.parse_args()
    input = options.input

    empirical_state_loc = '../inputs/level-6.csv'
    dominate_threshold = 1  # More than threshold
    zero_frequency_threshold = 20  # Less than threshold
    # target_otu = 0
    target_abund = 1
    theta = 0.1
    l_grid = 0.05
    iteration = 100000
    cv = 10
    train_len = 3
    t_range = 6
    output_dir = '/'.join([options.output, 'S-map'])
    output_dir_DD = '/'.join([options.output, 'DD_S-map'])
    uncontinous = True
    uncontinous_loc = [26,58,89,118,146,170]

    # Work in parallel
    num_cores = 40

    # Make ouput direction
    path_coefs = '/'.join([output_dir, 'coefs'])
    if not os.path.exists(path_coefs):
        os.makedirs('/'.join([output_dir, 'coefs']))
    path_fitresult = '/'.join([output_dir, 'fit_result'])
    if not os.path.exists(path_fitresult):
        os.makedirs('/'.join([output_dir, 'fit_result']))
    for direction in ['increase', 'decrease']:
        ## Make ouput direction
        path_coefs = '/'.join([output_dir_DD, direction, 'coefs'])
        if not os.path.exists(path_coefs):
            os.makedirs('/'.join([output_dir_DD, direction, 'coefs']))
        path_fitresult = '/'.join([output_dir_DD, direction, 'fit_result'])
        if not os.path.exists(path_fitresult):
            os.makedirs('/'.join([output_dir_DD, direction, 'fit_result']))

    # Select frequent and dominate genera
    abund = select_frequent_dominate_genera(input, dominate_threshold, zero_frequency_threshold, True)
    print('Output analyzed OTUs')
    abund.to_csv('/'.join([options.output,'abund.csv']))

    # Select the wanted genera in empirical states
    empirical_states = select_frequent_dominate_genera(empirical_state_loc, dominate_threshold,
                                                       zero_frequency_threshold, False)

    # Infer Jacobian matrices for each OTU and state
    Parallel(n_jobs=num_cores, backend='multiprocessing')(delayed(main)
                                                          (interest_otu) for interest_otu in range(abund.shape[1]))

