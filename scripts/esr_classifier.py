import sklearn.ensemble
import sklearn.cross_validation
import sklearn.preprocessing
from sklearn.cross_validation import KFold
from sklearn import cross_validation
from sklearn import datasets
from sklearn import svm
from sklearn import preprocessing
import numpy as np
import scipy as sp
import os
import pandas as pd
import h5py
import re
import scipy.optimize
from sklearn.naive_bayes import GaussianNB


classifier_globals = {'trained_scaler': None,
                      'forest_classifier' : None,
                      'svm_classifier' : None,
                      'nb_classifier' : None}



def esr_differential_evolution(data):
    bounds = [(np.min(data[:,1]),np.max(data[:,1])),
             (-2*(np.max(data[:,1]) - np.min(data[:,1])), 2*(np.max(data[:,1]) - np.min(data[:,1]))),
             (-2*(np.max(data[:,1]) - np.min(data[:,1])), 2*(np.max(data[:,1]) - np.min(data[:,1]))),
             (0.002,0.015),
             (np.min(data[:,0]),np.max(data[:,0])),
             (np.min(data[:,0]),np.max(data[:,0]))]
    args = (data[:,0], data[:,1])
    de_out = scipy.optimize.differential_evolution(spectrum_cost, bounds, args= args, popsize=50, strategy='rand2exp',
                                                   maxiter=250)
    return de_out

def spectrum_cost(x, *data):
    xd, yd = data
    C, A, B, gamma1, x1, x2 = x
    #data[:,0]
    #yd = data[:,1]
    m = np.size(yd)
    y_ideal = C-A*(gamma1**2/(gamma1**2+(xd-x1)**2)) - B*gamma1**2/(gamma1**2+(xd-x2)**2)
    J = 1.0/m * np.sum((yd-y_ideal)**2 * 1.0/(2.0*np.abs(yd)) )
    return J


def classify_spectrum(data):
    # check if the scaler and classifier are trained. if not, train them.
    if classifier_globals['trained_scaler'] == None:
        # train the scaler
        print 'Training scaler...'
        X_array, Y_array = import_training_data()
        X_norm = train_scaler(X_array)

    if classifier_globals['forest_classifier'] == None:
        # train the RF classifier
        print 'Training Forest classifier...'
        X_array, Y_array = import_training_data()
        X_norm = train_scaler(X_array)
        train_forest_classifier(X_norm, Y_array)

    if classifier_globals['svm_classifier'] == None:
        # train the RF classifier
        print 'Training SVM classifier...'
        X_array, Y_array = import_training_data()
        X_norm = train_scaler(X_array)
        train_svm_classifier(X_norm, Y_array)
    if classifier_globals['nb_classifier'] == None:
        # train the NB classifier
        print 'Training Naive Bayes classifier...'
        X_array, Y_array = import_training_data()
        X_norm = train_scaler(X_array)
        train_nb_classifier(X_norm,Y_array)

    # now fit the data using differential evolution followed by a sequentual quadratic programming (SQP) step
    X_sample = fit_esr_return_features(data)
    # normalize the sample using the scaler
    X_sample_normed = scale_sample(X_sample)
    # now predict the class using the classifier
    predicted_class = classify_data(X_sample_normed)
    return predicted_class


def fit_esr_return_features(data):
    de_out = esr_differential_evolution(data)
    # now use BFGS to really minimize
    args = (data[:,0], data[:,1])
    data = data
    bounds = [(np.min(data[:,1]),np.max(data[:,1])),
             (-1.2*(np.max(data[:,1]) - np.min(data[:,1])), 1.2*(np.max(data[:,1]) - np.min(data[:,1]))),
             (-1.2*(np.max(data[:,1]) - np.min(data[:,1])), 1.2*(np.max(data[:,1]) - np.min(data[:,1]))),
             (0.002,0.02),
             (np.min(data[:,0]),np.max(data[:,0])),
             (np.min(data[:,0]),np.max(data[:,0]))]
    # sqp step
    bf_out = scipy.optimize.minimize(spectrum_cost,de_out.x,args=args,bounds = bounds,
                                     method='SLSQP', options={'ftol':0.0001, 'maxiter':400})
    #print(bf_out)
    C, A, B, gamma1, x1, x2 = bf_out.x
    xd = data[:,0]
    y_fit = C-A*gamma1**2/(gamma1**2+(xd-x1)**2) - B*gamma1**2/(gamma1**2+(xd-x2)**2)
    #sns.lmplot(x="x", y="y", data=wm_df, fit_reg=False)
    #plt.plot(xd,.0*i + y_fit, color="k", alpha=1.0)
    J = spectrum_cost(bf_out.x, data[:,0], data[:,1])
    if x1 < x2:
        left_amp = A
        left_loc = x1
        right_amp = B
        right_loc = x2
    else:
        left_amp = B
        left_loc = x2
        right_amp = A
        right_loc = x1
    # compute standard deviation of first 10 points for signal-to-noise features
    std_dev = np.std(data[0:10,1])
    snrL = np.abs(A/std_dev)
    snrR = np.abs(B/std_dev)
    X_sample = np.array((left_amp/std_dev,right_amp/std_dev,np.log(gamma1),left_loc,right_loc))
    return X_sample

def import_training_data():
    filenameX = 'X_array_20160513.csv'
    filenameY = 'Y_array_20160513.csv'
    csvX = np.genfromtxt(filenameX, delimiter=",")
    csvY = np.genfromtxt(filenameY,delimiter=",")
    return csvX, csvY

def train_scaler(X_array):
    rbs = sklearn.preprocessing.StandardScaler()
    X_normalized = rbs.fit_transform(X_array)
    classifier_globals['trained_scaler'] = rbs
    return X_normalized

def scale_sample(X_sample):
    return classifier_globals['trained_scaler'].transform(X_sample)


def train_svm_classifier(X_normalized, Y_array):
    clf = svm.SVC(kernel='rbf', C=1.0, gamma=0.1).fit(X_normalized, Y_array)
    classifier_globals['svm_classifier'] = clf
    return clf

def train_nb_classifier(X_normalized, Y_array):
    gnb = GaussianNB()
    gnb.fit(X_normalized, Y_array)
    classifier_globals['nb_classifier'] = gnb
    return gnb


def train_forest_classifier(X_normalized, Y_array):
    # print(X_normalized.shape)
    # print('Length of XN is {}, length of yarray is {}'.format(np.size(X_normalized[:,0]), np.size(Y_array)))
    # n_est_list = [ 3, 8, 9, 10, 20, 30]
    # for n_est in n_est_list:
    #     # Create the random forest object which will include all the parameters
    #     # for the fit
    #     forest = sklearn.ensemble.RandomForestClassifier(n_estimators = n_est)
    #
    #
    #     kf = KFold(np.size(X_normalized[:,0]), n_folds=5)
    #     total_incorrect = 0
    #     total_tests = 0
    #     for train, test in kf:
    #         #print(train)
    #         #print('Starting next fold...')
    #         # Fit the training data to the Survived labels and create the decision trees
    #         forest = forest.fit(X_normalized[train,:],Y_array[train])
    #
    #         # Take the same decision trees and run it on the test data
    #         output = forest.predict(X_normalized[test,:])
    #         total_incorrect = total_incorrect + (Y_array[test] != output).sum()
    #         inc_idx = Y_array[test] != output
    #         #print('Predicted {}, actual was {}'.format(output[inc_idx], y_array[test][inc_idx]))
    #         total_tests = total_tests + np.size(Y_array[test])
    #     print('For n_est {}, total incorrect was {}/{}'.format(n_est, total_incorrect,total_tests))

    # create a final forest estimator with n_est = 9:
    final_forest = sklearn.ensemble.RandomForestClassifier(n_estimators = 9)
    final_forest.fit(X_normalized,Y_array)
    classifier_globals['forest_classifier'] = final_forest
    return final_forest

def classify_data(X_element):
    # SVM classifier
    #class_out = int(classifier_globals['svm_classifier'].predict(X_element)[0])
    if classifier_globals['forest_classifier'] == None:
        print 'Training RF classifier again...'
        X_array, Y_array = import_training_data()
        X_norm = train_scaler(X_array)
        final_forest = train_forest_classifier(X_norm, Y_array)
        classifier_globals['forest_classifier'] = final_forest
    else:
        final_forest = classifier_globals['forest_classifier']


    if classifier_globals['nb_classifier'] == None:
        print 'Training NB classifier again...'
        X_array, Y_array = import_training_data()
        X_norm = train_scaler(X_array)
        gnb = train_nb_classifier(X_norm, Y_array)
        classifier_globals['nb_classifier'] = gnb
    if classifier_globals['svm_classifier'] == None:
        print 'Training SVM classifier again...'
        X_array, Y_array = import_training_data()
        X_norm = train_scaler(X_array)
        svmo = train_svm_classifier(X_norm, Y_array)
        classifier_globals['svm_classifier'] = svmo
    # naive bayes
    classifier = classifier_globals['svm_classifier']
    class_out = int(classifier.predict(X_element)[0])
    return class_out

