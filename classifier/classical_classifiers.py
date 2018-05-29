__author__ = "Ehsaneddin Asgari"
__license__ = "Apache 2"
__version__ = "1.0.0"
__maintainer__ = "Ehsaneddin Asgari"
__email__ = "asgari@berkeley.edu"
__project__ = "LLP - DiTaxa"
__website__ = "https://llp.berkeley.edu/ditaxa/"


import sys
sys.path.append('../')
from sklearn.svm import LinearSVC, SVC
from classifier.cross_validation import KFoldCrossVal, PredefinedFoldCrossVal
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from utility.file_utility import FileUtility
import numpy as np
import codecs
import math
import operator


class SVM:
    '''
        Support vector machine classifier
    '''

    def __init__(self, X, Y, clf_model='LSVM'):
        if clf_model == 'LSVM':
            self.model = LinearSVC(C=1.0)  # , multi_class='ovr'
            self.type = 'linear'
        else:
            self.model = SVC(C=1.0, kernel='rbf')
            self.type = 'rbf'
        self.X = X
        self.Y = Y

    def tune_and_eval(self, results_file,
                      params=[{'C': [1000, 500, 200, 100, 50, 20, 10, 5, 2, 1, 0.2, 0.5, 0.01, 0.02, 0.05, 0.001],
                               'penalty': ['l1'], "tol": [1e-06, 1e-04], 'dual': [False, True], "fit_intercept": [True],
                               'loss': ['l2'], 'class_weight': ['balanced', None]}], njobs=50, kfold=10):
        '''
        K-fold cross-validation
        :param results_file: file to save the results
        :param params: parameters to be tuned
        :param njobs: number of cores
        :param kfold: number of folds
        :return:
        '''
        CV = KFoldCrossVal(self.X, self.Y, folds=kfold)
        CV.tune_and_evaluate(self.model, parameters=params, score='f1_macro', file_name=results_file + '_SVM',
                             n_jobs=njobs)

    def tune_and_eval_predefined(self, results_file, isolates, folds, params=[
        {'C': [1000, 500, 200, 100, 50, 20, 10, 5, 2, 1, 0.2, 0.5, 0.01, 0.02, 0.05, 0.001], 'penalty': ['l1'],
         "tol": [1e-06, 1e-04], 'dual': [False, True], "fit_intercept": [True], 'loss': ['l2'],
         'class_weight': ['balanced', None]}], njobs=50):
        '''
        :param results_file:
        :param isolates:
        :param folds:
        :param params:
        :param njobs:
        :return:
        '''
        self.CV = PredefinedFoldCrossVal(self.X, self.Y, isolates, folds)
        self.CV.tune_and_evaluate(self.model, parameters=params, score='f1_macro', file_name=results_file + '_SVM',
                                  n_jobs=njobs)


class LogRegression:
    '''
        LR classifier
    '''

    def __init__(self, X, Y):
        '''

        :param X:
        :param Y:
        '''
        self.model = LogisticRegression(C=1.0)
        self.X = X
        self.Y = Y

    def tune_and_eval(self, results_file,
                      params=[{'C': [1000, 500, 200, 100, 50, 20, 10, 5, 2, 1, 0.2, 0.5, 0.01, 0.02, 0.05, 0.001],
                               'penalty': ['l1'], "tol": [1e-06, 1e-04], 'dual': [False, True], "fit_intercept": [True],
                               'class_weight': ['balanced', None], 'solver': ['liblinear']}], njobs=50, kfold=10):
        '''
        :param results_file:
        :param params:
        :param njobs:
        :param kfold:
        :return:
        '''
        CV = KFoldCrossVal(self.X, self.Y, folds=kfold)
        CV.tune_and_evaluate(self.model, parameters=params, score='f1_macro', file_name=results_file + '_LR',
                             n_jobs=njobs)

    def tune_and_eval_predefined(self, results_file, isolates, folds, params=[
        {'C': [1000, 500, 200, 100, 50, 20, 10, 5, 2, 1, 0.2, 0.5, 0.01, 0.02, 0.05, 0.001], 'penalty': ['l1'],
         "tol": [1e-06, 1e-04], 'dual': [False, True], "fit_intercept": [True], 'class_weight': ['balanced', None],
         'solver': ['liblinear']}], njobs=50):
        '''
        :param results_file:
        :param isolates:
        :param folds:
        :param params:
        :param njobs:
        :return:
        '''
        self.CV = PredefinedFoldCrossVal(self.X, self.Y, isolates, folds)
        self.CV.tune_and_evaluate(self.model, parameters=params, score='f1_macro', file_name=results_file + '_LR',
                                  n_jobs=njobs)


class RFClassifier:
    '''
        Random forest classifier
    '''

    def __init__(self, X, Y):
        '''
        :param X:
        :param Y:
        '''
        self.model = RandomForestClassifier(bootstrap=True, criterion='gini',
                                            min_samples_split=2, max_features='auto', min_samples_leaf=1,
                                            n_estimators=1000)
        self.X = X
        self.Y = Y

    def tune_and_eval(self, results_file, params=None, feature_names=None, njobs=50, kfold=10):
        '''
        :param results_file:
        :param params:
        :param feature_names:
        :param njobs:
        :param kfold:
        :return:
        '''
        if params is None:
            params = [{"n_estimators": [100, 200, 500, 1000],
                       "criterion": ["entropy"],  # "gini",
                       'max_features': ['sqrt', 'auto'],  # 'auto',
                       'min_samples_split': [2, 5, 10],  # 2,5,10
                       'min_samples_leaf': [1], 'class_weight': ['balanced', None]}]
        self.CV = KFoldCrossVal(self.X, self.Y, folds=kfold)
        self.CV.tune_and_evaluate(self.model, parameters=params, score='f1_macro', file_name=results_file + '_RF',
                                  n_jobs=njobs)
        if feature_names is not None:
            try:
                [label_set, conf, best_score_, best_estimator_, cv_results_, best_params_,
                 (y_predicted, Y, label_set)] = FileUtility.load_obj(results_file + '_RF.pickle')
            except:
                [label_set, best_score_, best_estimator_, cv_results_, best_params_,
                 (Y, label_set)] = FileUtility.load_obj(results_file + '_RF.pickle')
            self.generate_RF_important_features(best_estimator_, feature_names, results_file, 500)

    def tune_and_eval_predefined(self, results_file, isolates, folds, params=None, feature_names=None, njobs=50):
        '''
        :param results_file:
        :param isolates:
        :param folds:
        :param params:
        :param feature_names:
        :param njobs:
        :return:
        '''
        if params is None:
            params = [{"n_estimators": [100, 200, 500, 1000],
                       "criterion": ["entropy"],  # "gini",
                       'max_features': ['sqrt', 'auto'],  # 'auto',
                       'min_samples_split': [2, 5, 10],  # 2,5,10
                       'min_samples_leaf': [1, 2], 'class_weight': ['balanced', None]}]
        self.CV = PredefinedFoldCrossVal(self.X, self.Y, isolates, folds)
        self.CV.tune_and_evaluate(self.model, parameters=params, score='f1_macro', file_name=results_file + '_RF',
                                  n_jobs=njobs)
        if feature_names is not None:
            try:
                [label_set, conf, best_score_, best_estimator_, cv_results_, best_params_,
                 (y_predicted, Y, label_set)] = FileUtility.load_obj(results_file + '_RF.pickle')
            except:
                [label_set, best_score_, best_estimator_, cv_results_, best_params_,
                 (Y, label_set)] = FileUtility.load_obj(results_file + '_RF.pickle')
            self.generate_RF_important_features(best_estimator_, feature_names, results_file, 1000)

    def generate_RF_important_features(self, clf_random_forest, feature_names, results_file, N):
        '''
        :param clf_random_forest:
        :param feature_names:
        :param results_file:
        :param N:
        :return:
        '''
        file_name = results_file + 'RF_features'
        clf_random_forest.fit(self.X, self.Y)
        std = np.std([tree.feature_importances_ for tree in clf_random_forest.estimators_], axis=0)

        scores = {feature_names[i]: (s, std[i]) for i, s in enumerate(list(clf_random_forest.feature_importances_)) if
                  not math.isnan(s)}
        scores = sorted(scores.items(), key=operator.itemgetter([1][0]), reverse=True)[0:N]
        f = codecs.open(file_name, 'w')
        f.write('\t'.join(['feature', 'score', 'std', '#I-out-of-' + str(np.sum(self.Y)),
                           '#O-out-of-' + str(len(self.Y) - np.sum(self.Y))]) + '\n')
        for w, score in scores:
            feature_array = self.X[:, feature_names.index(w)]
            pos = [feature_array[idx] for idx, x in enumerate(self.Y) if x == 1]
            neg = [feature_array[idx] for idx, x in enumerate(self.Y) if x == 0]
            f.write('\t'.join([str(w), str(score[0]), str(score[1]), str(np.sum(pos)), str(np.sum(neg))]) + '\n')
        f.close()


class KNN:
    '''
        K-nearest neighbor classifier
    '''

    def __init__(self, X, Y):
        '''
        :param X:
        :param Y:
        '''
        self.model = KNeighborsClassifier(n_neighbors=3)
        self.X = X
        self.Y = Y

    def tune_and_eval(self, results_file, params=None, njobs=50, kfold=10):
        '''
        :param results_file:
        :param params:
        :param njobs:
        :param kfold:
        :return:
        '''
        if params is None:
            params = [{"n_neighbors": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20], 'weights': ['uniform', 'distance']}]
        self.CV = KFoldCrossVal(self.X, self.Y, folds=kfold)
        self.CV.tune_and_evaluate(self.model, parameters=params, score='f1_macro', file_name=results_file + '_KNN',
                                  n_jobs=njobs)

    def tune_and_eval_predefined(self, results_file, isolates, folds, params=None, njobs=50):
        '''
        :param results_file:
        :param isolates:
        :param folds:
        :param params:
        :param njobs:
        :return:
        '''
        if params is None:
            params = [{"n_neighbors": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20], 'weights': ['uniform', 'distance']}]
        self.CV = PredefinedFoldCrossVal(self.X, self.Y, isolates, folds)
        self.CV.tune_and_evaluate(self.model, parameters=params, score='f1_macro', file_name=results_file + '_KNN',
                                  n_jobs=njobs)

