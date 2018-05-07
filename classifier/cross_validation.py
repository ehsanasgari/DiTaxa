__author__ = "Ehsaneddin Asgari"
__license__ = "Apache 2"
__version__ = "1.0.0"
__maintainer__ = "Ehsaneddin Asgari"
__email__ = "asgari@berkeley.edu"
__project__ = "LLP - BioCPE"
__website__ = "https://llp.berkeley.edu/16scpe/"


import sys

sys.path.append('../')
from sklearn.model_selection import GridSearchCV, StratifiedKFold, cross_val_predict, cross_val_score
from utility.file_utility import FileUtility
from sklearn.metrics import confusion_matrix, roc_auc_score
from sklearn.metrics.classification import precision_recall_fscore_support
from sklearn.metrics.scorer import make_scorer


class CrossValidator(object):
    '''
     The Abstract Cross-Validator
    '''

    def __init__(self, X, Y):
        self.X = X
        self.Y = Y
        self.scoring = {#'auc_score_macro': make_scorer(self.roc_auc_macro),
                        #'auc_score_micro': make_scorer(self.roc_auc_micro),
                        'accuracy': 'accuracy',
                        'scores_p_1': 'precision',
                        'scores_r_1': 'recall',
                        'scores_f_1_1': 'f1',
                        'scores_p_0': make_scorer(self.precision_0),
                        'scores_r_0': make_scorer(self.recall_0),
                        'scores_f_1_0': make_scorer(self.f1_0),
                        'tnr': make_scorer(self.TNR),
                        'precision_micro': 'precision_micro',
                        'precision_macro': 'precision_macro', 'recall_macro': 'recall_macro',
                        'recall_micro': 'recall_micro', 'f1_macro': 'f1_macro', 'f1_micro': 'f1_micro'}

    def TNR(self, y_true, y_pred):
        '''
        :param y_true:
        :param y_pred:
        :return: True-negative rate
        '''
        TP = 0
        FP = 0
        TN = 0
        FN = 0
        for i in range(len(y_pred)):
            if (y_true[i] == y_pred[i]) and y_pred[i] == 1:
                TP += 1
            if y_pred[i] == 1 and (y_true[i] != y_pred[i]):
                FP += 1
            if y_true[i] == y_pred[i] and y_pred[i] == 0:
                TN += 1
            if y_pred[i] == 0 and y_true[i] != y_pred[i]:
                FN += 1
        return float(TN / (TN + FP))

    def roc_auc_macro(self, y_true, y_score):
        return roc_auc_score(y_true, y_score, average="macro")

    def roc_auc_micro(self, y_true, y_score):
        return roc_auc_score(y_true, y_score, average="micro")

    def precision_0(self, y_true, y_pred, labels=None, average='binary', sample_weight=None):
        '''
        :param y_true:
        :param y_pred:
        :param labels:
        :param average:
        :param sample_weight:
        :return: calculate prec for neg class
        '''
        p, _, _, _ = precision_recall_fscore_support(y_true, y_pred,
                                                     beta=1,
                                                     labels=labels,
                                                     pos_label=0,
                                                     average=average,
                                                     warn_for=('f-score',),
                                                     sample_weight=sample_weight)
        return p

    def recall_0(self, y_true, y_pred, labels=None, average='binary', sample_weight=None):
        '''
        :param y_true:
        :param y_pred:
        :param labels:
        :param average:
        :param sample_weight:
        :return: calculate recall for neg class
        '''
        _, r, _, _ = precision_recall_fscore_support(y_true, y_pred,
                                                     beta=1,
                                                     labels=labels,
                                                     pos_label=0,
                                                     average=average,
                                                     warn_for=('f-score',),
                                                     sample_weight=sample_weight)
        return r

    def f1_0(self, y_true, y_pred, labels=None, average='binary', sample_weight=None):
        '''
        :param y_true:
        :param y_pred:
        :param labels:
        :param average:
        :param sample_weight:
        :return: calculate f1 for neg class
        '''
        _, _, f, _ = precision_recall_fscore_support(y_true, y_pred,
                                                     beta=1,
                                                     labels=labels,
                                                     pos_label=0,
                                                     average=average,
                                                     warn_for=('f-score',),
                                                     sample_weight=sample_weight)
        return  f


class KFoldCrossVal(CrossValidator):
    '''
        K-fold cross-validation tuning and evaluation
    '''

    def __init__(self, X, Y, folds=10, random_state=1):
        '''
        :param X:
        :param Y:
        :param folds:
        :param random_state:
        '''
        CrossValidator.__init__(self, X, Y)
        self.cv = StratifiedKFold(n_splits=folds, shuffle=True, random_state=random_state)
        self.X = X
        self.Y = Y

    def tune_and_evaluate(self, estimator, parameters, score='f1_macro', n_jobs=-1, file_name='results'):
        '''
        :param estimator:
        :param parameters:p
        :param score:
        :param n_jobs:
        :param file_name: directory/tuning/classifier/features/
        :return:
        '''
        # greed_search
        self.greed_search = GridSearchCV(estimator=estimator, param_grid=parameters, cv=self.cv, scoring=self.scoring,
                                         refit=score, error_score=0, n_jobs=n_jobs)

        label_set = list(set(self.Y))
        # fitting
        self.greed_search.fit(X=self.X, y=self.Y)
        y_predicted = cross_val_predict(self.greed_search.best_estimator_, self.X, self.Y,  cv=self.cv)
        conf = confusion_matrix(self.Y, y_predicted, labels=label_set)
        # save in file
        FileUtility.save_obj(file_name,
                             [label_set, conf, self.greed_search.best_score_, self.greed_search.best_estimator_,
                              self.greed_search.cv_results_, self.greed_search.best_params_,  (y_predicted, self.Y,label_set )])


class NestedCrossVal(CrossValidator):
    '''
    Nested cross-validation
    '''

    def __init__(self, X, Y, inner_folds=10, outer_folds=10, random_state=1, opt_f1_class=0):
        '''
        :param X:
        :param Y:
        :param inner_folds:
        :param outer_folds:
        :param random_state:
        :param opt_f1_class:
        '''
        CrossValidator.__init__(self, X, Y, opt_f1_class=opt_f1_class)
        self.inner_cv = StratifiedKFold(n_splits=inner_folds, shuffle=True, random_state=random_state)
        self.outer_cv = StratifiedKFold(n_splits=outer_folds, shuffle=True, random_state=random_state)

    def tune_and_evaluate(self, estimator, parameters, score='f1_macro', file_name='results'):
        '''
        :param estimator:
        :param parameters:
        :param score:
        :param file_name: directory/tuning/classifier/features/
        :return:
        '''
        # inner cross_validation
        self.greed_search = GridSearchCV(estimator=estimator, param_grid=parameters, cv=self.inner_cv,
                                         scoring=self.scoring, refit=score, error_score=0)
        # Nested CV with parameter optimization
        self.nested_score = cross_val_score(self.greed_search, X=self.X, y=self.Y, cv=self.outer_cv)

        # saving
        FileUtility.save_obj([self.greed_search, self.nested_score], file_name)


class PredefinedFoldCrossVal(CrossValidator):
    '''
        Predefined folds
    '''

    def __init__(self, X, Y, isolate_list, fold_file):
        '''
        :param X:
        :param Y:
        :param folds:
        :param random_state:
        '''
        CrossValidator.__init__(self, X, Y)

        map_to_idx = {isolate: idx for idx, isolate in enumerate(isolate_list)}
        splits = [[map_to_idx[item] for item in fold_list.split() if item in map_to_idx] for fold_list in
                  FileUtility.load_list(fold_file)]
        new_splits = []
        for i in range(len(splits)):
            train = [j for i in splits[:i] + splits[i + 1:] for j in i]
            test = splits[i]
            new_splits.append([train, test])
        self.cv = new_splits
        self.X = X
        self.Y = Y

    def tune_and_evaluate(self, estimator, parameters, score='f1_macro', n_jobs=-1, file_name='results'):
        '''
        :param estimator:
        :param parameters:p
        :param score:
        :param n_jobs:
        :param file_name: directory/tuning/classifier/features/
        :return:
        '''
        # greed_search
        self.greed_search = GridSearchCV(estimator=estimator, param_grid=parameters, cv=self.cv, scoring=self.scoring,
                                         refit=score, error_score=0, n_jobs=n_jobs)

        label_set = list(set(self.Y))
        # fitting
        self.greed_search.fit(X=self.X, y=self.Y)
        y_predicted = cross_val_predict(self.greed_search.best_estimator_, self.X, self.Y,  cv=self.cv)
        conf = confusion_matrix(self.Y, y_predicted, labels=label_set)
        # save in file
        FileUtility.save_obj(file_name,
                             [label_set, conf, self.greed_search.best_score_, self.greed_search.best_estimator_,
                              self.greed_search.cv_results_, self.greed_search.best_params_, (y_predicted, self.Y,label_set )])
