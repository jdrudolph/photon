import pandas as pd
import numpy as np
import scipy
from sklearn import metrics
from sklearn import cross_validation as cv
from sklearn.linear_model import LogisticRegression


def logistic_regression(df, features, N=10):
    """ run logistic regression model
    
    :param df: pandas.DataFrame: inlcudes features and label column 
    :param features: list of features on which to train
    :param N: number of cross-validations
    :returns (all_probas, coef, mean_auc, classifier, fig):
        prediction probabilities, mean coefficients, mean AUC of the ROC,
        matplotlib figure of the ROC
    """
    X = df[features].values
    y = df['label'].values
    classifier = LogisticRegression(class_weight='balanced')
    kfold = cv.StratifiedKFold(y, n_folds=N)
    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)
    all_probas = np.zeros_like(y, dtype=float)
    all_coef = []
    for i, (train, test) in enumerate(kfold):
        probas = classifier.fit(X[train], y[train]).predict_proba(X[test])
        all_probas[test] = probas[:,1]
        all_coef.append(classifier.coef_[0])
        fpr, tpr, _ = metrics.roc_curve(y[test], probas[:, 1])
        mean_tpr += scipy.interp(mean_fpr, fpr, tpr)
        mean_tpr[0] = 0.0
    mean_tpr = mean_tpr / len(kfold)
    mean_tpr[-1] = 1.0
    mean_auc = metrics.auc(mean_fpr, mean_tpr)
    coef = pd.DataFrame(all_coef, columns=features).mean()
    return all_probas, coef, mean_fpr, mean_tpr, mean_auc, classifier
