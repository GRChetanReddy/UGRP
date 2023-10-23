# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 12:33:06 2020

@author: Alexander van Teijlingen
"""

import pandas, numpy, h5py, json, time
import numpy as np
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.svm import LinearSVR, SVR
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

class prescreen:
    """
    As this takes the data from Judred it is all single letter encoded
    peptides, this will need to be addressed in larger datasets to avoid 
    constantly converting back and forth
    """
    def get_top_sample(self, N):
        prescreen_model = self.load_model()
        prescreen_model.fit(self.X_train, self.y_train.values.reshape(-1))
        predictions = prescreen_model.predict(self.X_test)
        prediction_indexes = predictions.argsort()[-N:][::-1]
        
        return list(self.X_test.index[prediction_indexes])
    def metrics(self, peptides, y_test):
        prescreen_model = self.load_model()
        prescreen_model.fit(self.X_train, self.y_train.values.reshape(-1))
        X_test = self.X_test.reindex(peptides)
        predictions = prescreen_model.predict(X_test)
        MSE = mean_squared_error(y_test, predictions)
        MAE = mean_absolute_error(y_test, predictions)
        r2 = r2_score(y_test, predictions)
        
        return MSE, MAE, r2
    
    def removeLogP(self, LogP):
        self.X_test = self.X_test[self.X_test["LogP WW"] < LogP]
        
    def test(self, peptides):
        prescreen_model = self.load_model()
        prescreen_model.fit(self.X_train, self.y_train.values.reshape(-1))
        tester = self.X_test.reindex(peptides)
        
        tester = tester.dropna(axis=0)
        p = prescreen_model.predict(tester)
        
        p = pandas.DataFrame(p, index=tester.index)
        return p
    
    def load_unknown(self, fname):
        h5_file = h5py.File(fname, 'r')
        names = h5_file["peptides"][()]
        data = h5_file["data"][()]  #[()]converts them to numpy.ndarray
        Features = h5_file["features"][()]
        h5_file.close()
        
        self.X_test = pandas.DataFrame(data, index=names, columns=Features)
        
    def add_data(self, peptide, AP):
        self.X_train.loc[peptide] = self.X_test.loc[peptide]
        self.y_train.loc[peptide] = AP
        self.X_test = self.X_test.drop(peptide)
    def move_data(self, peptides):
        self.X_train = self.X_test.reindex(peptides)
        self.X_test = self.X_test.drop(peptides)
    def remove_pep(self, peptides):
        self.X_test  = self.X_test.drop(peptides, errors = "ignore")
    def load_model(self):
        if self.model_type == "RF":
            self.params = {'bootstrap': False, 'criterion': 'mse', 'max_depth': None,
                       'max_features': 'sqrt', 'max_leaf_nodes': None, 'min_impurity_decrease': 0.0, 
                       'min_impurity_split': None, 'min_samples_leaf': 3, 'min_samples_split': 2, 
                       'min_weight_fraction_leaf': 0.0, 'n_estimators': 20, 'n_jobs': self.nJobs, 
                       'oob_score': False, 'random_state': 25975, 'verbose': False, 'warm_start': False}
            self.params['nJobs'] = self.nJobs
            model = RandomForestRegressor(**self.params)
        elif self.model_type == "rbf":
            with open("SVRrbf_best_HPO.json") as json_file:
                self.params = json.load(json_file)
            #self.params['nJobs'] = self.nJobs
            model = SVR(**self.params)
        elif self.model_type == "GBR":
            with open("GBR_best_HPO.json") as json_file:
                self.params = json.load(json_file)
            self.params['nJobs'] = self.nJobs
            model = GradientBoostingRegressor(**self.params)
        return model
    
    def __init__(self, ap_data, load_data=False, model_type="RF"):
        if load_data:
            h5_file = h5py.File(ap_data, 'r')
            names = h5_file["peptides"][()]
            data = h5_file["data"][()]  #[()]converts them to numpy.ndarray
            APs = h5_file["APs"][()]
            Features = h5_file["features"][()]
            h5_file.close()

            self.X_train = pandas.DataFrame(data, index=names, columns=Features)
            self.y_train = pandas.DataFrame(APs, index=names, columns=[1])
        else:
            Features = ["SPRatio", "NH2", "MW", "S", "LogP WW", "Z", "MaxASA", "RotRatio", "Bulkiness", "OH"]
            self.X_train = pandas.DataFrame(columns=Features)
            self.y_train = pandas.DataFrame(columns=[1])
        self.model_type = model_type
        self.nJobs = 1

    