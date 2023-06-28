import os
import json
import time
import logging
import pandas as pd
from shutil import move
from datetime import datetime

import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import IsolationForest
from sklearn.metrics import accuracy_score
import joblib

class DataPreProcessor:
    '''
    This class get the data, and prepare it for anomaly detection
    '''
    def __init__(self, output_dir,data):
        
        self.output_dir = output_dir
        self.data = data
        
    def remove_unnamed_columns(self):
        
        unnamed_columns = [col for col in self.data.columns if isinstance(col, str) and 'Unnamed' in col]
        self.data = self.data.drop(columns=unnamed_columns)
        return self.data
        
    def preprocess_time_columns(self):
        
        time_columns = [col for col in self.data.columns if isinstance(col, str) and 'time' in col.lower()]
        for col in time_columns:
            self.data[col] = pd.to_datetime(self.data[col])
        if time_columns:
            self.data.set_index(time_columns, inplace=True)
        return self.data
    
    def remove_null_columns(self, threshold=0.3):
        
        null_counts = self.data.isnull().sum()
        null_columns = null_counts[null_counts > threshold * self.data.shape[0]].index
        self.data = self.data.drop(columns=null_columns)
        return self.data
    
    def fill_nulls_with_mean(self):
        
        self.data.fillna(self.data.mean(), inplace=True)
        return self.data
    
    def all_preprocess(self):
        
        self.data = self.remove_unnamed_columns()
        self.data = self.preprocess_time_columns()
        self.data = self.remove_null_columns()
        self.data = self.fill_nulls_with_mean()
        return (self.data)
    


class ModelSaver:
    '''
    THis class is to fit a model, persist it, and evaluate the model
    '''
    def __init__(self, model_name, output_dir, DataPreProcessor):
        
        self.model_name = model_name # to persist the model
        self.output_dir = output_dir
        self.DataPreProcessor = DataPreProcessor
        self.model_file = os.path.join(self.output_dir, self.model_name)

    def split(self, train_start = '2018-04-01',train_stop = '2018-06-30', 
               test_start = '2018-07-01', test_stop ='2018-07-31', 
               valid_start = '2018-08-01', valid_stop = '2018-08-31'):
        
        # load preprocessed data
        data = self.DataPreProcessor.all_preprocess()
        
        # split: Training,testing,validation
        self.training_data = data.loc[train_start:train_stop]
        self.testing_data = data.loc[test_start:test_stop]
        self.validation_data = data.loc[valid_start:valid_stop]

    def transform(self, data): # data = self.training_data or self.testing_data or self.validation_data
        
        # transform the data on Standard Scaler
        m, n = data.shape
        X = data.iloc[:,:n-1] # ignore machine status columns
        scaler = StandardScaler()
        X = scaler.fit_transform(X)
        return X
    
    def modeling(self, algorithm_name, algorithm_method):
        
        X = self.transform(self.training_data)
        clf = algorithm_method.fit(X)
        self.model_file = os.path.join(self.output_dir, self.model_name)
        # Persist the model
        joblib.dump(clf, self.model_file)
        
    def prediction_evaluation(self):
        
        #load the model
        clf = joblib.load(self.model_file)
        
        # Convert categorical labels to numerical values
        label_mapping = {"NORMAL":1, "BROKEN": -1, "RECOVERING": -1}
        y_true_test = self.testing_data['machine_status'].map(label_mapping)
        y_true_valid = self.validation_data['machine_status'].map(label_mapping)
        
        # predict on test and valid
        X1 = self.transform(self.testing_data)
        X2 = self.transform(self.validation_data)
        
        y_pred_test = clf.predict(X1)
        y_pred_valid = clf.predict(X2)
        
        self.testing_data['pred'] = y_pred_test
        self.validation_data['pred'] = y_pred_valid      
        
        # accuracy score
        accuracy_test = accuracy_score(y_true_test, y_pred_test)
        accuracy_valid = accuracy_score(y_true_valid, y_pred_valid)
        
        return (accuracy_test, accuracy_valid)


class PlotSensor:
    '''
    This class get data and the given sensor name, plot and save it
    '''

    def __init__(self, sensor_name, data):
        
        self.sensor_name = sensor_name
        self.data = data
    
    def broken_recovery_normal(self):
        
        self.broken_rows = self.data[self.data['machine_status']=='BROKEN']
        self.recovery_rows = self.data[self.data['machine_status']=='RECOVERING']
        self.normal_rows = self.data[self.data['machine_status']=='NORMAL']
    
    def plot_sensor_anomalies(self):
        anomoly_rows = self.data[self.data['pred'] == -1]
        

        fig, ax = plt.subplots(figsize=(25, 3))
        ax.plot(self.data[self.sensor_name], color='grey')
        ax.plot(self.recovery_rows[self.sensor_name], linestyle='none', marker='o', color='yellow', markersize=5, label='recovering', alpha=0.5)
        ax.plot(self.broken_rows[self.sensor_name], linestyle='none', marker='X', color='red', markersize=20, label='broken')
        ax.plot(anomoly_rows[self.sensor_name], linestyle='none', marker='X', color='blue', markersize=4, label='anomaly predicted', alpha=0.1)
        ax.set_title(self.sensor_name)
        ax.legend()

        return fig