import os
import json
import time
import logging
import pandas as pd
from shutil import move
from datetime import datetime
from datetime import datetime

import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import IsolationForest
from sklearn.metrics import accuracy_score
import joblib

from sensor_related import DataPreProcessor
from sensor_related import ModelSaver
from sensor_related import PlotSensor

class DataProcessor:
    '''
    This class is used to check the input directory and do the following steps:
    1- If there is a new .csv file, do preprocess the file,
    2- Use the IsolationForest method to find anomalies,
    3- Evaluate the model on test and validation data,
    4- Plot with the given sensor names and save on image output directory.
    5- Move the file from input directory to output directory.
    Also, logging method is used to debugging. the logging file is saved on output directory.
    '''
    def __init__(self, config_file):
        '''
        parameter: config_file: all information regarding the directories and the sensors to draw are in it
        '''
        self.config = self.load_config(config_file)
        self.input_dir = self.config['input_directory']
        self.output_dir = self.config['output_directory']
        self.img_dir = self.config['img_directory']
        self.sensors_to_draw = self.config['sensors_to_draw']
        self.check_interval = self.config['check_interval']
        self.logger = self.setup_logger()

    def load_config(self, config_file):
        '''
        This method is used to load the config file
        '''
        with open(config_file) as f:
            config = json.load(f)
        return config

    def setup_logger(self):
        '''
        This method sets all options regarding the logger like name of file, direction, its level, format ...
        '''
        logger = logging.getLogger('DataProcessor')
        logger.setLevel(logging.INFO)

        log_file = os.path.join(self.output_dir, 'data_processor.log')
        file_handler = logging.FileHandler(log_file) # To create a log file just for this class
        file_handler.setLevel(logging.INFO)

        formatter = logging.Formatter('%(asctime)s -%(levelname)s - %(message)s ')
        file_handler.setFormatter(formatter)

        logger.addHandler(file_handler)

        return logger

    def get_new_data_files(self):
        '''
        It check if there are new data file and add in a list
        returns: The list of new files
        '''
        files = []
        for file in os.listdir(self.input_dir):
            if file.endswith('.csv'):
                files.append(file)
        return files

    def process_data_file(self, file):
        '''
        This method do all preprocess methods, modelling, and plotting of a given file and save the plot

        '''
        self.logger.info(f'Found new data file: {file}')

        try:
            self.logger.info(f'Loaded the file: {file}')
            file_path = os.path.join(self.input_dir, file)
            self.logger.info(f'data: {file_path}')
            data = pd.read_csv(file_path)
            
            preprocessor = DataPreProcessor(self.output_dir, data)
            # preprocessed_data = preprocessor.all_preprocess()

            model_saver = ModelSaver('model.pkl', self.output_dir, preprocessor)

            model_saver.split()
            model_saver.modeling(IsolationForest, IsolationForest())
            self.logger.info('model persisted in a .pkl file')

            predictions = model_saver.prediction_evaluation()
            self.logger.info('Received predictions')

            # self.save_predictions(file, preprocessed_data, predictions)
            self.generate_sensor_images(model_saver.testing_data)

            self.logger.info(f'Processed data file: {file}')
            self.remove_data_file(file)

        except Exception as e:
            self.logger.error(f'Error processing data file: {file}\n{str(e)}')


    def process_data_files(self):
        '''
        This method is to manage all new files, between every two files use time.sleep to stop for the given seconds.
        If in the input folder there is no new file, the runnig will be stopped by 'break'. 
        If we remove the 'break' command, after the time.sleep interval it again checkes the input directory! 
        '''
        self.logger.info('Data Processor started. Monitoring input directory...')

        while True:
            files = self.get_new_data_files()
            if files:
                for file in files:
                    self.process_data_file(file)
            else:
                self.logger.info('No new data files found. Resuming listening.')
                break

            time.sleep(self.check_interval)

    def generate_sensor_images(self, data):
        '''
        This method is to generate the image of the plots on the image directory
        '''
        for sensor_name in self.sensors_to_draw:
            plot_sensor= PlotSensor(sensor_name, data)
            plot_sensor.broken_recovery_normal()
            fig = plot_sensor.plot_sensor_anomalies()

            timestamp = datetime.now().strftime('%Y-%m')
            image_file = f'{timestamp}-{sensor_name}.png'
            image_path = os.path.join(self.img_dir, image_file)

            fig.savefig(image_path)
            plt.close(fig)

            self.logger.info(f'Saving image: {image_file}')

    def remove_data_file(self, file):
        '''
        This method is to move the processed file from input directory to output directory.
        It can be removed. But here we needed the data, so, hust moved
        '''
        try:
            file_path = os.path.join(self.input_dir, file)
            move(file_path, os.path.join(self.output_dir, file))
            self.logger.info(f'Removed data file: {file}')
        except Exception as e:
            self.logger.error(f'Error removing data file: {file}\n{str(e)}') 



if __name__ == '__main__':
    processor = DataProcessor('application.json')
    processor.process_data_files()