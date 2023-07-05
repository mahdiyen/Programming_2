import json
import linecache as lc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time

class CsvConverter:
    def __init__(self,file_name, header):
        self.file_name = file_name
        self.header = header
        self.keys = self.header.split(',')

    
    def convert_csv_to_json(self,lines):
        """
        Convert CSV lines to a JSON string.
        Args:
            lines (list): List of CSV lines.
        Returns:
            str: JSON string representing the converted data.
        """
        json_list=[]
        for line in lines:
            try:
                assert len(line) == len(self.keys)

                json_dict = {key: value for key, value in zip(self.keys, line)}
                json_list.append(json_dict)
            except AssertionError:
                print ("Number of items in line does not match the number of keys in header.")
            
        return json.dumps(json_list)
        



class Reader:
    def __init__(self, location, csvconverter, period):
        self.location = location
        self.csvconverter = csvconverter
        self.starter = 2
        self.period = period
        self.observers = set()
        self.is_reading = False
        
    def get_lines(self):
        """
        Read lines from the CSV file.
        Returns:
            str: JSON string representing the converted data.
        """
        lines_list = list()
        for i in range(self.starter, self.starter + self.period):
            line = lc.getline(self.csvconverter.file_name, i).strip()
            if line:
                lines_list.append(line.split(','))
            else:
                break
        self.starter = self.starter + self.period
        if not lines_list:
            return ""

        return self.csvconverter.convert_csv_to_json(lines_list)

    def add_observer(self, observer):
        """
        Add an observer to the Reader's list of observers.
        Args:
            observer: Observer instance to be added.
        """
        self.observers.add(observer)

    def remove_observer(self, observer):
        """
        Remove an observer from the Reader's list of observers if present.
        Args:
            observer: Observer instance to be removed.
        """
        self.observers.discard(observer)

    def notify_observers(self):
        """
        Notify all registered observers by calling their update() method.
        """
        for observer in self.observers:
            observer.update()

    def start_reading(self):
        """
        Start reading the CSV file, periodically notifying the observers.
        The method reads a chunk of lines from the CSV file, converts them to JSON,
        and notifies the observers. It repeats this process until all lines have been read.
        The method waits for 5 seconds between iterations.
        """
        self.is_reading = True
        while self.is_reading:
            self.data_json = self.get_lines()
            self.notify_observers()
            if not self.data_json:
                self.is_reading = False
            time.sleep(5)


class AverageYear:
    def __init__(self, reader):
        self.reader = reader
        self.average_year = None

    def calculate_average(self):
        if self.reader.data_json:
            data = pd.read_json(self.reader.data_json)
            self.average_year = data.mean()['J-D']
            return self.average_year
        else:
            return None

    def plot_average(self):
        if self.average_year is not None:
            if self.reader.data_json:
                data = pd.read_json(self.reader.data_json)
                years = data['Year']
                temperatures = data['J-D']

                plt.plot(years, temperatures)
                plt.xlabel('Year')
                plt.ylabel('Temperature Anomaly')
                plt.title('Average Yearly Temperature Anomaly')
                plt.axhline(self.average_year, color='red', linestyle='--', label='Average')
                plt.legend()
                plt.show()
            else:
                print("No data available.")
        else:
            print("Average not calculated.")


class AverageMonth:
    def __init__(self, reader):
        self.reader = reader
        self.average_month = None

    def calculate_average(self):
        if self.reader.data_json:
            data = pd.read_json(self.reader.data_json)
            monthly_data = data.iloc[:, 1:13]
            self.average_month = monthly_data.mean()
            return self.average_month
        else:
            return None

    def plot_average(self):
        if self.average_month is not None:
            if self.reader.data_json:
                data = pd.read_json(self.reader.data_json)
                years = data['Year']
                months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
                temperatures = self.average_month.values

                plt.plot(months, temperatures)
                plt.xlabel('Month')
                plt.ylabel('Temperature Anomaly')
                first_year = years.iloc[0]
                last_year = years.iloc[-1]
                title = f'Average Monthly Temperature Anomaly from year {first_year} to year {last_year}'
                plt.title(title)
                plt.show()
            else:
                print("No data available.")
        else:
            print("Average not calculated.")
            

class ReaderObserver:
    def __init__(self, reader):
        self.reader = reader

    def update(self):
         """
        Update method called when the observed Reader has new data.
        It creates instances of AverageYear and AverageMonth, calculates and plots the averages.
        """
        cons1 = AverageYear(self.reader)
        cons2 = AverageMonth(self.reader)
        cons1.calculate_average()
        cons1.plot_average()
        cons2.calculate_average()
        cons2.plot_average()


if __name__ == "__main__":
    prod = Reader('dSST.csv', CsvConverter('dSST.csv', header), 10)
    obs = ReaderObserver(prod)
    prod.add_observer(obs)
    prod.start_reading()