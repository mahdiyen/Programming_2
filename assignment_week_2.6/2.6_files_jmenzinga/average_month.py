import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


class AverageMonth:

    def __init__(self, reader) -> None:
        self.reader = reader

    def get_data(self):
        df = pd.read_json(self.reader.get_lines())

        while df.iloc[-1].isna()[1] == False:
            df = pd.concat([df, pd.read_json(self.reader.get_lines())],
                           ignore_index=True)

        self.df = df
        return df

    def get_averages(self):
        # prolly need to change this to inster. maybe combine it with get_data
        column = 0
        avg_month = pd.DataFrame()

        while column < len(self.df):
            avg_month[str(column)] = self.df[
                ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 0'Aug', 'Sep',
                 'Oct', 'Nov', 'Dec']
            ].iloc[column:column+5].mean(axis=0)

            column += 5

        self.avg_month = avg_month

    def plot(self):
        return self.avg_month.plot()
