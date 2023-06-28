import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class AverageYear:

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
        self.df['average'] = self.df[
            ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep',
             'Oct', 'Nov', 'Dec']].mean(axis=1)

    def plot(self):
        return self.df[['Year', 'Average']
                       ].plot(x='Year', ylabel='Average extremes (\u2103)')

    def all_in_one(self):
        df = pd.read_json(self.reader.get_lines())

        while df.iloc[-1].isna()[1] == False:
            df = pd.concat([df, pd.read_json(self.reader.get_lines())],
                           ignore_index=True)

            df['average'] = df[['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 
                                'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'
                                ]].mean(axis=1)

        plot = df[['Year', 'Average']
                  ].plot(x='Year', ylabel='Average extremes (\u2103)')

        plt.show()
