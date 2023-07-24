import socketserver
from http.server import SimpleHTTPRequestHandler as SimpleHandler
import json
import pandas as pd

PORT = 8080
FILENAME = '../../data/dSST.csv'

# Good realisation
# I just don't see the relation between this file and the dataHandler.py...
class DataProvider:
    '''
    The class is created to handle the data and generate the JSON response
    '''
    def __init__(self):
        FILENAME = '../../data/dSST.csv' # you already had this defined (7 lines up)
        self.data = pd.read_csv(FILENAME)

    def get_data(self, year_range=None):
        '''
        method takes an optional year_range parameter, which can be an integer representing a single year 
        or a list of two integers representing a range of years. It returns the filtered data as a JSON
        '''
        if year_range is None:
            return self.data.to_json(orient='records')

        if isinstance(year_range, int):
            filtered_data = self.data.loc[self.data['Year'] == year_range]
        elif isinstance(year_range, list) and len(year_range) == 2:
            filtered_data = self.data.loc[(self.data['Year'] >= year_range[0]) & (self.data['Year'] <= year_range[1])]
        else:
            raise ValueError("Invalid parameter")

        return filtered_data.to_json(orient='records')

class WeatherDataHandler(SimpleHandler):
    def __init__(self, *args, **kwargs):
        self.data_provider = DataProvider()
        super().__init__(*args, **kwargs)

    def do_GET(self):
        if self.path.startswith('/data'):
            self.handle_weather_data_request()
        else:
            self.send_error(404, 'Not found')

    def handle_weather_data_request(self):
        try:
            if self.path == '/data/all':
                json_data = self.data_provider.get_data()
            else:
                path_parts = self.path.split('/')
                if len(path_parts) == 4 and path_parts[2].isdigit() and path_parts[3].isdigit():
                    from_year, to_year = int(path_parts[2]), int(path_parts[3])
                    json_data = self.data_provider.get_data([from_year, to_year])
                elif len(path_parts) == 3 and path_parts[2].isdigit():
                    year = int(path_parts[2])
                    json_data = self.data_provider.get_data(year)
                else:
                    self.send_error(400, 'Bad Request')
                    return

            self.send_response(200)
            self.send_header('Content-type', 'application/json')
            self.end_headers()
            self.wfile.write(json_data.encode())
        except ValueError as e:
            self.send_error(400, str(e))
        except FileNotFoundError:
            self.send_error(404, 'Not found')


def main():
    with socketserver.TCPServer(("", PORT), WeatherDataHandler) as httpd:
        print("serving at port", PORT)
        httpd.serve_forever()


if __name__ == "__main__":
    main()