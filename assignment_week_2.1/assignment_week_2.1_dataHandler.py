import socketserver
from http.server import SimpleHTTPRequestHandler as SimpleHandler
import json
import pandas as pd

PORT = 8080
FILENAME = '../../data/dSST.csv'


class WeatherDataHandler(SimpleHandler):
    def do_GET(self):
        if self.path.startswith('/data'):
            self.handle_weather_data_request()
        else:
            self.send_error(404, 'Not found')

    def handle_weather_data_request(self):
        if self.path == '/data/all':
            self.send_weather_data()
        else:
            path_parts = self.path.split('/')
            if len(path_parts) == 4 and path_parts[2].isdigit() and path_parts[3].isdigit():
                self.send_weather_data_range(int(path_parts[2]), int(path_parts[3]))
            elif len(path_parts) == 3 and path_parts[2].isdigit():
                self.send_weather_data_year(int(path_parts[2]))
            else:
                self.send_error(400, 'Bad Request')

    def send_weather_data(self):
        try:
            with open(FILENAME, 'r') as file:
                header = file.readline().strip().split(',')
                lines = file.readlines()

            json_data = []
            for line in lines:
                values = line.strip().split(',')
                json_dict = {key: value for key, value in zip(header, values)}
                json_data.append(json_dict)

            self.send_response(200)
            self.send_header('Content-type', 'application/json')
            self.end_headers()
            self.wfile.write(json.dumps(json_data).encode())
        except FileNotFoundError:
            self.send_error(404, 'Not found')

    def send_weather_data_range(self, from_year, to_year):
        try:
            with open(FILENAME, 'r') as file:
                header = file.readline().strip().split(',')
                year_index = header.index('Year')
                lines = file.readlines()

            filtered_data = []
            for line in lines:
                values = line.strip().split(',')
                year = int(values[year_index])
                if from_year <= year <= to_year:
                    filtered_data.append(values)

            if not filtered_data:
                self.send_error(404, 'Not found')
                return

            json_data = []
            for values in filtered_data:
                json_dict = {key: value for key, value in zip(header, values)}
                json_data.append(json_dict)

            self.send_response(200)
            self.send_header('Content-type', 'application/json')
            self.end_headers()
            self.wfile.write(json.dumps(json_data).encode())
        except FileNotFoundError:
            self.send_error(404, 'Not found')

    def send_weather_data_year(self, year):
        try:
            with open(FILENAME, 'r') as file:
                header = file.readline().strip().split(',')
                year_index = header.index('Year')

                for line in file:
                    values = line.strip().split(',')
                    if int(values[year_index]) == year:
                        json_dict = {key: value for key, value in zip(header, values)}
                        self.send_response(200)
                        self.send_header('Content-type', 'application/json')
                        self.end_headers()
                        self.wfile.write(json.dumps(json_dict).encode())
                        return

            self.send_error(404, 'Not found')
        except FileNotFoundError:
            self.send_error(404, 'Not found')


def main():
    with socketserver.TCPServer(("", PORT), WeatherDataHandler) as httpd:
        print("serving at port", PORT)
        httpd.serve_forever()


if __name__ == "__main__":
    main()
