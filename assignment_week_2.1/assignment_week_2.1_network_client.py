import asyncio
import aiohttp
# Good

class NetworkClient:
    def __init__(self, base_url):
        self.base_url = base_url

    async def fetch_data(self, endpoint):
        """
        Fetches data from the given endpoint by sending an HTTP GET request.
        Returns the response data as JSON if the status is 200, otherwise returns None.
        """
        async with aiohttp.ClientSession() as session:
            async with session.get(self.base_url + endpoint) as response:
                if response.status == 200:
                    return await response.json()
                else:
                    return None

    async def process_data(self, endpoint, callback):
        """
        Fetches data from the given endpoint using fetch_data() method,
        and then calls the provided callback function on the retrieved data.
        Returns the result of the callback function.
        """
        data = await self.fetch_data(endpoint)
        if data:
            return await callback(data)
        else:
            return None


async def calculate_average_temperature(data):
    '''
    calculates the average temperature from the data.
    Receives data as input and returns a dictionary with the average temperature.
    '''
    temperatures = [float(entry['J-D']) for entry in data]
    return {'average_temperature': sum(temperatures) / len(temperatures)}

async def print_data(data):
    """
    Receives data as input and prints it.
    """    
    print(data)

async def main():

    base_url = 'http://localhost:8080/data/'
    client = NetworkClient(base_url)

    tasks = [
        client.process_data('all', calculate_average_temperature),
        client.process_data('1991', print_data),
        client.process_data('1991/2000', print_data),
    ]

    await asyncio.gather(*tasks)

asyncio.run(main())
