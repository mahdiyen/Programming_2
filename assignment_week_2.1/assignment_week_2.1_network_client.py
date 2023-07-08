import asyncio
import aiohttp

class NetworkClient:
    def __init__(self, base_url):
        self.base_url = base_url

    async def fetch_data(self, endpoint):
        async with aiohttp.ClientSession() as session:
            async with session.get(self.base_url + endpoint) as response:
                if response.status == 200:
                    return await response.json()
                else:
                    return None

    async def process_data(self, endpoint, callback):
        data = await self.fetch_data(endpoint)
        if data:
            return await callback(data)
        else:
            return None


async def calculate_average_temperature(data):
    # Example function to calculate the average temperature from the data
    temperatures = [float(entry['J-D']) for entry in data]
    return {'average_temperature': sum(temperatures) / len(temperatures)}

async def print_data(data):
    # Example function to print the received data
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
