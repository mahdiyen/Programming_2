import urllib.request, urllib.parse, urllib.error
from bs4 import BeautifulSoup
import ssl
import re

# I am missing the `__iter__` of `__next__` in this class (to change it into an iterable)
class Crawler:
    def __init__(self):
        self.url = "https://sport050.nl/sportaanbieders/alle-aanbieders/"

    # Good observeration that this method could be made static
    @staticmethod
    def hack_ssl():
        """ ignores the certificate errors"""
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        return ctx

    def open_url(self, url):
        """ reads url file as a big string and cleans the html file to make it
            more readable. input: url, output: soup object
        """
        ctx = self.hack_ssl()
        html = urllib.request.urlopen(url, context=ctx).read()
        soup = BeautifulSoup(html, 'html.parser')
        return soup

    def read_hrefs(self, soup):
        """ get from soup object a list of anchor tags,
            get the href keys and and prints them. Input: soup object
        """
        reflist = []
        tags = soup('a')
        for tag in tags:
            reflist.append(tag)
        return reflist

    def read_li(self, soup):
        lilist = []
        tags = soup('li')
        for tag in tags:
            lilist.append(tag)
        return lilist

    def get_phone(self, info):
        reg = r"(?:(?:00|\+)?[0-9]{4})?(?:[ .-][0-9]{3}){1,5}"
        phone = []
        for s in info:
            if 'Telefoon' in str(s):
                phone.append(str(s))
        try:
            phone = str(phone[0])
        except:
            phone = []
            for s in info:
                if re.findall(reg, str(s)):
                    phone.append(str(s))
            try:
                phone = str(phone[0])
            except:
                phone = ""   
        return phone.replace('Facebook', '').replace('Telefoon:', '')

    def get_email(self, soup):
        try: 
            email = []
            for s in soup:
                if '@' in str(s):
                    email.append(str(s))
            email = str(email[0])[4:-5]
            bs = BeautifulSoup(email, features="html.parser")
            email = bs.find('a').attrs['href'].replace('mailto:', '')
        except:
            email = ""
        return email

    def remove_html_tags(self, text):
        """Remove html tags from a string"""
        clean = re.compile('<.*?>')
        return re.sub(clean, '', text)

    def fetch_sidebar(self, soup):
        """ reads html file as a big string and cleans the html file to make it
            more readable. input: html, output: tables
        """
        sidebar = soup.findAll(attrs={'class': 'sidebar'})
        return sidebar[0]

    def extract(self, url):
        text = str(url)
        text = text[26:].split('"')[0] + "/"
        return text

    def crawl_site(self):
        print('fetch urls')
        s = self.open_url(self.url)
        reflist = self.read_hrefs(s)

        print('getting sub-urls')
        sub_urls = [s for s in reflist if '<a href="/sportaanbieders' in str(s)]
        sub_urls = sub_urls[3:]

        print('extracting the data')
        print(f'{len(sub_urls)} sub-urls')

        for sub in sub_urls:
            try:
                sub = self.extract(sub)
                site = self.url[:-16] + sub
                soup = self.open_url(site)    
                info = self.fetch_sidebar(soup)
                info = self.read_li(info)
                phone = self.get_phone(info)
                phone = self.remove_html_tags(phone).strip()
                email = self.get_email(info)
                email = self.remove_html_tags(email).replace("/", "")
                print(f'{site} ; {phone} ; {email}')
            except Exception as e:
                print(e)
                exit()

if __name__ == '__main__':
    crawler = Crawler()
    crawler.crawl_site()
