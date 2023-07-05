import socketserver
from http.server import SimpleHTTPRequestHandler as SimpleHandler
port = 9000
socketserver.TCPServer.allow_reuse_address = True
http = socketserver.TCPServer(('localhost',port), SimpleHandler)
print (f'serving on port{port}')
http.serve_forever()