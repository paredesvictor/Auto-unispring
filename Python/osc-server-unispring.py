import argparse
import sys
import os
from copy import deepcopy
import unispring as usp
from random import random
from pythonosc import dispatcher
from pythonosc import osc_server
from pythonosc import udp_client

class recursion_depth:
    def __init__(self, limit):
        self.limit = limit
        self.default_limit = sys.getrecursionlimit()
    def __enter__(self):
        sys.setrecursionlimit(self.limit)
    def __exit__(self, type, value, traceback):
        sys.setrecursionlimit(self.default_limit)

def update_unispring(addrs, args, *coord):
    print('updating...')
    n = 1
    while True:
        try:
            with recursion_depth(1000*n):
                temp_corpus = deepcopy(args[1]["corpus1"])
            break
        except:
            n += 1
    vertices = [(coord[i],1-coord[i+1]) for i in range(0,len(coord),2)]
    region = usp.RegionPolygon(vertices)
    temp_corpus.region = region
    temp_corpus.unispringUniform(1, 0.01, 0.02, limit=200*(len(vertices)/4))
    print('export')
    save_dir = args[1]['dir']
    directory = os.path.dirname(os.path.realpath(__file__))
    temp_corpus.exportJson(directory+'/remap.json')
    args[0].send_message("/unispring", save_dir)
    print('waiting...')

def init_unispring(addrs, args, max_path):
    directory = os.path.dirname(os.path.realpath(__file__))
    print('Creating corpus and region')
    vertices = ((0,0),(1,0),(1,1),(0,1))
    descX = "CentroidMean"
    descY = "PeriodicityMean"
    region = usp.RegionPolygon(vertices)
    corpus = usp.Corpus(directory + '/corpus.json',
    region, descX, descY, plot=False)
    print('uniformization...')
    limit = min(1000, max(200, len(corpus.getAllPoints())))
    corpus.unispringUniform(1, 0.01, 0.02, limit=limit)
    print('export')
    corpus.exportJson(directory + '/remap.json')
    args[1]['corpus1'] = corpus
    args[1]['dir'] = directory
    args[0].send_message("/unispring", directory)
    print('waiting...')

if __name__ == "__main__":
    parser_client = argparse.ArgumentParser()
    parser_client.add_argument("--ip", default="127.0.0.1")
    parser_client.add_argument("--port", type=int, default=8012)
    args_client = parser_client.parse_args()

    client = udp_client.SimpleUDPClient(args_client.ip, args_client.port)
    
    
    parser_server = argparse.ArgumentParser()
    parser_server.add_argument("--ip", default="127.0.0.1")
    parser_server.add_argument("--port", type=int, default=8011)
    args_server = parser_server.parse_args()
    
    corpus = {}
    dispatcher = dispatcher.Dispatcher()
    dispatcher.map("/region", update_unispring, client, corpus)
    dispatcher.map("/unispring", init_unispring, client, corpus)
    
    server = osc_server.ThreadingOSCUDPServer(
        (args_server.ip, args_server.port), dispatcher)
    print("Serving on {}".format(server.server_address))
    print('waiting...')
    server.serve_forever()