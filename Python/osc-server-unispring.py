
import argparse
import math
import random
import unispring as usp
from pythonosc import dispatcher
from pythonosc import osc_server
from pythonosc import udp_client

def apply_unispring(addrs, message):
    if message == "apply":
        descX = "CentroidMean"
        descY = "PeriodicityMean"
        region = RegionPolygon(vertices)
        corpus = Corpus('/Users/victorparedes/Documents/GitHub/Auto-unispring/corpus.json',
        region, descX, descY, plot=False)
        corpus.preUniformization(inSquareAuto=False)
        corpus.unispringUniform(1, 0.02, 0.01)
        corpus.exportJson('/Users/victorparedes/Documents/GitHub/Auto-unispring/remap.json')

def init_unispring(addrs, message):
    if message == "apply":
        print('Creating corpus and region')
        vertices = ((0,0),(1,0),(1,1),(0,1))
        descX = "CentroidMean"
        descY = "PeriodicityMean"
        region = usp.RegionPolygon(vertices)
        corpus = usp.Corpus('/Users/victorparedes/Documents/GitHub/Auto-unispring/corpus.json',
        region, descX, descY, plot=False)
        print('preUniformization')
        corpus.preUniformization(inSquareAuto=False)
        print('uniformization')
        corpus.unispringUniform(1, 0.02, 0.01)
        print('export')
        corpus.exportJson('/Users/victorparedes/Documents/GitHub/Auto-unispring/remap.json')

def new_region(addrs, *coord):
    vertices = [(coord[i],coord[i+1]) for i in range(0,len(coord),2)]
    region = usp.RegionPolygon(vertices)
    return region

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
    
    dispatcher = dispatcher.Dispatcher()
    dispatcher.map("/region", new_region)
    dispatcher.map("/unispring", init_unispring)
    
    
    server = osc_server.ThreadingOSCUDPServer(
        (args_server.ip, args_server.port), dispatcher)
    print("Serving on {}".format(server.server_address))
    server.serve_forever()