import argparse
import math
import random
import unispring as usp
from pythonosc import dispatcher
from pythonosc import osc_server
from pythonosc import udp_client

def update_unispring(addrs, args, coord):
    print('updating unispring...')
    vertices = [(coord[i],coord[i+1]) for i in range(0,len(coord),2)]
    region = usp.RegionPolygon(vertices)
    args[1]["corpus1"].region = region
    args[1]["corpus1"].unispringUniform(1, 0.01, 0.02)
    print('export')
    corpus.exportJson(dir + '/remap.json')
    args[1]['corpus1'] = corpus
    args[0].send_message("/unispring", 'done')

def init_unispring(addrs, args, dir):
    print('Creating corpus and region')
    vertices = ((0,0),(1,0),(1,1),(0,1))
    descX = "CentroidMean"
    descY = "PeriodicityMean"
    region = usp.RegionPolygon(vertices)
    corpus = usp.Corpus(dir + '/corpus.json',
    region, descX, descY, plot=False)
    print('preUniformization')
    corpus.preUniformization(inSquareAuto=False)
    print('uniformization...')
    corpus.unispringUniform(1, 0.01, 0.02)
    print('export')
    corpus.exportJson(dir + '/remap.json')
    args[1]['corpus1'] = corpus
    args[0].send_message("/unispring", 'done')

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
    server.serve_forever()