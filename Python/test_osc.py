import argparse
import numpy as np
import unispring as usp
from copy import deepcopy, copy
from pythonosc import dispatcher
from pythonosc import osc_server
from pythonosc import udp_client
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import json

def add_line(addrs, args, message):
    print(addrs)
    print(args)
    print(message)

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
    
    global_hash = {'buffer':{}}
    dispatcher = dispatcher.Dispatcher()
    dispatcher.map("/add_line", add_line, client, global_hash)

    
    server = osc_server.ThreadingOSCUDPServer(
        (args_server.ip, args_server.port), dispatcher)
    print("Serving on {}".format(server.server_address))
    print('waiting...')
    server.serve_forever()