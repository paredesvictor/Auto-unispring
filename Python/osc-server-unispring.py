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


def MinMaxScale(track):
    n_descr = len(track['1'][0])
    norm_track = {}
    list_min = [float('inf') for i in range(n_descr)]
    list_max = [float('-inf') for i in range(n_descr)]
    for key, table in track.items():
        for line in table:
            for i in range(1,n_descr):
                if line[i] < list_min[i]:
                    list_min[i] = line[i]
                elif line[i] > list_max[i]:
                    list_max[i] = line[i]
    for key, table in track.items():
        norm_track[key] = []
        for line in table:
            new_line = [line[0]]
            for i in range(1,n_descr):
                new_line.append((line[i]-list_min[i])/(list_max[i]-list_min[i]))
            norm_track[key].append(new_line)
    return norm_track


def add_line(addrs, args, *message):
    index = int(message[-2])
    buffer = str(message[-1])
    n_descr = len(message)-2
    descriptors = [message[i] for i in range(n_descr)]
    args[1]['buffer'][buffer][index] = descriptors


def add_buffer(addrs, args, *message):
    n_cols = int(message[2])
    n_rows = int(message[0])
    buffer = str(message[1])
    args[1]['buffer'][buffer] = [[0.0 for i in range(n_cols)] for j in range(n_rows)]


def eval_str(addrs, args, eval_string):
    eval(eval_string)


def reset_track(addrs, args, *unused):
    args[1]['buffer'] = {}


def dumpdone(addrs, args, buffer):
    n_lines = len(args[1]['buffer'][str(buffer)])
    print('done importing buffer ', buffer, ', ',n_lines, ' grains.' )


def create_norm_track(addrs, args, *unused):
    args[1]['norm_buffer'] = {}
    args[1]['norm_buffer'] = MinMaxScale(args[1]['buffer'])
    print('done normalizing')


def write_norm_track(addrs, args, *unused):
    for idx_buffer, track in args[1]['norm_buffer'].items():
        args[0].send_message('/buffer_index', int(idx_buffer))
        for i,line in enumerate(track):
            args[0].send_message('/append', line)
    print('done exporting normalized buffers')


def init_unispring(addrs, args, *descr):
    vertices = ((0,0),(1,0),(1,1),(0,1))
    region = usp.RegionPolygon(vertices)
    args[1]['corpus'] = usp.Corpus(args[1]['norm_buffer'], region, descr[0], descr[1])
    args[1]['corpus'].unispringUniform(1, 0.01, 0.02, exportPeriod=5, client=args[0])
    print('uniformization done')
    args[1]['corpus'].exportToMax(args[0])


def update_unispring(addrs, args, *coord):
    print('updating...')
    temp_corpus = deepcopy(args[1]["corpus"])
    vertices = [(coord[i],1-coord[i+1]) for i in range(0,len(coord),2)]
    print(vertices)
    region = usp.RegionPolygon(vertices)
    temp_corpus.region = region
    temp_corpus.unispringUniform(1, 0.01, 0.02, exportPeriod=5, client=args[0], limit=200*(len(vertices)/4))
    print('done updating')
    temp_corpus.exportToMax(args[0])


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
    dispatcher.map("/add_buffer", add_buffer, client, global_hash)
    dispatcher.map("/create_norm_track", create_norm_track, client, global_hash)
    dispatcher.map("/write_norm_track", write_norm_track, client, global_hash)
    dispatcher.map("/reset_track", reset_track, client, global_hash)
    dispatcher.map("/dumpdone", dumpdone, client, global_hash)
    dispatcher.map("/init_unispring", init_unispring, client, global_hash)
    dispatcher.map("/region", update_unispring, client, global_hash)
    dispatcher.map("/tsne", tSNE, client, global_hash)
    dispatcher.map("/print", print)
    dispatcher.map("/eval", eval_str, client, global_hash)
    
    server = osc_server.ThreadingOSCUDPServer(
        (args_server.ip, args_server.port), dispatcher)
    print("Serving on {}".format(server.server_address))
    print('waiting...')
    server.serve_forever()