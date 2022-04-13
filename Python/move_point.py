from pythonosc import dispatcher
from pythonosc import osc_server
from pythonosc import udp_client
import argparse
import numpy as np

def normalize(array):
    return 

def add_line(addrs, args, *message):
    index = int(message[-2])
    buffer = str(message[-1])
    n_descr = len(message)-3
    descriptors = [message[i+1] for i in range(n_descr)]
    args[1][buffer][index] = descriptors
    #args[0].send_message('/callback', truc)

def add_buffer(addrs, args, *message):
    print(message)
    n_cols = int(message[1])
    n_rows = int(message[0])
    buffer = str(message[2])
    args[1][buffer] = [[] for j in range(n_rows)]

def show(addrs, args, key):
    print(args[1][str(key)])

def eval_str(addrs, args, eval_string):
    eval(eval_string)

def reset_track(addrs, args, *unused):
    args[1]['from_max'] = []

def import_track(addrs, args, *unused):
    track = [np.asarray(buffer) for buffer in args[1]['from_max']]
    args[1]['track'] = [normalize(buf) for buf in track]
    track0 = args[1]['track'][0]
    print(track0[0,0])
    track1 = args[1]['track'][1]
    print(track1[0,0])

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
    
    global_hash = {'from_max':[]}
    dispatcher = dispatcher.Dispatcher()
    dispatcher.map("/add_line", add_line, client, global_hash)
    dispatcher.map("/add_buffer", add_buffer, client, global_hash)
    dispatcher.map("/import_track", import_track, client, global_hash)
    dispatcher.map("/print", print)
    dispatcher.map("/show", show, client, global_hash)
    dispatcher.map("/eval", eval_str, client, global_hash)
    dispatcher.map("/reset_track", reset_track, client, global_hash)
    
    server = osc_server.ThreadingOSCUDPServer(
        (args_server.ip, args_server.port), dispatcher)
    print("Serving on {}".format(server.server_address))
    print('waiting...')
    server.serve_forever()