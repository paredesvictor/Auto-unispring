from pythonosc import dispatcher
from pythonosc import osc_server
from pythonosc import udp_client
import argparse
import numpy as np

def normalize(table):
    n_descr = len(table[0])
    normalized_table = []
    list_min = [float('inf') for i in range(n_descr)]
    list_max = [float('-inf') for i in range(n_descr)]
    for line in table:
        for i in range(1,n_descr):
            if line[i] < list_min[i]:
                list_min[i] = line[i]
            elif line[i] > list_max[i]:
                list_max[i] = line[i]
    for line in table:
        new_line = [line[0]]
        for i in range(1,n_descr):
            new_line.append((line[i]-list_min[i])/(list_max[i]-list_min[i]))
        normalized_table.append(new_line)
    return normalized_table

def add_line(addrs, args, *message):
    index = int(message[-2])
    buffer = str(message[-1])
    n_descr = len(message)-2
    descriptors = [message[i] for i in range(n_descr)]
    args[1]['buffer'][buffer][index] = descriptors
    #args[0].send_message('/callback', truc)

def add_buffer(addrs, args, *message):
    n_cols = int(message[2])
    n_rows = int(message[0])
    buffer = str(message[1])
    args[1]['buffer'][buffer] = [[] for j in range(n_rows)]

def eval_str(addrs, args, eval_string):
    eval(eval_string)

def reset_track(addrs, args, *unused):
    args[1]['buffer'] = {}

def dumpdone(addrs, args, buffer):
    n_lines = len(args[1]['buffer'][str(buffer)])
    print('done importing buffer ', buffer, ', ',n_lines, ' grains.' )

def create_norm_track(addrs, args, *unused):
    args[1]['norm_buffer'] = {}
    for key, buffer in args[1]['buffer'].items():
        args[1]['norm_buffer'][key] = normalize(buffer)
    print('done normalizing')

def write_norm_track(addrs, args, *unused):
    tracks = args[1]['norm_buffer']
    for idx_buffer, track in tracks.items():
        args[0].send_message('/buffer_index', int(idx_buffer))
        for i,line in enumerate(track):
            args[0].send_message('/append', line)

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
    dispatcher.map("/print", print)
    dispatcher.map("/eval", eval_str, client, global_hash)
    
    server = osc_server.ThreadingOSCUDPServer(
        (args_server.ip, args_server.port), dispatcher)
    print("Serving on {}".format(server.server_address))
    print('waiting...')
    server.serve_forever()