import argparse
import math
import random
#import unispring as usp
from pythonosc import dispatcher
from pythonosc import osc_server
from pythonosc import udp_client


def callBackFunction(addrs, args, message):
    if random.random() < 0.5:
        args[0].send_message("/update", 'done')
    else:
        args[0].send_message("/update", 'failed')
        


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
    dispatcher.map("/message", print)
    dispatcher.map("/message", callBackFunction, client)
    
    
    server = osc_server.ThreadingOSCUDPServer(
        (args_server.ip, args_server.port), dispatcher)
    print("Serving on {}".format(server.server_address))
    server.serve_forever()
    