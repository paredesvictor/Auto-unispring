import argparse
import unispring as usp
from pythonosc import dispatcher
from pythonosc import osc_server
from pythonosc import udp_client
from region import RegionPolygon
import numpy as np
from sklearn import mixture

def gauss2D(x, y, mx, my, sigx, sigy, theta):
    a = np.cos(theta)**2/(2*sigx**2) + np.sin(theta)**2/(2*sigy**2)
    b = -np.sin(2*theta)/(4*sigx**2) + np.sin(2*theta)/(4*sigy**2)
    c = np.sin(theta)**2/(2*sigx**2) + np.cos(theta)**2/(2*sigy**2)
    gauss = np.exp(- a*(x-mx)**2 - 2*b*(x-mx)*(y-my) - c*(y-my)**2)
    return gauss

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
    args[1]['remaining_lines'][buffer] -= 1
    args[1]['nb_lines'][buffer] += 1
    if args[1]['remaining_lines'][buffer] == 0:
        print('buffer', buffer, ',',args[1]['nb_lines'][buffer], 'grains' )
        args[0].send_message('/next_buffer', 1)
        end_test = [item<=0 for i,item in args[1]['remaining_lines'].items()]
        len_test = args[1]['nb_buffer'] == len(end_test)
        if all(end_test) and len_test:
            args[0].send_message('/done_import', 1)
    if args[1]['nb_lines'][buffer] % args[1]['osc_batch_size'] == 0:
        args[0].send_message('/next_batch', 1)

def add_buffer(addrs, args, *message):
    n_cols = int(message[2])
    n_rows = int(message[0])
    buffer = str(message[1])
    args[1]['buffer'][buffer] = [[0.0 for i in range(n_cols)] for j in range(n_rows)]
    args[1]['remaining_lines'][buffer] = n_rows
    args[1]['nb_lines'][buffer] = 0
    args[0].send_message('/start_dump', 1)

def import_init(addrs, args, *message):
    print('Export from Max...')
    args[1]['buffer'] = {}
    args[1]['osc_batch_size'] = int(message[1])
    args[1]['nb_buffer'] = int(message[0])
    args[1]['nb_lines'] = {}
    args[1]['remaining_lines'] = {}
    args[0].send_message('/begin_import', 1)

def eval_str(addrs, args, eval_string):
    eval(eval_string)

def create_norm_track(addrs, args, *unused):
    args[1]['norm_buffer'] = {}
    args[1]['norm_buffer'] = MinMaxScale(args[1]['buffer'])
    args[0].send_message('/done_create', 1)

def write_norm_track(addrs, args, *unused):
    for idx_buffer, track in args[1]['norm_buffer'].items():
        args[0].send_message('/buffer_index', int(idx_buffer))
        for i,line in enumerate(track):
            args[0].send_message('/append', line)
    args[0].send_message('/done_norm', 1) 
    args[0].send_message('/update', 'update')
    print('----- Done')

def init_unispring(addrs, args, *descr):
    print('Uniformization...')
    vertices = ((0,0),(1,0),(1,1),(0,1))
    region = RegionPolygon(vertices)
    args[1]['corpus'] = usp.Corpus(
        args[1]['norm_buffer'], 
        region, descr[0]+1, 
        descr[1]+1
    )
    print('e : ',args[1]['corpus'].uniform(client=args[0]))
    args[1]['corpus'].exportToMax(args[0])
    args[0].send_message('/update', 'update')
    print('----- Done')

def update_unispring(addrs, args, *coord):
    print('Update...')
    vertices = [(coord[i],1-coord[i+1]) for i in range(0,len(coord),2)]
    region = RegionPolygon(vertices)
    args[1]["corpus"].region = region
    print('e : ',args[1]["corpus"].unispring(
        exportPeriod=1, 
        client=args[0], 
        limit=200*(len(vertices)/4)
        ))
    args[1]["corpus"].exportToMax(args[0])
    args[0].send_message('/update', 'update')
    print('----- Done')

def adapt_unispring(addrs, args, *unused):
    print('Adapt...')
    e = args[1]["corpus"].unispring(
        exportPeriod=1, 
        client=args[0], 
        limit=500, 
        hDist='from_table', 
        hTable=args[1]['expl_record']
        )
    print('e : ',e)
    args[1]["corpus"].exportToMax(args[0])
    args[0].send_message('/update', 'update')
    print('----- Done')

def scale_unispring(addrs, args, *coeff):
    args[1]['corpus'].exportToMax(args[0], coeff = coeff)

def corner_unispring(addrs, args, corner):
    if corner == 1:
        coeff = (0.5, 0, 0.5, 0)
        invert = (0, 0)
    elif corner == 2:
        coeff = (0.5, 0.5, 0.5, 0)
        invert = (1, 0)
    elif corner == 3:
        coeff = (0.5, 0, 0.5, 0.5)
        invert = (0, 1)
    elif corner == 4:
        coeff = (0.5, 0.5, 0.5, 0.5)
        invert = (1, 1)
    args[1]['corpus'].exportToMax(args[0], coeff = coeff, invert = invert)

def add_expl_point(addrs, args, *coord):
    i = int(args[1]['expl_record'].sum())
    args[1]['expl_points'][i, :] = np.array((coord[0], coord[1]))
    n = args[1]['expl_record'].shape[0]
    xIdx = int(coord[0] * n)
    yIdx = int(coord[1] * n)
    if 0 <= xIdx < n and 0 <= yIdx < n:
        args[1]['expl_record'][xIdx, yIdx] += 1

def init_expl(addrs, args, size):
    sr = 100
    dur = 300
    args[1]['expl_record'] = np.zeros((size, size))
    args[1]['expl_points'] = np.zeros((sr * dur, 2))

def print_expl(addrs, args, *unused):
    print(args[1]['expl_record'])

def interpolation(addrs, args, interp):
    args[1]['corpus'].exportToMax(args[0], interp)

def gaussian_attract(addrs, args, *mess):
    mx = (mess[0], )
    my = (mess[1], )
    sigx = (mess[2], )
    sigy = (mess[3], )
    theta = (mess[4], )
    args[1]['corpus'].simpleAttractor(
        mx, my, sigx, sigy, theta, client=args[0]
        )

def ident_gaussian(addrs, args, *mess):
    data = args[1]['expl_points']
    n_comp = 3
    gmm = mixture.GaussianMixture(
        n_components = n_comp, 
        covariance_type = 'full',
        init_params = 'kmeans'
    )
    gmm.fit(data)
    mx = []
    my = []
    sigx = []
    sigy = []
    theta = []
    for i in range(n_comp):
        v, w = np.linalg.eigh(gmm.covariances_[i])
        mean = gmm.means_[i]
        angle = np.arctan2(w[0][1], w[0][0])
        v = 2.0 * np.sqrt(2.0) * np.sqrt(v)
        mx.append(mean[0])
        my.append(mean[1])
        sigx.append(v[0])
        sigy.append(v[1])
        theta.append(0-angle)
    for i in range(3):
        print(i+1, '.')
        print('mean: ', mx[i], my[i])
        print('sig : ', sigx[i], sigy[i])
        print('angle: ',180*theta[i]/np.pi)
    args[1]['corpus'].simpleAttractor(
        mx, my, sigx, sigy, theta, client=args[0]
        )


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

    dispatcher.map("/export_init", import_init, client, global_hash)
    dispatcher.map("/add_line", add_line, client, global_hash)
    dispatcher.map("/add_buffer", add_buffer, client, global_hash)
    dispatcher.map("/create_norm_track", create_norm_track, client, global_hash)
    dispatcher.map("/write_norm_track", write_norm_track, client, global_hash)

    dispatcher.map("/init_unispring", init_unispring, client, global_hash)
    dispatcher.map("/adapt_unispring", adapt_unispring, client, global_hash)
    dispatcher.map("/scale_unispring", scale_unispring, client, global_hash)
    dispatcher.map("/corner_unispring", corner_unispring, client, global_hash)
    dispatcher.map("/region", update_unispring, client, global_hash)

    dispatcher.map("/init_expl", init_expl, client, global_hash)
    dispatcher.map("/add_expl_point", add_expl_point, client, global_hash)
    dispatcher.map("/print_expl", print_expl, client, global_hash)
    dispatcher.map("/interpolation", interpolation, client, global_hash)
    dispatcher.map("/gaussian_attract", gaussian_attract, client, global_hash)
    dispatcher.map("/ident_gaussian", ident_gaussian, client, global_hash)

    dispatcher.map("/print", print)
    dispatcher.map("/eval", eval_str, client, global_hash)

    server = osc_server.ThreadingOSCUDPServer(
        (args_server.ip, args_server.port), dispatcher)

    print("----- Serving on {}".format(server.server_address))
    server.serve_forever()