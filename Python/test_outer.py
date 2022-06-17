from scipy.spatial import ConvexHull
from unispring import Point
from numpy import mean, asarray, sqrt, pi, cos, sin, tan, sign
import matplotlib.pyplot as plt
import random as rd

#string = "0.533341 0.443199 0.461124 0.333393 0.235369 0.341284 0.184073 0.449201 0.366361 0.592812 0.381871 0.456201 0.355112 0.521218 0.321692 0.355933 0.335829 0.499081 0.501694 0.516926 0.224963 0.422581 0.016166 0.138031 0.084873 0.156853 0.276466 0.405155 0.32669 0.614651 0.250342 0.602289 0.171687 0.403733 0.476653 0.394199 0.330887 0.36904 0.445647 0.171542 0.293811 0.297335 0.10004 0.242068 0.101647 0.131399 0.33216 0.678671 0.426982 0.535679 0.38183 0.664461 0.150601 0.51489 0.164425 0.242656 0.35409 0.335979 0.372314 0.528419 0.201434 0.07032 0.56222 0.533036 0.736877 0.518619 0.290352 0.405068 0.1674 0.167554 0.269122 0.18855 0.559746 0.578843 0.704666 0.505914 0.932668 0.548123 0.670066 0.626163 0.741748 0.532013 0.834602 0.875518 0.431553 0.433747 0.345441 0.38521 0.419352 0.470731 0.669068 0.552542 0.479218 0.393032 0.585085 0.245646 0.579012 0.267689 0.327875 0.597135 0.125261 0. 0.377685 0.602101 0.767582 1. 0.648253 0.550314 0.694302 0.693364 0.845912 0.68392 0.64453 0.519549 0.801108 0.73224 0.434903 0.561222 0.332164 0.575129 0.083547 0.217152 0.514927 0.45382 0.697068 0.72071 0.453567 0.477648 0.396405 0.611031 0.477611 0.542573 0.438698 0.328845 0.52628 0.591951 0.45265 0.538925 0.686413 0.453525 0.281086 0.187226 0.579931 0.838417 0.342517 0.192873 0.41947 0.313064 0.386622 0.628791 0.424354 0.640153 0.231787 0.312107 0.772592 0.957769 0.784949 0.354123 0.226804 0.232295 0.579663 0.405865 0.417996 0.419061 0.3372 0.309536 0.573854 0.435242 0.52705 0.61138 1. 0.68059 0.827221 0.509968 0.224084 0.255841 0.245053 0.110856 0.611013 0.940906 0.291219 0.500104 0.248476 0.621241 0.25316 0.336038 0.161209 0.275135 0.312219 0.27193 0.250153 0.193189 0.318415 0.197308 0.496073 0.541758 0.122787 0.507507 0.255858 0.495058 0.223762 0.397947 0.415263 0.394128 0.225362 0.27094 0.20997 0.441178 0.269653 0.431944 0.414042 0.557085 0.105695 0.266894 0.332031 0.399626 0.22243 0.275773 0.221896 0.166882 0.154218 0.285558 0.170824 0.419713 0. 0.149428 0.142284 0.329449 0.10978 0.496257 0.290802 0.724902 0.507622 0.614905 0.299456 0.59643 0.134769 0.289035 0.169 0.268101 0.276389 0.480706 0.215989 0.360894 0.027138 0.113574 0.117958 0.284064 0.123688 0.28819 0.193348 0.324845 0.071452 0.165619 0.162508 0.304272"
#coord = [float(s) for s in string.split(' ')]
N=200
coord = [rd.random() for i in range(2*N)]
x = coord[0:len(coord):2]
y = coord[1:len(coord)+1:2]
pts = [Point(x[i], y[i]) for i in range(len(x))]

points = asarray([(p.x,p.y) for p in pts])
convex_hull = ConvexHull(points)   
#centroid = mean(points[convex_hull.vertices, :], axis=0)
centroid = mean(points, axis=0)
p_cent = Point(*centroid)
n_bord = int(1 + sqrt(len(pts))) * 4
n_bin = int(n_bord / 2)
angle_bin = 2 * pi / n_bin
sorted_pts = sorted(pts, key=lambda p : Point.vecOrientation(p, p_cent))
bins = [[p for p in sorted_pts if -pi+i*angle_bin <= p.vecOrientation(p_cent) < -pi+(i+1)*angle_bin] 
            for i in range(n_bin)]
hull = []
for bin in bins:
    if len(bin) == 0:
        continue
    elif len(bin) == 1:
        p1 = bin[0]
        hull.append(p1)
    else:
        p1 = bin[0]
        p2 = bin[1]
        for p in bin[2:]:
            if p.distTo(p_cent) > p1.distTo(p_cent) or p.distTo(p_cent) > p2.distTo(p_cent):
                p1 = max([p1,p2], key=lambda p:Point.distTo(p,p_cent))
                p2 = p
        hull.append([p1, p2])

n = 29
x = []
y = []
for i in range(n):
    theta = i * 2 * pi / n
    if theta == 0 or theta == pi:
        x = cos(theta)
    else:
        x = sign(cos(theta))*min(1,abs(1/tan(theta)))
    if theta == pi/2 or theta == 3 * pi / 2:
        y = sin(theta)
    else:
        y = sign(sin(theta))*min(1,abs(tan(theta)))
# plt.plot(x,y,'.')
# for i in range(n_bin):
#     xp = [centroid[0], 0.5*cos(angle_bin*i)+centroid[0]]
#     yp = [centroid[1], 0.5*sin(angle_bin*i)+centroid[1]]
#     plt.plot(xp,yp,color='black')
# for i in range(len(hull)-1):
#     plt.plot((hull[i].x,hull[i+1].x),(hull[i].y,hull[i+1].y),color='orange')
# plt.plot((hull[-1].x,hull[0].x),(hull[-1].y,hull[0].y),color='orange')
# plt.show()