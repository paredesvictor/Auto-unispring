from scipy.spatial import Delaunay
from numpy import arctan2, sqrt, sin, cos, asarray, clip, exp
from scipy.stats import multivariate_normal as norm
from math import ceil

def hFuncFromTable(x, y, table):
    n = table.shape[0]
    xIdx = int(clip(x * n, 0, n - 1))
    yIdx = int(clip(y * n, 0, n - 1))
    val = table[xIdx, yIdx]
    return val / table.max() + 0.25

def hFunction(type, table):
    if type == 'uniform':
        return lambda x, y : 1
    elif type == 'gaussian':
        mx = 0.5
        my = 0.5
        s = 0.02
        f_norm = norm([mx,my],[[s,0],[0,s]])
        p = f_norm.pdf([0.5, 0.5])
        return lambda x, y : 1 + f_norm.pdf([x, y])/p
    elif type == 'progressive':
        return lambda x, y : 1 + x
    elif type == 'from_table':
        return lambda x, y : hFuncFromTable(x, y, table)


class Corpus():
    
    def __init__(self, track, region, descrX, descrY):
        # import json
        # descriptors extraction
        self.buffers = []
        for key,buffer in track.items():
            self.buffers.append(Buffer(buffer, descrX, descrY, int(key)))
        # initialize normalize bool
        self.is_norm = False
        self.region = region
        self.normalize()
        self.preUniformization()
    
    def getAllPoints(self):
        allPoints = []
        for buffer in self.buffers:
            allPoints += buffer.points
        return allPoints

    def getScalingFactor(self, l0):
        targetArea = 0
        nPair = 0
        for point in self.getAllPoints():
            for near in point.near:
                nPair += 1
                midX ,midY = point.midTo(near)
                targetArea += 1/self.hDist(midX, midY)**2
        return sqrt(nPair*l0**2/targetArea)

    def normalize(self):
        pointsX = []
        pointsY = []
        for buffer in self.buffers:
            pointsX += [point.x for point in buffer.points]
            pointsY += [point.y for point in buffer.points]
        upperX = max(pointsX)
        lowerX = min(pointsX)
        upperY = max(pointsY)
        lowerY = min(pointsY)
        for buffer in self.buffers:
            buffer.normalize(upperX, lowerX, upperY, lowerY)
        self.is_norm = True
            
    def preUniformization(self, resize=False, og=(0,0), s = 1, inSquareAuto=False):
        if inSquareAuto:
            p1, p2 = self.region.getBoundingBox()
            sideX, sideY, origin = p2.x-p1.x, p2.y-p1.y, p1
        else:
            sideX, sideY, origin = s, s, Point(og[0],og[1])
        allPoints = self.getAllPoints()
        nbPoints = len(allPoints)
        if not resize:
            allPoints.sort(key=Point.getX)
            for i in range(nbPoints):
                allPoints[i].x = (i / (nbPoints - 1)) * sideX + origin.x
            allPoints.sort(key=Point.getY)
            for i in range(nbPoints):
                allPoints[i].y = (i / (nbPoints - 1)) * sideY + origin.y
        else:
            for i in range(nbPoints):
                allPoints[i].x = allPoints[i].x * sideX + origin.x
            for i in range(nbPoints):
                allPoints[i].y = allPoints[i].y * sideY + origin.y
        
    def delaunayTriangulation(self):
        allCoord = []
        for buffer in self.buffers:
            allCoord += [[pt.x, pt.y] for pt in buffer.points]
        allCoord = asarray(allCoord)
        triangulation = Delaunay(allCoord)
        self.updateNearPoints(triangulation)
        return triangulation
    
    def updateNearPoints(self, triangulation):
        allPoints = self.getAllPoints()
        for point in allPoints:
            point.resetNear()
            point.updateOrigin()
        for tri in triangulation.simplices:
            p1 = allPoints[tri[0]]
            p2 = allPoints[tri[1]]
            p3 = allPoints[tri[2]]
            if p1 not in p2.near:
                p2.near.append(p1)
                p1.near.append(p2)
            if p1 not in p3.near:
                p3.near.append(p1)
                p1.near.append(p3)
            if p2 not in p3.near:
                p3.near.append(p2)
                p2.near.append(p3)
    
    def uniform(self, client=None, store=True):
        c1, c2 = self.unispring(uniform=True, exportPeriod=int(bool(client)), client=client)
        if store:
            for p in self.getAllPoints():
                p.storeUni()
        return c1, c2
        
    def unispring(self, uniform=False ,exportPeriod=0, client=None, limit=0, hDist='uniform', hTable=None):
        allPoints = self.getAllPoints()
        # change hDist function
        self.hDist = hFunction(hDist, hTable)
        # pre-uniformization
        if uniform:
            self.preUniformization(resize=False, inSquareAuto=False)
        else:
            self.preUniformization(resize=True, inSquareAuto=False)
        # first triangulation
        self.delaunayTriangulation()
        # l0 is calculated for a uniform target distribution with all 
        # delaunay triangles becoming equilateral, does not account for points 
        # on the borders
        nbPoints = len(allPoints)
        uniform_density = nbPoints / self.region.getArea()
        l0 = sqrt(2/(sqrt(3)*uniform_density))
        print('l0 : ',l0)
        exit = False
        count = 0
        c = 0
        while not exit:
            exit = True
            count += 1
            updateTri = False
            hScale = self.getScalingFactor(l0)
            #  fScale = 1
            # 0.4 for 121 points
            fScale = 0.35 * max(0,(1 - (3*count/nbPoints)**2)) + 0.8
            for point in allPoints:
                for near in point.near:
                    midX ,midY = point.midTo(near)
                    f = hScale / self.hDist(midX, midY) - point.distTo(near)
                    if f > 0:
                        near.repulsiveForce(f * fScale, point)
            for point in allPoints:
                isInside, closestPoint = self.region.isInside(point)
                if not isInside:
                    point.moveTo(closestPoint)
                if point.moveDist() > l0/3:
                    exit = False
                point.update()
                if point.distFromOrigin() > l0:
                    updateTri = True
            if updateTri :
                try :
                    c += 1
                    self.delaunayTriangulation()
                except :
                    self.plot(tri=True)
                    print('error')
                    exit = True
            if exportPeriod != 0  and count%exportPeriod == 0:
                self.exportToMax(client)
            if limit != 0 and count > limit:
                print('forced exit')
                exit = True
        for point in allPoints:
            point.resetNear()
        return count, c

    def simpleAttractorOld(self, G, coeff, recallUni=True, client=None):
        mx = [0.4]
        my = [0.5]
        for i in range(len(mx)):
            s_x = coeff
            s_y = 2 * coeff
            theta = 0
            a = cos(theta)**2 / (2 * s_x**2) + sin(theta)**2 / (2 * s_y**2)
            b = sin(2*theta) / (4 * s_x**2) - sin(2*theta)**2 / (4 * s_y**2)
            c = sin(theta)**2 / (2 * s_x**2) + cos(theta)**2 / (2 * s_y**2)
            f_norm = norm([mx[i],my[i]],[[a,b],[b,c]])
            center_gaussian = Point(mx[i],my[i])
            max_gaussian = f_norm.pdf((mx[i],my[i]))
            allPoints = self.getAllPoints()
            for p in allPoints:
                if recallUni:
                    p.recallUni()
                fg = G * f_norm.pdf((p.x, p.y)) / max_gaussian
                fa = G / p.distTo(center_gaussian)**2
                fs = G * p.distTo(center_gaussian)
                f = fs
                k = 1
                l = f / k
                l = min(l, p.distTo(center_gaussian))
                p.attractiveForce(l, center_gaussian)
        for p in allPoints:
            p.update()
        if client:
            self.exportToMax(client)

    def simpleAttractor(self, coeff, recallUni=True, client=None):
        mx = [0.2, 0.5, 0.8]
        my = [0.5, 0.2, 0.8]
        # mx = [0.5]
        # my = [0.5]
        for i in range(len(mx)):
            def dg(x, y): 
                sx = 0.25
                sy = 0.25
                deriv_gauss = exp(-1 * ((x - mx[i])**2 / (2 * sx**2) + (y - my[i])**2 / (2 * sy**2)))
                deriv_gauss *= sqrt((x - mx[i])**2 / sx**2 + (y - my[i])**2 / sy**2)
                return deriv_gauss
            center_gaussian = Point(mx[i],my[i])
            max_gaussian = 0.86
            allPoints = self.getAllPoints()
            for p in allPoints:
                if recallUni:
                    p.recallUni()
                # fg = G * f_norm.pdf((p.x, p.y)) / max_gaussian
                # fa = G / p.distTo(center_gaussian)**2
                # fs = G * p.distTo(center_gaussian)
                f = 0.184 * dg(p.x, p.y) / max_gaussian
                c = 0.55 # 0.55: unused
                border = 0.5 - c
                if (border < p.x < 1-border and border < p.y < 1-border):
                    k = min(9 * (p.distTo(center_gaussian)/coeff)**7 + 1, 10)
                else:
                    k = (10 - 9*c**7+1)/(sqrt(2)/2 - c) * p.distTo(center_gaussian) + 9*c**7+1
                    # k = 10
                l = f / k
                l = min(l, p.distTo(center_gaussian))
                p.attractiveForce(l, center_gaussian)
        for p in allPoints:
            p.update()
        if client:
            self.exportToMax(client)
    
    def exportToMax(self, client, itrp=0):
        for buffer in self.buffers:
            client.send_message('/buffer_index', buffer.id)
            uniX = [p.x*(1-itrp) + p.uni_x*itrp for p in buffer.points]
            uniY = [p.y*(1-itrp) + p.uni_y*itrp for p in buffer.points]
            n_rows = len(uniX)
            steps = ceil(n_rows/200)
            client.send_message('/matrixcol', 7)
            for i in range(steps):
                if i != steps-1:
                    client.send_message('/set_matrix', [i*200] + uniX[i*200:(i+1)*200])
                else :
                    client.send_message('/set_matrix', [i*200] + uniX[i*200:])
            client.send_message('/matrixcol', 8)
            for i in range(steps):
                if i != steps-1:
                    client.send_message('/set_matrix', [i*200] + uniY[i*200:(i+1)*200])
                else :
                    client.send_message('/set_matrix', [i*200] + uniY[i*200:])
            client.send_message('/refresh', 1)


class Buffer():
    
    def __init__(self, descr, idxX, idxY, nbId):
        self.id = nbId
        self.points = [Point(grain[idxX],grain[idxY]) for grain in descr]
        
    def normalize(self, upperX, lowerX, upperY, lowerY):
        widthX = upperX - lowerX
        widthY = upperY - lowerY
        for point in self.points:
            point.x = (point.x - lowerX) / widthX
            point.y = (point.y - lowerY) / widthY


class Point():
    
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.uni_x = x
        self.uni_y = y
        self.ogX = x
        self.ogY = y
        self.near = []
        self.pushX = 0.0
        self.pushY = 0.0
    
    def midTo(self, point):
        midX = (self.x + point.x)/2
        midY = (self.y + point.y)/2
        return midX, midY

    def getX(self):
        return self.x
    
    def getY(self):
        return self.y
    
    def distTo(self, point):
        return sqrt((self.x-point.x)**2 + (self.y-point.y)**2)

    def repulsiveForce(self, f, point):
        angle = arctan2(self.y - point.y, self.x - point.x)
        self.pushX += f * cos(angle)
        self.pushY += f * sin(angle)
        
    def attractiveForce(self, f, point):
        angle = arctan2(self.y - point.y, self.x - point.x)
        self.pushX -= f * cos(angle)
        self.pushY -= f * sin(angle)

    def update(self):
        self.x += self.pushX
        self.y += self.pushY
        self.pushX = 0.0
        self.pushY = 0.0
    
    def resetToLast(self):
        self.pushX = 0
        self.pushY = 0
        
    def updateOrigin(self):
        self.ogX = self.x
        self.ogY = self.y
        
    def distFromOrigin(self):
        return sqrt((self.x-self.ogX)**2 + (self.y-self.ogY)**2)
    
    def resetNear(self):
        self.near = []
        
    def moveTo(self, pt):
        nextX = pt.x
        nextY = pt.y
        self.pushX = nextX - self.x
        self.pushY = nextY - self.y

    def moveDist(self):
        return sqrt(self.pushX ** 2 + self.pushY ** 2)

    def vecOrientation(self, point):
        x = point.x - self.x
        y = point.y - self.y
        return arctan2(y,x)

    def storeUni(self):
        self.uni_x = self.x
        self.uni_y = self.y
    
    def recallUni(self):
        self.x = self.uni_x
        self.y = self.uni_y

    def __str__(self):
        return str(self.x) + ' ' + str(self.y)

    def __repr__(self):
        return str(self.x) + ' ' + str(self.y)

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import random as rd
    from region import RegionPolygon
    from time import sleep

    class CorpusExt(Corpus):

        def __init__(self, track, region, descrX, descrY):
            Corpus.__init__(self, track, region, descrX, descrY)

        def exportToMax(self, client=None):
            plt.clf()
            print('plotting')
            for buffer in self.buffers:
                for p in buffer.points:
                    for pt in p.near:
                        plt.plot([p.x,pt.x],[p.y,pt.y],'--',color='grey')
                    plt.arrow(p.x, p.y, p.pushX, p.pushY, color='red', lw=1)
                    plt.plot(p.x,p.y,'.',color='orange')
            plt.pause(0.05)

    string = "0.533341 0.443199 0.461124 0.333393 0.235369 0.341284 0.184073 0.449201 0.366361 0.592812 0.381871 0.456201 0.355112 0.521218 0.321692 0.355933 0.335829 0.499081 0.501694 0.516926 0.224963 0.422581 0.016166 0.138031 0.084873 0.156853 0.276466 0.405155 0.32669 0.614651 0.250342 0.602289 0.171687 0.403733 0.476653 0.394199 0.330887 0.36904 0.445647 0.171542 0.293811 0.297335 0.10004 0.242068 0.101647 0.131399 0.33216 0.678671 0.426982 0.535679 0.38183 0.664461 0.150601 0.51489 0.164425 0.242656 0.35409 0.335979 0.372314 0.528419 0.201434 0.07032 0.56222 0.533036 0.736877 0.518619 0.290352 0.405068 0.1674 0.167554 0.269122 0.18855 0.559746 0.578843 0.704666 0.505914 0.932668 0.548123 0.670066 0.626163 0.741748 0.532013 0.834602 0.875518 0.431553 0.433747 0.345441 0.38521 0.419352 0.470731 0.669068 0.552542 0.479218 0.393032 0.585085 0.245646 0.579012 0.267689 0.327875 0.597135 0.125261 0. 0.377685 0.602101 0.767582 1. 0.648253 0.550314 0.694302 0.693364 0.845912 0.68392 0.64453 0.519549 0.801108 0.73224 0.434903 0.561222 0.332164 0.575129 0.083547 0.217152 0.514927 0.45382 0.697068 0.72071 0.453567 0.477648 0.396405 0.611031 0.477611 0.542573 0.438698 0.328845 0.52628 0.591951 0.45265 0.538925 0.686413 0.453525 0.281086 0.187226 0.579931 0.838417 0.342517 0.192873 0.41947 0.313064 0.386622 0.628791 0.424354 0.640153 0.231787 0.312107 0.772592 0.957769 0.784949 0.354123 0.226804 0.232295 0.579663 0.405865 0.417996 0.419061 0.3372 0.309536 0.573854 0.435242 0.52705 0.61138 1. 0.68059 0.827221 0.509968 0.224084 0.255841 0.245053 0.110856 0.611013 0.940906 0.291219 0.500104 0.248476 0.621241 0.25316 0.336038 0.161209 0.275135 0.312219 0.27193 0.250153 0.193189 0.318415 0.197308 0.496073 0.541758 0.122787 0.507507 0.255858 0.495058 0.223762 0.397947 0.415263 0.394128 0.225362 0.27094 0.20997 0.441178 0.269653 0.431944 0.414042 0.557085 0.105695 0.266894 0.332031 0.399626 0.22243 0.275773 0.221896 0.166882 0.154218 0.285558 0.170824 0.419713 0. 0.149428 0.142284 0.329449 0.10978 0.496257 0.290802 0.724902 0.507622 0.614905 0.299456 0.59643 0.134769 0.289035 0.169 0.268101 0.276389 0.480706 0.215989 0.360894 0.027138 0.113574 0.117958 0.284064 0.123688 0.28819 0.193348 0.324845 0.071452 0.165619 0.162508 0.304272"
    coord = [float(s) for s in string.split(' ')]
    # N=200
    # coord = [rd.random() for i in range(2*N)]
    x = coord[0:len(coord):2]
    y = coord[1:len(coord)+1:2]
    track = {'0':[[x[i], y[i]] for i in range(len(x))]}
    vertices = ((0,0),(1,0),(1,1),(0,1))
    region = RegionPolygon(vertices)
    corpus = CorpusExt(track, region, 0, 1)
    with plt.ion():
        corpus.unispring(0.02, 0.06, exportPeriod=1, client=None, limit=500)
    corpus.exportToMax()
    plt.show()