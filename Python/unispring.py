"""
Created on Wed Feb 23 13:58:22 2022

@author: victorparedes
"""
# Import
from scipy.spatial import Delaunay
from numpy import arctan2, sqrt, sin, cos, asarray
from scipy.stats import multivariate_normal as norm
from math import ceil

def hFunction(type):
    if type == 'uniform':
        f = lambda x, y : 1
    if type == 'gaussian':
        mx = 0.5
        my = 0.5
        s = 0.02
        f_norm = norm([mx,my],[[s,0],[0,s]])
        p = f_norm.pdf([0.5, 0.5])
        f = lambda x, y : 1 + f_norm.pdf([x, y])/p
    if type == 'progressive':
        f = lambda x, y : 1 + x
    return f

class Corpus():
    
    def __init__(self, track, region, descrX, descrY, hDist = 'uniform'):
        # import json
        # descriptors extraction
        self.buffers = []
        for key,buffer in track.items():
            self.buffers.append(Buffer(buffer, descrX, descrY, int(key)))
        # initialize normalize bool
        self.is_norm = False
        self.region = region
        self.hDist = hFunction(hDist)
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
        
    def unispringUniform(self, minDist, maxDist, exportPeriod=0, client=None, limit=0):
        # first triangulation
        self.delaunayTriangulation()
        self.preUniformization(resize=True, inSquareAuto=True)
        allPoints = self.getAllPoints()
        # l0 is calculated for a uniform target distribution with all 
        # delaunay triangles becoming equilateral, does not account for points 
        # on the borders
        nbPoints = len(allPoints)
        uniform_density = nbPoints / self.region.getArea()
        l0 = sqrt(2/(sqrt(3)*uniform_density))
        print('l0 : ',l0)
        exit = False
        count = 0
        while not exit:
            exit = True
            count += 1
            updateTri = False
            hScale = self.getScalingFactor(l0)
            fScale = 1
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
                if sqrt(point.pushX**2 + point.pushY**2) > minDist:
                    exit = False
                point.update()
                if point.distFromOrigin() > maxDist:
                    updateTri = True
            if updateTri :
                try :
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
        return count
    
    def exportToMax(self, client):
        for buffer in self.buffers:
            client.send_message('/buffer_index', buffer.id)
            uniX = [point.x for point in buffer.points]
            uniY = [point.y for point in buffer.points]
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
        self.ogX = x
        self.ogY = y
        self.near = []
        self.pushX = 0.0
        self.pushY = 0.0
        self.l0 = 0
        self.f = 0
    
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
        
    def update(self):
        self.x += self.pushX
        self.y += self.pushY
        self.pushX = 0.0
        self.pushY = 0.0
        
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
    
    def distToCenter(self):
        return self.distTo(Point(0.5, 0.5))

    def vecOrientation(self, point):
        x = point.x - self.x
        y = point.y - self.y
        return arctan2(y,x)

if __name__ == '__main__':
    print('no main process')