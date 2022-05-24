"""
Created on Wed Feb 23 13:58:22 2022

@author: victorparedes
"""
# Import
from scipy.spatial import Delaunay, KDTree
from numpy import arctan2, sqrt, sin, cos, asarray, degrees
from math import pi, ceil


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
        
    def meanDistance(self):
        d = 0
        for point in self.getAllPoints():
            for near in point.near:
                d += point.distTo(near)
        return d / (2 * len(self.getAllPoints()))
            
    def preUniformization(self, resize=False, og=(0,0), s = 1, inSquareAuto=False):
        if inSquareAuto:
            p1, p2 = self.region.getBoundingBox()
            sideX, sideY, origin = p2.x-p1.x, p2.y-p1.y, p1
        else:
            sideX, sideY, origin = s, s, Point(og[0],og[1])
        
        allPoints = []
        for buffer in self.buffers:
            allPoints += buffer.points
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
        allPoints = []
        for buffer in self.buffers:
            for point in buffer.points:
                allPoints.append(point)
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
        
    def unispringUniform(self, k, minDist, maxDist, exportPeriod=0, client=None, limit=0):
        # first triangulation
        self.delaunayTriangulation()
        self.preUniformization(resize=True, inSquareAuto=True)
        allPoints = []
        for buffer in self.buffers:
            allPoints += buffer.points
        nbPoints = len(allPoints)
        # l0 is calculated for a uniform target distribution with all 
        # delaunay triangles becoming equilateral, does not account for points 
        # on the borders
        uniform_density = nbPoints / self.region.getArea()
        l0 = sqrt(2/(sqrt(3)*uniform_density))
        exit = False
        count= 0
        while not exit:
            exit = True
            count += 1
            updateTri = False
            for point in allPoints:
                for nearPoint in point.near:
                    f = k * (l0 - point.distTo(nearPoint))
                    if f > 0:
                        nearPoint.repulsiveForce(f, point)
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


class BorderPoint(Point):
    
    def __init__(self, x, y, edge):
        Point.__init__(self, x, y)
        self.edge = edge
    
    def __eq__(self, pt):
        return self.x == pt.x and self.y == pt.y


class RegionPolygon():
    
    def __init__(self, lstVertices, density=50):
        self.vertices = []
        self.edges = []
        self.points = []
        self.tree = None
        self.addVertices(lstVertices)
        self.calculateBorders(density)
        self.generateTree()
    
    def calculateBorders(self, density):
        self.points = []
        for edge in self.edges:
            self.points += edge.segment(density)
        
    def addVertices(self, lstVertices):
        for pt in lstVertices:
            self.vertices.append(Point(pt[0], pt[1]))
        for i in range(1, len(self.vertices)):
            self.edges.append(Edge(self.vertices[i-1], self.vertices[i]))
        self.edges.append(Edge(self.vertices[-1], self.vertices[0]))
    
    def generateTree(self):
        points = asarray([[pt.x, pt.y] for pt in self.points])
        self.tree = KDTree(points)
    
    def getBarycenter(self):
        barycenter = Point(0,0)
        for point in self.vertices:
            barycenter.x += point.x
            barycenter.y += point.y
        barycenter.x /= len(self.vertices)
        barycenter.y /= len(self.vertices)
        return barycenter
    
    def inSquare(self):
        center = self.getBarycenter()
        distToClosest = self.tree.query((center.x,center.y))[0]
        side = 2 * distToClosest / sqrt(2)
        origin = Point(center.x - side/2, center.y - side/2)
        return side, origin
        
    def getArea(self):
        area = 0
        for i in range(len(self.vertices)-1):
            p1 = self.vertices[i]
            p2 = self.vertices[i+1]
            area += p1.x * p2.y - p2.x * p1.y
        p1 = self.vertices[-1]
        p2 = self.vertices[0]
        area += p1.x * p2.y - p2.x * p1.y
        return area / 2

    def getBoundingBox(self):
        ver1 = self.vertices[0]
        ver2 = self.vertices[1]
        xmin = min(ver1.x, ver2.x)
        xmax = max(ver1.x, ver2.x)
        ymin = min(ver1.y, ver2.y)
        ymax = max(ver1.y, ver2.y)
        for vertice in self.vertices:
            if vertice.x < xmin:
                xmin = vertice.x
            elif vertice.x > xmax:
                xmax = vertice.x
            if vertice.y < ymin:
                ymin = vertice.y
            elif vertice.y > ymax:
                ymax = vertice.y
        return Point(xmin, ymin), Point(xmax, ymax)
    
    def isInside(self, pt):
        nextX = pt.x + pt.pushX
        nextY = pt.y + pt.pushY
        closestPoints = [self.points[i] 
                         for i in self.tree.query((nextX, nextY), k=2)[1]]
        isIn1 = closestPoints[0].edge.isRightSide(pt)
        isIn2 = closestPoints[1].edge.isRightSide(pt)  
        return isIn1 and isIn2, closestPoints[0]


class Edge():
    
    def __init__(self, p1, p2):
        self.start = p1
        self.end = p2
    
    def segment(self, density):
        p1 = self.start
        p2 = self.end
        dist = p1.distTo(p2)
        nbPoints = int(dist * density)
        xSeg = [i*(p2.x-p1.x)/nbPoints + p1.x for i in range(nbPoints)]
        ySeg = [i*(p2.y-p1.y)/nbPoints + p1.y for i in range(nbPoints)]
        return [BorderPoint(xSeg[i], ySeg[i], self) for i in range(nbPoints)]
    
    def isRightSide(self, pt):
        og = self.start
        ref = self.end
        refAngle = arctan2(ref.y - og.y, ref.x - og.x)
        ptAngle = arctan2(pt.y - og.y, pt.x - og.x)
        return sin(refAngle - ptAngle) < 0.00001


class RegionCircle():
    
    def __init__(self, centerX, centerY, radius):
        self.center = Point(centerX, centerY)
        self.radius = radius
        
    def getArea(self):
        return pi * self.radius**2
    
    def getBoundingBox(self):
        xmin = self.center.x - self.radius
        ymin = self.center.y - self.radius
        xmax = self.center.x + self.radius
        ymax = self.center.y + self.radius
        return Point(xmin, ymin), Point(xmax, ymax)
    
    def isInside(self, pt):
        nextX = pt.x + pt.pushX
        nextY = pt.y + pt.pushY
        angle = arctan2(nextY - self.center.y, nextX - self.center.x)
        distFromCenter = self.center.distTo(pt)
        closestX = self.radius * cos(angle) + self.center.x
        closestY = self.radius * sin(angle) + self.center.y
        return distFromCenter <= self.radius, Point(closestX, closestY)
    
    def inSquare(self):
        origin = Point(0,0)
        angle = arctan2(-1,-1)
        origin.x = self.radius * cos(angle) + self.center.x
        origin.y = self.radius * sin(angle) + self.center.y
        side = 2 * self.radius / sqrt(2)
        return side, origin


if __name__ == '__main__':
    print('AH')