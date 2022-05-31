from unispring import Point
from scipy.spatial import KDTree
from numpy import arctan2, sqrt, sin, cos, asarray
from math import pi

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
    print('no main process')