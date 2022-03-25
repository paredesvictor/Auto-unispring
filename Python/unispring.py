#%%
"""
Created on Wed Feb 23 13:58:22 2022

@author: victorparedes
"""
# Import
import json
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay, ConvexHull, KDTree
from numpy import arctan2, sqrt, sin, cos, asarray, degrees
from math import pi

plt.rcParams['figure.dpi'] = 100
    
class Corpus():
    
    def __init__(self, fileName, region, descrX, descrY, plot=False):
        '''
        Take a json file from mubu and initialise a corpus by extracting 
        descriptors from the second track (can manage multiple buffers).

        Parameters
        ----------
        fileName : string
            path to the json file created by mubu, the file must have at least
            2 tracks with the second one containing sound descriptors.
        region : Region object
            the user-defined region in which the points will be scattered.
        '''
        # import json
        with open(fileName, 'r') as dataFile:
            rawData = dataFile.read()
        self.data = json.loads(rawData)
        descr = self.data['tracks'][1]
        descrNames = descr['mxColNames']
        # descriptors extraction
        self.buffers = []
        for i in range(len(descr['buffers'])):
            bufferDescr = {}
            nbPts = len(descr['buffers'][i]['mxData'])
            for j,name in enumerate(descrNames):
                descrData = descr['buffers'][i]['mxData']
                bufferDescr[name] = descrData[j:nbPts:len(descrNames)]
            self.buffers.append(Buffer(bufferDescr, i))
        # initialize normalize bool
        self.is_norm = False
        self.region = region
        # point extractions from choice of axis
        self.extractPoints(descrX, descrY)
        if plot:
            self.plot()
        self.normalize()
        self.preUniformization()
    
    def extractPoints(self, descrNameX, descrNameY):
        '''
        Create points by extracting 2 descriptors from the descriptor matrix,
        one for X, one for Y. Does this for each buffer.

        Parameters
        ----------
        descrNameX : int
            index of the descriptor for X.
        descrNameY : int
            index of the descriptor for Y.
        '''
        for buffer in self.buffers:
            buffer.extractPoints(descrNameX, descrNameY)
    
    def getAllPoints(self):
        allPoints = []
        for buffer in self.buffers:
            allPoints += buffer.points
        return allPoints
        
    def normalize(self):
        '''
        Normalize the points to ([0,1],[0,1])
        '''
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
        '''
        Distribute points in the bounding box of the region.

        Parameters
        ----------
        resize : bool, optional
            just move the points in the destination box, no sort. The default is False
        og : tuple, optional
            origin of the square (bottom-left). The default is (0,0).
        s : float, optional
            length of the side of the square. The default is 1.
        inSquareAuto : bool, optional
            autocalculate the bounding box of the region. The default is False.
        '''
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
        '''
        Perform a Delaunay triangulation on all points using
        scipy.spatial toolkit.
        '''
        allCoord = []
        for buffer in self.buffers:
            allCoord += [[pt.x, pt.y] for pt in buffer.points]
        allCoord = asarray(allCoord)
        triangulation = Delaunay(allCoord)
        self.updateNearPoints(triangulation)
        return triangulation
    
    def convexHull(self):
        allCoord = [[pt.x, pt.y] for pt in self.getAllPoints()]
        allCoord = asarray(allCoord)
        convexHull = ConvexHull(allCoord)
        return convexHull
    
    def updateNearPoints(self, triangulation):
        '''
        Update neighboors of each point in the corpus after triangulation

        Parameters
        ----------
        triangulation : triangulation
            Triangulation obtained using scipy.spatial toolkit.
        '''
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
        
    def unispringUniform(self, k, minDist, maxDist, plotPeriod=0, limit=0):
        '''
        Perform a distribution of the corpus points in the user-defined region
        using a spring-mass physical model.
        

        Parameters
        ----------
        k : float
            stiffness of the spring mass assuming a uniform mass of 1.
        minDist : float
            The target minimal displacement of all the points used as an exit 
            condition.
        maxDist : float
            The maximal displacement above which the triangulation is updated.
        plotPeriod : bool, optional
            Number of steps before plotting a graph. The default is 0 (no plot).

        Returns
        -------
        count : int
            The total number of steps it took to reach minimal displacement
            condtion.
        '''
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
            if plotPeriod != 0  and count%plotPeriod == 0:
                self.plot(tri=False)
            if limit != 0 and count > limit:
                print('forced exit')
                exit = True
        return count
            
    def plot(self, tri=False, show=True):
        '''
        Plot a graph of the corpus points.

        Parameters
        ----------
        tri : bool, optional
            display Delaunay trinagulation. The default is False.
        show : bool, optional
            show immediately. The default is True.
        '''
        for buffer in self.buffers :
            pointsX = [point.x for point in buffer.points]
            pointsY = [point.y for point in buffer.points]
            fig = plt.scatter(pointsX, pointsY, s = 8)
        if self.is_norm:
            if tri:
                allPoints = []
                for buffer in self.buffers:
                    allPoints += buffer.points
                for point in allPoints:
                    for nearPoint in point.near:
                        segX = [point.x, nearPoint.x]
                        segY = [point.y, nearPoint.y]
                        plt.plot(segX, segY, color='black',
                                 linestyle='dashed', linewidth = 1,
                                 alpha = 0.3)
            fig.set_cmap('hot')
            fig.axes.get_xaxis().set_visible(False)
            fig.axes.get_yaxis().set_visible(False)
            fig.axes.set_aspect('equal')
            fig.axes.set_xlim([-0.1,1.1])
            fig.axes.set_ylim([-0.1,1.1])
        else:
            fig.axes.set_aspect('auto')
        if show:
            plt.show()
    
    def exportJson(self, name):
        '''
        Export the corpus into a new json file with added descriptors 
        corresponding to the new distribution.
        '''
        descr = self.data['tracks'][1]
        descr['mxColNames'].append('unispringX')
        descr['mxColNames'].append('unispringY')
        for buffer in self.buffers:
            buffer.updateJson(self.data, descr['mxCols'])
        descr['mxCols'] += 2
        with open(name, 'w') as file:
            json.dump(self.data, file, indent=2)
        
        
class Buffer():
    
    def __init__(self, descr, nbId):
        self.id = nbId
        self.descr = descr
        self.points = []
    
    def extractPoints(self, descrNameX, descrNameY):
        descrX = self.descr[descrNameX]
        descrY = self.descr[descrNameY]
        self.points = [Point(descrX[i], descrY[i]) for i in range(len(descrX))]
        
    def normalize(self, upperX, lowerX, upperY, lowerY):
        widthX = upperX - lowerX
        widthY = upperY - lowerY
        for point in self.points:
            point.x = (point.x - lowerX) / widthX
            point.y = (point.y - lowerY) / widthY
    
    def updateJson(self, data, nbDescr):
        allCoord = []
        for point in self.points:
            allCoord.append([point.x, point.y])
        vals = data['tracks'][1]['buffers'][self.id]['mxData']
        
        updateVals = []
        for i in range(len(vals)):
            updateVals.append(vals[i])
            if i%nbDescr == 6:
                updateVals.append(float(allCoord[i//nbDescr][0]))
                updateVals.append(float(allCoord[i//nbDescr][1]))
        data['tracks'][1]['buffers'][self.id]['mxData'] = updateVals
        data['tracks'][1]['buffers'][self.id]['mxCols'] += 2
        
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
        # add a random displacement around bounds to avoid points overlap
        # when moved to bounds
        nextX = pt.x #+ random() * 0.0001
        nextY = pt.y #+ random() * 0.0001
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
    
    def __init__(self, lstVertices):
        self.vertices = []
        self.edges = []
        self.points = []
        self.tree = None
        self.addVertices(lstVertices)
        self.calculateBorders()
        self.generateTree()
    
    def calculateBorders(self, density=50):
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

    def plot(self, show=True, pltBbox=False):
        pointsX = [point.x for point in self.points]
        pointsY = [point.y for point in self.points]
        fig = plt.scatter(pointsX, pointsY, s = 8)
        fig.set_cmap('hot')
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        fig.axes.set_aspect('equal')
        fig.axes.set_xlim([-0.1,1.1])
        fig.axes.set_ylim([-0.1,1.1])
        if pltBbox:
            p1, p2 = self.getBoundingBox()
            rectangle = plt.Rectangle((p1.x,p1.y), p2.x-p1.x, p2.y-p1.y, fc='none', ec="red")
            fig.axes.add_patch(rectangle)
        if show:
            plt.show()
            
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
    
    def plot(self, show=True, pltBbox=False):
        circle = plt.Circle((self.center.x,self.center.y),
                            self.radius, fill=False)
        fig = plt.scatter(self.center.x,self.center.y,marker='x')
        fig.set_cmap('hot')
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        fig.axes.set_aspect('equal')
        fig.axes.set_xlim([-0.1,1.1])
        fig.axes.set_ylim([-0.1,1.1])
        fig.axes.add_patch(circle)
        if pltBbox:
            p1, p2 = self.getBoundingBox()
            rectangle = plt.Rectangle((p1.x,p1.y), p2.x-p1.x, p2.y-p1.y, fc='none', ec="red")
            fig.axes.add_patch(rectangle)
        if show:
            plt.show()

if __name__ == '__main__':
    ## region building --> trigonometric rotation !
    vertices = ((0,0),(1,0),(1,1),(0,1))
    coord = "0.2919999957084656 0.7099999785423279 0.671999990940094 0.7139999866485596 0.6399999856948853 0.2939999997615814 0.3160000145435333 0.2919999957084656 0.7239999771118164 0.7279999852180481 0.7059999704360962 0.25600001215934753 0.2540000081062317 0.272000014781951".split(' ')
    vertices2 = [(int(coord[i]),1-int(coord[i+1])) for i in range(0,len(coord),2)]
    regionSquare = RegionPolygon(vertices)
    regionPoly = RegionPolygon(vertices2)
    regionCircle = RegionCircle(0.5,0.5,0.5)
    
    region = regionSquare
    
    ## corpus creation
    descX = 'CentroidMean'
    descY = 'PeriodicityMean'
    corpus = Corpus('/Users/victorparedes/Documents/GitHub/Auto-unispring/corpus.json',
    region, descX, descY, plot=True)

    ## pre-uniformization
    corpus.plot()
    corpus.unispringUniform(1, 0.02, 0.01, plotPeriod=0)
    corpus.plot()
    regionPoly.plot()
#%%
    corpus.region = regionPoly
    corpus.unispringUniform(1, 0.02, 0.01, plotPeriod=0)
    corpus.plot()


    ## export final corpus to json
    #corpus.exportJson('data_sympoiesis_uni.json')