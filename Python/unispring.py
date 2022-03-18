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
        
    def getBorderPoints(self):
        allPoints = self.getAllPoints()
        allCoord = asarray([(pt.x, pt.y) for pt in self.getAllPoints()])
        hull = [allPoints[idx] for idx in ConvexHull(allCoord).vertices]
        tree = KDTree(allCoord)
        
        maxDist = self.meanDistance()
        
        current = hull[0]
        prev = hull[-1]
        border = [current]
        print(int(degrees(current.vecOrientation(prev))))
        while True:
            nearest = tree.query_ball_point((current.x,current.y),
                                            r=maxDist)
            maxAngle = 2 * pi
            for pt in [allPoints[idx] for idx in nearest]:
                isCloseToPrev = prev.distTo(pt) < maxDist
                    
                if pt != prev and pt != current:
                    anglePrev = current.vecOrientation(prev)
                    angleNext = current.vecOrientation(pt)
                    angle = angleNext - anglePrev
                    if angle < 0:
                        angle += 2 * pi
                    if angle < maxAngle:
                        print(round(pt.x, 2), round(pt.y, 2),
                              int(degrees(anglePrev)), int(degrees(angleNext)),
                              int(degrees(angle)), ' x'
                              )
                        maxAngle = angle
                        nextPoint = pt
                    else:
                        print(round(pt.x, 2), round(pt.y, 2),
                              int(degrees(anglePrev)), int(degrees(angleNext)),
                              int(degrees(angle))
                              )
                        
            if nextPoint==border[0]:
                print('done')
                return border
            if len(border) > 100:
                print('aie')
                return border
            
            border.append(nextPoint)
    
            self.plot(tri=False, show=False)
            plt.plot(border[-2].x, border[-2].y, 'x')
            for i in range(1,len(border)):
                p1 = border[i-1]
                p2 = border[i]
                x = (p1.x, p2.x)
                y = (p1.y, p2.y)
                plt.plot(x,y,color='blue')
            plt.show()
            input((round(nextPoint.x, 2), round(nextPoint.y, 2),
                   int(degrees(anglePrev)), int(degrees(angleNext))
                  ))
            
            prev, current = current, nextPoint
    
    def meanDistance(self):
        d = 0
        for point in self.getAllPoints():
            for near in point.near:
                d += point.distTo(near)
        return d / (2 * len(self.getAllPoints()))
            
    def preUniformization(self, og=(0,0), s = 1, inSquareAuto=False):
        '''
        Distribute points in a square before using unispring.

        Parameters
        ----------
        og : tuple, optional
            origin of the square (bottom-left). The default is (0,0).
        s : float, optional
            length of the side of the square. The default is 1.
        inSquareAuto : bool, optional
            autocalculate the insquare of the region. The default is False.
        '''
        if inSquareAuto:
            side, origin = self.region.inSquare()
        else:
            side, origin = (s, Point(og[0],og[1]))
        
        allPoints = []
        for buffer in self.buffers:
            allPoints += buffer.points
        nbPoints = len(allPoints)
        allPoints.sort(key=Point.getX)
        for i in range(nbPoints):
            allPoints[i].x = (i / (nbPoints - 1)) * side + origin.x
        allPoints.sort(key=Point.getY)
        for i in range(nbPoints):
            allPoints[i].y = (i / (nbPoints - 1)) * side + origin.y
        
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
        
    def unispringUniform(self, k, minDist, maxDist, plotPeriod=0):
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
        self.data['tracks'].append(dict(self.data['tracks'][1]))
        self.data['numTracks'] += 1
        descr = self.data['tracks'][2]
        descr['name'] = "descr-audio-uni"
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
    
    def isInside(self, pt):
        nextX = pt.x + pt.pushX
        nextY = pt.y + pt.pushY
        closestPoints = [self.points[i] 
                         for i in self.tree.query((nextX, nextY), k=2)[1]]
        isIn1 = closestPoints[0].edge.isRightSide(pt)
        isIn2 = closestPoints[1].edge.isRightSide(pt)  
        return isIn1 and isIn2, closestPoints[0]

    def plot(self, show=True):
        pointsX = [point.x for point in self.points]
        pointsY = [point.y for point in self.points]
        fig = plt.scatter(pointsX, pointsY, s = 8)
        fig.set_cmap('hot')
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        fig.axes.set_aspect('equal')
        fig.axes.set_xlim([-0.1,1.1])
        fig.axes.set_ylim([-0.1,1.1])
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
    
    def plot(self, show=True):
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
        if show:
            plt.show()

if __name__ == '__main__':
    ## region building --> trigonometric rotation !
    vertices = ((0,0),(1,0),(1,1),(0,1))
    vertices2 = ((0.1,0.1),(0.8,0.2),(0.8,0.6),(0.6,0.5),(0.05,0.68))
    regionSquare = RegionPolygon(vertices)
    regionPoly = RegionPolygon(vertices2)
    regionCircle = RegionCircle(0.5,0.5,0.5)
    
    region = regionSquare
    
    ## corpus creation
    descX = 'CentroidMean'
    descY = 'PeriodicityMean'
    corpus = Corpus('/Users/victorparedes/Documents/Max 8/Projects/Sympoiesis/corpus/data_sympoiesis_cut.json',
    region, descX, descY, plot=False)

    ## pre-uniformization
    corpus.preUniformization(inSquareAuto=False)
    #tri = corpus.delaunayTriangulation()
    #corpus.plot(tri=False, show=False)
    corpus.unispringUniform(1, 0.02, 0.01, plotPeriod=100)
    #border = corpus.getBorderPoints()
    corpus.plot()

    ## export final corpus to json
    #corpus.exportJson('data_sympoiesis_uni.json')