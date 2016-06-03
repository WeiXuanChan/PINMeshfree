'''
File: discretizer.py
Description: function definition
History:
    Date    Programmer SAR# - Description
    ---------- ---------- ----------------------------
  Author: w. x. chan 29Apr2016           - Created

'''

import numpy as np
from . import pinm as pinm
from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot
from matplotlib.widgets import Button
import random 
    
def stlImport(filePath,coords):
    # Create a new plot
    figure = pyplot.figure()
    pyplot.subplots_adjust(bottom=0.2)
    axes = mplot3d.Axes3D(figure)
    # Load the STL files and add the vectors to the plot
    modelMesh = mesh.Mesh.from_file(filePath)
    indexedTri=[]
    for n in range(len(modelMesh.vectors)):
        indexedTri.append(mplot3d.art3d.Poly3DCollection([modelMesh.vectors[n]],facecolors='b'))
        axes.add_collection3d(indexedTri[n])
    indexedTri[0].set_facecolor('k')
    
    scale = modelMesh.points.flatten(-1)
    axes.auto_scale_xyz(scale, scale, scale)
            
    callback = DomainSelector(indexedTri)
    
    axprev = pyplot.axes([0.7, 0.05, 0.1, 0.075])
    axnext = pyplot.axes([0.81, 0.05, 0.1, 0.075])
    axselect = pyplot.axes([0.05, 0.05, 0.15, 0.075])
    axaddToDomain = pyplot.axes([0.05, 0.85, 0.15, 0.075])
    axswapSelected = pyplot.axes([0.8, 0.85, 0.15, 0.075])
    bnext = Button(axnext, 'Next')
    bnext.on_clicked(callback.next)
    bprev = Button(axprev, 'Previous')
    bprev.on_clicked(callback.prev)
    bselect = Button(axselect, '(un)Select')
    bselect.on_clicked(callback.select)
    baddToDomain = Button(axaddToDomain, 'Add Domain')
    baddToDomain.on_clicked(callback.addToDomain)
    bswapSelected = Button(axswapSelected, 'Swap Selected')
    bswapSelected.on_clicked(callback.swapSelected)
    
    # Show the plot to the screen
    #pyplot.connect('key_press_event', callback.keyPressed)
    pyplot.show()
    
    maindomain=pinm.Domain('')
    subdomain=[]
    for domainNumber in range(callback.domainCount):
        subdomain.append(pinm.Domain(''))
        maindomain.addNode(subdomain[domainNumber])
    normalVector=[]
    vertices=[]
    minvert={}
    maxvert={}
    for n in range(len(modelMesh.normals)):
        normalVector.append({})
        vertices.append([])
        for keyIndex in range(len(coords)):
            normalVector[n][coords[keyIndex]]=modelMesh.normals[n][keyIndex]
        for m in range(3):
            temp_vert={}
            for keyIndex in range(len(coords)):
                if coords[keyIndex] not in minvert:
                    minvert[coords[keyIndex]]=modelMesh.vectors[n][m][keyIndex]
                else:
                    minvert[coords[keyIndex]]=min(minvert[coords[keyIndex]],modelMesh.vectors[n][m][keyIndex])
                if coords[keyIndex] not in maxvert:
                    maxvert[coords[keyIndex]]=modelMesh.vectors[n][m][keyIndex]
                else:
                    maxvert[coords[keyIndex]]=max(maxvert[coords[keyIndex]],modelMesh.vectors[n][m][keyIndex])
                temp_vert[coords[keyIndex]]=modelMesh.vectors[n][m][keyIndex]
            vertices[n].append(temp_vert)
    domainVertices=[]
    for n in range(8):
        temp_domainVertices={}
        for key in range(len(coords)):
            if (key==0 and (n in [1,2,5,6])) or (key==1 and (n in [2,3,6,7])) or (key==2 and (n in [4,5,6,7])):
                temp_domainVertices[coords[key]]=maxvert[coords[key]]
            else:
                temp_domainVertices[coords[key]]=minvert[coords[key]]
        domainVertices.append(temp_domainVertices)
    for n in range(len(callback.domainInfo)):
        temp_sub2domain=pinm.Domain('',norm=normalVector[n])
        temp_sub2domain.setCentroid(vertices[n])
        subdomain[callback.domainInfo[n]].addNode(temp_sub2domain)
    maindomain.setCentroid(domainVertices)
    return (maindomain,subdomain)


def filterNodes(domainList,nodalDistribution,closeness=0.2): #first in list is prioritized to keep
    nodes=[]
    for domain in domainList:
        for node in domain.nodes():
            nodes.append(node)
    for n in range(len(nodes)):
        if nodes[n].domain!=None:
            nodalSpacing=nodalDistribution(nodes[n].pos)
            closeNodalSpacing=multiplyDictionary(nodalSpacing,closeness)
            linkedNodes=findNodes(nodes[n].pos,domainList,distance=closeNodalSpacing,searchDepth=-1.)
            for temp_linkNode in linkedNodes:
                if temp_linkNode is not nodes[n]:
                    temp_domain=temp_linkNode.domain
                    temp_domain.removeNode(temp_linkNode)
                    temp_linkNode.domain=None
                    while len(temp_domain.subDomain)==0:
                        if temp_domain.superDomain==None:
                            break
                        else:
                            temp2_domain=temp_domain.superDomain
                            temp2_domain.removeNode(temp_domain)
                            temp_domain=temp2_domain
                
    return; 
def secondaryLinkNode(targetNode,primarylinkIdentifier,secondarylinkIdentifier='secondary'):
    nodes=[]
    targetNode.addLink(secondarylinkIdentifier,targetNode)
    for node in targetNode.link[primarylinkIdentifier]:
        for temp_node in node.link[primarylinkIdentifier]:
            targetNode.addLink(secondarylinkIdentifier,node)
    return;
def primaryLinkNodes(domainList,nodalDistribution,linkIdentifier='primary',closeness=1.5):#influence is function with dictionary input and output
    nodes=[]
    for domain in domainList:
        for node in domain.nodes():
            nodes.append(node)
    for n in range(len(nodes)):
        nodalSpacing=nodalDistribution(nodes[n].pos)
        expandedNodalSpacing=multiplyDictionary(nodalSpacing,closeness)
        linkedNodes=findNodes(nodes[n].pos,domainList,distance=expandedNodalSpacing,searchDepth=-1.)
        addNodesToLink=[]
        for temp_linkNode in linkedNodes:
            if temp_linkNode is not nodes[n]:
                addNodesToLink.append(temp_linkNode)
        addNodesToLink.insert(0,nodes[n])
        nodes[n].addLink(linkIdentifier,addNodesToLink)
    return;
def duplicateNode(coordinateIdentifier,value,nodePorting,newNodes,domainList,targetDomain):
    nodeInDomain=[]
    for domain in domainList:
        if type(domain) is pinm.Node:
            new_pos=node.pos.copy()
            new_pos[coordinateIdentifier]=value
            newNode=pinm.Node(new_pos)
            tempCopy=domain.variable.copy()
            for key in tempCopy:
                newNode.addvariable(key,tempCopy[key])
            newNode.addLink('copied from',domain)
            tempCopy=node.linkBasis.copy()
            for key in tempCopy:
                newNode.setLinkBasis(key,tempCopy[key])
            newNode.setNorm(node.norm.copy())
            newNode.setNormLink(node.normLink)
            for n in range(len(node.material)):
                newNode.addMaterial(n,node.material[n])
            tempCopy=node.variableLink.copy()
            for key in tempCopy:
                newNode.setVariableLink(key,tempCopy[key])
            nodePorting[domain]=newNode
            newNodes.append(newNode)
            nodeInDomain.append(newNode)
        else:
            newDomain=pinm.Domain('')
            newDomain.pos=domain.pos.copy()
            newDomain.maxDistance=domain.maxDistance.copy()
            newDomain.pos[coordinateIdentifier]=value
            newDomain.maxDistance[coordinateIdentifier]=0.
            nodeInDomain.append(newDomain)
            duplicateNode(coordinateIdentifier,value,nodePorting,newNodes,domain.subDomain,newDomain)
    targetDomain.addNode(nodeInDomain)
    return nodeInDomain
def extrudeDimension(domainList,coordinateIdentifier,valueList,prevLinkIdentifier='',nextLinkIdentifier=''):
    newDomain=[]
    prevNodes=[]
    for m in range(len(valueList)):
        newNodes=[]
        nodePorting={}
        newDomain.append([])
        for domain in domainList:
            tempDomain=pinm.Domain('')
            tempDomain.pos=domain.pos.copy()
            tempDomain.maxDistance=domain.maxDistance.copy()
            tempDomain.pos[coordinateIdentifier]=value
            tempDomain.maxDistance[coordinateIdentifier]=0.
            newDomain[-1].append(tempDomain)
            duplicateNode(coordinateIdentifier,value,nodePorting,newNodes,domain.subDomain,tempDomain)
        for new_node in newNodes:
            for temp_linkIdentifier in new_node.link['copied from'][0].link:
                if temp_linkIdentifier!='copied from':
                    tempList=[]
                    for linkNode in new_node.link['copied from'][0].link[temp_linkIdentifier]:
                        tempList.append(nodePorting[linkNode])
                    new_node.addLink(temp_linkIdentifier,tempList)
        if (prevLinkIdentifier!='' or nextLinkIdentifier!='') and len(prevNodes)!=0:
            for n in range(len(newNodes)):
                if prevLinkIdentifier!='':
                    prevNodes[n].addLink(linkIdentifier,newNodes[n])
                if nextLinkIdentifier!='':    
                    newNodes[n].addLink(linkIdentifier,prevNodes[n])
        prevNodes=newNodes[:]
    return newDomain
def arrangeExtrudeDimension(domainList,coordinateIdentifier,valueList,prevLinkIdentifier='',nextLinkIdentifier='',newDomainNameAddOn=' new',firstDomainNameAddOn='',lastDomainNameAddOn=''): 
    nameList=[]
    for domain in domainList:
        nameList.append(domain.name) 
    subDomain=extrudeDimension(domainList,coordinateIdentifier,valueList,prevLinkIdentifier=prevLinkIdentifier,nextLinkIdentifier=nextLinkIdentifier)
    newDomain=[]
    startDomain=[]
    endDomain=[]
    for n in range(len(subDomain[0])):
        startCount=0
        endCountReduce=0
        if firstDomainNameAddOn!='':
            if firstDomainNameAddOn!=lastDomainNameAddOn:
                subDomain[0][n].setDomainName(nameList[n]+firstDomainNameAddOn)
            startDomain.append(subDomain[0][n])
            startCount=1
        if lastDomainNameAddOn!='':
            if firstDomainNameAddOn!=lastDomainNameAddOn:
                subDomain[-1][n].setDomainName(nameList[n]+lastDomainNameAddOn)
            else:
                domainGroup=pinm.Domain(nameList[n]+firstDomainNameAddOn)
                domainGroup.addNode([subDomain[0][n],subDomain[-1][n]])
            endDomain.append(subDomain[-1][n])
            endCountReduce=1
        leftOverDomain=[]
        for m in range(startCount,len(subDomain)-endCountReduce):
            leftOverDomain.append(subDomain[m][n])
        if len(leftOverDomain)!=0:
            tempDomain=pinm.Domain(nameList[n]+newDomainNameAddOn)
            tempDomain.addNode(leftOverDomain)
            newDomain.append(tempDomain)
    return (newDomain,startDomain,endDomain)
    
def meshSurfaceDomainTriangle(subDomain,nodalDistribution):
    for domain in subDomain:
        for sub2domain in domain.subDomain:
            toBeFurtherMeshed=meshTriangleSpliting(sub2domain,nodalDistribution)
            while len(toBeFurtherMeshed)>0:
                copy_toBeFurtherMeshed=toBeFurtherMeshed
                toBeFurtherMeshed=[]
                for new_domain in copy_toBeFurtherMeshed:
                    for temp_domain in meshTriangleSpliting(new_domain,nodalDistribution):
                        toBeFurtherMeshed.append(temp_domain)
    return;
def meshMainDomain(mainDomain,boundaryDomainList,nodalDistribution,meshOuterNode=False):
    innerNodesDomain=pinm.Domain('')
    innerNodesDomain.setCentroid(mainDomain.vertices)
    mainDomain.addNode(innerNodesDomain)
    toBeFurtherMeshed=meshVolume(innerNodesDomain,boundaryDomainList,nodalDistribution,meshOuter=meshOuterNode)
    while len(toBeFurtherMeshed)>0:
        copy_toBeFurtherMeshed=toBeFurtherMeshed
        toBeFurtherMeshed=[]
        for new_domain in copy_toBeFurtherMeshed:
            for temp_domain in meshVolume(new_domain,boundaryDomainList,nodalDistribution,meshOuter=meshOuterNode):
                toBeFurtherMeshed.append(temp_domain)
    return innerNodesDomain

def meshTriangleSpliting(domain,nodalDistribution):  #nodalDistribution is a function with both i/o dictionary objects     
    toBeFurtherMeshed=[]
    subDomain=[]
    #check for odd triangle
    sidelength=[0.,0.,0.]
    maxSideLength=0.
    minSideLength=float('inf')
    maxSideIndex=-1
    minSideIndex=-1
    for n in range(len(domain.vertices)):
        for coord in domain.vertices[0]:
            sidelength[n]+=(domain.vertices[n][coord]-domain.vertices[n-1][coord])**2.
        sidelength[n]=np.sqrt(sidelength[n])
        if sidelength[n]>maxSideLength:
            maxSideLength=sidelength[n]
            maxSideIndex=n
        if sidelength[n]<minSideLength:
            minSideLength=sidelength[n]
            minSideIndex=n
    NodeSpacing=nodalDistribution(domain.pos)
    newPoint=multiplyDictionary(addDictionary([domain.vertices[maxSideIndex],domain.vertices[maxSideIndex-1]]),0.5)
    tri1Domain=pinm.Domain('',norm=domain.normalVector)
    tri1Domain.setCentroid([newPoint,domain.vertices[maxSideIndex],domain.vertices[maxSideIndex-2]])
    tri2Domain=pinm.Domain('',norm=domain.normalVector)
    tri2Domain.setCentroid([newPoint,domain.vertices[maxSideIndex-2],domain.vertices[maxSideIndex-1]])
    temp_total=0.
    for coord in NodeSpacing:
        temp_total+=NodeSpacing[coord]**2.
    nodeDis=np.sqrt(temp_total)
    if nodeDis<(sum(sidelength)/3.):
        subDomain.append(tri1Domain)   
        subDomain.append(tri2Domain)
        toBeFurtherMeshed.append(tri1Domain)
        toBeFurtherMeshed.append(tri2Domain)
    else:
        subDomain.append(pinm.Node(tri1Domain.pos,norm=domain.normalVector))
        subDomain.append(pinm.Node(tri2Domain.pos,norm=domain.normalVector))
    domain.addNode(subDomain)
    return toBeFurtherMeshed
def meshVolume(domain,boundaryDomainList,nodalDistribution,meshOuter=False):  #nodalDistribution is a function with both i/o dictionary objects     
    if meshOuter:
        meshOuterCoef=-1.
    else:
        meshOuterCoef=1.
    NodeSpacing=nodalDistribution(domain.pos)
    addNodeInstead=1
    for coord in domain.maxDistance:
        if domain.maxDistance[coord]>NodeSpacing[coord]:
            addNodeInstead=0
    centerPlane=[]
    centerPlaneMidPoints=[]
    for n in range(4):
        centerPlane.append(multiplyDictionary(addDictionary([domain.vertices[n],domain.vertices[4+n]]),0.5))
    for n in range(3):
        centerPlaneMidPoints.append(multiplyDictionary(addDictionary([centerPlane[n],centerPlane[n+1]]),0.5))
    centerPlaneMidPoints.append(multiplyDictionary(addDictionary([centerPlane[3],centerPlane[0]]),0.5))
    planeCentroid=[]
    midPoints=[]
    for m in range(2):
        midPoints.append([])
        for n in range(3):
            midPoints[m].append(multiplyDictionary(addDictionary([domain.vertices[m*4+n],domain.vertices[m*4+n+1]]),0.5))
        midPoints[m].append(multiplyDictionary(addDictionary([domain.vertices[m*4+3],domain.vertices[m*4]]),0.5))
    for m in range(2):
        planeCentroid.append(multiplyDictionary(addDictionary([midPoints[m][0],midPoints[m][2]]),0.5))
    subDomain=[]
    toBeFurtherMeshed=[]
    for m in range(2):
        for n in range(4):
            temp_subdomain=pinm.Domain('')
            temp_vertices=[midPoints[m][n-1],domain.vertices[4*m+n],midPoints[m][n],planeCentroid[m],
                           centerPlaneMidPoints[n-1],centerPlane[n],centerPlaneMidPoints[n],domain.pos]
            temp_subdomain.setCentroid(temp_vertices)
            temp_boundaryNode=findNodes(temp_subdomain.pos,boundaryDomainList)
            distancebetween={}
            for coord in temp_boundaryNode.pos:
                distancebetween[coord]=np.absolute(temp_boundaryNode.pos[coord]-temp_subdomain.pos[coord])
            boundaryNodes=findNodes(temp_subdomain.pos,boundaryDomainList,distance=distancebetween)
            innerNode=True
            for boundaryNode in boundaryNodes:
                boundaryNodeCentroid=boundaryNode.pos
                boundaryNodeNorm=boundaryNode.norm
                dotProduct=0.
                normamplitude=0.
                for coords in temp_subdomain.pos:
                    dotProduct+= (temp_subdomain.pos[coords]-boundaryNodeCentroid[coords])*boundaryNodeNorm[coords]
                    normamplitude+=boundaryNodeNorm[coords]**2.
                dotProduct=dotProduct/np.sqrt(normamplitude)
                for coords in temp_subdomain.maxDistance:
                    if (temp_subdomain.maxDistance[coords]*(1-addNodeInstead))<(meshOuterCoef*dotProduct):
                        innerNode=False
                        break
                if innerNode==False:
                    break
            if innerNode:
                if addNodeInstead==1:
                    temp_node=pinm.Node(temp_subdomain.pos,norm=domain.normalVector)
                    subDomain.append(temp_node)
                else:
                    toBeFurtherMeshed.append(temp_subdomain)
                    subDomain.append(temp_subdomain)
    domain.addNode(subDomain)
    return toBeFurtherMeshed;
def findNodes(position,domainList,distance=None,searchDepth=-1.):#assign search depth to -1 for nodes
    temp_searchDepth=searchDepth
    if distance==None:
        findNearest=True
    else:
        findNearest=False
    if findNearest:
        referenceDomain=None
        minDistanceSq=float("inf")
        otherDomain=[]
        for domain in domainList:
            temp_distSq=0.
            if bool(domain.pos):
                for coords in position:
                    temp_distSq+=(position[coords]-domain.pos[coords])**2.
                if minDistanceSq>temp_distSq:
                    minDistanceSq=temp_distSq
                    referenceDomain=domain
            else:
                for allDomain in domain.subDomain:
                    otherDomain.append(allDomain)
        if len(otherDomain)!=0:
            if type(referenceDomain) is pinm.Domain:
                for includeDomain in referenceDomain.subDomain:
                    otherDomain.append(includeDomain)
            elif type(referenceDomain) is pinm.Node:
                otherDomain.append(referenceDomain)
            nodes=findNodes(position,otherDomain,searchDepth=temp_searchDepth)
        elif (type(referenceDomain) is not pinm.Node) and searchDepth!=0:
            nodes=findNodes(position,referenceDomain.subDomain,searchDepth=(temp_searchDepth-1))
        else:
            nodes=referenceDomain
        return nodes
    else:
        nodes=[]
        for domain in domainList:
            toAdd=True
            if bool(domain.pos):
                if type(domain) is not pinm.Node:
                    maxDistance=domain.maxDistance
                else:
                    maxDistance={}
                    for coords in position:
                        maxDistance[coords]=0.
                for coords in position:
                    if np.absolute(position[coords]-domain.pos[coords])>(maxDistance[coords]+distance[coords]):
                            toAdd=False
            if toAdd:
                if type(domain) is not pinm.Node:
                    for temp_nodes in findNodes(position,domain.subDomain,distance):
                        nodes.append(temp_nodes)
                else:
                    nodes.append(domain)
        return nodes    
def addDictionary(a):
    result={}
    for dicts in a:
        for key in dicts:
            if key in result:
                result[key]+=dicts[key]
            else:
                result[key]=dicts[key]
    return result
def multiplyDictionary(a,b):
    result={}
    for key in a:
        result[key]=a[key]*b
    return result
def randonNumber(a,b):
    #excludes boundary
    result==a
    while result==a or result==b:
        result=random.uniform(a,b)
    return result

def plotNodes(nodes,coordinate=['x','y','z']):
    figure = pyplot.figure()
    axes = mplot3d.Axes3D(figure)
    coordinateKey=[]
    numOfNodes=len(nodes)
    coords=np.zeros((3,numOfNodes))
    for n in range(numOfNodes):
        for m in range(len(coords)):
            coords[m][n]=nodes[n].pos[coordinate[m]]
            
    axes.scatter(coords[0], coords[1], coords[2])
    pyplot.show()

          
class DomainSelector:
    def __init__(self,collectionList):
        self.ind = 0
        self.collectionList=collectionList
        self.selectedIndex=[]
        self.domainInfo=[]
        self.domainCount=1
        self.end=False
        self.keyFunc={'l':self.nextFunc,
                      'k':self.prevFunc,
                      's':self.selectFunc,
                      'a':self.addToDomainFunc}
        for n in collectionList:
            self.selectedIndex.append(False)
            self.domainInfo.append(0)
        self.maxIndex=len(collectionList)-1
    def next(self, event):
        self.nextFunc()
    def prev(self, event):
        self.prevFunc()
    def select(self, event):
        self.selectFunc()
        self.nextFunc()
    def addToDomain(self, event):
        self.addToDomainFunc()
    def swapSelected(self, event):
        self.swapSelectedFunc()
 #   def keyPressed(self,event):
 #       self.keyFunc[event.key]() #find code error
    def nextFunc(self):
        if not(self.end):
            if self.selectedIndex[self.ind]:
                self.collectionList[self.ind].set_facecolor('g')
            else:
                self.collectionList[self.ind].set_facecolor('b')
            self.ind += 1
            if self.ind>self.maxIndex:
                    self.ind = 0
            while self.domainInfo[self.ind]!=0:
                self.ind += 1
                if self.ind>self.maxIndex:
                    self.ind = 0
            if self.selectedIndex[self.ind]:
                self.collectionList[self.ind].set_facecolor('r')
            else:
                self.collectionList[self.ind].set_facecolor('k')
            pyplot.draw()
    def prevFunc(self):
        if not(self.end):
            if self.selectedIndex[self.ind]:
                self.collectionList[self.ind].set_facecolor('g')
            else:
                self.collectionList[self.ind].set_facecolor('b')
            self.ind -= 1
            if self.ind<0:
                    self.ind = self.maxIndex
            while self.domainInfo[self.ind]!=0:
                self.ind -= 1
                if self.ind<0:
                    self.ind = self.maxIndex
            if self.selectedIndex[self.ind]:
                self.collectionList[self.ind].set_facecolor('r')
            else:
                self.collectionList[self.ind].set_facecolor('k')
            pyplot.draw()
    def selectFunc(self):
        if not(self.end):
            if self.selectedIndex[self.ind]:
                self.collectionList[self.ind].set_facecolor('k')
                self.selectedIndex[self.ind]=False
            else:
                self.collectionList[self.ind].set_facecolor('r')
                self.selectedIndex[self.ind]=True
            pyplot.draw()
    def addToDomainFunc(self):
        for n in range(len(self.selectedIndex)):
            if self.selectedIndex[n]:
                self.selectedIndex[n]=False
                self.domainInfo[n]=self.domainCount
                self.collectionList[n].set_facecolor('none')
        self.domainCount +=1
        self.end=True
        for n in range(len(self.domainInfo)):
            if self.domainInfo[n]==0:
                self.ind = n
                self.collectionList[self.ind].set_facecolor('k')
                self.end=False
                break
        pyplot.draw()
    def swapSelectedFunc(self):
        for n in range(len(self.selectedIndex)):
            if self.domainInfo[n]==0:
                if self.selectedIndex[n]:
                    self.selectedIndex[n]=False
                    if n==self.ind:
                        self.collectionList[n].set_facecolor('k')
                    else:
                        self.collectionList[n].set_facecolor('b')
                else:
                    self.selectedIndex[n]=True
                    if n==self.ind:
                        self.collectionList[n].set_facecolor('r')
                    else:
                        self.collectionList[n].set_facecolor('g')
        pyplot.draw()





'''
#old code

def meshTriangle(domain,nodalDistribution):  #nodalDistribution is a function with both i/o dictionary objects     
    toBeFurtherMeshed=[]
    subDomain=[]
    #check for odd triangle
    sidelength=[0.,0.,0.]
    maxSideLength=0.
    minSideLength=float('inf')
    maxSideIndex=-1
    minSideIndex=-1
    for n in range(len(domain.vertices)):
        for coord in domain.vertices[0]:
            sidelength[n]+=(domain.vertices[n][coord]-domain.vertices[n-1][coord])**2.
        sidelength[n]=np.sqrt(sidelength[n])
        if sidelength[n]>maxSideLength:
            maxSideLength=sidelength[n]
            maxSideIndex=n
        if sidelength[n]<minSideLength:
            minSideLength=sidelength[n]
            minSideIndex=n
    if maxSideLength>(minSideLength*2):#two long one short
        newPoints=[]
        newPoints.append(addDictionary([multiplyDictionary(domain.vertices[minSideIndex-2],minSideLength/(minSideLength+sidelength[minSideIndex-2])),
                                        multiplyDictionary(domain.vertices[minSideIndex],sidelength[minSideIndex-2]/(minSideLength+sidelength[minSideIndex-2]))]))
        newPoints.append(addDictionary([multiplyDictionary(domain.vertices[minSideIndex-2],minSideLength/(minSideLength+sidelength[minSideIndex-1])),
                                        multiplyDictionary(domain.vertices[minSideIndex-1],sidelength[minSideIndex-1]/(minSideLength+sidelength[minSideIndex-1]))]))
        triDomain=pinm.Domain('')
        triDomain.setCentroid([newPoints[0],newPoints[1],domain.vertices[minSideIndex-2]])
        continueMesh=meshTriangle(triDomain,nodalDistribution)
        subDomain.append(triDomain)
        for domain in continueMesh:
            toBeFurtherMeshed.append(domain)
        quadDomain=pinm.Domain('')
        quadDomain.setCentroid([newPoints[0],newPoints[1],domain.vertices[minSideIndex-1],domain.vertices[minSideIndex]])
        continueMesh=meshQuadrilateral(triDomain,nodalDistribution)
        subDomain.append(quadDomain)
        for domain in continueMesh:
            toBeFurtherMeshed.append(domain)
    elif np.absolute(sidelength[maxSideIndex-1]/sidelength[maxSideIndex-2]+sidelength[maxSideIndex-2]/sidelength[maxSideIndex-1]-sidelength[maxSideIndex]**2./sidelength[maxSideIndex-1]/sidelength[maxSideIndex-2])<1: #more than 60 degree
        newPoint=multiplyDictionary(addDictionary([domain.vertices[maxSideIndex],domain.vertices[maxSideIndex-1]]),0.5)
        tri1Domain=pinm.Domain('')
        tri1Domain.setCentroid([newPoint,domain.vertices[maxSideIndex],domain.vertices[maxSideIndex-2]])
        continueMesh=meshTriangle(tri1Domain,nodalDistribution)
        subDomain.append(tri1Domain)
        for domain in continueMesh:
            toBeFurtherMeshed.append(domain)
        tri2Domain=pinm.Domain('')
        tri2Domain.setCentroid([newPoint,domain.vertices[maxSideIndex-2],domain.vertices[maxSideIndex-1]])
        continueMesh=meshTriangle(tr2Domain,nodalDistribution)
        subDomain.append(tri2Domain)
        for domain in continueMesh:
            toBeFurtherMeshed.append(domain)
    else:
        midPoints=[]
        for n in range(2):
            midPoints.append(multiplyDictionary(addDictionary([domain.vertices[n],domain.vertices[n+1]]),0.5))
        midPoints.append(multiplyDictionary(addDictionary([domain.vertices[2],domain.vertices[0]]),0.5))
        for n in range(3):
            temp_subdomain=pinm.Domain('',norm=domain.normalVector)
            temp_vertices=[midPoints[n-1],domain.vertices[n],midPoints[n],domain.pos]
            temp_subdomain.setCentroid(temp_vertices)
            for toBeFurtherMeshedDomain in meshQuadrilateral(temp_subdomain,nodalDistribution):
                toBeFurtherMeshed.append(toBeFurtherMeshedDomain)
            subDomain.append(temp_subdomain)
    domain.addNode(subDomain)
    return toBeFurtherMeshed
def meshQuadrilateral(domain,nodalDistribution):  #nodalDistribution is a function with both i/o dictionary objects     
    toBeFurtherMeshed=[]
    NodeSpacing=nodalDistribution(domain.pos)
    addNodeInstead=True
    for coord in domain.maxDistance:
        if domain.maxDistance[coord]>NodeSpacing[coord]:
            addNodeInstead=False
    midPoints=[]
    for n in range(3):
        midPoints.append(multiplyDictionary(addDictionary([domain.vertices[n],domain.vertices[n+1]]),0.5))
    midPoints.append(multiplyDictionary(addDictionary([domain.vertices[3],domain.vertices[0]]),0.5))
    subDomain=[]
    for n in range(4):
        temp_subdomain=pinm.Domain('',norm=domain.normalVector)
        temp_vertices=[midPoints[n-1],domain.vertices[n],midPoints[n],domain.pos]
        temp_subdomain.setCentroid(temp_vertices)
        if addNodeInstead:
            temp_node=pinm.Node(temp_subdomain.pos,norm=domain.normalVector)
            subDomain.append(temp_node)
        else:
            toBeFurtherMeshed.append(temp_subdomain)
            subDomain.append(temp_subdomain)
    domain.addNode(subDomain)
    return toBeFurtherMeshed
def meshSurfaceDomainQuad(subDomain,nodalDistribution):
    for domain in subDomain:
        for sub2domain in domain.subDomain:
            toBeFurtherMeshed=meshTriangle(sub2domain,nodalDistribution)
            while len(toBeFurtherMeshed)>0:
                copy_toBeFurtherMeshed=toBeFurtherMeshed
                toBeFurtherMeshed=[]
                for new_domain in copy_toBeFurtherMeshed:
                    for temp_domain in meshQuadrilateral(new_domain,nodalDistribution):
                        toBeFurtherMeshed.append(temp_domain)
    return;
'''