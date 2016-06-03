'''
File: pinm.py
Description: Class definition
History:
    Date    Programmer SAR# - Description
    ---------- ---------- ----------------------------
  Author: w. x. chan 29Apr2016           - Created

'''
'''

'''

import numpy as np
import autoD as ad
import sys
from matplotlib import pyplot
from matplotlib.widgets import Button
from matplotlib.collections import LineCollection
import time
from scipy import sparse
from scipy.sparse import linalg
import scipy
global lastOuputTime
lastOuputTime=0.
'''
--------------------Enable save and load session-------------------
Current code unable to pickle nodes after link and maximum recursion depth exceeded while pickling
'''
try:
   import cPickle as pickle
except:
   import pickle
#
class Session:
    def __init__(self,fileName):
        self.nodes=None
        self.saveToFile=fileName
        self.objectDomain={}
        self.objectBasis={}
        self.objectMaterial={}
        self.objectTrack={}
        self.objectSolver={}
    def addDomain(self,domainName,domainObject):
        self.objectDomain[domainName]=domainObject
        return;
    def addEquation(self,eqnName,eqnObject):
        self.objectEquation[eqnName]=eqnObject
        return;
    def addBasis(self,name,object):
        self.objectBasis[name]=object
        return;
    def addMaterial(self,name,object):
        self.objectMaterial[name]=object
        return;
    def addTrack(self,name,object):
        self.objectTrack[name]=object
        return;
    def addSolver(self,name,object):
        self.objectSolver[name]=object
        return;
    def saveTo(self,fileName):
        self.saveToFile=fileName
global currentSession
currentSession=Session('')
def saveSessionTo(fileName):
    global currentSession
    currentSession.saveToFile=fileName
def saveSession(fileName='',fast=False):
    global currentSession
    if fileName!='':
        saveSessionTo(fileName)
    #sanitize linker to refrain from recursion depth limit
    
    if fast:
        currentSession.nodes=None
    else:
        nodes=[]
        count=0
        for domain in currentSession.objectDomain.values():
            for node in domain.nodes():
                node.setIndex(count)
                nodes.append(node)
                count+=1
        print("Saving progress : 'converting nodes'")
        numOfNodes=len(nodes)
        for m in range(numOfNodes):
            updateProgress(float(m+1)/numOfNodes)
            for linkIdentifier in nodes[m].link:
                for n in range(len(nodes[m].link[linkIdentifier])):
                    if type(nodes[m].link[linkIdentifier][n]) is Node:
                        ind=nodes[m].link[linkIdentifier][n].ind
                        nodes[m].link[linkIdentifier][n]=ind
        currentSession.nodes=nodes
    with open(currentSession.saveToFile, "wb") as file:
        pickle.dump(currentSession, file)
    if currentSession.nodes!=None:
        for domain in currentSession.objectDomain.values():
            for node in domain.nodes():
                for linkIdentifier in node.link:
                    for n in range(len(node.link[linkIdentifier])):
                        if type(node.link[linkIdentifier][n]) is int:
                            node.link[linkIdentifier][n]=currentSession.nodes[node.link[linkIdentifier][n]]
    currentSession.nodes=None
    
def loadSession(file):
    global currentSession
    with open(file, "rb") as f:
        result = pickle.load(f)
    currentSession=result
    if currentSession.nodes!=None:
        for domain in currentSession.objectDomain.values():
            for node in domain.nodes():
                for linkIdentifier in node.link:
                    for n in range(len(node.link[linkIdentifier])):
                        if type(node.link[linkIdentifier][n]) is int:
                            node.link[linkIdentifier][n]=currentSession.nodes[node.link[linkIdentifier][n]]
    currentSession.nodes=None
    return result
def setSessionAsCurrent(sessionObject):
    global currentSession
    currentSession=sessionObject
    return;
def updateProgress(progress):
    global lastOuputTime
    if ((time.time()-lastOuputTime)>30.) or (progress >= 1):
        barLength = 10 # Modify this to change the length of the progress bar
        status = ""
        if isinstance(progress, int):
            progress = float(progress)
        if not isinstance(progress, float):
            progress = 0
            status = "error: progress var must be float\r\n"
        if progress >= 1:
            progress = 1
            status = "Done...\r\n"
        if progress < 0:
            text = str(-progress)
        else:
            block = int(round(barLength*progress))
            text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), progress*100, status)
        sys.stdout.write(text)
        sys.stdout.flush()
        lastOuputTime=time.time()
'''
--------------------End of Session class and functions-------------------
'''

class Node:
    def __init__(self,pos,norm={}): #pos and norm are type dict pos={'x':1,'y':2}
        self.type='node'
        self.pos=pos
        self.variable={}
        self.variableSolveToggle={}
        self.deltaVariable={}
        self.variableCal={}
        self.eqn={}
        self.eqnSolveToggle={}
        self.link={} #include self if basis used requires
        self.linkBasis={}
        self.momentMatrix={}
        self.transformation={}
        self.choleskyDecomposition={}
        self.norm=norm #normal, vector, for use of surface equation
        self.normLink=''
        self.material=[]
        self.variableLink={}
        self.domain=None
    def setIndex(self,ind):#temporary
        self.ind=ind
    def setDomain(self,domain):
        self.domain=domain
        return;
    def setDeltaVariable(self,variableIdentifier,value):
        self.deltaVariable[variableIdentifier]=value
        return;
    def resetDeltaVariable(self):
        for variableIdentifier in self.deltaVariable:
            self.deltaVariable[variableIdentifier]=0.
        return;
    def addVariable(self,variableIdentifier,init_value):
        #variable_name,init_value are lists
        if type(variableIdentifier) is list:
            for n in range(len(variableIdentifier)):
                if variableIdentifier[n] not in self.variable:
                    self.variableCal[variableIdentifier[n]]=ad.Function(pointInterpolationMethod,self,variableIdentifier[n])
                    self.deltaVariable[variableIdentifier[n]]=0.
                    self.variableSolveToggle[variableIdentifier[n]]=True
                self.variable[variableIdentifier[n]]=init_value[n]
                
        else:
            if variableIdentifier not in self.variable:
                self.variableCal[variableIdentifier]=ad.Function(pointInterpolationMethod,self,variableIdentifier)
                self.deltaVariable[variableIdentifier]=0.
                self.variableSolveToggle[variableIdentifier]=True
            self.variable[variableIdentifier]=init_value
            
        return;
    def removeVariable(self,variableIdentifier):
        if type(variableIdentifier) is list:
            for n in range(len(variableIdentifier)):
                del self.variable[variableIdentifier[n]]
                del self.deltaVariable[variableIdentifier[n]]
                del self.variableLink[variableIdentifier[n]]
                del self.variableCal[variableIdentifier[n]]
                del self.variableSolveToggle[variableIdentifier[n]]
        else:
            del self.variable[variableIdentifier]
            del self.variableLink[variableIdentifier]
            del self.variableCal[variableIdentifier]
            del self.variableSolveToggle[variableIdentifier]
        return;
    def toggleVariable(self,variableIdentifier,switch):
        if type(variableIdentifier) is list:
            for n in range(len(variableIdentifier)):
                if switch:
                    self.variableSolveToggle[variableIdentifier[n]]=True
                else:
                    self.variableSolveToggle[variableIdentifier[n]]=False
        else:
            if switch:
                self.variableSolveToggle[variableIdentifier]=True
            else:
                self.variableSolveToggle[variableIdentifier]=False
    def updateNewVariable(self):
        for variableIdentifier in self.variable:
            if self.variableSolveToggle:
                if float('-inf')<self.deltaVariable[variableIdentifier]<float('inf'):
                    self.variable[variableIdentifier]=self.variable[variableIdentifier]+self.deltaVariable[variableIdentifier]
                else:
                    print('Waring! ',self.deltaVariable[variableIdentifier],' encountered in deltaVariable')
        return;
    def setNorm(self,norm):
        self.norm=norm
        return;
    def setNormLink(self,linkIdentifier):
        self.normLink=linkIdentifier
        return;
    def setLinkBasis(self,linkIdentifier,basis):
        self.linkBasis[linkIdentifier]=basis
        return;
    def addMaterial(self,materialIndex,material):
        while len(self.material)<=materialIndex:
            self.material.append(None)
        self.material[materialIndex]=Material(material.name) #materialIndex is use for differentiating internal and external of surface etc.
        #fixed properties to node
        for propertyIdentifier in material.properties:
            if callable(material.properties[propertyIdentifier]):
                tempcall=material.properties[propertyIdentifier](self)
                self.material[materialIndex].setProperty(propertyIdentifier,tempcall)
            else:
                self.material[materialIndex].setProperty(propertyIdentifier,ad.Constant(material.properties[propertyIdentifier]))
        return;
    def addEquation(self,eqnClass):
        tempcall=eqnClass(self)
        self.eqn[tempcall.name]=tempcall
        self.eqnSolveToggle[tempcall.name]=True
        return;
    def toggleEquation(self,eqnIdentifier,switch):
        if type(eqnIdentifier) is list:
            for n in range(len(eqnIdentifier)):
                if switch:
                    self.eqnSolveToggle[eqnIdentifier[n]]=True
                else:
                    self.eqnSolveToggle[eqnIdentifier[n]]=False
        else:
            if switch:
                self.eqnSolveToggle[eqnIdentifier]=True
            else:
                self.eqnSolveToggle[eqnIdentifier]=False
    def setPos(self,value):
        self.pos=pos #pos is a dict pos={'x':1,'y':2}
        return;
    def addLink(self,linkIdentifier,nodalObject):
        if linkIdentifier not in self.link:
            self.link[linkIdentifier]=[]
            self.momentMatrix[linkIdentifier]=np.zeros(0)
            self.transformation[linkIdentifier]=np.zeros(0)
            self.choleskyDecomposition[linkIdentifier]=np.zeros(0)
        if type(nodalObject) is list:
            for n in range(len(nodalObject)):
                if nodalObject[n] not in self.link[linkIdentifier]:
                    self.link[linkIdentifier].append(nodalObject[n])
        else:
            if nodalObject not in self.link[linkIdentifier]:
                self.link[linkIdentifier].append(nodalObject)
        return;
    def removeLink(self,linkIdentifier,nodalObject=None):
        if linkIdentifier in self.link:
            if type(nodalObject) is list:
                for n in range(len(nodalObject)):
                    self.link[linkIdentifier].remove(nodalObject[n])
            elif nodalObject==None:
                del self.link[linkIdentifier]
                del self.momentMatrix[linkIdentifier]
                del self.transformation[linkIdentifier]
                del self.choleskyDecomposition[linkIdentifier]
                if linkIdentifier in self.linkBasis:
                    del self.linkBasis[linkIdentifier]
            else:
                self.link[linkIdentifier].remove(nodalObject)
        else:
            print('"',linkIdentifier,'" does not exists.')
        return;
    def setVariableLink(self,variableIdentifier,linkIdentifier):
        if type(variableIdentifier) is list:
            for n in range(len(variableIdentifier)):
                self.variableLink[variableIdentifier[n]]=linkIdentifier[n]
        else:
            self.variableLink[variableIdentifier]=linkIdentifier
        return;
    def updateShapeFunction(self,linkIdentifier):
        emptydOrder={}
        newMomentMatrix_temp=[]
        self.linkBasis[linkIdentifier].setNode(self)
        for n in self.link[linkIdentifier]:
            newMomentMatrix_temp.append(self.linkBasis[linkIdentifier].cal(n.pos,emptydOrder))
        newMomentMatrix=np.vstack(tuple(newMomentMatrix_temp))
        if not(np.array_equal(self.momentMatrix[linkIdentifier],newMomentMatrix)):
            self.momentMatrix[linkIdentifier]=newMomentMatrix
            shape=newMomentMatrix.shape
            if shape[0]==shape[1]:
                self.transformation[linkIdentifier]=np.transpose(self.momentMatrix[linkIdentifier])
                leastSQ=np.dot(self.transformation[linkIdentifier],self.momentMatrix[linkIdentifier])
                self.choleskyDecomposition[linkIdentifier]=scipy.linalg.cholesky(leastSQ)
                #self.shapeFunctionMatrix[linkIdentifier]=np.linalg.inv(self.momentMatrix[linkIdentifier])
            elif shape[0]>shape[1]:
                maxDistance={}
                distance=[]
                for n in range(len(self.link[linkIdentifier])):
                    distance.append({})
                    for coord in self.pos:
                        distance[n][coord]=np.absolute(self.link[linkIdentifier][n].pos[coord]-self.pos[coord])
                        if coord not in maxDistance:
                            maxDistance[coord]=distance[-1][coord]
                        else:
                            maxDistance[coord]=max(maxDistance[coord],distance[n][coord])
                weight=[]
                for n in range(len(self.link[linkIdentifier])):
                    tempWeight=0.
                    for coord in self.pos:
                        tempWeight+=(distance[n][coord]/maxDistance[coord])**2.
                    weight.append(np.exp(-tempWeight/0.2))
                weightMatrix=np.diag(weight)
                self.transformation[linkIdentifier]=np.dot(np.transpose(self.momentMatrix[linkIdentifier]),weightMatrix)
                leastSQ=np.dot(self.transformation[linkIdentifier],self.momentMatrix[linkIdentifier])
                self.choleskyDecomposition[linkIdentifier]=scipy.linalg.cholesky(leastSQ)
                #self.shapeFunctionMatrix[linkIdentifier]=np.dot(invLeastSQ,transformation)
            else:
                print('Not enough support nodes')
            
        return;
  
class Basis(ad.Function):
    def __init__(self,basis,basisName,*specificArgs):
        self.type='basis'
        self.name=basisName
        self.specificArgs=specificArgs
        if basisName not in currentSession.objectBasis:
            currentSession.addBasis(basisName,self)
        ad.Function.__init__(self,basis,None,*specificArgs)
    def changeSpecificArgs(self,*specificArgs):
        self.specificArgs=specificArgs
    def setNode(self,nodalObject):
        specificArgs=self.specificArgs
        self.changeArgs(nodalObject,*specificArgs)
        return;
        
class Domain:
    def __init__(self,domainName,norm={}):
        self.type='domain'
        self.name=domainName
        if domainName!='':
            if domainName not in currentSession.objectDomain:
                currentSession.addDomain(domainName,self)
        self.subDomain=[] #list of subdomains or nodes
        self.superDomain=None
        self.normalVector=norm
        self.maxDistance={}
        self.vertices=None
        self.pos={}
    def setCentroid(self,vertices): #note this is just the average of the vertices not the true centroid
        self.vertices=vertices
        maxDistance={}
        centroid={}
        numOfVertices=len(vertices)
        for coord in vertices[0]:
            total=0.
            for vert in vertices:
                total+=vert[coord]
            centroid[coord]=total/numOfVertices
        self.pos=centroid
        for coord in vertices[0]:
            maxDistance[coord]=0.
            for vert in vertices:
                temp_distance=np.absolute(vert[coord]-centroid[coord])
                if temp_distance>maxDistance[coord]:
                    maxDistance[coord]=temp_distance
        self.maxDistance=maxDistance
        return;
    def addNode(self,nodalObjects):
        if type(nodalObjects) is list:
            for n in range(len(nodalObjects)):
                if nodalObjects[n] not in self.subDomain:
                    self.subDomain.append(nodalObjects[n])
                    if type(nodalObjects[n]) is Node:
                        nodalObjects[n].setDomain(self)
                    else:
                        nodalObjects[n].setSuperDomain(self)
        else:
            if nodalObjects not in self.subDomain:
                self.subDomain.append(nodalObjects)
                if type(nodalObjects) is Node:
                    nodalObjects.setDomain(self)
                else:
                    nodalObjects.setSuperDomain(self)
        return;
    def removeNode(self,nodes):
        if type(nodes) is list:
            for n in range(len(nodes)):
                self.subDomain.remove(nodes[n])
        else:
            self.subDomain.remove(nodes)
        return;
    def setDomainName(self,domainName):
        self.name=domainName
        if domainName not in currentSession.objectDomain:
            currentSession.addDomain(domainName,self)
        return;
    def setSuperDomain(self,superDomain):
        if self.superDomain !=None:
            print('Warning! Overwriting superDomain.')
            self.superDomain.removeNode(self)
        self.superDomain=superDomain
        return;
    def nodes(self):
        nodes=[]
        for subDomain in self.subDomain:
            if type(subDomain) is Node:
                nodes.append(subDomain)
            else:
                for temp_nodes in subDomain.nodes():
                    nodes.append(temp_nodes)
        return nodes
    def setMaterial(self,materialIndex,materialObject):
        for node in self.nodes():
            node.addMaterial(materialIndex,materialObject)
        return;
    def addEquation(self,eqnClass):
        if type(eqnClass) is list:
            for n in range(len(eqnClass)):
                for node in self.nodes():
                    node.addEquation(eqnClass[n])
        else:
            for node in self.nodes():
                node.addEquation(eqnClass)
    def toggleEquation(self,eqnIdentifier,switch):
        for node in self.nodes():
            node.toggleEquation(eqnIdentifier,switch)
    def setBasis(self,linkIdentifier,basisObject):
        for node in self.nodes():
            node.setLinkBasis(linkIdentifier,basisObject)
    def setVariable(self,variableIdentifier,init_value):
        for node in self.nodes():
            node.addVariable(variableIdentifier,init_value)
    def setVariableLink(self,variableIdentifier,linkIdentifier):
        for node in self.nodes():
            node.setVariableLink(variableIdentifier,linkIdentifier)
    def toggleVariable(self,variableIdentifier,switch):
        for node in self.nodes():
            node.toggleVariable(variableIdentifier,switch)
class Material:
    def __init__(self,materialName):
        self.type='material'
        self.name=materialName
        if materialName not in currentSession.objectMaterial:
            currentSession.addMaterial(materialName,self)
        self.properties={}
    def setProperty(self,propertyIdentifier,classOrValue):
        self.properties[propertyIdentifier]=classOrValue
        return;
    def removeProperty(self,propertyIdentifier):
        del self.properties[propertyIdentifier]
        return;

class Track:
    def __init__(self,trackName):
        self.type='track'
        self.tracker={}
        self.trackerColor={}
        self.trackerToggle={}
        self.name=trackName
        if trackName not in currentSession.objectTrack:
            currentSession.addTrack(trackName,self)
        self.value={}
    def addTracker(self,trackerIdentifier,color,func): #func=ad.Function()
        self.tracker[trackerIdentifier]=func
        self.value[trackerIdentifier]=0.
        self.trackerColor[trackerIdentifier]=color
        self.trackerToggle[trackerIdentifier]=True
        return;
    def removeTracker(self,trackerIdentifier):
        del self.tracker[trackerIdentifier]
        del self.value[trackerIdentifier]
        del self.trackerColor[trackerIdentifier]
        del self.trackerToggle[trackerIdentifier]
        return;
    def toggleTracker(self,trackerIdentifier,switch):
        if type(trackerIdentifier) is list:
            for n in range(len(trackerIdentifier)):
                if switch:
                    self.trackerToggle[trackerIdentifier[n]]=True
                else:
                    self.trackerToggle[trackerIdentifier[n]]=False
        else:
            if switch:
                self.trackerToggle[trackerIdentifier]=True
            else:
                self.trackerToggle[trackerIdentifier]=False
    def update(self):
        for trackerIdentifier in self.tracker:
            self.value[trackerIdentifier]=self.tracker[trackerIdentifier].cal({},{})

class Solver:
    def __init__(self,solverName):
        self.type='solver'
        self.name=solverName
        if solverName not in currentSession.objectSolver:
            currentSession.addSolver(solverName,self)
        self.domain={}
        self.track=None
        self.stop=False
        self.errorMonintorTracker=''
        self.errorMonintorErrorValue=0
        self.errorMonintorInvert=False
        self.nodes=[]
        self.jMatrixRow=[]
        self.jMatrixCol=[]
        self.jMatrixData=[]
        self.jMatrixShape=None
        self.varIndex=[]
        self.indexVar=[]
        self.equationIndex=[]
        self.equationNodeIndex=[]
        self.fxValue=[]
    def reset(self):
        self.stop=False
    def setTrack(self,trackObject): #input Track object
        self.track=trackObject
        return;
    def addDomain(self,domainObject): #input Domain object
        self.domain[domainObject.name]=domainObject
        return;
    def removeDomain(self,domainName):
        del self.domain[domainName]
        return;
    def addEquationIntoDomains(self,equationClass,domainName):
        if type(domainName) is list:
            for n in range(len(domainName)):
                self.domain[domainName[n]].addEquation(self,eqnClass)
        else:
            self.domain[domainName].addEquation(self,eqnClass)
        return;
    def setErrorMonintor(self,trackerIdentifier,errorValue,invert=False):
        self.errorMonintorTracker=trackerIdentifier
        self.errorMonintorErrorValue=errorValue
        self.errorMonintorInvert=invert
    def stopCheck(self):
        if self.errorMonintorInvert==False:
            if self.track.value[self.errorMonintorTracker]<self.errorMonintorErrorValue:
                self.stop=True
        elif self.errorMonintorInvert==True:
            if self.track.value[self.errorMonintorTracker]>self.errorMonintorErrorValue:
                self.stop=True
    '''
    def iterate_stationary(self):
        #Newton-Raphson iteration to find variable values
        #Gauss-Seidel Iteration Methods to find delta of variable values in Newton-Raphson iteration
        emptydOrder={}
        
        #determine Gauss-Seidel Iteration coefficients ie the Jacobian matrix
        ###########################its not always converging!!!!
        jMatrix={}
        fxValue={}
        for domainIdentifier in self.domain:
            for node in self.domain[domainIdentifier].nodes():
                node.resetDeltaVariable()
        print("Solving Progress : 'Jacobian coefficients'")
        for domainIdentifier in self.domain:
            jMatrix[domainIdentifier]=[]
            fxValue[domainIdentifier]=[]
            print("Solving Progress : 'Domain ",domainIdentifier,"'")
            nodes=self.domain[domainIdentifier].nodes()
            nodeNum=len(nodes)
            for n in range(len(nodes)):
                jMatrix[domainIdentifier].append({})
                fxValue[domainIdentifier].append({})
                for equationName in nodes[n].eqn:
                    if nodes[n].eqnSolveToggle[equationName]:
                        tempfxValue=nodes[n].eqn[equationName].cal(nodes[n].pos,emptydOrder)
                        fxValue[domainName][n][equationName]=tempfxValue
                        if not(float('-inf')<tempfxValue<float('inf')):
                            print('Waring! ',tempfxValue,' encountered for fxValue')
                        jMatrix[domainIdentifier][n][equationName]={}
                        for variableIdentifier in nodes[n].variable:
                            for link in nodes[n].link[nodes[n].variableLink[variableIdentifier]]:
                                if link.variableSolveToggle[variableIdentifier]:
                                    var=link.variableCal[variableIdentifier]
                                    tempjMatrixValue=nodes[n].eqn[equationName].cal(nodes[n].pos,{var:1})
                                    jMatrix[domainIdentifier][n][equationName][var]=tempjMatrixValue
                                    if not(float('-inf')<tempjMatrixValue<float('inf')):
                                        print('Waring! ',tempjMatrixValue,' encountered in calculating Jacobian')
                updateProgress(float(n+1)/nodeNum)
        #Gauss-Seidel Iteration Method with alteration
        converged=False
        iterCount=0
        print('Solving Progress : Gauss-Seidel Iteration Method')
        while not(converged):
            converged=True
            for domainIdentifier in self.domain:
                nodes=self.domain[domainIdentifier].nodes()
                nodeNum=len(nodes)
                for n in range(len(nodes)):
                    for equationName in nodes[n].eqn:
                        if nodes[n].eqnSolveToggle[equationName]:
                            update={}
                            for variableIdentifier in nodes.variable:
                                if node.variableSolveToggle[variableIdentifier]:
                                    var=node.variableCal[variableIdentifier]
                                    if jMatrix[domainIdentifier][n][equationName][var]!=0.:
                                        sumtotal=0.
                                        for variable in jMatrix[domainIdentifier][n][equationName]:
                                            linkNode,linkVariableIdentifier=variable.checkArgs()
                                            if variable!=var and linkNode.variableSolveToggle[linkVariableIdentifier]:
                                                sumtotal+=jMatrix[domainIdentifier][n][equationName][variable]*linkNode.deltaVariable[linkVariableIdentifier]
                                        update[var]=(fxValue[domainName][n][equationName]-sumtotal)/jMatrix[domainIdentifier][n][equationName][var]
                            for var in update:
                                oldVar=var.checkArgs()[0].deltaVariable[var.checkArgs()[1]]
                                var.checkArgs()[0].setDeltaVariable(var.checkArgs()[1],update[var])
                                if update[var]!=0:
                                    if np.absolute(np.absolute(oldVar/update[var]-1.))>self.errorMonintorErrorValue:
                                        converged=False
                                elif np.absolute(np.absolute(oldVar))>self.errorMonintorErrorValue:
                                    converged=False
            iterCount-=1
            updateProgress(iterCount)
        print(' -------Converged')
        #update variables
        for domain in self.domain.values():
            for node in domain.nodes():
                node.updateNewVariable()  
        return;
    '''
    def calFx(self,start=0,stop=0):
        #calculate the Jacobian matrix
        if stop<=start:
            stop=self.jMatrixShape[0]
        tempNum=stop-start
        print('Solving Progress : Updating Function value')
        result=[]
        for n in range(self.jMatrixShape[0]):
            tempfxValue=self.equationIndex[n].cal(self.nodes[self.equationNodeIndex[n]].pos,{})
            self.fxValue[n]=-tempfxValue
            result.append(-tempfxValue)
            if not(float('-inf')<tempfxValue<float('inf')):
                print('Waring! ',tempfxValue,' encountered for fxValue')
            updateProgress(float(n+1)/tempNum)
        return result
    def calJacobianMatrix(self,initialRun=False,start=0,stop=0):
        if stop<=start:
            stop=len(self.jMatrixRow)
        #calculate the Jacobian matrix
        tempNum=stop-start
        print('Solving Progress : Updating Jacobian Matrix')
        result=[]
        for n in range(start,stop):
            if self.equationIndex[self.jMatrixRow[n]].nonLinear or initialRun:
                tempjMatrixValue=self.equationIndex[self.jMatrixRow[n]].cal(self.nodes[self.equationNodeIndex[self.jMatrixRow[n]]].pos,{self.indexVar[self.jMatrixCol[n]]:1})
                self.jMatrixData[n]=tempjMatrixValue
                result.append(tempjMatrixValue)
                if not(float('-inf')<tempjMatrixValue<float('inf')):
                    print('Waring! ',tempjMatrixValue,' encountered in calculating Jacobian')
            updateProgress(float(n+1)/tempNum)
        return result
    def iterate(self):
        jMatrixSparse=sparse.csr_matrix((self.jMatrixData, (self.jMatrixRow, self.jMatrixCol)), shape=self.jMatrixShape)
        fxValueMat=np.array(self.fxValue)
        sumAll=0.
        for allValue in self.fxValue:
            sumAll+=np.absolute(np.absolute(allValue))**2.
        print(sumAll)
        tempResult=linalg.lsqr(jMatrixSparse,fxValueMat)
        newDeltaVar=tempResult[0]
        #update variables
        for n in range(len(self.varIndex)):
            for variableIdentifier in self.varIndex[n]:
                self.nodes[n].setDeltaVariable(variableIdentifier,newDeltaVar[self.varIndex[n][variableIdentifier]])
            self.nodes[n].updateNewVariable()  
        return;
    def resetJMatrix(self):
        self.nodes=[]
        indCount=0
        for domainIdentifier in self.domain:
            for node in self.domain[domainIdentifier].nodes():
                node.ind=indCount
                self.nodes.append(node)
                indCount+=1
        self.varIndex=[]
        for node in self.nodes:
            node.resetDeltaVariable()
            self.varIndex.append({})
        self.jMatrixRow=[]
        self.jMatrixCol=[]
        self.jMatrixData=[]
        self.jMatrixShape=None
        self.indexVar=[]
        self.equationIndex=[]
        self.equationNodeIndex=[]
        self.fxValue=[]
        if self.errorMonintorTracker=='' or self.errorMonintorTracker=='default error track':
            if self.track==None:
                newTrack=Track('error')
                self.setTrack(newTrack)
            defaultErrorTrack=ad.Function(defaultError,self.nodes)
            self.track.addTracker('default error track','k',defaultErrorTrack)
            self.setErrorMonintor('default error track',0.001)
        equationCount=0
        varCount=0
        nodeNum=len(self.nodes)
        for n in range(len(self.nodes)):
            for equationName in self.nodes[n].eqn:
                if self.nodes[n].eqnSolveToggle[equationName]:
                    self.equationIndex.append(self.nodes[n].eqn[equationName])
                    self.equationNodeIndex.append(n)
                    self.fxValue.append(0.)
                    for variableIdentifier in  self.nodes[n].variable:
                        for link in self.nodes[n].link[self.nodes[n].variableLink[variableIdentifier]]:
                            if link.variableSolveToggle[variableIdentifier]:
                                var=link.variableCal[variableIdentifier]
                                if variableIdentifier not in self.varIndex[link.ind]:
                                    self.varIndex[link.ind][variableIdentifier]=varCount
                                    self.indexVar.append(link.variableCal[variableIdentifier])
                                    varCount+=1
                                self.jMatrixCol.append(self.varIndex[link.ind][variableIdentifier])
                                self.jMatrixRow.append(equationCount)
                                self.jMatrixData.append(0.)
                    equationCount+=1
        self.jMatrixShape=(equationCount, varCount)
    def solve(self,fullSolve=True):
        countplot={}
        trackplot={}
        if len(self.nodes)==0:
            initRun=True
            self.resetJMatrix()
        else:
            initRun=False
        for trackIdentifier in self.track.value:
            countplot[trackIdentifier]=[]
            trackplot[trackIdentifier]=[]
        if not(fullSolve):
            self.iterate()
            self.track.update()
            return self.track.value[self.errorMonintorTracker]
        figure = pyplot.figure()
        pyplot.subplots_adjust(bottom=0.2)
        pyplot.xlabel('Iteration')
        pyplot.ylabel(self.errorMonintorTracker)
        callback = solverControl(self)
        axpause = pyplot.axes([0.7, 0.05, 0.1, 0.075])
        axcont = pyplot.axes([0.81, 0.05, 0.1, 0.075])
        axstop = pyplot.axes([0.05, 0.05, 0.15, 0.075])
        bstop = Button(axstop, 'Stop')
        bstop.on_clicked(callback.stop)
        bpause = Button(axpause, 'Pause')
        bpause.on_clicked(callback.pause)
        bcont = Button(axcont, 'Continue')
        bcont.on_clicked(callback.cont)
    
        count=0
        ymin=float('inf')
        ymax=10.**-30.
        pyplot.ion()
        ax = pyplot.axes()
        ax.set_yscale('log')
        while self.stop==False:
            while callback.pausing==False and self.stop==False:
                if count==0 and initRun:
                    self.calJacobianMatrix(True)
                else:
                    self.calJacobianMatrix()
                self.calFx()
                self.iterate()
                self.track.update()
                for trackIdentifier in self.track.value:
                    countplot[trackIdentifier].append(count)
                    trackplot[trackIdentifier].append(self.track.value[trackIdentifier])
                    ymin=min(ymin,self.track.value[trackIdentifier])
                    ymax=max(ymax,self.track.value[trackIdentifier])
                    if count>0:
                        line_segment = LineCollection([[(countplot[trackIdentifier][-2],trackplot[trackIdentifier][-2]),(countplot[trackIdentifier][-1],trackplot[trackIdentifier][-1])]],
                                       linewidths=(0.5, 1, 1.5, 2),
                                       linestyles='solid',
                                       colors=self.track.trackerColor[trackIdentifier])
                        ax.add_collection(line_segment)
                if count>0:
                    ax.set_xlim(0, count)
                    if ymin<(10.**-30.):
                        ymin=10.**-30.
                    ax.set_ylim(ymin, ymax)
                    pyplot.draw()
                    pyplot.pause(0.05)
                self.stopCheck()
                count+=1
            pyplot.pause(0.05)
        pyplot.show()
def defaultError(x,dOrder,nodeList):
    error=0.
    maxVar={}
    for node in nodeList:
        for var in node.variable:
            if var in maxVar:
                maxVar[var]=max(maxVar[var],np.absolute(np.absolute(node.variable[var])))
            else:
                maxVar[var]=np.absolute(node.variable[var])
    for node in nodeList:
        for var in node.variable:
            if maxVar[var]!=0:
                error=max(error,np.absolute(np.absolute(node.deltaVariable[var]))/maxVar[var])
            else:
                error=1.
    return error
class solverControl:
    def __init__(self,solverObject):
        self.solver=solverObject
        self.pausing=False
    def stop(self, event):
        self.stopFunc()
    def pause(self, event):
        self.pauseFunc()
    def cont(self, event):
        self.contFunc()
    def stopFunc(self):
        self.solver.stop=True
    def pauseFunc(self):
        self.pausing=True
    def contFunc(self):
        self.pausing=False

'''
------------------------Method-----------------------------------------
'''
def pointInterpolationMethod(x,dOrder,nodalObject,variableIdentifier):
    #check basis variable
    linkIdentifier=nodalObject.variableLink[variableIdentifier]
    nodalObject.linkBasis[linkIdentifier].changeArgs(nodalObject)
    nodalObject.updateShapeFunction(linkIdentifier)
    
    for n in range(len(nodalObject.link[linkIdentifier])):
        if nodalObject.link[linkIdentifier][n].variableCal[variableIdentifier] in dOrder:
            valMatrix_temp=np.zeros(len(nodalObject.link[linkIdentifier]))
            valMatrix_temp[n]=1.
            valMatrix=np.transpose(valMatrix_temp)
            new_dOrder=dOrder.copy()
            new_dOrder[nodalObject.link[linkIdentifier][n].variableCal[variableIdentifier]]=dOrder[nodalObject.link[linkIdentifier][n].variableCal[variableIdentifier]]-1
            break
    else:
        valMatrix_temp=[]
        for m in nodalObject.link[linkIdentifier]:
            valMatrix_temp.append(m.variable[variableIdentifier])
        valMatrix=np.vstack(tuple(valMatrix_temp))
        new_dOrder=dOrder.copy()     
    basisVal=nodalObject.linkBasis[linkIdentifier].cal(x,new_dOrder)
    tempValMatrix=np.dot(nodalObject.transformation[linkIdentifier],valMatrix)
    basisCoef=scipy.linalg.cho_solve((nodalObject.choleskyDecomposition[linkIdentifier],False),tempValMatrix)
    result=np.dot(basisVal,basisCoef)
    return np.sum(result) #ensure return float

