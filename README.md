# Point Interpolated Newton-Rapshon Meshfree Method (General Solver)

Keywords: Meshfree, Finite Element Analysis, Computational fluid dynamics

###Requirement

Numpy

AutoD - https://github.com/WeiXuanChan/autoD

#####Basic Requirement (standard module)
-pinm module: sys, matplotlib, time, scipy

-discretizer module: stl, mpl_toolkits, matplotlib

###Description

General numerical solver using automatic differentiation and Point Interpolation Method with Meshfree Method.

Enables user defined governing equation (non-linear accepted), material properties and basis.

###Work flow
Import -> Mesh -> Filter and Link Nodes -> Setup Solver (variables, equations, material) -> see results (detailed postprocessing currently in progress)
#####Import
######mainDomain,subDomainList=discretizer.stlImport(FILENAME)
Import *.stl files. mainDomain includes all triangular surfaces and subdomain is selected with graphical aid. (press "add domain" button to add selected triangles to domain subDomainList[1],subDomainList[2],subDomainList[3].. etc, subDomainList[0] is reserved for unselected domain)
#####Mesh
######discretizer.meshSurfaceDomainTriangle(domainList,nodalDistribution)
Mesh domains in domainList according to nodalDistribution (See Glossary).
######innerNodesDomain=discretizer.meshMainDomain(mainDomain,boundaryDomainList,nodalDistribution,meshOuterNode=False)
Generate a domain containing nodes boundaed by list of domains in boundaryDomainList with distribution defined in the function nodalDistribution (See Glossary). Out put from discretizer.stlImport(FILENAME) can be used as mainDomain. Alternatively, it can be generated using mainDomain=createMainDomain(minvert,maxvert,coords), with minvert, maxvert as list of minimun and maximum values of coordinates in coords (list).
#####Filter and Link Nodes
######discretizer.filterNodes(domainList,nodalDistribution,closeness=0.2)
Remove nodes in domainList closer than closeness*(spacing result generated from nodalDistribution function)
######discretizer.primaryLinkNodes(domainList,nodalDistribution,linkIdentifier='primary',closeness=1.5)
Link nodes in domainList closer than closeness*(spacing result generated from nodalDistribution function) and label the linkage as linkIdentifier
#####Setup Solver
######materialObject=pinm.Material(materialName)
Create materialObject with materialName and set property with .setProperty(propertyIdentifier,classOrValue). propertyIdentifier is the name of the property (e.g 'density' etc) and classOrvalue can be a value for constant property or a user defined nodalObjectCallableClass (see glossary).
######pinm.Domain.setMaterial(self,materialIndex,materialObject)
Set materialObject to domain as an integer index (materialIndex).
######basisObject=pinm.Basis(self,basis,basisName,\*specificArgs)
Can be used to turn functions into basisObject. (see basis library for example of function (basis) to be used as input in pinm.Basis(self,basis,basisName,\*specificArgs), specificArgs are inputs after nodalObject)
######pinm.Domain.setBasis(linkIdentifier,basisObject)
set an basisObject as the basis for the link with linkIdentifier., (see basislibrary.py for example)
######pinm.Domain.setVariable(variableIdentifier,init_value)
In the specified domain, set or initialize a variable of the name variableIdentifier and an initial value (init_value).
######pinm.Domain.setVariableLink(variableIdentifier,linkIdentifier)
In the specified domain, set the link to be used for the variable
######pinm.Domain.toggleVariable(variableIdentifier,switch)
Toggle on/off the variables to be solved. Off variables will not change in value (can be used as Dirichlet boundary conditions)
######pinm.Domain.addEquation(eqnClass)
In the specified domain, add a nodalObjectCallableClass as the governing equation (see equationlibrary.py for examples, user-defined subclass (as nodalObjectCallableClass) with inputs (except nodalObject) of classes in equationlibrary.py defined can be used).
######pinm.Domain.toggleEquation(eqnIdentifier,switch)
Toggle on/off the equaitions to be solved.
######solver=pinm.Solver(solverName)
Create solver object
######pinm.Solver.addDomain(domain)
Add domain to solver.
#####Solve
######pinm.Solver.solve()
Solve variables in domains of the solver with the equations.
###Additional control
#####Tracking
######pinm.Track(trackName)
Create trackObject
######pinm.Track.addTracker(self,trackerIdentifier,color,func):
Add a tracker to trackObject. The value calculated with func (autoDObject) is ploted in color during solving iterations.
######pinm.Solver.setTrack(trackObject)
Set a trackObject to a solverObject
######pinm.Solver.setErrorMonintor(self,trackerIdentifier,errorValue,invert=False)
The solver will look at its track object with trackerIdentifier and set it as the error monintor (once the value of the trackerIdentifier is <(when invert=False) errorValue, it is determined as converged and the solver stops)
#####Saving and Loading
######pinm.saveSessionTo(FILENAME)
Save current session to FILENAME. In order to save objects (track,solver,domain etc), they must have a non-empty string as name.
######pinm.saveSession(FILENAME='')
Save current session to the last FILENAME used.
######session=pinm.loadSession(FILENAME)
Load session from FILENAME to current session.
######pinm.setSessionAsCurrent(sessionObject)
set any sessionObject as current session
#####Glossary
######nodalDistribution
a user defined function: outputDict=nodalDistribution(inputDict), where the spacing of the node (e.g {'x':0.1,'y':0.1,'z':0.2}) is calculated with a position (e.g {'x':1,'y':3,'z',2}) input.
######nodalObjectCallableClass
a callable class which take a single input (self,pinm.Node) during initiation and has the function .cal(x,dOrder) (see autoD documentation)
######autoDObject
a class object with function .cal(x,dOrder) (see autoD documentation)

###Function
######pinm.Node
useful properties include:

.pos : return the position of node (dictionary)

.variable[variableIdentifier] : return variable nodal value

.variableCal[variableIdentifier] : return autoDObject (for calculating true variable value with local support domain). The differential of the variable can also be calculated, i.e autoD.Differentiate(nodalObject.variableCal[variableIdentifier],{'x':1}) see equationlibrary for examples.

.material[materialIndex][propertyIdentifier]: return autoDObject of property in material index of node


