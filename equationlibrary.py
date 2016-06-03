import numpy as np
import autoD as ad

class Variable:
    def __init__(self,eqnName,nodalObject,variableIdentifier,add_dOrder):
        self.name=eqnName
        self.dependent=['ALL']
        self.nonLinear=False
        self.result=ad.Differentiate(nodalObject.variableCal[variableIdentifier],add_dOrder)
    def cal(self,x,dOrder):
        return self.result.cal(x,dOrder)
class Incompressible4d:
    def __init__(self,nodalObject):
        self.name='incompressible4d'
        self.dependent=['ALL']
        self.nonLinear=False
        desityChange=ad.Differentiate(nodalObject.material[0].properties['density'],{'t':1})
        flowx=ad.Differentiate(nodalObject.material[0].properties['density']*nodalObject.variableCal['U'],{'x':1})
        flowy=ad.Differentiate(nodalObject.material[0].properties['density']*nodalObject.variableCal['V'],{'y':1})
        flowz=ad.Differentiate(nodalObject.material[0].properties['density']*nodalObject.variableCal['W'],{'z':1})
        self.result=desityChange+flowx+flowy+flowz
    def cal(self,x,dOrder):
        return self.result.cal(x,dOrder)
class Incompressible4dGeneral:
    def __init__(self,nodalObject,variableList,coordList):
        self.name='incompressible4d'
        self.dependent=['ALL']
        self.nonLinear=False
        desityChange=ad.Differentiate(nodalObject.material[0].properties['density'],{'t':1})
        flow=[]
        for n in range(len(coordList)):
            flow.append(ad.Differentiate(nodalObject.material[0].properties['density']*nodalObject.variableCal[variableList[n]],{coordList[n]:1}))
        flow.append(desityChange)
        self.result=ad.Addition(flow)
    def cal(self,x,dOrder):
        return self.result.cal(x,dOrder)
class ConvectionDiffusion4d:
    def __init__(self,eqnName,nodalObject,solvingVariableIdentifier,diffusivityIdentifier,internalSourceEquation=None,externalSourceEquation=None,nonLinear=False):
        self.name=eqnName
        self.dependent=['ALL']
        self.nonLinear=nonLinear
        variation=ad.Differentiate(nodalObject.variableCal[solvingVariableIdentifier],{'t':1})
        convectionU=nodalObject.variableCal['U']*ad.Differentiate(nodalObject.variableCal[solvingVariableIdentifier],{'x':1})
        convectionV=nodalObject.variableCal['V']*ad.Differentiate(nodalObject.variableCal[solvingVariableIdentifier],{'y':1})
        convectionW=nodalObject.variableCal['W']*ad.Differentiate(nodalObject.variableCal[solvingVariableIdentifier],{'z':1})
        diffusion=ad.Multiply([ad.Addition([ad.Differentiate(nodalObject.variableCal[solvingVariableIdentifier],{'x':2}),
                                            ad.Differentiate(nodalObject.variableCal[solvingVariableIdentifier],{'y':2}),
                                            ad.Differentiate(nodalObject.variableCal[solvingVariableIdentifier],{'z':2})]),
                               -1.,
                               nodalObject.material[0].properties[diffusivityIdentifier]])
        if internalSourceEquation!=None:
            internalSource=internalSourceEquation(nodalObject)
        else:
            internalSource=0.
        if externalSourceEquation!=None:
            externalSource=ad.Multiply([-1.,externalSourceEquation(nodalObject)])
        else:
            externalSource=0.
        self.result=ad.Addition([variation,convectionU,convectionV,convectionW,diffusion,internalSource,externalSource])
    def cal(self,x,dOrder):
        return self.result.cal(x,dOrder)
class Incompressible3d:
    def __init__(self,nodalObject):
        self.name='incompressible3d'
        self.dependent=['ALL']
        self.nonLinear=False
        flowx=ad.Differentiate(ad.Multiply([nodalObject.material[0].properties['density'],nodalObject.variableCal['U']]),{'x':1})
        flowy=ad.Differentiate(ad.Multiply([nodalObject.material[0].properties['density'],nodalObject.variableCal['V']]),{'y':1})
        flowz=ad.Differentiate(ad.Multiply([nodalObject.material[0].properties['density'],nodalObject.variableCal['W']]),{'z':1})
        self.result=ad.Addition([flowx,flowy,flowz])
    def cal(self,x,dOrder):
        return self.result.cal(x,dOrder)
class IncompressibleLinearizedNavierStokes4d:
    def __init__(self,eqnName,nodalObject,solvingVariableIdentifier,diffusivityIdentifier,internalSourceEquation=None,externalSourceEquation=None):
        self.name=eqnName
        self.dependent=['ALL']
        self.nonLinear=False
        variation=ad.Differentiate(nodalObject.variableCal[solvingVariableIdentifier],{'t':1})
        diffusion=ad.Multiply([ad.Addition([ad.Differentiate(nodalObject.variableCal[solvingVariableIdentifier],{'x':2}),
                                            ad.Differentiate(nodalObject.variableCal[solvingVariableIdentifier],{'y':2}),
                                            ad.Differentiate(nodalObject.variableCal[solvingVariableIdentifier],{'z':2})]),
                               -1.,
                               nodalObject.material[0].properties[diffusivityIdentifier]])
        if internalSourceEquation!=None:
            internalSource=internalSourceEquation(nodalObject)
        else:
            internalSource=0.
        if externalSourceEquation!=None:
            externalSource=ad.Multiply([-1.,externalSourceEquation(nodalObject)])
        else:
            externalSource=0.
        self.result=ad.Addition([variation,diffusion,internalSource,externalSource])
    def cal(self,x,dOrder):
        return self.result.cal(x,dOrder)
class microStreamingConvection:
    def __init__(self,eqnName,nodalObject,convectionVariableIdentifier):
        convectionU=ad.Multiply([nodalObject.variableCal['U'],ad.Conjugate(ad.Differentiate(nodalObject.variableCal[convectionVariableIdentifier],{'x':1}))])
        convectionV=ad.Multiply([nodalObject.variableCal['V'],ad.Conjugate(ad.Differentiate(nodalObject.variableCal[convectionVariableIdentifier],{'y':1}))])
        convectionW=ad.Multiply([nodalObject.variableCal['W'],ad.Conjugate(ad.Differentiate(nodalObject.variableCal[convectionVariableIdentifier],{'z':1}))])
        convection=ad.Addition([convectionU,convectionV,convectionW])
    def cal(self,x,dOrder):
        return self.convection.cal(x,dOrder)    
        
        
        
        
        