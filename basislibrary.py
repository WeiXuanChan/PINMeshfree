import numpy as np
import autoD as ad

def poly4dBasis(x,dOrder,nodalObject):
    const=ad.Constant(1.)
    scalar_x=ad.Scalar('x')
    scalar_y=ad.Scalar('y')
    scalar_z=ad.Scalar('z')
    scalar_t=ad.Scalar('t')
    xx=ad.Power(scalar_x,2.)
    yy=ad.Power(scalar_y,2.)
    zz=ad.Power(scalar_z,2.)
    xy=ad.Multiply([scalar_x,scalar_y])
    xz=ad.Multiply([scalar_x,scalar_z])
    yz=ad.Multiply([scalar_y,scalar_z])
    result=np.zeros(11)
    result[0]=const.cal(x,dOrder)
    result[1]=scalar_x.cal(x,dOrder)
    result[2]=scalar_y.cal(x,dOrder)
    result[3]=scalar_z.cal(x,dOrder)
    result[4]=scalar_t.cal(x,dOrder)
    result[5]=xx.cal(x,dOrder)
    result[6]=yy.cal(x,dOrder)
    result[7]=zz.cal(x,dOrder)
    result[8]=xy.cal(x,dOrder)
    result[9]=xz.cal(x,dOrder)
    result[10]=yz.cal(x,dOrder) 
    return result

def poly3dBasis(x,dOrder,nodalObject,omega):
    new_dOrder=dOrder.copy()
    coef=1.
    if 't' in new_dOrder:
        while new_dOrder['t']>0:
            coef=coef*1j*omega
            new_dOrder['t']-=1
    const=ad.Constant(1.)
    scalar_x=ad.Scalar('x')
    scalar_y=ad.Scalar('y')
    scalar_z=ad.Scalar('z')
    xx=ad.Power(scalar_x,2.)
    yy=ad.Power(scalar_y,2.)
    zz=ad.Power(scalar_z,2.)
    xy=ad.Multiply([scalar_x,scalar_y])
    xz=ad.Multiply([scalar_x,scalar_z])
    yz=ad.Multiply([scalar_y,scalar_z])
    result=np.zeros(10)
    result[0]=const.cal(x,new_dOrder)
    result[1]=scalar_x.cal(x,new_dOrder)
    result[2]=scalar_y.cal(x,new_dOrder)
    result[3]=scalar_z.cal(x,new_dOrder)
    result[4]=xx.cal(x,new_dOrder)
    result[5]=yy.cal(x,new_dOrder)
    result[6]=zz.cal(x,new_dOrder)
    result[7]=xy.cal(x,new_dOrder)
    result[8]=xz.cal(x,new_dOrder)
    result[9]=yz.cal(x,new_dOrder) 
    return result*coef


