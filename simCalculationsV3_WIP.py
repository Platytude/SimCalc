import math
from math import *
import os
import sys
import numpy as np

#this is meant to be a refactoring of simcalc using standard libraries and numpy in order to clean things up and make them a bit more maintainable

#a note about the geometry used in these calculations:
#The co-ordinate system is assumed to be as follows:
#the Z axis is the vertical axis of the sim, perpendicular to the ground
#the y axis is the horizontal axis of the sim running from front to back
#the x axis is the horizontal axis of the sim running from right to left
#PHI is the angle measured from the z axis to a line
#THETA is the angle measured from the x axis to the projection of the line
#The origin of the coordinate system centers on the axis of rotation of the pivot joint of the sim

#CONSTANTS:
DOTZERO_PRECISION = 0.000001
POINTZERO_PRECISION = 0.001

#because somehow this isn't a default function....
def normalize(vect):
    vect = vect / np.linalg.norm(vect)
    return(vect)

def vAngle(v1, v2):
    angle = acos(np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)))
    return(angle)

#calculate linear speed based on CTC length.
#ctcLength in mm, mSpeed in rpm, mTorque in nm , currentAngle in degrees, unit is metric or imperial
#find out what default theta is for simtools/smc3
def getCTCParams(ctcLength, thetaMax=55, mSpeed = None, mTorque = None, currentAngle = 0, unit = None):
    currentAngle = currentAngle + 90
    
    #default for return travel variables
    travel = None
    vMax = None
    fMax = None
    fCurrent = None
    vCurrent = None

    if (ctcLength != None) and (thetaMax != None): 
        travel = ctcLength*cos(radians(90-thetaMax)) #returns maximum linear travel IN ONE DIRECTION FROM THE ZERO POINT
        print('Maximum linear Travel in is +/-',travel)

    if mSpeed != None:
        vMax = ctcLength * (mSpeed*0.10472)  #mSpeed in RPM so needs to be scaled to rad/s, calculates theoretical max velocity
        print('Maximum linear velocity is ',vMax)

    if mTorque != None:
        fMax = mTorque/(ctcLength)
        print('Maximum linear force is ',fMax);

    if currentAngle != None:
        fCurrent = abs( fMax * cos(radians(90 - currentAngle)))
        print('Current linear force is ',fCurrent,' at a ctc angle of ',currentAngle,' degrees')
        vCurrent = abs(vMax * cos(radians(90- currentAngle)))
        print('Current linear velocity is ',vCurrent,' at a ctc angle of ',currentAngle,' degrees')

    return(travel, vMax, fMax, fCurrent, vCurrent)


#returns a vector object containing the cartesian position of the ctc lever in 3d space, with the pivot point set as the origin of the system
#length is the length of the lever
#restAngle is the angle the CTC lever makes with the "x axis" if we're looking straight at it from the outside, measured counterclockwise.
#ctcAngle is the angle of rotation away from the rest position of the lever
#rotMt is the xyz coordinate vector of the rod attachment point to the frame
#mtrPos is the xyz coordinate vector of the point where the motor axis connects to the ctc lever. (the centre of rotation of the CTC lever)
#mtrAngle is the angle that the axis of the motor makes with the y-axis of our coordinate system

def getCTCGeometry3d(length, restAngle, ctcAngle, mtrPos, mtrAngle=0):
    #the position of the lever will follow a circular arc through space, centered about the axis of rotation of the lever
    #theta is the absolute angle of rotation around the axis, measured starting from negative y.
    restAngle = 180-restAngle #dear lord forgive me I have sinned but I can't remember why the f*** this has to be here
    thetaCTC = restAngle - ctcAngle
    print('Absolute ThetaCTC is ',thetaCTC)
    
    #use polar coordinate conversion to find the cartesian points
    posx =  sin(radians(mtrAngle)) * length * cos(radians(thetaCTC))
    print('posx',posx)
    posy =  cos(radians(mtrAngle)) * length * cos(radians(thetaCTC))
    print('posy',posy)
    posz =  length*sin(radians(thetaCTC))
    print('posz',posz)

    #place the output into an array, add to the original motor position:
    dist = np.array([posx, posy, posz])
    position = mtrPos + dist
    
    #now get the direction vector of the force being output by the lever:
    #find phi and theta in spherical coordinates
    phi = radians(90- thetaCTC)
    print('phi is ',degrees(phi))
    theta = radians(90 + mtrAngle)
    print('theta is ',degrees(theta))
    
    dfx = sin(phi+radians(90))*cos(theta) 
    dfy = sin(phi+radians(90))*sin(theta) 
    dfz = cos(phi+radians(90)) 

    f = np.array([dfx, dfy, dfz])
    fHat = normalize(f)
    #print('Force Direction Unit Vector is:')
    return(position,  fHat)


#returns optimal resting angle for CTC lever
#length is CTC length
def getCTCRestAngle(length, mtrPos, rodMt, mtrAngle=0):
    #see diagram of triangle and geometry to get an explanationn for this.
    
    mtrPos = np.array([0, mtrPos[1], mtrPos[2]])
    rodMt = np.array([0, rodMt[1], rodMt[2]])

    a = rodMt - mtrPos
    a = np.linalg.norm(a)

    b = abs(rodMt[2] - mtrPos[2])

    theta = asin(b/a)
    phi = acos(length/a)
    restAngle = 180 - degrees(theta + phi)

    if (mtrPos[1] > 0):
        restAngle = 180 - restAngle
        
    return(restAngle)

#returns the pitch and roll travel parameters of the simulator based on two CTC positions, the rod mount position, and rod length
def simTravel(ctcPosi, ctcPosf, rodMt, rodLen, pivLenUser=None):
    #shorter variable names so that the equations are at least legible
    pivPos = np.array([0,0,0])
    pivLen = np.linalg.norm(rodMt)
    
    #calculate total travel and angle by finding points of travel at min and max ctc travel.
    #sphereSphereIntersect(S1, r1, S2, r2, constraintx=None, constrainty=None):
    #return(xConstraintIntersect1, xConstraintIntersect2, yConstraintIntersect1, yConstraintIntersect2)
    minPosVects = sphereSphereIntersect(pivPos, pivLen, ctcPosi, rodLen, rodMt[0], rodMt[1])
    minPosPitch1 = minPosVects[0]
    minPosPitch2 = minPosVects[1]
    minPosRoll1 = minPosVects[2]
    minPosRoll2 = minPosVects[3]

    #determine which returned vector is the correct vector (this class of geometry problem will have two solutions, one supurious):
    minPitchDist1 = rodMt-minPosPitch1
    minPitchDist2 = rodMt-minPosPitch2
    if (np.linalg.norm(minPitchDist1) < np.linalg.norm(minPitchDist2)):
        minPitchPos = minPitchDist1
    else:
        minPitchPos = minPitchDist2

    minRollDist1 = rodMt-minPosRoll1
    minRollDist2 = rodMt-minPosRoll2
    if (np.linalg.norm(minRollDist1) < np.linalg.norm(minRollDist2)):
        minRollPos = minRollDist1
    else:
        minRollPos = minRollDist2
    
    #do the same for the lever at maximum position
    maxPosVects = sphereSphereIntersect(pivPos, pivLen, ctcPosf, rodLen, rodMt[0], rodMt[1])
    maxPosPitch1 = maxPosVects[0]
    maxPosPitch2 = maxPosVects[1]
    maxPosRoll1 = maxPosVects[2]
    maxPosRoll2 = maxPosVects[3]

    maxPitchDist1 = rodMt-maxPosPitch1
    maxPitchDist2 = rodMt-maxPosPitch2
    if (np.linalg.norm(maxPitchDist1) < np.linalg.norm(maxPitchDist2)):
        maxPitchPos = maxPitchDist1
    else:
        maxPitchPos = maxPitchDist2

    maxRollDist1 = rodMt-maxPosRoll1
    maxRollDist2 = rodMt-maxPosRoll2
    if (np.linalg.norm(maxRollDist1) < np.linalg.norm(maxRollDist2)):
        maxRollPos = maxRollDist1
    else:
        maxRollPos = maxRollDist2

    #we're now left with maxRollPos, minRollPos, maxPitchPos, maxRollPos, find angles and distances:
    #linear roll travel distance
    rollTrav = maxRollPos-minRollPos
    rollTrav = np.linalg.norm(rollTrav)
    
    #linear pitch travel
    pitchTrav = maxPitchPos-minPitchPos
    pitchTrav = np.linalg.norm(pitchTrav)
    
    #angular pitch & roll travel
    pitchAngle = degrees(2*atan((pitchTrav/2)/pivLen))
    rollAngle = degrees(2*atan((rollTrav/2)/pivLen))

    #package
    pitchParams = [pitchTrav, pitchAngle]
    rollParams = [rollTrav, rollAngle]
    
    return(pitchParams, rollParams)
    
#The following function will output the force values of pitch and roll for a 2DOF simulator with CTC levers.
#CTCpos is the set of XYZ coordinates of the end of the CTC lever
#rodMt is the set of XYZ coordinates of the rod mount point on the frame
#ctcForce is the magnitude of the applied force at the CTC lever.
#forceDir is the direction vector of the applied force of the CTC lever.
#motorAngle in degrees
def forceCalcs3(ctcPos, rodMt, ctcForce, forceDir, motorAngle, ctcAngle):
    
    motorAngle = radians(motorAngle)
    ctcAngle = radians(ctcAngle)
    
    #the easiest way to do this will be to do a change of basis of the system to make it in-line with the CTC forces.
    x = np.array([1,0,0])
    y = np.array([0,1,0])
    z = np.array([0,0,1])
    
    #the corresponding vectors for the applied forces will be the axial, tangent, and applied forces at the CTC joint
    #the applied force direction is easy (forceDir), to find the axial one:
    #since the motor is tilted inwards by motorAngle (in the XY plane):
    
    appliedForceDir = normalize(forceDir)
    print("applied force direction is",appliedForceDir)
    
    axialForceDir = np.array([cos(motorAngle), sin(motorAngle), 0])
    axialForceDir = normalize(axialForceDir)
    print("axial force direction is",axialForceDir)
    
    tangentForceDir = np.cross(axialForceDir, appliedForceDir)
    tangentForceDir = normalize(tangentForceDir)
    print("tangent force direction is",tangentForceDir)
    
    #now for the transformed coordinate system unit vectors:
    Xt = tangentForceDir
    Yt = axialForceDir
    Zt = appliedForceDir
    
    #now to transform the ctcPos and rodMt points from XYZ coodrinates to our new ones, we perform a change of basis:
    #the following is the basis matrix to transform from our "force basis" to R2, 
    basisMatrixInverse = np.array([[Xt[0] , Yt[0] , Zt[0] ],
                                   [Xt[1] , Yt[1] , Zt[1] ],
                                   [Xt[2] , Yt[2] , Zt[2] ]])
    basisMatrix = np.linalg.inv(basisMatrixInverse) #we invert to get the matrix to transform from R2 to the force matrix:
    print("basis Matrix is",basisMatrix)
    
    #to get the coordinates of the old points in the new basis, we simply multiply.
    ctcPosTransformed = np.dot(basisMatrix, ctcPos)
    print("transformed CTCpos is",ctcPosTransformed)
    rodMtTransformed = np.dot(basisMatrix, rodMt)
    print("transformed rodMt is",rodMtTransformed)
    
    #now for the vector going from the CTC to the rod mount:
    r = rodMtTransformed - ctcPosTransformed
    rHat = normalize(r)
    print("rHat (transformed) is",rHat)
    
    #Scale R by the magnitude of the applied force to get the total (transformed) force vector Ft:
    Ft = ctcForce*rHat
    print("Transformed force vector is",Ft)
    
    #Now un-transform the force vector:
    Force = np.dot(basisMatrixInverse,Ft)
    print("Force vector in regular coordinate system is",Force,"with magnitude",np.linalg.norm(Force))
    
    #now we do the projection stuff to find the pitch/roll forces:
    #pitch first:
    rodMtyz = np.copy(rodMt)
    rodMtyz[0]=0 #we want only the y and z components
    print("rodMtyz is",rodMtyz)
    Fyz = np.copy(Force) #numpy doesn't copy variable data on reassingment, only creates a new pointer
    Fyz[0]=0 #we want only the y and z components
    print("Fyz is",Fyz)
    thetaPitch = vAngle(rodMtyz, Fyz)
    print("pitch diff angle is",degrees(thetaPitch))
    pitchForce = np.linalg.norm(Fyz)*sin(thetaPitch)
    print("pitch force is",pitchForce)
    pitchTorque = pitchForce*np.linalg.norm(rodMtyz)
    print("pitch torque is",pitchTorque)
    
    #now for roll:
    print("rodMt",rodMt)
    rodMtxz = np.copy(rodMt)
    rodMtxz[1] = 0
    print("rodMtxz is",rodMtxz)
    Fxz = np.copy(Force)
    Fxz[1] = 0
    print("Fxz is",Fxz)
    thetaRoll = vAngle(rodMtxz, Fxz)
    print("Roll diff angle is",degrees(thetaRoll))
    rollForce = np.linalg.norm(Fxz)*sin(thetaRoll)
    print("roll force is",rollForce)
    rollTorque = rollForce*np.linalg.norm(rodMtxz)
    print("roll torque is",rollTorque)
    
    return(pitchForce, pitchTorque, rollForce, rollTorque)
    
    
#this function returns two points located on the circle of intersection of two spheres.
#S1, R1 and S2, R2 are the positions and radii of two spheres located in 3d space
#the function will then calculate the circle of intersection of the two spheres
#after this, the circle will examine constraintx and constrainty, and attempt to find
#the points on the circle that satisfy the constraints.
#ie, the point on the circle where x=1 and y=3
#see diagrams in documentation (TODO)
def sphereSphereIntersect(S1, r1, S2, r2, constraintx=None, constrainty=None):
    yconstraintIntersect = None
    xconstraintIntersect = None

    a = S2 - S1
    print("vect a is ", a)
    aHat = normalize(a)
    
    if (np.linalg.norm(a) > (r1 + r2)):
        print("Error, spheres do not intersect")

    thetaS = acos((np.linalg.norm(a)**2 + r1**2 - r2**2)/(2*np.linalg.norm(a)*r1));
    print('theta s', degrees(thetaS))
    s = r1*cos(thetaS)
    print('s',s)
    h = r1 * sin(thetaS)
    print('h',h)
    M = S1 + s*aHat
    print('point M is ', M)

    someVect = np.array([-aHat[0], aHat[1], 0])
    u = np.cross(someVect, aHat)
    uHat = normalize(u)
    v = np.cross(uHat, aHat)
    vHat = normalize(v)

    print('uHat:', uHat)
    print('vHat:', vHat)

    if ((np.dot(vHat, uHat)>DOTZERO_PRECISION) or (np.dot(vHat, aHat)>DOTZERO_PRECISION) or (np.dot(aHat, uHat)>DOTZERO_PRECISION)):
        print("error, basis is not orthonormal to a")

    #Oh boy, here we go
    kx = (constraintx - M[0]) / h
    print('kx',kx)
    ky = (constrainty - M[1]) / h
    print('ky',ky)
    divx = sqrt(uHat[0]**2 + vHat[0]**2)
    print('divx',divx)
    divy = sqrt(uHat[1]**2 + vHat[1]**2)
    print('divy',divy)

    Ux = uHat[0]/divx
    print('Ux', Ux)
    Vx = vHat[0]/divx
    print('Vx', Vx)

    Uy = uHat[1]/divy
    print('Uy', Uy)
    Vy = vHat[1]/divy
    print('Vy', Vy)

    Kx = kx/divx
    print('Kx', kx)
    Ky = ky/divy
    print('Ky', ky)

    betaX = atan2(Ux, Vx)
    print('Beta X', betaX)
    betaY = atan2(Uy, Vy)
    print('Beta Y', betaY)
        
    tx = asin(Kx) - betaX
    print('tx',tx)
    ty = asin(Ky) - betaY
    print('tx',ty)
    
    tx2 = pi + asin(Kx) + betaX;
    print('tx2',tx2)
    ty2 = pi + asin(Ky) + betaY;
    print('ty2',ty2)
 
    xConstraintIntersect1 = M + h*cos(tx)*uHat + h*sin(tx)*vHat
    yConstraintIntersect1 = M + h*cos(ty)*uHat + h*sin(ty)*vHat
    
    xConstraintIntersect2 = M + h*cos(tx2)*uHat - h*sin(tx2)*vHat
    yConstraintIntersect2 = M + h*cos(ty2)*uHat - h*sin(ty2)*vHat

    print('\n\n')
    print('x constraint is: ',constraintx)
    print('x constraint intersect 1:',xConstraintIntersect1)
    print('x constraint intersect 2:', xConstraintIntersect2)

    print('y constraint is: ',constrainty)
    print('y constraint intersect 1:',yConstraintIntersect1)
    print('y constraint intersect 2:', yConstraintIntersect2)
    #I'm so sorry
    
    return(xConstraintIntersect1, xConstraintIntersect2, yConstraintIntersect1, yConstraintIntersect2)



 
 
if (__name__ == "__main__"):
    #testing stuff
    length = 1
    ctcLength = length
    restAngle = 45
    ctcAngle = 0
    origin = np.array([0,0,0])
    mtrPos = np.array([22.4,-45.26,-7.9])
    mtrAngle = 15
    rodMt = np.array([5.37, -21.99, 27.68])
    ctcPosi = np.array([23.15, -48.1, -7.7])
    ctcPosf = np.array([22.62, -46.18, -5.32])
    rodLen = 47.44
    mSpeed = 120
    mTorque = 1
    thetaMax = 55
    currentAngle = 0
    
    ctcParams = getCTCParams(ctcLength, thetaMax, mSpeed, mTorque, currentAngle = 0, unit = None)
    ctcForce = ctcParams[2]
    
    ctcGeometry3d = getCTCGeometry3d(length, restAngle, ctcAngle, mtrPos, mtrAngle)
    ctcPos = ctcGeometry3d[0]
    forceDir = ctcGeometry3d[1]
    
    #x = getCTCGeometry3d(length, restAngle, ctcAngle, mtrPos, mtrAngle)
    #print(x)
    
    #y = getCTCRestAngle(length, mtrPos, rodMt, mtrAngle)
    #print(y)
    
    #tst = sphereSphereIntersect(origin, np.linalg.norm(rodMt), ctcPosi, rodLen, rodMt[0], rodMt[1])
    #print(tst)
    
    #tst = simTravel(ctcPosi, ctcPosf, rodMt, rodLen)
    #print(tst)
    


    
    tst = forceCalcs3(ctcPos, rodMt, ctcForce, forceDir, mtrAngle, restAngle)
    
