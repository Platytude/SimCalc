from math import *
import math
import os
import sys
import numpy

#a note about the geometry used in these calculations:
#The co-ordinate system is assumed to be as follows:
#the Z axis is the vertical axis of the sim, perpendicular to the ground
#the y axis is the horizontal axis of the sim running from front to back
#the x axis is the horizontal axis of the sim running from side to side
#PHI is the angle measured from the z axis to a line
#THETA is the angle measured from the x axis to the projection of the line
#The origin of the coordinate system centers on the axis of rotation of the pivot joint of the sim

#CONSTANTS:
##################################################################
#Since we have floating point math, checking dot products will always have error
#sets the threshold for what is considered a "zero" dot product for two orthogonal vectors
DOTZERO_PRECISION = 0.000001
POINTZERO_PRECISION = 0.001 #Sets the "zero" threshold for if two points are the same when comparing after fpoint caclculations

#Classes:
###################################################################

class Point():
    def __init__(self, x, y, z, master=None, *args, **kwargs):
        #print("x is ",x," y is ",y," z is ",z)
        self.x=x
        self.y=y
        selx.z=z
        
        
#DO NOT CHANGE THE VARIABLES UNLESS YOU EXPECT THE LENGTH NOT TO BE USED        
class Vector():
    def __init__(self, x, y, z, master=None, *args, **kwargs):
        #print("x is ",x," y is ",y," z is ",z)
        self.x=x
        self.y=y
        self.z=z
        #self.length = math.sqrt(self.x*self.x + self.y*self.y + self.z*self.z) #uncomment this is shit breaks, and comment out the @property method
    @property
    def length(self):
        return math.sqrt(self.x*self.x + self.y*self.y + self.z*self.z)

    #returns a vector object that points in the same direction as the parent, but with length 1.
    def unitVect(self):
        xu = self.x/self.length
        yu = self.y/self.length
        zu = self.z/self.length
        return(Vector(xu, yu, zu))

    def print(self):
        print("(x=",self.x,", y=",self.y,", z=",self.z,", length=",self.length,')')

class Vector2():
    def __init__(self, x, y, z, master=None, *args, **kwargs):
        #print("x is ",x," y is ",y," z is ",z)
        self.x=x
        self.y=y
        self.z=z

    @property
    def length(self):
        return math.sqrt(self.x*self.x + self.y*self.y + self.z*self.z)

#CLASS RELATED FUNCTIONS

#returns the dot product of the two vectors. Technically not a class, but should be near here
def dot(v1, v2):
    dotProduct = (v1.x*v2.x) + (v1.y*v2.y) + (v1.z*v2.z)
    return(dotProduct)

def cross(v1, v2):
    crossProduct = Vector(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x)
    return(crossProduct)

def vSubtract(v1, v2):
    v3 = Vector(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z)
    return(v3)

def vAdd(v1, v2):
    v3 = Vector(v1.x+v2.x, v1.y+v2.y, v1.z+v2.z)
    return(v3)

def vScale(scalar, vector):
    v2 = Vector(scalar*vector.x, scalar*vector.y, scalar*vector.z)
    return(v2)


def vAngle(v1, v2):
    angle = acos(dot(v1,v2)/(v1.length*v2.length))
    return(angle)
        

#FUNCTIONS:
###################################################################

#calculate linear speed based on CTC length.
#ctcLength in mm, mSpeed in rpm, mTorque in nm , currentAngle in degrees, unit is metric or imperial
#find out what default theta is for simtools/smc3
def getCTCParams(ctcLength, thetaMax=55, mSpeed = None, mTorque = None, currentAngle = 0, unit = None):
    if unit == 'metric':
        pass
    if unit == 'imperial':
        pass

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
        vMax = ctcLength * mSpeed*0.10472  #mSpeed in RPM so needs to be scaled to rad/s, calculates theoretical max velocity
        print('Maximum linear velocity is ',vMax)

    if mTorque != None:
        fMax = mTorque/(ctcLength)
        print('Maximum linear force is ',fMax);

    if currentAngle != None:
        fCurrent = abs( fMax * cos(radians(90 - currentAngle)))
        print('Current linear force is ',fCurrent,' at a ctc angle of ',currentAngle,' degrees')
        vCurrent = abs(vMax * cos(radians(90- currentAngle)))
        print('Current linear velocity is ',vCurrent,' at a ctc angle of ',currentAngle,' degrees');

    return(travel, vMax, fMax, fCurrent, vCurrent)

#returns a vector object containing the cartesian position of the ctc lever in 3d space, with the pivot point set as the origin of the system
#length is the length of the lever
#restAngle is the angle the CTC lever makes with the "x axis" if we're looking straight at it, measured counterclockwise.
#ctcAngle is the angle of rotation away from the rest position of the lever
#rotMt is the xyz coordinate vector of the rod attachment point to the frame
#mtrPos is the xyz coordinate vector of the point where the motor axis connects to the ctc lever. (the centre of rotation of the CTC lever)
#mtrAngle is the angle that the axis of the motor makes with the y-axis of our coordinate system
def getCTCGeometry3d(length, restAngle, ctcAngle, mtrPos, mtrAngle=0):
    #the position of the lever will follow a circular arc through space, centered about the axis of rotation of the lever
    #theta is the absolute angle of rotation around the axis, measured starting from negative y.
    restAngle = 180-restAngle
    thetaCTC = restAngle - ctcAngle
    print('Absolute ThetaCTC is ',thetaCTC)
    
    posx =  sin(radians(mtrAngle)) * length * cos(radians(thetaCTC))
    print('posx',posx)
    posy =  cos(radians(mtrAngle)) * length * cos(radians(thetaCTC))
    print('posy',posy)
    posz =  length*sin(radians(thetaCTC))
    print('posz',posz)

    dist = Vector(posx, posy, posz)
    #print('\nposition vector for angle',ctcAngle,'is:')
    #dist.print()

    ctcposx = posx + mtrPos.x
    ctcposy = posy + mtrPos.y
    ctcposz = posz + mtrPos.z

    position = Vector(ctcposx, ctcposy, ctcposz)


    #now get the direction vector of the force being output by the lever:
    #find phi and theta in spherical coordinates
    phi = radians(90- thetaCTC)
    print('phi is ',degrees(phi))
    theta = radians(90 + mtrAngle)
    print('theta is ',degrees(theta))
    
    dfx = sin(phi+radians(90))*cos(theta) 
    dfy = sin(phi+radians(90))*sin(theta) 
    dfz = cos(phi+radians(90)) 

    f = Vector(dfx, dfy, dfz)
    fHat = f.unitVect()
    print('Force Direction Unit Vector is:')
    fHat.print()
    
    return(position,  fHat)

def getCTCRestAngle(length, mtrPos, rodMt, mtrAngle=0):
    #see diagram of triangle and geometry to get an explanationn for this.
    
    mtrPos = Vector(0, mtrPos.y, mtrPos.z)
    rodMt = Vector(0, rodMt.y, rodMt.z)

    a = vSubtract(rodMt, mtrPos)
    a = a.length

    b = abs(rodMt.z - mtrPos.z)

    theta = asin(b/a)
    phi = acos(length/a)
    restAngle = 180 - degrees(theta + phi)

    
    
    if (mtrPos.y > 0):
        restAngle = 180 - restAngle
        
    return restAngle


#returns  angular and linear travel of a sim based on the parameters below:
#ctcPosi, ctcPosf, rodMt, are vector objects (see classes above)
#ctcPosi is the initial position of the CTC rose joint
#ctcPosf is the final position of the CTC rose joint
#rodMt is the point where the rod mount connects to the frame
#rodLen is the rod length
#pivLen is the pivot lever length, for verification purposes only.
#the platform pivot point always sits at the origin
#x runs horizontally, parallel to the side-to-side axis of the sim
#y runs horizontally, parallel to the front-to-back axis of the sim
#z runs vertically
#will return maximum pitch and roll travels, as well as the corresponding angles of pitch and roll

#this turned out to be far harder than I initially thought
def simTravel(ctcPosi, ctcPosf, rodMt, rodLen, pivLenUser=None):
    #shorter variable names so that the equations are at least legible
    pivPos = Vector(0,0,0)
    pivLen = rodMt.length
    
    #calculate total travel and angle by finding points of travel at min and max ctc travel.
    #sphereSphereIntersect(S1, r1, S2, r2, constraintx=None, constrainty=None):
    #return(xConstraintIntersect1, xConstraintIntersect2, yConstraintIntersect1, yConstraintIntersect2)
    
    minPosVects = sphereSphereIntersect(pivPos, pivLen, ctcPosi, rodLen, rodMt.x, rodMt.y)
    minPosPitch1 = minPosVects[0]
    minPosPitch2 = minPosVects[1]
    minPosRoll1 = minPosVects[2]
    minPosRoll2 = minPosVects[3]

    #determine which vector is the correct vector:
    minPitchDist1 = vSubtract(rodMt, minPosPitch1)
    minPitchDist2 = vSubtract(rodMt, minPosPitch2)
    if (minPitchDist1.length < minPitchDist2.length):
        minPitchPos = minPitchDist1
    else:
        minPitchPos = minPitchDist2

    minRollDist1 = vSubtract(rodMt, minPosRoll1)
    minRollDist2 = vSubtract(rodMt, minPosRoll2)
    if (minRollDist1.length < minRollDist2.length):
        minRollPos = minRollDist1
    else:
        minRollPos = minRollDist2
    
    maxPosVects = sphereSphereIntersect(pivPos, pivLen, ctcPosf, rodLen, rodMt.x, rodMt.y)
    maxPosPitch1 = maxPosVects[0]
    maxPosPitch2 = maxPosVects[1]
    maxPosRoll1 = maxPosVects[2]
    maxPosRoll2 = maxPosVects[3]

    maxPitchDist1 = vSubtract(rodMt, maxPosPitch1)
    maxPitchDist2 = vSubtract(rodMt, maxPosPitch2)
    if (maxPitchDist1.length < maxPitchDist2.length):
        maxPitchPos = maxPitchDist1
    else:
        maxPitchPos = maxPitchDist2

    maxRollDist1 = vSubtract(rodMt, maxPosRoll1)
    maxRollDist2 = vSubtract(rodMt, maxPosRoll2)
    if (maxRollDist1.length < maxRollDist2.length):
        maxRollPos = maxRollDist1
    else:
        maxRollPos = maxRollDist2

    #we're now left with maxRollPos, minRollPos, maxPitchPos, maxRollPos, find angles and distances:
    rollTrav = vSubtract(maxRollPos, minRollPos)
    rollTrav = rollTrav.length

    pitchTrav = vSubtract(maxPitchPos, minPitchPos)
    pitchTrav = pitchTrav.length

    pitchAngle = degrees(2*atan((pitchTrav/2)/pivLen))
    rollAngle = degrees(2*atan((rollTrav/2)/pivLen))

    pitchParams = [pitchTrav, pitchAngle]
    rollParams = [rollTrav, rollAngle]
    
    return(pitchParams, rollParams)

    
                


#returns the 3d coordinates of the intersections between a sphere and a 2D circle
#cCirc and cSphere are Vector objects for the center of the sphere and the circle respectively
#rCirc and rSphere are the radii of the circle and the sphere
#circNormal is a vector object describing the normal to the plane that contains the circle
def sphereCircIntersects(cCirc, rCirc, circNormal, cSphere, rSphere):
    #first find an orthonormal basis for the plane containing the circle , comprised of vectors u,v.
    #generate first vector arbitrarily.

    #generate first vector by crossing arbitrary starting point with the normal
    #So long as the normal has no z component, uHat will always be the "y" component of the circle, if we look at it from the normal, it will be the x component
    vect = Vector(0,0,1)
    u = cross(vect,circNormal)
    uHat = u.unitVect()
    print('uHat')
    uHat.print()
                  
                  

    #generate second by taking the cross product with the normal, vHat will always be the "z" component of the circle, if we look at it from the normal of the circle it will be the y component
    v = cross(circNormal, uHat)
    vHat = v.unitVect()
    print('\nvHat')
    vHat.print()

    #w should be the same as the normal to the plane.
    wHat = circNormal.unitVect()
    print('\nwHat')
    wHat.print()

    print('\nPosition vector of circle centre')
    cCirc.print()

    print('\nPosition vector of sphere centre')
    cSphere.print()

    #vector with tip at the the circle centre and tail at the sphere centre
    spheretoPlane = vSubtract(cCirc, cSphere)
    print('\nvector from sphere centre to a point on the plane (circle centre) is')
    spheretoPlane.print()

    #find the distance from the centre of the sphere to the centre of the circle-sphere by dotting the normal with this vector
    dist = dot(spheretoPlane, wHat)/wHat.length
    print('\ndistance from sphere centre to circlesphere centre is ',dist)

    #distance vector from the sphere centre to the circlesphere centre
    distVect = vScale(dist, wHat)
    print('\nvector from sphere centre to circleSphere centre is ')
    distVect.print()

    #position vector of the circle-sphere centre
    cCircSphere = vAdd(distVect, cSphere)
    print('\nPosition vector of circle-sphere centre')
    cCircSphere.print()

    #now we need to parametrize the circle and circleSphere in terms of u and v, instead of x and y.
    #for any vector V in the plane, we can represent it as v = auHat + bvHat
    #circle first:

    #come back to this. This is proving to be a royal pain in the ass without matrix operators


#DOUBLE CHECK THIS SHIT
#returns force and force efficiency values for the sim in the pitch and roll directions
#ctcPos is the position of the ctc lever
#rodMt is the position of the rose joint where it attaches to the frame
#force is the amount of perpenicular force that the CTC lever is applying to the rod.
#mtrAngle is the angle that the motor makes with the Y-axis
#ctcAngle is the absolute angle of the CTC lever as measured from the "X-axis" if you look directly at the lever.
#add torque at u joint
def forceCalcs(ctcPos, rodMt, force, forceDir):

    ##THIS IS WRONG I THINK
    #subtract start from end to get direction vector
    #vector going along the linkage rod of the system
    rodVect = vSubtract(rodMt, ctcPos)
    #rodVect.print()

    pivVectYZ = Vector(0, rodMt.y, rodMt.z)#pitch
    pivVectXZ = Vector(rodMt.x, 0, rodMt.z)#roll

    print('forceDir value is ',forceDir.length)
    forceDir.print()
    
    forceVect = vScale(force, forceDir)
    forceVectYZ = Vector(0,forceVect.y, forceVect.z)
    forceVectXZ = Vector(forceVect.x, 0, forceVect.z)

    #the pitch direction vector will be perpendicular to the YZ projection of the pivot vector    
    pitchVect = Vector(0, pivVectYZ.z, -pivVectYZ.y)
    print('Pitch force direction vector')
    pitchVect.print()
    
    if (dot(pitchVect,pivVectYZ) > DOTZERO_PRECISION):
        print('error, pitch vector is not perpendicular to pivot vector')
    else:
        pass
        #print('Pitch vector is perpendicular to projection of pivot')
        

    #same for the roll vector
    rollVect = Vector(-pivVectXZ.z, 0, pivVectXZ.x)
    print('Roll force direction vector')
    rollVect.print()
    print('CHECK THAT THIS IS RIGHT BEFORE V1.0')
    if (dot(rollVect,pivVectXZ) > DOTZERO_PRECISION):
        print('error, roll vector is not perpendicular to pivot vector')
    else:
        pass
        #print('Roll vector is perpendicular to projection of pivot')

    #this is the angle difference between the applied force direction and the ideal direction of pitch force application (90 degrees from the pivot vector)
    pitchAngle = vAngle(pitchVect, forceVectYZ)
    print('pitchforce diff angle is', degrees(pitchAngle))

    pitchEff = abs( cos(pitchAngle))
    pitchForce = forceVectYZ.length * pitchEff
    pitchTorque = pitchForce * pivVectYZ.length

    #print('Pitch force is ',pitchForce,'with a directional efficiency of',pitchEff)

    rollAngle = vAngle(rollVect, forceVectXZ)
    print('\nrollforce diff angle is', degrees(rollAngle))

    rollEff = abs( cos(rollAngle))
    rollForce = forceVectXZ.length*rollEff
    rollTorque = rollForce * pivVectXZ.length
    
    #print('Roll force is ',rollForce,'with a directional efficiency of',rollEff)

    return(pitchEff, pitchForce, rollEff, rollForce, pitchTorque, rollTorque)


#takes in two vectors being the centres of each sphere, and two scalars being the lengths r1 and r2
#calculates the intersections between each sphere, given a constraint in either x or y
#returns vectors based on these constraints.
def sphereSphereIntersect(S1, r1, S2, r2, constraintx=None, constrainty=None):
    yconstraintIntersect = None
    xconstraintIntersect = None

    a = vSubtract(S2, S1)
    print("vect a is ")
    a.print()
    aHat = a.unitVect()
    
    if (a.length > (r1 + r2)):
        print("Error, spheres do not intersect")

    thetaS = acos((a.length**2 + r1**2 - r2**2)/(2*a.length*r1)); #this uses an alternate method, instead of finding a, we find the included angle.
    print('theta s', degrees(thetaS))
    s = r1*cos(thetaS)
    print('s',s)
    h = r1 * sin(thetaS)
    print('h',h)
    M = vAdd(S1, vScale(s, aHat))
    print('point M is ')
    M.print()

    someVect = Vector(-aHat.x, aHat.y, 0)
    u = cross(someVect, aHat)
    uHat = u.unitVect()
    v = cross(uHat, aHat)
    vHat = v.unitVect()

    print('uHat:')
    uHat.print()
    print('vHat:')
    vHat.print()

    if ((dot(vHat, uHat)>DOTZERO_PRECISION) or (dot(vHat, aHat)>DOTZERO_PRECISION) or (dot(aHat, uHat)>DOTZERO_PRECISION)):
        print("error, basis is not orthonormal to a")

    kx = (constraintx - M.x) / h
    print('kx',kx)
    ky = (constrainty - M.y) / h
    print('ky',ky)
    divx = sqrt(uHat.x**2 + vHat.x**2)
    print('divx',divx)
    divy = sqrt(uHat.y**2 + vHat.y**2)
    print('divy',divy)

    Ux = uHat.x/divx
    print('Ux', Ux)
    Vx = vHat.x/divx
    print('Vx', Vx)

    Uy = uHat.y/divy
    print('Uy', Uy)
    Vy = vHat.y/divy
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

    C_xconstraintX = M.x + h*cos(tx)*uHat.x + h*sin(tx)*vHat.x;
    C_xconstraintY = M.y + h*cos(tx)*uHat.y + h*sin(tx)*vHat.y;
    C_xconstraintz = M.z + h*cos(tx)*uHat.z + h*sin(tx)*vHat.z;

    xConstraintIntersect1 = Vector(C_xconstraintX, C_xconstraintY, C_xconstraintz)
    

    C_yconstraintX = M.x + h*cos(ty)*uHat.x + h*sin(ty)*vHat.x;
    C_yconstraintY = M.y + h*cos(ty)*uHat.y + h*sin(ty)*vHat.y;
    C_yconstraintz = M.z + h*cos(ty)*uHat.z + h*sin(ty)*vHat.z;

    

    yConstraintIntersect1 = Vector(C_yconstraintX, C_yconstraintY, C_yconstraintz)

    C2_xconstraintX = M.x + h*cos(tx2)*uHat.x - h*sin(tx2)*vHat.x;
    C2_xconstraintY = M.y + h*cos(tx2)*uHat.y - h*sin(tx2)*vHat.y;
    C2_xconstraintz = M.z + h*cos(tx2)*uHat.z - h*sin(tx2)*vHat.z;

    xConstraintIntersect2 = Vector(C2_xconstraintX, C2_xconstraintY, C2_xconstraintz)

    C2_yconstraintX = M.x + h*cos(ty2)*uHat.x - h*sin(ty2)*vHat.x;
    C2_yconstraintY = M.y + h*cos(ty2)*uHat.y - h*sin(ty2)*vHat.y;
    C2_yconstraintz = M.z + h*cos(ty2)*uHat.z - h*sin(ty2)*vHat.z;

    yConstraintIntersect2 = Vector(C2_yconstraintX, C2_yconstraintY, C2_yconstraintz)

    print('\n\n')
    
    
    print('x constraint is: ',constraintx)
    print('x constraint intersect 1:')
    xConstraintIntersect1.print()
    print('x constraint intersect 2:')
    xConstraintIntersect2.print()

    print('y constraint is: ',constrainty)
    print('y constraint intersect 1:')
    yConstraintIntersect1.print()
    print('y constraint intersect 2:')
    yConstraintIntersect2.print()

    return(xConstraintIntersect1, xConstraintIntersect2, yConstraintIntersect1, yConstraintIntersect2)


        
        
if (__name__ == "__main__"):
    #TENTATIVE VALUES FOR MY SIM, with a "short" lever (mostly rough values, take with a sizeable grain of salt)
    rodMt = Vector(5.37, -21.99, 27.68)
    ctcPosi = Vector(23.15, -48.1, -7.7)
    ctcPosf = Vector(22.62, -46.18, -5.32)
    rodLen = 47.44

    #simTravel(ctcPosi, ctcPosf, rodMt, rodLen, 'roll')



    origin = Vector(0,0,0)
    cCirc = Vector(23.15, -48.1, -7.7)
    cSphere = Vector(5.37, -21.99, 27.68)
    circNormalVect = Vector(1,0.2,0)
    circNormal = circNormalVect.unitVect()
    rCirc = 2
    rSphere = 47
    #sphereCircIntersects(cCirc, rCirc, circNormal, cSphere, rSphere)

    length = 2
    restAngle = 45
    ctcAngle = 0
    mtrPos = Vector(22.4, -45.26, -7.9)
    mtrAngle = 25

    #tst = sphereSphereIntersect(origin, rodMt.length, ctcPosi, rodLen, rodMt.x, rodMt.y)

    
    
    
    #print(x)
    
    #p = getCTCGeometry3d(length, restAngle, 60, mtrPos, mtrAngle=0)
    #ctcPosi = p[0]
    #print('ctcPosi is ')
    #ctcPosi.print()
    
    
    #v = getCTCGeometry3d(length, restAngle, -60, mtrPos, mtrAngle=0)
    #ctcPosf = v[0]
    #print('ctcPosf is ')
    #ctcPosf.print()
    
    #x = simTravel(ctcPosi, ctcPosf, rodMt, rodLen)
    #print(x)
        
    #print('\nchange in position is')
    #c = vSubtract(v,p)
    #c.print()

    #TESTING FOR FORCE AND EFFICIENCY VALUES
    #initial rod position is "rodMt"
    ctcParams = getCTCGeometry3d(length, restAngle, ctcAngle, mtrPos, mtrAngle)
    ctcPos = ctcParams[0]
    ctcPos.print()
    forceDir = ctcParams[1]
    forceDir.print()

    tst = forceCalcs(ctcPos, rodMt, 1, forceDir)
    print(tst)


    #def simTravel(ctcPosi, ctcPosF, rodMt, rodLen, constraint=None, pivLen=None):
    #def getCTCParams(ctcLength, thetaMax=65, mSpeed = None, mTorque = None, currentAngle = 90, unit = None):
    #def getCTCGeometry3d(length, restAngle, ctcAngle, mtrPos, mtrAngle=0):
    #def forceCalcs(ctcPos, rodMt, force, forceDir):



















    
    
    
