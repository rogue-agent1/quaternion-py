#!/usr/bin/env python3
"""Quaternion math for 3D rotations."""
import sys,math

class Quat:
    def __init__(self,w=1,x=0,y=0,z=0):self.w,self.x,self.y,self.z=w,x,y,z
    def __mul__(self,o):
        return Quat(self.w*o.w-self.x*o.x-self.y*o.y-self.z*o.z,
                    self.w*o.x+self.x*o.w+self.y*o.z-self.z*o.y,
                    self.w*o.y-self.x*o.z+self.y*o.w+self.z*o.x,
                    self.w*o.z+self.x*o.y-self.y*o.x+self.z*o.w)
    def conj(self):return Quat(self.w,-self.x,-self.y,-self.z)
    def norm(self):return math.sqrt(self.w**2+self.x**2+self.y**2+self.z**2)
    def normalize(self):n=self.norm();return Quat(self.w/n,self.x/n,self.y/n,self.z/n)
    def rotate(self,v):
        p=Quat(0,*v);r=self*p*self.conj()
        return(r.x,r.y,r.z)
    @staticmethod
    def from_axis_angle(axis,angle):
        s=math.sin(angle/2);c=math.cos(angle/2)
        n=math.sqrt(sum(a**2 for a in axis))
        return Quat(c,axis[0]/n*s,axis[1]/n*s,axis[2]/n*s)
    def to_euler(self):
        sy=2*(self.w*self.y-self.z*self.x)
        sy=max(-1,min(1,sy))
        return(math.atan2(2*(self.w*self.x+self.y*self.z),1-2*(self.x**2+self.y**2)),
               math.asin(sy),
               math.atan2(2*(self.w*self.z+self.x*self.y),1-2*(self.y**2+self.z**2)))
    def __repr__(self):return f"Quat({self.w:.4f},{self.x:.4f},{self.y:.4f},{self.z:.4f})"

def slerp(q1,q2,t):
    dot=q1.w*q2.w+q1.x*q2.x+q1.y*q2.y+q1.z*q2.z
    if dot<0:q2=Quat(-q2.w,-q2.x,-q2.y,-q2.z);dot=-dot
    if dot>0.9995:
        return Quat(q1.w+t*(q2.w-q1.w),q1.x+t*(q2.x-q1.x),q1.y+t*(q2.y-q1.y),q1.z+t*(q2.z-q1.z)).normalize()
    theta=math.acos(dot);s=math.sin(theta)
    w1=math.sin((1-t)*theta)/s;w2=math.sin(t*theta)/s
    return Quat(w1*q1.w+w2*q2.w,w1*q1.x+w2*q2.x,w1*q1.y+w2*q2.y,w1*q1.z+w2*q2.z)

def main():
    if len(sys.argv)>1 and sys.argv[1]=="--test":
        # 90° rotation around Z axis
        q=Quat.from_axis_angle((0,0,1),math.pi/2)
        v=q.rotate((1,0,0))
        assert abs(v[0])<1e-10 and abs(v[1]-1)<1e-10 and abs(v[2])<1e-10
        # Identity
        qi=Quat(1,0,0,0);v2=qi.rotate((1,2,3))
        assert all(abs(a-b)<1e-10 for a,b in zip(v2,(1,2,3)))
        # Norm
        assert abs(q.norm()-1)<1e-10
        # Composition: two 90° = 180°
        q2=q*q;v3=q2.rotate((1,0,0))
        assert abs(v3[0]+1)<1e-10 and abs(v3[1])<1e-10
        # Slerp
        q_mid=slerp(Quat(1,0,0,0),q,0.5)
        assert abs(q_mid.norm()-1)<1e-10
        print("All tests passed!")
    else:
        q=Quat.from_axis_angle((0,1,0),math.pi/4)
        print(f"45° around Y: {q}")
        print(f"Rotate (1,0,0): {q.rotate((1,0,0))}")
if __name__=="__main__":main()
