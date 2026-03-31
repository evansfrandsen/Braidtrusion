# -*- coding: utf-8 -*-
# title : SplineBraidModel_V4_functionFile.py
# contributors : Pierre-Richard Junior Doffonsou (pierre-richard-junior.doffonsou@polymtl.ca) 
# date : 21-august-2023

#-------------------------------------------------------------------#
#                        Function file                              #
#-------------------------------------------------------------------#

# modules 
import numpy as np
import os
from math import *
from prettytable import PrettyTable
# import trimesh 



def clear_terminal():
    """Function that clears the terminal"""
    system = os.name
    if system == 'posix':
        os.system('clear')
    elif system == 'nt':
        os.system('cls')

def ellipse(aby,bby,npts):
    """Function that creates and cloud of points to plot an ellipse"""
    t = np.linspace(0, 2*np.pi,npts)
    EllipseCloud = np.zeros((npts,3))
    EllipseCloud[:,0] = bby * np.cos(t)
    EllipseCloud[:,1] = aby * np.sin(t)
    return EllipseCloud   

def rotate(points,vdir,vini,axe):
    """Function that rotates points around an axis"""
    # Normalisation des vecteurs 
    vdir = vdir / np.linalg.norm(vdir) 
    vini = vini / np.linalg.norm(vini)
    theta = np.arccos(np.dot(vdir,vini))
    if axe.lower() == 'x': 
        # rotation autour de l'axe x
        Rx = np.array([
            [1, 0, 0], 
            [0, np.cos(theta), np.sin(theta)], 
            [0, -np.sin(theta), np.cos(theta)]
        ])
        pts = np.dot(points,Rx)
    elif axe.lower() == 'y': 
        # rotation autour de l'axe y
        Ry = np.array([
            [np.cos(theta), 0, -np.sin(theta)], 
            [0, 1, 0], 
            [np.sin(theta), 0, np.cos(theta)]
        ])
        pts = np.dot(points,Ry)
    elif axe.lower() == 'z': 
        # rotation autour de l'axe z
        theta = np.arctan2(vdir[1], vdir[0])
        Rz = np.array([
            [np.cos(theta), np.sin(theta), 0], 
            [-np.sin(theta), np.cos(theta), 0], 
            [0, 0, 1]
        ])
        pts = np.dot(points,Rz)    
    return pts

def rotation_arbitrary(angle, axis, A):
    # On s'assure que l'axe est unitaire 
    axis = axis / np.linalg.norm(axis)
    wx, wy, wz = axis[0], axis[1], axis[2]
    
    R = np.array([
        [cos(angle) + wx**2*(1-cos(angle)), wx*wy*(1-cos(angle))-wz*sin(angle), wy*sin(angle)+wx*wz*(1-cos(angle))],
        [wz*sin(angle)+wx*wy*(1-cos(angle)), cos(angle)+wy**2*(1-cos(angle)), -wx*sin(angle)+wy*wz*(1-cos(angle))],
        [-wy*sin(angle)+wx*wz*(1-cos(angle)), wx*sin(angle)+wy*wz*(1-cos(angle)), cos(angle)+wz**2*(1-cos(angle))]
    ])
    A_ = np.zeros(A.shape)
    for i in range(A_.shape[0]): 
        A_[i,:] = np.dot(R,A[i,:])

    return A_

def director_vector(x1,x2):
    return (x2-x1)/np.linalg.norm(x2-x1)

# def rotate_arbitrary(points,vdir):
#     """Allow the rotation of points around an arbitrary axis"""
#     vnorm = np.array([0.,0.,1.])
#     vdir = vdir / np.linalg.norm(vdir)
#     axis = np.cross(vnorm,vdir)
#     axis = axis / np.linalg.norm(axis)
#     theta = - np.arccos(np.dot(vnorm,vdir))

#     a_skew = np.array([
#         [0, -axis[2], axis[1]], 
#         [axis[2], 0, -axis[0]], 
#         [-axis[1], axis[0], 0]
#     ])
#     R = np.cos(theta) * np.eye(3) + (1 - np.cos(theta)) * axis.reshape(3,1) * axis + np.sin(theta) * a_skew
#     return np.dot(points, R)

def rotation_matrix(phi,theta,alpha):
    """Function that returns a 3D rotation matrix
       phi : angle of rotation about Z axis 
       theta : angle of rotation about the new x1' axis (braiding angle)
       alpha : angle of rotation about the y axis (crimp angle)
       The angles are in radians
    """
    Rx = np.array([
        [1,          0,           0],
        [0, cos(theta), -sin(theta)], 
        [0, sin(theta),  cos(theta)]
    ])
    Ry = np.array([
        [ cos(phi), 0, sin(phi)],
        [        0, 1,        0],
        [-sin(phi), 0, cos(phi)]
    ])
    Rz = np.array([
        [cos(alpha), -sin(alpha), 0],
        [sin(alpha),  cos(alpha), 0],
        [         0,           0, 1]
    ])
    return np.dot(Ry,Rx,Rz)

def translate(points,vtrans):
    """Function that allows the translation of points in a particular direction"""
    u = np.zeros(points.shape)
    for i in range(points.shape[1]):
        u[:,i] = points[:,i] + vtrans[i]
    return u

def transformation(vector, theta_x, theta_y, theta_z, trans):
    """
    vector : array containing the coordinates of the point or vector to transform
    theta_x : angle of rotation around x axis
    theta_y : angle of rotation around y axis 
    theta_z : angle of rotation around z axis 
    trans : array containing the translation coordinates 
    """

    # parameter : boolean value --> 0 for a vector, 1 for a point
    parameter = vector[-1,0]
    
    # Rotation matrix around z axis (Rz)
    Rz = np.array([
        [cos(theta_z), -sin(theta_z),             0,           0],
        [sin(theta_z),  cos(theta_z),             0,           0],
        [           0,             0,             1,           0], 
        [           0,             0,             0,   parameter]
    ])
    
    # Rotation matrix around x axis (Rx)
    Rx = np.array([
        [           1,             0,             0,    trans[0]],
        [           0,  cos(theta_x),  sin(theta_x),    trans[1]],
        [           0, -sin(theta_x),  cos(theta_x),    trans[2]], 
        [           0,             0,             0,   parameter]
    ])

    # Rotation matrix around y axis (Ry) combined with a translation array
    Ry = np.array([
        [ cos(theta_y),             0, -sin(theta_y),           0],
        [            0,             1,             0,           0],
        [ sin(theta_y),             0,  cos(theta_y),           0], 
        [            0,             0,             0,   parameter]
    ])

    # y = Mx where Y is the transformed point and M the transformation matrix 
    y = np.dot(Rz,vector)
    y = np.dot(Ry,y)
    y = np.dot(Rx,y)
    return y[0:3].reshape((1,3))[0]


def computeYarnsLength(braid):
    """Compute the length of yarns"""
    pass
