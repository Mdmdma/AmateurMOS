import numpy as np
import matplotlib as plt

focal_length = 3.434
aperture = 0.508
numeric_aperture = 6.8


def degree_to_dist(degree, focal_length=focal_length):
    return np.asin(np.deg2rad(degree / 3600)) * focal_length 

def dist_to_degree(dist, focal_length = focal_length):
    return np.rad2deg(np.sin(dist/focal_length)) * 3600

dist_to_degree(0.000012) * 1000