#!/usr/bin/env python
# Copyright (c) 2017 Zi Wang
from push_world import *
import sys


if __name__ == '__main__':
    rx = float(sys.argv[1])
    ry = float(sys.argv[2])
    gx = float(sys.argv[4])
    gy = float(sys.argv[5])
    simu_steps = int(float(sys.argv[3]) * 10)
    
    # set it to False if no gui needed
    world = b2WorldInterface(False)
    oshape, osize, ofriction, odensity, bfriction, hand_shape, hand_size  = 'circle', 1, 0.01, 0.05, 0.01, 'rectangle', (0.3,1) 
    thing,base = make_thing(500, 500, world, oshape, osize, ofriction, odensity, bfriction, (0,0))

    init_angle = np.arctan(ry/rx)
    robot = end_effector(world, (rx,ry), base, init_angle, hand_shape, hand_size)
    ret = simu_push(world, thing, robot, base, simu_steps)
    ret = np.linalg.norm(np.array([gx, gy]) - ret)
    sys.stdout.write(str(ret))