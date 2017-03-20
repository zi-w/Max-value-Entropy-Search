#!/usr/bin/env python
# Copyright (c) 2017 Zi Wang
from push_world import *
import sys

# difference to generate_simudata is an input that control angle of push
if __name__ == '__main__':
    rx = float(sys.argv[1])
    ry = float(sys.argv[2])
    gx = float(sys.argv[4])
    gy = float(sys.argv[5])
    init_angle = float(sys.argv[6])
    simu_steps = int(float(sys.argv[3]) * 10)
    # Set the parameter to True if need gui
    world = b2WorldInterface(False)
    oshape, osize, ofriction, odensity, bfriction, hand_shape, hand_size  = 'circle', 1, 0.01, 0.05, 0.01, 'rectangle', (0.3,1) 
    thing,base = make_thing(500, 500, world, oshape, osize, ofriction, odensity, bfriction, (0,0))
    xvel = -rx;
    yvel = -ry;
    regu = np.linalg.norm([xvel,yvel])
    xvel = xvel / regu * 10;
    yvel = yvel / regu * 10;
    robot = end_effector(world, (rx,ry), base, init_angle, hand_shape, hand_size)
    ret = simu_push2(world, thing, robot, base, xvel, yvel, simu_steps)
    ret = np.linalg.norm(np.array([gx, gy]) - ret)
    sys.stdout.write(str(ret))