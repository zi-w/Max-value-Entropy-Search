#!/usr/bin/env python
# Author: Zi Wang
from push_world import *
import sys

# difference to generate_simudata is an input that control angle of push
if __name__ == '__main__':
    rx = float(sys.argv[1])
    ry = float(sys.argv[2])
    xvel = float(sys.argv[3])
    yvel = float(sys.argv[4])
    simu_steps = int(float(sys.argv[5]) * 10)
    init_angle = float(sys.argv[6])
    rx2 = float(sys.argv[7])
    ry2 = float(sys.argv[8])
    xvel2 = float(sys.argv[9])
    yvel2 = float(sys.argv[10])
    simu_steps2 = int(float(sys.argv[11]) * 10)
    init_angle2 = float(sys.argv[12])
    rtor = float(sys.argv[13])
    rtor2 = float(sys.argv[14])
    gx = float(sys.argv[15])
    gy = float(sys.argv[16])
    gx2 = float(sys.argv[17])
    gy2 = float(sys.argv[18])
    
    
    world = b2WorldInterface(False)
    oshape, osize, ofriction, odensity, bfriction, hand_shape, hand_size  = 'circle', 1, 0.01, 0.05, 0.01, 'rectangle', (1,0.3) #'circle', 0.3#
    #thing,base = make_thing(500, 500, world, oshape, osize, ofriction, odensity, bfriction, (0,0))
    base = make_base(500, 500, world)
    thing = make_1thing(base, world, 'rectangle', (0.5,0.5), ofriction, odensity, (0, 2))
    thing2 = make_1thing(base, world, 'circle', 1, ofriction, odensity, (0,-2))
    #xvel = np.cos(init_angle)*5;
    #yvel = np.sin(init_angle)*5;
    robot = end_effector(world, (rx,ry), base, init_angle, hand_shape, hand_size)
    robot2 = end_effector(world, (rx2,ry2), base, init_angle2, hand_shape, hand_size)
    (ret1, ret2) = simu_push_2robot2thing(world, thing, thing2, robot, robot2, base, xvel, yvel, xvel2, yvel2, rtor, rtor2, simu_steps, simu_steps2)
    #print ret1, ret2
    ret1 = np.linalg.norm(np.array([gx, gy]) - ret1)
    ret2 = np.linalg.norm(np.array([gx2, gy2]) - ret2)
    sys.stdout.write(str(ret1+ret2))