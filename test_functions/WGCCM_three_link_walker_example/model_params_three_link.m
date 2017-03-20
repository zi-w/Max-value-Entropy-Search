function [r,m,Mh,Mt,L,g]=model_params_three_link
%model_params_three_link.m
%
%   Model parameters for three-link legged biped.
%
% This file is associated with the book Feedback Control of Dynamic 
% Bipedal Robot Locomotion by Eric R. Westervelt, Jessy W. Grizzle, 
% Christine Chevallereau, Jun-Ho Choi, and Benjamin Morris published 
% by Taylor & Francis/CRC Press in 2007.
% 
% Copyright (c) 2007 by Eric R. Westervelt, Jessy W. Grizzle, Christine
% Chevallereau, Jun-Ho Choi, and Benjamin Morris.  This code may be
% freely used for noncommercial ends.  If use of this code in part or in
% whole results in publication, proper citation must be included in that
% publication.  This code comes with no guarantees or support.
% 
% Eric Westervelt
% 20 February 2007

r=1;   % length of a leg
m=5;   % mass of a leg
Mh=15; % mass of hips
Mt=10; % mass of torso
L=0.5; % distance between hips and torso
g=9.8; % acceleration due to gravity