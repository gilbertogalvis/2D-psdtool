clc; close all; clear all;

polyshape = {'5', [30,15]};
imsize = [550, 550];
distparams = [55, 25];
N = 3;

im = polygon_shape_placement(imsize, polyshape, distparams, N, 'on');

imshow(im)
size(im)