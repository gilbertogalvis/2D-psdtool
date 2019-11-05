#################################################
2D-psdtool - 2D Polygonal Shape Distribution Tool
#################################################

:Authors: John Kay, Gilberto Galvis
:Email: john.kays2020@gmail.com, galvisgilberto@gmail.com
:Version: $revision: 0.1.1 $

This project provides a tool to generate 2D images that have randomly distributed polygonal shapes along it.

Requirements
------------

- MatLab: The interface tool is based on MatLab, so it is required to have MatLab installed

Installation
------------

- Clone this repository on your machine either using http or ssh

Usage
-----

The tool is basically composed of a function called ``polygon_shape_placement`` . Then, just execute this function with the desired input parameters. In this way, an array of data with the image content is returned as output. This matrix can then be manipulated either to display, save or manipulate the image. 

Usage Example
=============

Here is an example that can be executed using a MatLab script (such as the test.m script placed in this repository) or also directly in the MatLab Command Window.

.. code:: shell
	
	clc; close all; clear all;

	polyshape = {'5', [30,15]};
	imsize = [550, 550];
	distparams = [55, 25];
	N = 3;

	im = polygon_shape_placement(imsize, polyshape, distparams, N, 'on');

	imshow(im)

Please feel free to change the parameters to the values you want

Notes
+++++

- If you need to know in detail about the input parameters of the function, please open the file ``polygon_shape_placement.m``. In that file you can see the source code as well as the explanatory documentation of the function

- Please feel free to review the report.pdf document where we show more examples of this tool.

