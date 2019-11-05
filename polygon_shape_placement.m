%------------------------------------------------------------------------------%
%---------------------- polygon shape placement function ----------------------%
%                                                                              %
% im = polygon_shape_placement (imsize, polyshape, distparams, N, polyrotation)% 
% creates the 2D binary image im of size imsize with the polinogonal shapes    %
% specified by polyshape placed along it. The size of the polygonal shapes va- %
% ries according to a normal distribution defined by the distribution parame-  %
% ters given by the polyshape input argument.The center-to-center distance     %
% along the x and y directions between the polygonal shapes follows a pseudo-  %
% normal distribution defined by distparams. The im image has the option to be %
% an image composed of N binary layers. Additionally, with the argument poly-  %
% rotation you can make random rotations on the polygonal shapes.              %
%                                                                              %
%  -imsize:       Vector of two elements specifying the dimensions of the 2D   %
%                 output image, im                                             %
%  >> imsize = [m, n]  % im of size m x n                                      %
%                                                                              %
%  -polyshape:    Array of two cells defining the pattern of the polygonal     %
%                 shape. The first element is a string specifying the type of  %
%                 shape desired (see the table below to know the types of sha- %
%                 pes available). The second element is a vector of two ele-   %
%                 ments specifying the mean & standar deviation of a pseudo-   %
%                 normal distribution used to define the size of each polygon. %
%  >> polyshape = {'square', [meansize, stdsize]}                              %
%                                      squares of a x a where 'a'              %
%                                      follows a ~ N(meansize, stdsize)        %
%                                                                              %
%  -distparams:   Vector of two elements specifying the mean & standar devia-  %
%                 tion error of a pseudo-normal distribution used to define    %
%                 the center-to-center distance                                %
%  >> distparams = [mean, std]                                                 %
%                                                                              %
%  -N:            Scalar defining the number of layers of the binary im image. %
%                 Each layer represents a subsample of the resulting image     %
%  >> N = 2  % Two layers                                                      %
%                                                                              %
%  -polyrotation  By default this variable is 'on'. You can set it to 'off' to %
%                 deactivate the rotations.                                    %
%  >> polyrotation = 'off'  % deactivate the rotation                          %
%                                                                              %
%  Table: polyshape options                                                    %
%  +------------+-------------+--------------+---------+-----------------+     %
%  | shape      |  long form  |  short form  |  size   |  description    |     %
%  +------------+-------------+--------------+---------+-----------------+     %
%  | general    |             |  '4', '5'... | [d]     | d: diameter*    |     %
%  | polygon    |             |  '6', '7'... | [n]     | n: corners**    |     %
%  +------------+-------------+--------------+---------+-----------------+     %
%  | square     | 'square'    |  's'         | [a]     | a: side         |     %
%  +------------+-------------+--------------+---------+-----------------+     %
%  | rectangle  | 'rectangle' |  'r'         | [a, b]  | a: base         |     %
%  |            |             |              |         | b: lateral side |     %
%  +------------+-------------+--------------+---------+-----------------+     %
%  | circle     | 'circle'    |  'c'         | [d]     | r: diameter     |     %
%  +------------+-------------+--------------+---------+-----------------+     %
%  | triangle   | 'triangle'  |  '^'         | [b, h]  | b: base         |     %
%  |            |             |              |         | h: height       |     %
%  +------------+-------------+--------------+---------+-----------------+     %
%  | triangle   |             |  'v'         | [b, h]  | invertided      |     %
%  +------------+-------------+--------------+---------+-----------------+     %
%  | triangle   |             |  '<'         | [b, h]  | left orient..   |     %
%  +------------+-------------+--------------+---------+-----------------+     %
%  | triangle   |             |  '>'         | [b, h]  | right orient..  |     %
%  +------------+-------------+--------------+---------+-----------------+     %
%                                                                              %
%  *  The general polygons are created by means of a circumference, in which   %
%        the polygons are circumscribed. So, the diameter d is the diameter of %
%        the circumference                                                     %
%  ** n is the number of corners you want the polygon to have. This value is   %
%        taken directly from the numeric character that is placed as input to  %
%        the general polygon option                                            %
%                                                                              %
%                                                                              %
%  NOTE: Now all the size parameters of the polygonal shapes (a, b, d and h)   %
%        will be defined by the distribution parameters specified in polyshape %
%                                                                              %
%------------------------------------------------------------------------------%

function im = polygon_shape_placement(imsize, polyshape, distparams, N, polyrotation)

%--- Initializations
	
	% Getting polyrotation
	if ~exist('polyrotation')
		polyrotation = 'on';
	end
	
	% Getting the distribution parameters
	pmean = distparams(1);
	pstd = floor(distparams(2) / 2);

	% Estimating the number of shapes per dimension
	nshapes = [ceil(imsize(1)/pmean), ceil(imsize(2)/pmean)];

%--- Get the piths for all the shapes

	% Getting it from a pseudo normally distribution defined by
	%	the pmean & pstd
	pitchs = randi([-pstd, pstd], [nshapes, 2]);

%--- Get the size for all the shapes

	% Getting it from a pseudo normally distribution defined by
	%	the meansize & stdsize ditribution parameters
	sizeparams = polyshape{2};
	sizemax = sizeparams(1) + sizeparams(2);
	sizemin = sizeparams(1) - sizeparams(2);
	sizes = randi([sizemin, sizemax], [nshapes, 2]);

%--- Get the angles for the rotation

	% Getting uniform random rotation
	if strcmpi(polyrotation, 'on')
		rotations = 360*rand([nshapes]);
	end

%--- Placement radonming the shape at the image im

	% Creating the output binary image
	im = zeros( [(nshapes+2)*pmean, N] );

	% Creating a variable for counting the shapes
	ishape = 1;

	% Sweep over whole image
	%	loop extern to control the sweep at x direction
	for x = 1:nshapes(2)

		% Variable for counting the shapes at the current column
		term = 1;
		ishape = x;

		%	loop intern to control the sweep at y direction
		for y = 1:nshapes(1)

			% Getting the dynaminc size mask
			polymask = getshape({polyshape{1}, sizes(y,x,:)});

			% Making the polygonal rotation
			if strcmpi(polyrotation, 'on')
				polymask = imrotate(polymask, rotations(x,y), 'bilinear');
			end

			% Defining dynamic variables to control the sweep
			xmask = size(polymask, 2) - 1;
			ymask = size(polymask, 1) - 1;
			xi = floor( (pmean - xmask) / 2 );
			yi = floor( (pmean - ymask) / 2 );

			dx = pitchs(y, x, 1);
			dy = pitchs(y, x, 2);

			xs = x * pmean + 1 + xi + dx;
			ys = y * pmean + 1 + yi + dy;

			xe = xs + xmask;
			ye = ys + ymask;

			% Getting the current layer
			ilayer = mod(ishape-1, N) + 1;

			% Placement the shape at the corresponding layer
			im(ys:ye, xs:xe, ilayer) = im(ys:ye, xs:xe, ilayer) + polymask;

			% Counting the ishape variable
			ishape = ishape + term;
		end
	end
	offset = pmean + 1;
	im = im(offset:offset+imsize(1)-1, offset:offset+imsize(2)-1, :);
	im = min(im, 1);
end

%------------------------------------------------------------------------------%
%----------------------------- getshape function ------------------------------%
%                                                                              %
% polymask = getshape(polyshape) creates a polygonal shape (refered here as a  %
% polygonal mask) with features given by polyshape cell array (input argument).%
% Here the polyshape have the same specification as the corresponding input    %
% parameter of main function.                                                  %
%                                                                              %
%------------------------------------------------------------------------------%
function polymask = getshape(polyshape)
	shape = polyshape{1};
	masksize = polyshape{2};

	if strcmpi(shape, 's') || strcmpi(shape, 'square')
		a = masksize(1);
		polymask = squaremask(a);
	elseif strcmpi(shape, 'r') || strcmpi(shape, 'rectangle')
		a = masksize(1);
		b = masksize(2);
		polymask = rectanglemask(a, b);
	elseif strcmpi(shape, 'c') || strcmpi(shape, 'circle')
		d = masksize(1);
		polymask = circlemask(d);
	elseif str2num(shape)
		d = masksize(1);
		n = str2num(shape);
		polymask = generalmask(d, n);
	elseif strcmpi(shape, '^') || strcmpi(shape, 'triangle') || ...
		   strcmpi(shape, 'v') || strcmpi(shape, '<') || strcmpi(shape, '>')
		b = masksize(1);
		h = masksize(2);
		if strcmpi(shape, 'triangle')
			orientation = '^';
		else
			orientation = shape;
		end
		polymask = trianglemask(b, h, orientation);
	end
end

%------------------------------------------------------------------------------%
%----------------------------- squaremask function ----------------------------%
%                                                                              %
% polymask = squaremask(a) creates a square polygonal shape (refered here as a %
% square mask) with side 'a' according with the description of main function   %
%                                                                              %
%------------------------------------------------------------------------------%
function polymask = squaremask(a)
	polymask = ones(a + 1);
end


%------------------------------------------------------------------------------%
%--------------------------- rectanglemask function ---------------------------%
%                                                                              %
% polymask = rectanglemask(a, b) creates a rectangle polygonal shape (refered  %
% here as a rectangle mask) with sides 'a' and 'b' according with the descrip- %
% tion of main function                                                        %
%                                                                              %
%------------------------------------------------------------------------------%
function polymask = rectanglemask(a, b)
	polymask = ones(b + 1, a + 1);
end


%------------------------------------------------------------------------------%
%----------------------------- circlemask function ----------------------------%
%                                                                              %
% polymask = circlemask(d)) creates a circle polygonal shape (refered here as  %
% a circle mask) with diameter 'd' according with the description of main      %
% function                                                                     %
%                                                                              %
%------------------------------------------------------------------------------%
function polymask = circlemask(d)
	r = floor(d/2);
	[maskcols, maskrows] = meshgrid(1:d+1, 1:d+1);
	x0 = r + 1;
	y0 = x0;
	polymask = (maskrows - y0).^2 + (maskcols - x0).^2 <= (r).^2;
end


%------------------------------------------------------------------------------%
%--------------------------- trianglemask function ----------------------------%
%                                                                              %
% polymask = trianglemask(b, h, orientation) creates a triangle polygonal sha- %
% pe  (refered here as a triangle mask) with base 'b', height 'h' and orienta- %
% tion according with the description of main function                         %
%                                                                              %
%------------------------------------------------------------------------------%
function polymask = trianglemask(b, h, orientation)
	switch orientation
		case '^'
			[maskcols, maskrows] = meshgrid(1:b+1, h+1:-1:1);
			x0 = floor(b/2);
			y0 = h;
			m = 2*h/b;
			polymask = (maskrows-y0) <= -abs( m*(maskcols - x0) );
		case 'v'
			[maskcols, maskrows] = meshgrid(1:b+1, 1:h+1);
			x0 = floor(b/2);
			y0 = h;
			m = 2*h/b;
			polymask = (maskrows-y0) <= -abs( m*(maskcols - x0) );
		case '>'
			[maskcols, maskrows] = meshgrid(1:h+1, 1:b+1);
			y0 = floor(b/2);
			x0 = h;
			m = 2*h/b;
			polymask = (maskcols-x0) <= -abs( m*(maskrows - y0) );
		case '<'
			[maskcols, maskrows] = meshgrid(h+1:-1:1, 1:b+1);
			y0 = floor(b/2);
			x0 = h;
			m = 2*h/b;
			polymask = (maskcols-x0) <= -abs( m*(maskrows - y0) );
	end	
end


%------------------------------------------------------------------------------%
%--------------------------- generalmask function -----------------------------%
%                                                                              %
% polymask = generalmask(d, n) creates a general polygonal shape (refered here %
% as a general poligonal mask) where the polygon is ciscunscript in a circle   %
% of diameter 'd' and it have 'n' corners,  according with the description of  %
% main function                                                                %
%                                                                              %
%------------------------------------------------------------------------------%
function polymask = generalmask(d, n)
	%--- Initializations
	tita = [0:(2*pi/n):2*pi]+pi/2;
	r = floor(d/2);
	x = r*cos(tita);
	y = r*sin(tita);
	[X, Y] = meshgrid(1:d, d:-1:1);
	x0 = r; y0 = x0;
	polymask = ones(size(X));

	Xr = X - x0;
	Yr = Y - y0;

	%--- Generate the mask
	for i=1:n
		m = (y(i+1)-y(i))/(x(i+1)-x(i));
		if (y(i)+y(i+1))/2 >= 0
			p = Yr <= m*( Xr - x(i) ) + y(i);
		else
			p = Yr >= m*( Xr - x(i) ) + y(i);
		end
		polymask = polymask.*p;
	end
end