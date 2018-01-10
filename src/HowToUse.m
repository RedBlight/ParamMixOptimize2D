% Code files insidefolders can be accessible by adding the folder path:
addpath( 'FuncWithD' );
addpath( 'ParamPolContour' );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% B_FuncWithD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "B_FuncWithD" is an abstract class. Classes that hold mathematical
% functions and their derivatives can be derived from this class. If a
% function is to be used multiple times with same parameters, it has a
% cache mechanism inside it and can be locked to return only cached values
% instead of doing calculations again with same parameters. Lock and cache
% mechanism is defined in the abstract class and derived classes should
% only define mathematical functions and their derivatives. "SuperFormula"
% is an example derived class.

% Example usage of "SuperFormula", it has 4 derivatives defined in it:
% a, b,
% m, n1, n2, n3,
% alpha, beta,
% center, rotation
superFormula = SuperFormula( ...
    5.08, 6.90, ...
    10.0, 0.139, 0.139, 0.972, ...
    4.29, 0.10, ...
    [ 0; 0 ], 0 ...
);

% An example closed parameter vector:
N = 256;
paramStart = -pi;
paramEnd = pi;
paramRange = paramEnd - paramStart;
paramDiff = paramRange / N;
param = ( 0 : N ) * paramDiff + paramStart + paramDiff / 2;
param( end ) = param( 1 );

% All derivatives:
[ RD0, RD1, RD2, RD3, RD4 ] = superFormula.RadiusFullD( param, 4 );

% 0, 1, and 2nd derivatives:
[ RD0, RD1, RD2 ] = superFormula.RadiusFullD( param, 2 );

% Only 4th derivative:
RD4 = superFormula.RadiusFullD( param, 4 );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% B_ParamPolContour
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "B_ParamPolContour" is an abstract class for parametrized polar contours.
% It has all methods for arc length, curvature, cartesian points etc.
% defined in it. These methods are automaically inherited by derived
% classes and is accessible by them. Derived classes should only define
% parameter differential function and parameter derivatives with respect
% to polar parameter, and the rest is done in the abstract class.

% Currently, we have 8 different parametrizations:

% Good old polar parametrization:
% PolParamPolContour( funcWithD, splineN = 2^16 )
polContour = PolParamPolContour( superFormula, 2^16 );

% Arc length parametrization:
% ArcParamPolContour( funcWithD, splineN = 2^16 )
arcContour = ArcParamPolContour( superFormula, 2^16 );

% Tangent angle parametrization:
% TanParamPolContour( funcWithD, splineN = 2^16 )
tanContour = TanParamPolContour( superFormula, 2^16 );

% Smooth tangent angle parametrization:
% TasParamPolContour( funcWithD, splineN = 2^16, smoother = 1 )
tasContour = TasParamPolContour( superFormula, 2^16, 1 );

% Curvature parametrization:
% CurParamPolContour( funcWithD, splineN = 2^16 )
curContour = CurParamPolContour( superFormula, 2^16 );

% Smooth curvature parametrization:
% CusParamPolContour( funcWithD, splineN = 2^16, smoother = 1 )
cusContour = CusParamPolContour( superFormula, 2^16, 1 );

% Superposition of parameters parametrization:
% SupParamPolContour( funcWithD, contourList, coefList, splineN = 2^16 )
supContour = SupParamPolContour( ...
    superFormula, ...
    { tasContour, arcContour }, ...
    [ 1, 1 ], ...
    2^16 ...
);

% Geometric mean of parameters parametrization:
% GemParamPolContour( funcWithD, contourList, powerList, splineN = 2^16 )
gemContour = SupParamPolContour( ...
    superFormula, ...
    { tasContour, arcContour }, ...
    [ 1, 1 ], ...
    2^16 ...
);

% splineN is the amount of points used to solve parameter differential and
% build the parameter function. Parameter functions and their differentials
% behave like a polynomial, solving the differential with Newton-Cotes
% formulas converges very quickly. In most situations, error is 0 with
% 2^16 points. After solving the differential, newParam-polParam
% correspondence is formed with spline interpolation. Spline interpolation
% works much, MUCH MUCH, better than linear interpolation for this.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% GemContour Optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% After some observations, I found that geometric mean of smooth tangent
% angle parametrization and arc lenght parametrization generates the
% smallest derivatives for the cartesian contour and approximate it's area
% much better with smaller numbers of samples.

% As an example, 3rd derivatives from the gemContour below are ~100 times
% smaller than arcContour or polContour.
tasContour = TasParamPolContour( superFormula, 2^16, 0.739122405275166 );
arcContour = ArcParamPolContour( superFormula, 2^16 );
gemContour = GemParamPolContour( ...
    superFormula, ...
    { tasContour, arcContour }, ...
    [ 1, 2.552762845719015 ] ...
);

% But, we don't have any formulation for the best tangent angle smoother
% and the best arc length power. "SearchOptimisedValues" function exists to
% find those values with a crude, simple "genetic algorithm", or maybe
% "stochastic hill climbing". 

% SearchOptimisedValues( funcWithD, generationN, populationN, splineN, sampleN )
[ tanSmoother, arcPower ] = SearchOptimisedValues( superFormula, 1000, 8, 2^9, 2^11 );

% Genetic algorithm starts with a value set in the middle of the search
% space. It generates "populationN" amount of new value sets, with some
% random mutations. Then selects the most optimised value set in the
% population, and starts a new generation again from it. Repeats this
% process for "generationN" times.
% This process is slow and requires a few minutes to settle on a good
% optimized value.
