clc;
clear all;
format long;

addpath( 'FuncWithD' );
addpath( 'ParamPolContour' );
SuperFormulaList;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% superFormula = superFormula_5Star;
% 
% tasContour = TasParamPolContour( superFormula, 2^16, 0.744234858869572 );
% arcContour = ArcParamPolContour( superFormula, 2^16 );
% gemContour = GemParamPolContour( ...
%     superFormula, ...
%     { tasContour, arcContour }, ...
%     [ 1, 2.553424164133645 ] ...
% );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

superFormula = superFormula_Spiky;

tasContour = TasParamPolContour( superFormula, 2^16, 4.111607467940276 );
arcContour = ArcParamPolContour( superFormula, 2^16 );
gemContour = GemParamPolContour( ...
    superFormula, ...
    { tasContour, arcContour }, ...
    [ 1, 2.472519747712406 ] ...
);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

testContour = arcContour;

N = 2^6;
paramStart = -pi;
paramEnd = pi;
paramRange = paramEnd - paramStart;
paramDiff = paramRange / N;
param = ( 0 : N ) * paramDiff + paramStart + paramDiff / 2;
param( end ) = param( 1 );

newToPol = testContour.PolParam( param );

nVec = 1 : N+1;
nVecSpline = 1 : 0.1 : N+1;

[ GCD0, GCD1, GCD2, GCD3 ] = testContour.PointFullDnew( param, 3 );
[ PCD0, PCD1, PCD2, PCD3 ] = testContour.PointFullDpol( param, 3 );

Gx = spline( nVec, GCD0(1,:), nVecSpline );
Gy = spline( nVec, GCD0(2,:), nVecSpline );

Px = spline( nVec, PCD0(1,:), nVecSpline );
Py = spline( nVec, PCD0(2,:), nVecSpline );

figure( 51 );
plot( GCD0(1,:), GCD0(2,:), 'red-o', Gx, Gy, 'blue', 'LineSmoothing', 'on' );

figure( 52 );
plot( PCD0(1,:), PCD0(2,:), 'red-o', Px, Py, 'blue', 'LineSmoothing', 'on' );

% figure( 53 );
% plot( GCD2(1,:), GCD2(2,:), 'red', PCD2(1,:), PCD2(2,:), 'blue', 'LineSmoothing', 'on' );
% 
% figure( 54 );
% plot( GCD3(1,:), GCD3(2,:), 'red', PCD3(1,:), PCD3(2,:), 'blue', 'LineSmoothing', 'on' );



























