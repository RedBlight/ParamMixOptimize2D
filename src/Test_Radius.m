clc;
clear all;
format long;

addpath( 'FuncWithD' );
addpath( 'ParamPolContour' );
SuperFormulaList;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

superFormula = superFormula_5Star;

tasContour = TasParamPolContour( superFormula, 2^16, 0.744234858869572 );
arcContour = ArcParamPolContour( superFormula, 2^16 );
gemContour = GemParamPolContour( ...
    superFormula, ...
    { tasContour, arcContour }, ...
    [ 1, 2.553424164133645 ] ...
);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

testContour = gemContour;

N = 2^16;
paramStart = -pi;
paramEnd = pi;
paramRange = paramEnd - paramStart;
paramDiff = paramRange / N;
param = ( 0 : N ) * paramDiff + paramStart + paramDiff / 2;
param( end ) = param( 1 );

[ GRD0, GRD1, GRD2, GRD3 ] = testContour.RadiusFullDnew( param, 3 );
[ PRD0, PRD1, PRD2, PRD3 ] = testContour.RadiusFullDpol( param, 3 );

figure( 1 );
plot( param, GRD0, 'red', param, PRD0, 'blue', 'LineSmoothing', 'on' );

figure( 2 );
plot( param, GRD1, 'red', param, PRD1, 'blue', 'LineSmoothing', 'on' );

figure( 3 );
plot( param, GRD2, 'red', param, PRD2, 'blue', 'LineSmoothing', 'on' );

figure( 4 );
plot( param, GRD3, 'red', param, PRD3, 'blue', 'LineSmoothing', 'on' );



























