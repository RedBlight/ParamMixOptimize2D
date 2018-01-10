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

N = 2^10;
paramStart = -pi;
paramEnd = pi;
paramRange = paramEnd - paramStart;
paramDiff = paramRange / N;
param = ( 0 : N ) * paramDiff + paramStart + paramDiff / 2;
param( end ) = param( 1 );

[ GALD1, GALD2, GALD3 ] = testContour.ArcLengthDnew( param, 3 );
[ PALD1, PALD2, PALD3 ] = testContour.ArcLengthDpol( param, 3 );

figure( 1 );
plot( param, GALD1, 'red', param, PALD1, 'blue', 'LineSmoothing', 'on' );

figure( 2 );
plot( param, GALD2, 'red', param, PALD2, 'blue', 'LineSmoothing', 'on' );

figure( 3 );
plot( param, GALD3, 'red', param, PALD3, 'blue', 'LineSmoothing', 'on' );































