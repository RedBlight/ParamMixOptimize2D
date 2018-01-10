clc;
clear all;
format long;

addpath( 'FuncWithD' );
addpath( 'ParamPolContour' );
SuperFormulaList;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

superFormula = superFormula_5Star;

[ tanSmoother, arcPower ] = SearchOptimisedValues( superFormula, 1000, 8, 2^9, 2^11 );

disp( 'Final tanSmoother:' );
disp( tanSmoother );
disp( 'Final arcPower:' );
disp( arcPower );































