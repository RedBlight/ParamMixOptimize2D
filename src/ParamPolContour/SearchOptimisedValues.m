function [ tanSmoother, arcPower ] = SearchOptimisedValues( funcWithD, generationN, populationN, splineN, sampleN )
    
    % tan smoother
    smtSearchStart = 0;
    smtSearchEnd = 50;

    % arc power
    alpSearchStart = -50;
    alpSearchEnd = 50;

    smtRange = smtSearchEnd - smtSearchStart;
    alpRange = alpSearchEnd - alpSearchStart;
    smtDiff = smtRange / 2;
    alpDiff = alpRange / 2;
    smtDiffStart = smtDiff;
    alpDiffStart = alpDiff;
    smt = smtSearchStart + smtDiff;
    alp = alpSearchStart + alpDiff;
    
    paramLength = 2 * pi;
    paramStart = -pi;
    paramDiff = paramLength / ( sampleN - 1 );
    paramVec = ( 0 : ( sampleN - 1 ) ) * paramDiff + paramStart + paramDiff / 2;

    maxVal = 1e200;
    
    %%%%%%%%%%%%%%%%
    smtLog = zeros( 1, generationN );
    alpLog = zeros( 1, generationN );
    maxLog = zeros( 1, generationN );
    minIt = 1;
    %%%%%%%%%%%%%%%%

    for k = 1 : generationN
        angleVec = ( 2 * pi / populationN ) .* ( ( 0 : ( populationN - 1 ) ) + rand() );
        cosVec = cos( angleVec ) .* 2 .* rand( 1, populationN );
        sinVec = sin( angleVec ) .* 2 .* rand( 1, populationN );
        smtDiffVec = smtDiff .* cosVec;
        alpDiffVec = alpDiff .* sinVec;
        smtVec = smt + smtDiffVec;
        smtVec( smtVec < smtSearchStart ) = smtSearchStart;
        smtVec( smtVec > smtSearchEnd ) = smtSearchEnd;
        alpVec = alp + alpDiffVec;
        alpVec( alpVec < alpSearchStart ) = alpSearchStart;
        alpVec( alpVec > alpSearchEnd ) = alpSearchEnd;

        maxCur = zeros( 1, populationN );
        parfor m = 1 : populationN
            tasCon = TasParamPolContour( funcWithD, splineN, smtVec( m ) );
            arcCon = ArcParamPolContour( funcWithD, splineN );
            gemCon = GemParamPolContour( ...
                funcWithD, ...
                { tasCon, arcCon }, ...
                [ 1, alpVec( m ) ], ...
                splineN ...
            );
            GD3 = gemCon.PointFullDgem( paramVec, 3 );
            maxCurX = max( GD3( 1, : ) );
            maxCurY = max( GD3( 2, : ) );
            maxCur( m ) = sqrt( maxCurX .* maxCurX + maxCurY .* maxCurY );
        end
        
        [ maxItVal, maxItIndex ] = min( maxCur );
        
        if maxItVal < maxVal
            maxVal = maxItVal;
            smt = smtVec( maxItIndex );
            alp = alpVec( maxItIndex );
            smtDiff = smtDiffStart;
            alpDiff = alpDiffStart;
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            minIt = minIt + 1;
            smtLog( minIt ) = smt;
            alpLog( minIt ) = alp;
            %%%%%%%%%%%%%%%%%%%%%%%%
        else
            smtDiff = smtDiff ./ 2;
            alpDiff = alpDiff ./ 2;
        end
        
        if smtDiff < 1e-8 && alpDiff < 1e-8
            smtDiff = smtDiffStart;
            alpDiff = alpDiffStart;
        end

        %%%%%%%%%%%%%%%%%%%%%%
        disp( smt );
        disp( alp );
        disp( maxVal );
        disp( k );
        maxLog( k ) = maxVal;
        %%%%%%%%%%%%%%%%%%%%%%
        
        
    end
    
    %%%%%%%%%%%%%%%%
    smtLog = smtLog( 1 : minIt );
    alpLog = alpLog( 1 : minIt );
    figure( 701 );
    plot( smtLog, alpLog, 'red-x', 'LineSmoothing', 'on' );
    figure( 702 );
    plot( 1 : generationN, maxLog, 'red', 'LineSmoothing', 'on' );
    %%%%%%%%%%%%%%%%
    
    tanSmoother = smt;
    arcPower = alp;
end
