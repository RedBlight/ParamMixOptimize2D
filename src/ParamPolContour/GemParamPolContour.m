classdef GemParamPolContour < B_ParamPolContour
    
    properties
        funcWithD_;
        contourList_;
        powerList_;
        powerListT_;
        powerListNormal_;
        powerListNormalT_;
        contourN_;
        powerSum_;
    end
    
    methods
        
        % ( funcWithD, contourList, powerList )
        % ( funcWithD, contourList, powerList, splineN )
        function obj = GemParamPolContour( varargin )
            obj.funcWithD_ = varargin{ 1 };
            obj.contourList_ = varargin{ 2 };
            obj.powerList_ = varargin{ 3 };
            if nargin > 3
                obj.splineN_ = varargin{ 4 };
            end
            obj.powerListT_ = transpose( obj.powerList_ );
            obj.contourN_ = length( obj.powerList_ );
            obj.powerSum_ = sum( abs( obj.powerList_ ) );
            obj.powerListNormal_ = obj.powerList_ ./ obj.powerSum_;
            obj.powerListNormalT_ = transpose( obj.powerListNormal_ );
            obj.SolveParamDifferential();
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Abstract functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function newParamDifferential = NewParamDifferential_( obj, polParam )
            polParamN = length( polParam );
            subParamDMat = zeros( obj.contourN_, polParamN );
            for k = 1 : obj.contourN_
                subParamDMat( k, : ) = obj.contourList_{ k }.NewParamD_( polParam, 1 );
            end
            powerMat = repmat( obj.powerListNormalT_, 1, polParamN );
            newParamDifferential = prod( subParamDMat .^ powerMat )  ;
        end
        
        function varargout = NewParamD_( obj, polParam, drvN )
            singleMode = nargout == 1;
            polParamN = length( polParam );
            
            subParamMat = cell( obj.contourN_, drvN );
            for k = 1 : obj.contourN_
                [ subParamMat{ k, 1 : drvN } ] = obj.contourList_{ k }.NewParamD_( polParam, drvN );
            end
            
            drvIndex = 1;
            if( drvN >= drvIndex )
                subParamD1Mat = zeros( obj.contourN_, polParamN );
                for k = 1 : obj.contourN_
                    subParamD1Mat( k, : ) = subParamMat{ k, 1 };
                end
                powerMat = repmat( obj.powerListNormalT_, 1, polParamN );
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = obj.newParamFactor_ .* ( ...
                        prod( subParamD1Mat .^ powerMat ) ...
                    );
                end
            end
            
            drvIndex = 2;
            if( drvN >= drvIndex )
                subParamD2Mat = zeros( obj.contourN_, polParamN );
                for k = 1 : obj.contourN_
                    subParamD2Mat( k, : ) = subParamMat{ k, 2 };
                end
                if( ~singleMode || drvIndex == drvN )                    
                    sumN = zeros( 1, polParamN );
                    for n = 1 : obj.contourN_
                        prodK = ones( 1, polParamN );
                        for k = 1 : obj.contourN_
                            if k == n
                                continue;
                            end
                            prodK = prodK .* ( subParamD1Mat( k, : ) .^ powerMat( k, : ) );
                        end
                        sumN = sumN + powerMat( n, : ) .* ( subParamD1Mat( n, : ) .^ ( powerMat( n, : ) - 1 ) ) .* subParamD2Mat( n, : ) .* prodK; 
                    end
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = obj.newParamFactor_ .* ( ...
                        sumN ...
                    );
                end
            end
            
            drvIndex = 3;
            if( drvN >= drvIndex )
                subParamD3Mat = zeros( obj.contourN_, polParamN );
                for k = 1 : obj.contourN_
                    subParamD3Mat( k, : ) = subParamMat{ k, 3 };
                end
                if( ~singleMode || drvIndex == drvN )
                    sumN = zeros( 1, polParamN );
                    for n = 1 : obj.contourN_
                        prodK = ones( 1, polParamN );
                        for k = 1 : obj.contourN_
                            if k == n
                                continue;
                            end
                            prodK = prodK .* ( subParamD1Mat( k, : ) .^ powerMat( k, : ) );
                        end
                        sumK = zeros( 1, polParamN );
                        for k = 1 : obj.contourN_
                            if k == n
                                continue;
                            end
                            prodM = ones( 1, polParamN );
                            for m = 1 : obj.contourN_
                                if m == n || m == k
                                    continue;
                                end
                                prodM = prodM .* ( subParamD1Mat( m, : ) .^ powerMat( m, : ) );
                            end
                            sumK = sumK + powerMat( k, : ) .* ( subParamD1Mat( k, : ) .^ ( powerMat( k, : ) - 1 ) ) .* subParamD2Mat( k, : ) .* prodM;
                        end    
                        sumN = sumN + powerMat( n, : ) .* ( ...
                            ( ( powerMat( n, : ) - 1 ) .* subParamD1Mat( n, : ) .^ ( powerMat( n, : ) - 2 ) .* subParamD2Mat( n, : ) .* subParamD2Mat( n, : ) + ...
                            subParamD1Mat( n, : ) .^ ( powerMat( n, : ) - 1 ) .* subParamD3Mat( n, : ) ) .* prodK + ...
                            subParamD1Mat( n, : ) .^ ( powerMat( n, : ) - 1 ) .* subParamD2Mat( n, : ) .* sumK ...
                        );
                    end
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = obj.newParamFactor_ .* ( ...
                        sumN ...
                    );
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Renamed functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function varargout = GemParamD( obj, polParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.NewParamD_( polParam, drvN );
        end
        
        function varargout = ArcLengthDgem( obj, gemParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.ArcLengthDnew( gemParam, drvN );
        end
        
        function varargout = TanAngleDgem( obj, gemParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.TanAngleDnew( gemParam, drvN );
        end
        
        function varargout = PointFullDgem( obj, gemParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.PointFullDnew( gemParam, drvN );
        end
        
        
    end
    
end

