classdef SupParamPolContour < B_ParamPolContour
    
    properties
        funcWithD_;
        contourList_;
        coefList_;
        coefListT_;
        coefListNormal_;
        coefListNormalT_;
        contourN_;
        coefSum_;
    end
    
    methods
        
        % ( funcWithD, contourList, coefList )
        % ( funcWithD, contourList, coefList, splineN )
        function obj = SupParamPolContour( varargin )
            obj.funcWithD_ = varargin{ 1 };
            obj.contourList_ = varargin{ 2 };
            obj.coefList = varargin{ 3 };
            if nargin > 3
                obj.splineN_ = varargin{ 4 };
            end
            obj.coefListT_ = transpose( coefList );
            obj.contourN_ = length( coefList );
            obj.coefSum_ = sum( coefList );
            obj.coefListNormal_ = coefList ./ obj.coefSum_;
            obj.coefListNormalT_ = transpose( obj.coefListNormal_ );
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
            coefMat = repmat( obj.coefListNormalT_, 1, polParamN );
            newParamDifferential = sum( coefMat .* subParamDMat ) ;
        end
        
        function varargout = NewParamD_( obj, polParam, drvN )
            singleMode = nargout == 1;
            polParamN = length( polParam );
            
            subParamMat = cell( obj.contourN_, drvN );
            for k = 1 : obj.contourN_
                [ subParamMat{ k, 1 : drvN } ] = obj.contourList_{ k }.NewParamD_( polParam, drvN );
            end
            
            subParamDMat = zeros( obj.contourN_, polParamN );
            
            drvIndex = 1;
            if( drvN >= drvIndex )
                if( ~singleMode || drvIndex == drvN )
                    for k = 1 : obj.contourN_
                        subParamDMat( k, : ) = subParamMat{ k, drvIndex };
                    end
                    coefMat = repmat( obj.coefListNormalT_, 1, polParamN );
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = obj.newParamFactor_ .* ( ...
                        sum( coefMat .* subParamDMat ) ...
                    );
                end
            end
            
            drvIndex = 2;
            if( drvN >= drvIndex )
                if( ~singleMode || drvIndex == drvN )
                    for k = 1 : obj.contourN_
                        subParamDMat( k, : ) = subParamMat{ k, drvIndex };
                    end
                    coefMat = repmat( obj.coefListNormalT_, 1, polParamN );
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = obj.newParamFactor_ .* ( ...
                        sum( coefMat .* subParamDMat ) ...
                    );
                end
            end
            
            drvIndex = 3;
            if( drvN >= drvIndex )
                if( ~singleMode || drvIndex == drvN )
                    for k = 1 : obj.contourN_
                        subParamDMat( k, : ) = subParamMat{ k, drvIndex };
                    end
                    coefMat = repmat( obj.coefListNormalT_, 1, polParamN );
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = obj.newParamFactor_ .* ( ...
                        sum( coefMat .* subParamDMat ) ...
                    );
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Renamed functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function varargout = SupParamD( obj, polParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.NewParamD_( polParam, drvN );
        end
        
        function varargout = ArcLengthDsup( obj, supParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.ArcLengthDnew( supParam, drvN );
        end
        
        function varargout = TanAngleDsup( obj, supParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.TanAngleDnew( supParam, drvN );
        end
        
        function varargout = PointFullDsup( obj, supParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.PointFullDnew( supParam, drvN );
        end
        
        
    end
    
end

