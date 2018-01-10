classdef ArcParamPolContour < B_ParamPolContour
    
    properties
        funcWithD_
    end
    
    methods
        
        % ( funcWithD )
        % ( funcWithD, splineN )
        function obj = ArcParamPolContour( varargin )
            obj.funcWithD_ = varargin{ 1 };
            if nargin > 1
                obj.splineN_ = varargin{ 2 };
            end
            obj.SolveParamDifferential();
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Abstract functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function newParamDifferential = NewParamDifferential_( obj, polParam )
            newParamDifferential = obj.ArcLengthDpol( polParam, 1 );
        end
        
        function varargout = NewParamD_( obj, polParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.ArcLengthDpol( polParam, drvN );
            for k = 1 : nargout
                varargout{ k } =  varargout{ k } .* obj.newParamFactor_;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Renamed functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function varargout = ArcParamD( obj, polParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.NewParamD_( polParam, drvN );
        end
        
        function varargout = ArcLengthDarc( obj, arcParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.ArcLengthDnew( arcParam, drvN );
        end
        
        function varargout = TanAngleDarc( obj, arcParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.TanAngleDnew( arcParam, drvN );
        end
        
        function varargout = PointFullDarc( obj, arcParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.PointFullDnew( arcParam, drvN );
        end
        
        
    end
    
end

