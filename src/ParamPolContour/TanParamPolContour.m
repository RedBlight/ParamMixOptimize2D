classdef TanParamPolContour < B_ParamPolContour
    
    properties
        funcWithD_
    end
    
    methods
        
        % ( funcWithD )
        % ( funcWithD, splineN )
        function obj = TanParamPolContour( varargin )
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
            newParamDifferential = obj.TanAngleDpol( polParam, 1 );
        end
        
        function varargout = NewParamD_( obj, polParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.TanAngleDpol( polParam, drvN );
            for k = 1 : nargout
                varargout{ k } =  varargout{ k } .* obj.newParamFactor_;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Renamed functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function varargout = TanParamD( obj, polParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.NewParamD_( polParam, drvN );
        end
        
        function varargout = ArcLengthDtan( obj, tanParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.ArcLengthDnew( tanParam, drvN );
        end
        
        function varargout = TanAngleDtan( obj, tanParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.TanAngleDnew( tanParam, drvN );
        end
        
        function varargout = PointFullDtan( obj, tanParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.PointFullDnew( tanParam, drvN );
        end
        
        
    end
    
end

