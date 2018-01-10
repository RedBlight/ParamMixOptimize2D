classdef CurParamPolContour < B_ParamPolContour
    
    properties
        funcWithD_
    end
    
    methods
        
        % ( funcWithD )
        % ( funcWithD, splineN )
        function obj = CurParamPolContour( varargin )
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
            newParamDifferential = obj.CurvatureDpol( polParam, 1 );
        end
        
        function varargout = NewParamD_( obj, polParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.CurvatureDpol( polParam, drvN );
            for k = 1 : nargout
                varargout{ k } =  varargout{ k } .* obj.newParamFactor_;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Renamed functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function varargout = CurParamD( obj, polParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.NewParamD_( polParam, drvN );
        end
        
        function varargout = ArcLengthDcur( obj, curParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.ArcLengthDnew( curParam, drvN );
        end
        
        function varargout = TanAngleDcur( obj, curParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.TanAngleDnew( curParam, drvN );
        end
        
        function varargout = PointFullDcur( obj, curParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.PointFullDnew( curParam, drvN );
        end
        
        
    end
    
end

