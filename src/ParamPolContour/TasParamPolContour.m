classdef TasParamPolContour < B_ParamPolContour
    
    properties
        funcWithD_
        smoother_ = 1;
    end
    
    methods
        
        % ( funcWithD )
        % ( funcWithD, splineN )
        % ( funcWithD, splineN, smoother )
        function obj = TasParamPolContour( varargin )
            obj.funcWithD_ = varargin{ 1 };
            if nargin > 1
                obj.splineN_ = varargin{ 2 };
            end
            if nargin > 2
                obj.smoother_ = varargin{ 3 };
            end
            obj.SolveParamDifferential();
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Abstract functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function newParamDifferential = NewParamDifferential_( obj, polParam )
            [ RD{ 1 : 3 } ] = obj.funcWithD_.RadiusFullD( polParam, 2 );
            RD0p2 = RD{ 1 } .* RD{ 1 };
            RD1p2 = RD{ 2 } .* RD{ 2 };
            RD0tRD2 = RD{ 1 } .* RD{ 3 };

            AD1 = RD1p2 - RD0tRD2;
            BD1 = RD0p2 + RD1p2;
            CD1 = 1 + AD1 ./ BD1;
            DD1 = CD1 .* CD1 + obj.smoother_ .^ 2;

            DD1r2 = sqrt( DD1 );
            
            newParamDifferential = DD1r2;
        end
        
        function varargout = NewParamD_( obj, polParam, drvN )
            singleMode = nargout == 1;
            indexVec = 1 : ( drvN + 2 );
            [ RD{ indexVec } ] = obj.funcWithD_.RadiusFullD( polParam, drvN + 1 );
            
            drvIndex = 1;
            if( drvN >= drvIndex )
                RD0p2 = RD{ 1 } .* RD{ 1 };
                RD1p2 = RD{ 2 } .* RD{ 2 };
                RD0tRD2 = RD{ 1 } .* RD{ 3 };
            
                AD1 = RD1p2 - RD0tRD2;
                BD1 = RD0p2 + RD1p2;
                CD1 = 1 + AD1 ./ BD1;
                DD1 = CD1 .* CD1 + obj.smoother_ .^ 2;
                
                DD1r2 = sqrt( DD1 );
                
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = obj.newParamFactor_ .* ( ...
                        DD1r2 ...
                    );
                end
            end
            
            drvIndex = 2;
            if( drvN >= drvIndex )
                
                RD1tRD2 = RD{ 2 } .* RD{ 3 };
                BD1p2 = BD1 .* BD1;
                
                AD2 = RD1tRD2 - RD{ 1 } .* RD{ 4 };
                BD2 = 2 .* ( RD{ 1 } .* RD{ 2 } + RD1tRD2 );
                AD2tBD1mAD1tBD2 = AD2 .* BD1 - AD1 .* BD2;
                CD2 = ( AD2tBD1mAD1tBD2 ) ./ ( BD1p2 );
                DD2 = 2 .* CD1 .* CD2;
                
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = obj.newParamFactor_ .* ( ...
                        DD2 ./ ( 2 .* DD1r2 ) ...
                    );
                end
            end
            
            drvIndex = 3;
            if( drvN >= drvIndex )
                RD2p2 = RD{ 3 } .* RD{ 3 };
                
                AD3 = RD2p2 - RD{ 1 } .* RD{ 5 };
                BD3 = 2 .* ( RD1p2 + RD0tRD2 + RD2p2 + RD{ 2 } .* RD{ 4 } );
                CD3 = ( ( AD3 .* BD1 - AD1 .* BD3 ) .* BD1 - 2 .* AD2tBD1mAD1tBD2 .* BD2 ) ./ ( BD1p2 .* BD1 );
                DD3 = 2 .* ( CD2 .* CD2 + CD1 .* CD3 );

                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = obj.newParamFactor_ .* ( ...
                        ( DD3 ./ DD1r2 - DD2 .* DD2 ./ ( 2 .* DD1r2 .* DD1r2 .* DD1r2 ) ) ./ 2 ...
                    );
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Renamed functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function varargout = TasParamD( obj, polParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.NewParamD_( polParam, drvN );
        end
        
        function varargout = ArcLengthDtas( obj, tasParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.ArcLengthDnew( tasParam, drvN );
        end
        
        function varargout = TanAngleDtas( obj, tasParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.TanAngleDnew( tasParam, drvN );
        end
        
        function varargout = PointFullDtas( obj, tasParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.PointFullDnew( tasParam, drvN );
        end
        
        
    end
    
end

