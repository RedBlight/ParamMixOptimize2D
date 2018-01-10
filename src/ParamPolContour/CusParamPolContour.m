classdef CusParamPolContour < B_ParamPolContour
    
    properties
        funcWithD_
        smoother_ = 1;
    end
    
    methods
        
        % ( funcWithD )
        % ( funcWithD, splineN )
        % ( funcWithD, splineN, smoother )
        function obj = CusParamPolContour( varargin )
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
            AL = sqrt( RD0p2 + RD1p2 );
            ALp3 = AL .* AL .* AL;

            AD1 = RD0p2 + 2 .* RD1p2 - RD0tRD2;
            BD1 = ALp3;
            AD1rBD1 = AD1 ./ BD1;
            CD1 = AD1rBD1 .* AD1rBD1 + obj.smoother_ .^ 2;
            
            newParamDifferential = sqrt( CD1 );
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
                AL = sqrt( RD0p2 + RD1p2 );
                ALp3 = AL .* AL .* AL;
            
                AD1 = RD0p2 + 2.* RD1p2 - RD0tRD2;
                BD1 = ALp3;
                AD1rBD1 = AD1 ./ BD1;
                CD1 = AD1rBD1 .* AD1rBD1 + obj.smoother_ .^ 2;
                sqCD1 = sqrt( CD1 );
                
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = obj.newParamFactor_ .* ( ...
                        sqCD1 ...
                    );
                end
            end
            
            drvIndex = 2;
            if( drvN >= drvIndex )
                RD0tRD1 = RD{ 1 } .* RD{ 2 };
                RD0tRD3 = RD{ 1 } .* RD{ 4 };
                RD1tRD2 = RD{ 2 } .* RD{ 3 };
                BD1p2 = BD1 .* BD1;
                K1 = ( RD0tRD1 + RD1tRD2 );
                
                AD2 = 2 .* RD0tRD1 + 3 .* RD1tRD2 - RD0tRD3;
                BD2 = 3 .* AL .* K1;
                AD2BD1mAD1BD2 = AD2 .* BD1 - AD1 .* BD2;
                K2 = AD2BD1mAD1BD2 ./ ( BD1p2 );
                CD2 = 2 .* AD1rBD1 .* ( K2 );
                
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = obj.newParamFactor_ .* ( ...
                        CD2 ./ ( 2 .* sqCD1 ) ...
                    );
                end
            end
            
            drvIndex = 3;
            if( drvN >= drvIndex )
                RD2p2 = RD{ 3 } .* RD{ 3 };
                RD0tRD4 = RD{ 1 } .* RD{ 5 };
                RD1tRD3 = RD{ 2 } .* RD{ 4 };
                
                AD3 = 2 .* ( RD1p2 + RD0tRD2 ) + 3 .* RD2p2 + 2 .* RD1tRD3 - RD0tRD4;
                BD3 = 3 .* ( K1 .* K1 ./ AL + AL .* ( RD1p2 + RD0tRD2 + RD2p2 + RD1tRD3 ) );
                AD3BD1mAD1BD3 = AD3 .* BD1 - AD1 .* BD3;
                CD3 = 2 .* ( K2 .* K2 + AD1rBD1 .* ( AD3BD1mAD1BD3 .* BD1 - 2 .* AD2BD1mAD1BD2 .* BD2 ) ./ ( BD1p2 .* BD1 ) );
                
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = obj.newParamFactor_ .* ( ...
                        ( 1 / 2 ) .* ( CD3 ./ sqCD1 - CD2 .* CD2 ./ ( 2 .* sqCD1 .* sqCD1 .* sqCD1 ) ) ...
                    );
                end
            end
%             if( drvN == 3 );
%                 paramD = circshift( polParam,[0 1]) - circshift( polParam,[0 -1]);
%                 AD1D = circshift( varargout{1},[0 1]) - circshift( varargout{1},[0 -1]);
%                 AD2N = AD1D ./ paramD;
%                 AD2D = circshift( varargout{2},[0 1]) - circshift( varargout{2},[0 -1]);
%                 AD3N = AD2D ./ paramD;
%                 figure(601);
%                 plot( polParam, varargout{2}, 'blue', polParam, AD2N, 'red', 'LineSmoothing', 'on' );
%                 figure(602);
%                 plot( polParam, varargout{3}, 'blue', polParam, AD3N, 'red', 'LineSmoothing', 'on' );
%             end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Renamed functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function varargout = CusParamD( obj, polParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.NewParamD_( polParam, drvN );
        end
        
        function varargout = ArcLengthDcus( obj, cusParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.ArcLengthDnew( cusParam, drvN );
        end
        
        function varargout = TanAngleDcus( obj, cusParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.TanAngleDnew( cusParam, drvN );
        end
        
        function varargout = PointFullDcus( obj, cusParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.PointFullDnew( cusParam, drvN );
        end
        
        
    end
    
end

