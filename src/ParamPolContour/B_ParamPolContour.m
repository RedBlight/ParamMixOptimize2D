classdef( Abstract = true ) B_ParamPolContour < handle
    
    properties( Abstract = true )
        funcWithD_;
    end
    
    properties
        drvN_ = 3;
        paramStart_ = -pi;
        paramEnd_ = pi;
        paramLength_ = 2 * pi;
        paramCenter_ = 0;
        splineN_ = 2^16;
        polParamSpline_;
        newParamSpline_;
        newParamFactor_;
    end
    
    methods( Abstract = true )
        newParamDifferential = NewParamDifferential_( obj, polParam );
        varargout = NewParamD_( obj, polParam, drvN );
    end
    
    methods( Access = protected )
        
        function SolveParamDifferential( obj )
            polParamDiff = obj.paramLength_ ./ ( obj.splineN_ - 1 );
            obj.polParamSpline_ = ( 0 : ( obj.splineN_ - 1 ) ) .* polParamDiff + obj.paramStart_;
            
            polParamStart = obj.polParamSpline_( 1 : ( end - 1 ) );
            polParamEnd = obj.polParamSpline_( 2 : end );
            polParamDiffM = polParamDiff ./ 4;
            polParamM1 = polParamStart + polParamDiffM;
            polParamM2 = polParamM1 + polParamDiffM;
            polParamM3 = polParamM2 + polParamDiffM;
            
            newParamDfrStart = obj.NewParamDifferential_( polParamStart );
            newParamDfrM1 = obj.NewParamDifferential_( polParamM1 );
            newParamDfrM2 = obj.NewParamDifferential_( polParamM2 );
            newParamDfrM3 = obj.NewParamDifferential_( polParamM3 );
            newParamDfrEnd = obj.NewParamDifferential_( polParamEnd );
            
            integralCells = ( polParamDiff ./ 90 ) .* ( ...
                7 .* newParamDfrStart + 32 .* newParamDfrM1 + 12 .* newParamDfrM2 + 32 .* newParamDfrM3 + 7 .* newParamDfrEnd ...
            );
            
            cellN = length( integralCells );
            obj.newParamSpline_ = zeros( 1, cellN + 1 );
            cError = 0;
            for k = 1 : cellN
                Y = integralCells( k ) - cError;
                T = obj.newParamSpline_( k ) + Y;
                cError = T - obj.newParamSpline_( k );
                cError = cError - Y;
                obj.newParamSpline_( k + 1 ) = T;
            end

            obj.newParamFactor_ = obj.paramLength_ / obj.newParamSpline_( end );
            obj.newParamSpline_ = obj.newParamSpline_ .* obj.newParamFactor_ + obj.paramStart_;
        end
        
    end
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% New param as a function of polar param
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function newParam = NewParam( obj, polParam )
            newParam = spline( obj.polParamSpline_, obj.newParamSpline_, polParam );
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Polar param as a function of new param
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function polParam = PolParam( obj, newParam )
            polParam = spline( obj.newParamSpline_, obj.polParamSpline_, newParam );
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Pol param D/ new param
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function varargout = PolParamD( obj, polParam, drvN )
            singleMode = nargout == 1;
            indexVec = 1 : drvN;
            [ newParamD{ indexVec } ] = obj.NewParamD_( polParam, drvN );
            
            drvIndex = 1;
            if( drvN >= drvIndex )
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = ...
                        1 ./ newParamD{ 1 };
                end
            end
            
            drvIndex = 2;
            if( drvN >= drvIndex )
                newParamD1p3 = newParamD{ 1 } .* newParamD{ 1 } .* newParamD{ 1 };
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = ...
                        - newParamD{ 2 } ./ newParamD1p3;
                end
            end
            
            drvIndex = 3;
            if( drvN >= drvIndex )
                newParamD1p5 = newParamD1p3 .* newParamD{ 1 } .* newParamD{ 1 };
                newParamD2p2 = newParamD{ 2 } .* newParamD{ 2 };
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = ...
                        ( 3 .* newParamD2p2 - newParamD{ 3 } .* newParamD{ 1 } ) ./ newParamD1p5;
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Arc length D/ pol param
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function varargout = ArcLengthDpol( obj, polParam, drvN )
            singleMode = nargout == 1;
            indexVec = 1 : ( drvN + 1 );
            [ RD{ indexVec } ] = obj.funcWithD_.RadiusFullD( polParam, drvN );
            
            drvIndex = 1;
            if( drvN >= drvIndex )
                RD0p2 = RD{ 1 } .* RD{ 1 };
                RD1p2 = RD{ 2 } .* RD{ 2 };
                APD = sqrt( RD0p2 + RD1p2 );
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = ...
                        APD;
                end
            end
            
            drvIndex = 2;
            if( drvN >= drvIndex )
                RD0aRD2 = RD{ 1 } + RD{ 3 };
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = ...
                        RD{ 2 } .* RD0aRD2 ./ APD;
                end
            end
            
            drvIndex = 3;
            if( drvN >= drvIndex )
                RD1p3 = RD1p2 .* RD{ 2 };
                RD1p4 = RD1p3 .* RD{ 2 };
                APDp3 = APD .* APD .* APD;
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = ...
                        ( RD1p4 - RD{ 3 } .* ( RD{ 1 } .* RD1p2 - RD0p2 .* RD0aRD2 ) + ( RD0p2 .* RD{ 2 } + RD1p3 ) .* RD{ 4 } ) ./ APDp3;
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Arc length D/ new param
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function varargout = ArcLengthDnew( obj, newParam, drvN )
            singleMode = nargout == 1;
            polParam = obj.PolParam( newParam );
            indexVec = 1 : drvN;
            [ arcLengthDpol{ indexVec } ] = obj.ArcLengthDpol( polParam, drvN );
            [ polParamD{ indexVec } ] = obj.PolParamD( polParam, drvN );
            
            drvIndex = 1;
            if( drvN >= drvIndex )
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = ...
                        arcLengthDpol{ 1 } .* polParamD{ 1 };
                end
            end
            
            drvIndex = 2;
            if( drvN >= drvIndex )
                polParamD1p2 = polParamD{ 1 } .* polParamD{ 1 };
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = ...
                        arcLengthDpol{ 2 } .* polParamD1p2 + arcLengthDpol{ 1 } .* polParamD{ 2 };
                end
            end
            
            drvIndex = 3;
            if( drvN >= drvIndex )
                polParamD1p3 = polParamD1p2 .* polParamD{ 1 };
                polParamD1tD2t3 = 3 .* polParamD{ 1 } .* polParamD{ 2 };
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = ...
                        arcLengthDpol{ 3 } .* polParamD1p3 + arcLengthDpol{ 2 } .* polParamD1tD2t3 + arcLengthDpol{ 1 } .* polParamD{ 3 };
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Tangent angle D/ pol param
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function varargout = TanAngleDpol( obj, polParam, drvN )
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
                DD1 = CD1 .* CD1;
                
                DD1r2 = sqrt( DD1 );
                
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = ...
                        DD1r2;
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
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = ...
                        DD2 ./ ( 2 .* DD1r2 );
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
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = ...
                        ( DD3 ./ DD1r2 - DD2 .* DD2 ./ ( 2 .* DD1r2 .* DD1r2 .* DD1r2 ) ) ./ 2;
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Tangent angle D/ new param
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function varargout = TanAngleDnew( obj, newParam, drvN )
            singleMode = nargout == 1;
            polParam = obj.PolParam( newParam );
            indexVec = 1 : drvN;
            [ tanAngleDpol{ indexVec } ] = obj.TanAngleDpol( polParam, drvN );
            [ polParamD{ indexVec } ] = obj.PolParamD( polParam, drvN );
            
            drvIndex = 1;
            if( drvN >= drvIndex )
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = ...
                        tanAngleDpol{ 1 } .* polParamD{ 1 };
                end
            end
            
            drvIndex = 2;
            if( drvN >= drvIndex )
                polParamD1p2 = polParamD{ 1 } .* polParamD{ 1 };
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = ...
                        tanAngleDpol{ 2 } .* polParamD1p2 + tanAngleDpol{ 1 } .* polParamD{ 2 };
                end
            end
            
            drvIndex = 3;
            if( drvN >= drvIndex )
                polParamD1p3 = polParamD1p2 .* polParamD{ 1 };
                polParamD1tD2t3 = 3 .* polParamD{ 1 } .* polParamD{ 2 };
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = ...
                        tanAngleDpol{ 3 } .* polParamD1p3 + tanAngleDpol{ 2 } .* polParamD1tD2t3 + tanAngleDpol{ 1 } .* polParamD{ 3 };
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Curvature D/ pol param
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function varargout = CurvatureDpol( obj, polParam, drvN )
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
            
                AD1 = RD0p2 + RD1p2 - RD0tRD2;
                BD1 = ALp3;
                AD1rBD1 = AD1 ./ BD1;
                CD1 = AD1rBD1 .* AD1rBD1;
                sqCD1 = sqrt( CD1 );
                
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = ...
                        sqCD1;
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
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = ...
                        CD2 ./ ( 2 .* sqCD1 );
                end
            end
            
            drvIndex = 3;
            if( drvN >= drvIndex )
                RD2p2 = RD{ 3 } .* RD{ 3 };
                RD0tRD4 = RD{ 1 } .* RD{ 5 };
                RD1tRD3 = RD{ 2 } .* RD{ 4 };
                
                AD3 = 2 .* ( RD1p2 .* RD0tRD2 ) + 3 .* RD2p2 - ( RD1tRD3 + RD0tRD4 );
                BD3 = 3 .* ( K1 .* K1 ./ AL + AL .* ( RD1p2 + RD0tRD2 + RD2p2 + RD1tRD3 ) );
                AD3BD1mAD1BD3 = AD3 .* BD1 - AD1 .* BD3;
                CD3 = 2 .* ( K2 .* K2 + AD1rBD1 .* ( AD3BD1mAD1BD3 .* BD1 - 2 .* AD2BD1mAD1BD2 .* BD2 ) ./ ( BD1p2 .* BD1 ) );
                
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = ...
                        ( 1 / 2 ) .* ( CD3 ./ sqCD1 - CD2 ./ ( 2 .* sqCD1 .* sqCD1 .* sqCD1 ) );
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Curvature D/ new param
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function varargout = CurvatureDnew( obj, newParam, drvN )
            singleMode = nargout == 1;
            polParam = obj.PolParam( newParam );
            indexVec = 1 : drvN;
            [ curvatureDpol{ indexVec } ] = obj.CurvatureDpol( polParam, drvN );
            [ polParamD{ indexVec } ] = obj.PolParamD( polParam, drvN );
            
            drvIndex = 1;
            if( drvN >= drvIndex )
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = ...
                        curvatureDpol{ 1 } .* polParamD{ 1 };
                end
            end
            
            drvIndex = 2;
            if( drvN >= drvIndex )
                polParamD1p2 = polParamD{ 1 } .* polParamD{ 1 };
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = ...
                        curvatureDpol{ 2 } .* polParamD1p2 + curvatureDpol{ 1 } .* polParamD{ 2 };
                end
            end
            
            drvIndex = 3;
            if( drvN >= drvIndex )
                polParamD1p3 = polParamD1p2 .* polParamD{ 1 };
                polParamD1tD2t3 = 3 .* polParamD{ 1 } .* polParamD{ 2 };
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * drvIndex } = ...
                        curvatureDpol{ 3 } .* polParamD1p3 + curvatureDpol{ 2 } .* polParamD1tD2t3 + curvatureDpol{ 1 } .* polParamD{ 3 };
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Cartesian points D/ pol param
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function varargout = PointFullDpol( obj, polParam, drvN )
            singleMode = nargout == 1;
            indexVec = 1 : ( drvN + 1 );
            [ RD{ indexVec } ] = obj.funcWithD_.RadiusFullD( polParam, drvN );
            
            drvIndex = 0;
            if( drvN >= drvIndex )
                center = obj.funcWithD_.center_;
                ct = cos( polParam );
                st = sin( polParam );
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * ( drvIndex + 1 ) } = [ ...
                        RD{ 1 } .* ct + center( 1 ); ...
                        RD{ 1 } .* st + center( 2 ) ...
                    ];
                end
            end
            
            drvIndex = 1;
            if( drvN >= drvIndex )
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * ( drvIndex + 1 ) } = [ ...
                        RD{ 2 } .* ct - RD{ 1 } .* st; ...
                        RD{ 2 } .* st + RD{ 1 } .* ct ...
                    ];
                end
            end
            
            drvIndex = 2;
            if( drvN >= drvIndex )
                rD2x = RD{ 3 } - RD{ 1 };
                rD2y = 2 .* RD{ 2 };
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * ( drvIndex + 1 ) } = [ ...
                        rD2x .* ct - rD2y .* st; ...
                        rD2x .* st + rD2y .* ct ...
                    ];
                end
            end
            
            drvIndex = 3;
            if( drvN >= drvIndex )
                rD3x = RD{ 4 } - 3.* RD{ 2 };
                rD3y = 3 .* RD{ 3 } - RD{ 1 };
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * ( drvIndex + 1 ) } = [ ...
                        rD3x .* ct - rD3y .* st; ...
                        rD3x .* st + rD3y .* ct ...
                    ];
                end
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Cartesian points D/ new param
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function varargout = PointFullDnew( obj, newParam, drvN )
            singleMode = nargout == 1;
            polParam = obj.PolParam( newParam );
            indexVec = 1 : ( drvN + 1 );
            [ pointFullDpol{ indexVec } ] = obj.PointFullDpol( polParam, drvN );
            [ polParamD{ 1 : drvN } ] = obj.PolParamD( polParam, drvN );
            
            drvIndex = 0;
            if( drvN >= drvIndex )
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * ( drvIndex + 1 ) } = ...
                        pointFullDpol{ 1 };
                end
            end
            
            drvIndex = 1;
            if( drvN >= drvIndex )
                polParamD1 = repmat( polParamD{ 1 }, 2, 1 );
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * ( drvIndex + 1 ) } = ...
                        pointFullDpol{ 2 } .* polParamD1;
                end
            end
            
            drvIndex = 2;
            if( drvN >= drvIndex )
                polParamD1p2 = polParamD1 .* polParamD1;
                polParamD2 = repmat( polParamD{ 2 }, 2, 1 );
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * ( drvIndex + 1 ) } = ...
                        pointFullDpol{ 3 } .* polParamD1p2 + pointFullDpol{ 2 } .* polParamD2;
                end
            end
            
            drvIndex = 3;
            if( drvN >= drvIndex )
                polParamD1p3 = polParamD1p2 .* polParamD1;
                polParamD1tD2t3 = 3 .* polParamD1 .* polParamD2;
                polParamD3 = repmat( polParamD{ 3 }, 2, 1 );
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * ( drvIndex + 1 ) } = ...
                        pointFullDpol{ 4 } .* polParamD1p3 + pointFullDpol{ 3 } .* polParamD1tD2t3 + pointFullDpol{ 2 } .* polParamD3;
                end
            end
            
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Radius D/ pol param
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function varargout = RadiusFullDpol( obj, polParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.funcWithD_.RadiusFullD( polParam, drvN );
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Radius D/ new param
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function varargout = RadiusFullDnew( obj, newParam, drvN )
            singleMode = nargout == 1;
            polParam = obj.PolParam( newParam );
            indexVec = 1 : ( drvN + 1 );
            [ RD{ indexVec } ] = obj.funcWithD_.RadiusFullD( polParam, drvN );
            [ polParamD{ 1 : drvN } ] = obj.PolParamD( polParam, drvN );
            
            drvIndex = 0;
            if( drvN >= drvIndex )
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * ( drvIndex + 1 ) } = ...
                        RD{ 1 };
                end
            end
            
            drvIndex = 1;
            if( drvN >= drvIndex )
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * ( drvIndex + 1 ) } = ...
                        RD{ 2 } .* polParamD{ 1 };
                end
            end
            
            drvIndex = 2;
            if( drvN >= drvIndex )
                polParamD1p2 = polParamD{ 1 } .* polParamD{ 1 };
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * ( drvIndex + 1 ) } = ...
                        RD{ 3 } .* polParamD1p2 + RD{ 2 } .* polParamD{ 2 };
                end
            end
            
            drvIndex = 3;
            if( drvN >= drvIndex )
                polParamD1p3 = polParamD1p2 .* polParamD{ 1 };
                polParamD1tD2t3 = 3 .* polParamD{ 1 } .* polParamD{ 2 };
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * ( drvIndex + 1 ) } = ...
                        RD{ 4 } .* polParamD1p3 + RD{ 3 } .* polParamD1tD2t3 + RD{ 2 } .* polParamD{ 3 };
                end
            end
            
        end
        
    end
    
end

