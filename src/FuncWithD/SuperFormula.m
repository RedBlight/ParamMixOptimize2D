classdef SuperFormula < B_FuncWithD
    
    properties
        drvN_;
        
        a_;
        b_;
        m_;
        n1_;
        n2_;
        n3_;
        alpha_
        beta_;
        center_;
        rotation_;
    end
    
    methods
        
        function obj = SuperFormula( ...
            a, b, m, n1, n2, n3, alpha, beta, center, rotation ...
        )
            % Base class parameter
            obj.drvN_ = 4;
            
            % SuperFormula parameters
            obj.a_ = a;
            obj.b_ = b;
            obj.m_ = m;
            obj.n1_ = n1;
            obj.n2_ = n2;
            obj.n3_ = n3;
            obj.alpha_ = alpha;
            obj.beta_ = beta;
            obj.center_ = center;
            obj.rotation_ = rotation;
        end
        
        function varargout = RadiusFullD( obj, polParam, drvN )
            [ varargout{ 1 : nargout } ] = obj.FuncValD( polParam, drvN );
        end
        
    end
    
    methods( Access = protected )
        
        function varargout = FuncValD_( obj, param, drvN )
            polParam = param + obj.rotation_;
            m = obj.m_;
            n1 = obj.n1_;
            n2 = obj.n2_;
            n3 = obj.n3_;
            alpha = obj.alpha_;
            beta = obj.beta_;
            a = obj.a_;
            b = obj.b_;
            
            singleMode = nargout == 1;
            
            drvIndex = 0;
            if( drvN >= drvIndex )
                mq1 = m / 4;
                alphap2 = alpha * alpha;
                betap2 = beta * beta;
                tq1 = polParam .* mq1;
                sq1 = sin( tq1 );
                cq1 = cos( tq1 );
                sq1p2 = sq1 .* sq1;
                cq1p2 = cq1 .* cq1;

                % E
                EV1 = betap2 + sq1p2;
                EV2 = sqrt( EV1 );
                ED0 = EV2 ./ b;
                
                % D
                DV1 = alphap2 + cq1p2;
                DV2 = sqrt( DV1 );
                DD0 = DV2 ./ a;

                % C
                ED0n0 = ED0 .^ n3;
                CD0 = ED0n0;

                % B
                DD0n0 = DD0 .^ n2;
                BD0 = DD0n0;

                % A
                AD0 = BD0 + CD0;

                % R

                n1i = 1 / n1;
                AD0n0 = AD0 .^ ( - n1i );
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * ( drvIndex + 1 ) } = ...
                        AD0n0;
                end
            end
            
            drvIndex = 1;
            if( drvN >= drvIndex )
                mq2 = mq1 * 2;
                tq2 = polParam .* mq2;
                sq2 = sin( tq2 );

                % E
                EK1 = m / ( 8 * b );
                ED1 = EK1 .* sq2 ./ EV2;

                % D
                DK1 = - m / ( 8 * a );
                DD1 = DK1 .* sq2 ./ DV2;

                % C
                ED0n1 = ED0n0 ./ ED0;
                CD1 = n3 .* ED0n1 .* ED1;

                % B
                DD0n1 = DD0n0 ./ DD0;
                BD1 = n2 .* DD0n1 .* DD1;

                % A
                AD1 = BD1 + CD1;

                % R
                AD0n1 = AD0n0 ./ AD0;
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * ( drvIndex + 1 ) } = ...
                        - n1i .* AD0n1 .* AD1;
                end
            end
            
            drvIndex = 2;
            if( drvN >= drvIndex )
                mq4 = mq1 * 4;
                mp2 = m * m;
                tq4 = polParam .* mq4;
                cq2 = cos( tq2 );
                cq4 = cos( tq4 );

                % E
                EK2 = - mp2 / ( 128 * b );
                EK7 = - 4 * ( 2 * betap2 + 1 );
                EV3 = EV1 .* EV2;
                ED2 = EK2 .* ( EK7 .* cq2 + cq4 + 3 ) ./ EV3;

                % D
                DK2 = - mp2 / ( 128 * a );
                DK7 = 4 * ( 2 * alphap2 + 1 );
                DV3 = DV1 .* DV2;
                DD2 = DK2 .* ( DK7 .* cq2 + cq4 + 3 ) ./ DV3;

                % C
                ED0n2 = ED0n1 ./ ED0;
                ED1p2 = ED1 .* ED1;
                CK1 = n3 - 1;
                CV1 = ED0 .* ED2;
                CD2 = n3 .* ED0n2 .* ( CV1 + CK1 .* ED1p2 );

                % B
                DD0n2 = DD0n1 ./ DD0;
                DD1p2 = DD1 .* DD1;
                BK1 = n2 - 1;
                BV1 = DD0 .* DD2;
                BD2 = n2 .* DD0n2 .* ( BV1 + BK1 .* DD1p2 );

                % A
                AD2 = BD2 + CD2;

                % R
                n1ip2 = n1i * n1i;
                AD0n2 = AD0n1 ./ AD0;
                AD1p2 = AD1 .* AD1;
                RK1 = n1 + 1;
                RV1 = AD0 .* AD2;
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * ( drvIndex + 1 ) } = ...
                        n1ip2 .* AD0n2 .* ( RK1 .* AD1p2 - n1 * RV1 );
                end
            end
            
            drvIndex = 3;
            if( drvN >= drvIndex )
                mp3 = mp2 * m;
                alphap4 = alphap2 * alphap2;
                betap4 = betap2 * betap2;

                % E
                EK3 = - mp3 / ( 1024 * b );
                EK5 = betap4 + betap2;
                EK6 = 32 * EK5;
                EV4 = EV1 .* EV3;
                ED3 = EK3 .* sq2 .* ( EK6 + EK7 .* cq2 + cq4 + 3 ) ./ EV4;

                % D
                DK3 = mp3 / ( 1024 * a );
                DK5 = alphap4 + alphap2;
                DK6 = 32 * DK5;
                DV4 = DV1 .* DV3;
                DD3 = DK3 .* sq2 .* ( DK6 + DK7 .* cq2 + cq4 + 3 ) ./ DV4;

                % C
                ED0p2 = ED0 .* ED0;
                ED0n3 = ED0n2 ./ ED0;
                CK2 = n3 - 2;
                CD3 = n3 .* ED0n3 .* ( ED3 .* ED0p2 + CK1 .* ED1 .* ...
                    ( 3 .* CV1 + CK2 .* ED1p2 ) );

                % B
                DD0p2 = DD0 .* DD0;
                DD0n3 = DD0n2 ./ DD0;
                BK2 = n2 - 2;
                BD3 = n2 .* DD0n3 .* ( DD3 .* DD0p2 + BK1 .* DD1 .* ...
                    ( 3 .* BV1 + BK2 .* DD1p2 ) );

                % A
                AD3 = BD3 + CD3;

                % R
                n1p2 = n1 * n1;
                n1ip3 = n1ip2 * n1i;
                AD0p2 = AD0 .* AD0;
                AD0n3 = AD0n2 ./ AD0;
                RK2 = 2 * n1 + 1;
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * ( drvIndex + 1 ) } = ...
                        - n1ip3 .* AD0n3 .* ( n1p2 .* AD0p2 .* AD3 + ...
                        RK1 .* AD1 .* ( RK2 .* AD1p2 - 3 .* n1 .* RV1 ) );
                end
            end
            
            drvIndex = 4;
            if( drvN >= drvIndex )
                mq6 = mq1 * 6;
                mq8 = mq1 * 8;
                mp4 = mp3 * m;
                alphap6 = alphap4 * alphap2;
                betap6 = betap4 * betap2;
                tq6 = polParam .* mq6;
                tq8 = polParam .* mq8;
                cq6 = cos( tq6 );
                cq8 = cos( tq8 );

                % E
                EK4 = mp4 / ( 2048 * 16 * b );
                EK8 = 7 * ( 64 * EK5 + 5 );
                EK9 = - 8 * ( 2 * betap2 + 1 );
                EKA = - 4 * ( 16 * EK5 - 7 );
                EKB = - 8 * ( 64 * betap6 + 96 * betap4 + 46 * betap2 + 7 );
                EV5 = EV1 .* EV4;
                ED4 = EK4 .* ( EK8 + EK9 * cq6 + EKA * cq4 + EKB * cq2 + cq8 ) ./ EV5;
                
                % D
                DK4 = mp4 / ( 2048 * 16 * a );
                DK8 = 7 * ( 64 * DK5 + 5 );
                DK9 = 8 * ( 2 * alphap2 + 1 );
                DKA = - 4 * ( 16 * DK5 - 7 );
                DKB = 8 * ( 64 * alphap6 + 96 * alphap4 + 46 * alphap2 + 7 );
                DV5 = DV1 .* DV4;
                DD4 = DK4 .* ( DK8 + DK9 * cq6 + DKA * cq4 + DKB * cq2 + cq8 ) ./ DV5;

                % C
                ED0p3 = ED0p2 .* ED0;
                ED0n4 = ED0n3 ./ ED0;
                ED1p4 = ED1p2 .* ED1p2;
                CK3 = n3 - 3;
                CD4 = n3 .* ED0n4 .* ( ED4 .* ED0p3 + CK1 .* ( ...
                    3 .* CV1 .* CV1 + CK2 .* CK3 .* ED1p4 + ...
                    4 .* ED0p2 .* ED3 .* ED1 + ...
                    6 .* CK2 .* CV1 .* ED1p2 ) );

                % B
                DD0p3 = DD0p2 .* DD0;
                DD0n4 = DD0n3 ./ DD0;
                DD1p4 = DD1p2 .* DD1p2;
                BK3 = n2 - 3;
                BD4 = n2 .* DD0n4 .* ( DD4 .* DD0p3 + BK1 .* ( ...
                    3 .* BV1 .* BV1 + BK2 .* BK3 .* DD1p4 + ...
                    4 .* DD0p2 .* DD3 .* DD1 + ...
                    6 .* BK2 .* BV1 .* DD1p2 ) );

                % A

                AD4 = BD4 + CD4;

                % R
                n1p3 = n1p2 * n1;
                n1ip4 = n1ip3 * n1i;
                AD0p3 = AD0p2 .* AD0;
                AD0n4 = AD0n3 ./ AD0;
                AD1p4 = AD1p2 .* AD1p2;
                RK3 = n1 * ( 6 * n1 + 5 ) + 1;
                if( ~singleMode || drvIndex == drvN )
                    varargout{ singleMode * 1 + ~singleMode * ( drvIndex + 1 ) } = ...
                        n1ip4 .* AD0n4 .* ( RK1 .* ( 3 .* n1p2 .* RV1 .* RV1 + ...
                        RK3 .* AD1p4 + 4 .* n1p2 .* AD0p2 .* AD3 .* AD1 - ...
                        6 .* n1 .* RK2 .* RV1 .* AD1p2 ) - n1p3 .* AD0p3 .* AD4 );
                end
            end
            
        end
        
    end
    
end