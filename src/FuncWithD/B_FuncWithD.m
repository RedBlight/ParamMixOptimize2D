classdef( Abstract = true ) B_FuncWithD < handle
    
    properties
        lockState__;
        storedLockState__;
        
        % cached values
        param__;
        funcValD__;
    end
    
    properties( Abstract = true )
        drvN_;
    end
    
    methods( Abstract = true, Access = protected )
        varargout = FuncValD_( obj, param, drvN );
    end
    
    methods( Access = protected )
        
        function varargout = FuncValD( obj, param, drvN )
            if nargout == 1
                if obj.lockState__
                    varargout{ 1 } = obj.funcValD__{ drvN + 1 };
                else
                    varargout{ 1 } = obj.FuncValD_( param, drvN );
                end
            else
                indexVec = 1 : ( drvN + 1 );
                if obj.lockState__
                    varargout( indexVec ) = obj.funcValD__( indexVec );
                else
                    [ varargout{ indexVec } ] = obj.FuncValD_( param, drvN );
                end
            end
        end
        
    end
    
    methods
        
        function Cache( obj, param, drvN )
            obj.param__ = param;
            obj.funcValD__ = cell( 1, obj.drvN_ + 1 );
            [ obj.funcValD__{ 1 : ( drvN + 1 ) } ] = obj.FuncValD_( param, drvN );
        end
        
        function Lock( obj )
            obj.lockState__ = true;
        end
        
        function Unlock( obj )
            obj.lockState__ = false;
        end
        
        function StoreLockState( obj )
            obj.storedLockState__ = obj.lockState__;
        end
        
        function RestoreLockState( obj )
            obj.lockState__ = obj.storedLockState__;
        end
        
    end
    
end

