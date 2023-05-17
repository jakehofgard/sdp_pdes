% Returns solver parameters as a Yalmip sdpsettings structure
% Parameter selected only for MOSEK and SeDuMi, otherwise using Yalmip's
% default ones

function options = getSolverParams( SDPsolver, sosmodel )

if(strcmpi(SDPsolver,'mosek'))
    if(~exist('mosekopt.m'))
        disp('***MOSEK not found, leaving the solver selection to Yalmip***')
        SDPsolver = '';
    end
end

if(strcmpi(SDPsolver,'sedumi'))
    if(~exist('sedumi.m'))
        disp('***SeDuMi not found, leaving the solver selection to Yalmip***')
        SDPsolver = '';
    end
end

switch SDPsolver
    case 'mosek'
        disp('SDP Solver: MOSEK')
       options = sdpsettings('solver','mosek-sdp','sos.model',sosmodel,...
            'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',1e-12,...
            'mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS',1e-12,...
            'mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP', 1e-12,...
            'mosek.MSK_DPAR_INTPNT_TOL_PFEAS', 1e-12,...
            'mosek.MSK_DPAR_INTPNT_TOL_DFEAS', 1e-12,...
            'mosek.MSK_DPAR_INTPNT_TOL_REL_GAP', 1e-12,...
            'mosek.MSK_DPAR_INTPNT_CO_TOL_INFEAS', 1e-12,...
            'mosek.MSK_DPAR_INTPNT_CO_TOL_MU_RED', 1e-12,...
            'mosek.MSK_DPAR_INTPNT_CO_TOL_NEAR_REL', 1.0,...
            'mosek.MSK_IPAR_INTPNT_MAX_ITERATIONS', 100,...
            'mosek.MSK_DPAR_INTPNT_TOL_STEP_SIZE', 1e-18);
    case 'sedumi'
        disp('SDP Solver: SeDuMi')
        options = sdpsettings('solver','sedumi','sedumi.maxiter',100,'sedumi.eps',1e-12,'sos.model',sosmodel);
    case ''
        options = sdpsettings('sos.model',sosmodel);
    otherwise
        options = sdpsettings('sos.model',sosmodel);
end


