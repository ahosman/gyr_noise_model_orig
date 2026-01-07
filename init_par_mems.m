%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Gyro MEMS  - Generic
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Updated 02/20/2024

% Channel definition:
% Ra = X
% Rs = Y
% Z = Z

% Without Interchange
% Mapping
% PAD          |  MEMS         | ASIC              | ASIC            
%              |               | (before mapping)  | (after mapping) 
% -------------------------------------------------------------------
% CGP1|CGN1    |  C (MEMS y)   | ch3               | ch2             
% CGP2|CGN2    |  S (MEMS x)   | ch1               | ch3             
% CGP3|CGN3    |  Z (MEMS z)   | ch2               | ch1
%
% gyr_axis_ex = "zxy"

function par = init_par_mems(parTbl)

    car = @(x) x(1);
    cadr = @(x) x(2);
    ratio = 1.0;

    %% MEMS Parameters
    par.AP_CGM   = sum(parTbl.AP(find(string(parTbl.PAD)=="CGM")))*1E-15;
    par.AN_CGM   = sum(parTbl.AN(find(string(parTbl.PAD)=="CGM")))*1E-15;

    par.AP_SUB   = parTbl.AP(find(string(parTbl.PAD)=="SUBL"))*1E-15 + ...
        parTbl.AP(find(string(parTbl.PAD)=="SUB"))*1E-15;
    par.AN_SUB   = parTbl.AN(find(string(parTbl.PAD)=="SUBL"))*1E-15 + ...
        parTbl.AN(find(string(parTbl.PAD)=="SUB"))*1E-15;

    par.DP_CGM   = sum(parTbl.DP(find(string(parTbl.PAD)=="CGM")))*1E-15;
    par.DN_CGM   = sum(parTbl.DN(find(string(parTbl.PAD)=="CGM")))*1E-15;

    par.DP_SUB   = parTbl.DP(find(string(parTbl.PAD)=="SUBL"))*1E-15 + ...
        parTbl.DP(find(string(parTbl.PAD)=="SUB"))*1E-15;
    par.DN_SUB   = parTbl.DN(find(string(parTbl.PAD)=="SUBL"))*1E-15 + ...
        parTbl.DN(find(string(parTbl.PAD)=="SUB"))*1E-15;


    par.CGP1_SUB = parTbl.CGP1(find(string(parTbl.PAD)=="SUBL"))*1E-15 ...
        + parTbl.CGP1(find(string(parTbl.PAD)=="SUB"))*1E-15 ...
        + sum(parTbl.CGP1(find(string(parTbl.PAD)=="QP")))*1E-15 ...
        + sum(parTbl.CGP1(find(string(parTbl.PAD)=="QN")))*1E-15;
    par.CGN1_SUB = parTbl.CGN1(find(string(parTbl.PAD)=="SUBL"))*1E-15 ...
        + parTbl.CGN1(find(string(parTbl.PAD)=="SUB"))*1E-15 ...
        + sum(parTbl.CGN1(find(string(parTbl.PAD)=="QP")))*1E-15 ...
        + sum(parTbl.CGN1(find(string(parTbl.PAD)=="QN")))*1E-15;
    par.CGP1_CGM = sum(parTbl.CGP1(find(string(parTbl.PAD)=="CGM")))*1E-15;
    par.CGN1_CGM = sum(parTbl.CGN1(find(string(parTbl.PAD)=="CGM")))*1E-15;
    par.CGP1_ALL = parTbl.CGP1(find(string(parTbl.PAD)=="cap sum with wires"))*1E-15;
    par.CGN1_ALL = parTbl.CGN1(find(string(parTbl.PAD)=="cap sum with wires"))*1E-15;



    par.CGP2_SUB = parTbl.CGP2(find(string(parTbl.PAD)=="SUBL"))*1E-15 ...
        + parTbl.CGP2(find(string(parTbl.PAD)=="SUB"))*1E-15 ...
        + sum(parTbl.CGP2(find(string(parTbl.PAD)=="QP")))*1E-15 ...
        + sum(parTbl.CGP2(find(string(parTbl.PAD)=="QN")))*1E-15;
    par.CGN2_SUB = parTbl.CGN2(find(string(parTbl.PAD)=="SUBL"))*1E-15 ...
        + parTbl.CGN2(find(string(parTbl.PAD)=="SUB"))*1E-15 ...
        + sum(parTbl.CGN2(find(string(parTbl.PAD)=="QP")))*1E-15 ...
        + sum(parTbl.CGN2(find(string(parTbl.PAD)=="QN")))*1E-15;
    par.CGP2_CGM = sum(parTbl.CGP2(find(string(parTbl.PAD)=="CGM")))*1E-15;
    par.CGN2_CGM = sum(parTbl.CGN2(find(string(parTbl.PAD)=="CGM")))*1E-15;
    par.CGP2_ALL = parTbl.CGP2(find(string(parTbl.PAD)=="cap sum with wires"))*1E-15;
    par.CGN2_ALL = parTbl.CGN2(find(string(parTbl.PAD)=="cap sum with wires"))*1E-15;

    par.CGP3_SUB = parTbl.CGP3(find(string(parTbl.PAD)=="SUBL"))*1E-15*ratio ...
        + parTbl.CGP3(find(string(parTbl.PAD)=="SUB"))*1E-15 ...
        + sum(parTbl.CGP3(find(string(parTbl.PAD)=="QP")))*1E-15 ...
        + sum(parTbl.CGP3(find(string(parTbl.PAD)=="QN")))*1E-15;
    par.CGN3_SUB = parTbl.CGN3(find(string(parTbl.PAD)=="SUBL"))*1E-15*ratio ...
        + parTbl.CGN3(find(string(parTbl.PAD)=="SUB"))*1E-15 ...
        + sum(parTbl.CGN3(find(string(parTbl.PAD)=="QP")))*1E-15 ...
        + sum(parTbl.CGN3(find(string(parTbl.PAD)=="QN")))*1E-15;
    par.CGP3_CGM = sum(parTbl.CGP3(find(string(parTbl.PAD)=="CGM")))*1E-15;
    par.CGN3_CGM = sum(parTbl.CGN3(find(string(parTbl.PAD)=="CGM")))*1E-15;
    par.CGP3_ALL = parTbl.CGP3(find(string(parTbl.PAD)=="cap sum with wires"))*1E-15 ...
        - (1-ratio)*parTbl.CGP3(find(string(parTbl.PAD)=="SUBL"))*1E-15;
    par.CGN3_ALL = parTbl.CGN3(find(string(parTbl.PAD)=="cap sum with wires"))*1E-15 ...
        - (1-ratio)*parTbl.CGN3(find(string(parTbl.PAD)=="SUBL"))*1E-15;

    par.CGP1_CGN1 = par.CGP1_CGM-par.CGN1_CGM;
    par.CGP2_CGN2 = par.CGP2_CGM-par.CGN2_CGM;
    par.CGP3_CGN3 = par.CGP3_CGM-par.CGN3_CGM;

    par.Cpar_DP  = parTbl.DP(find(string(parTbl.PAD)=="cap sum with wires"))*1E-15;
    par.Cpar_DN  = parTbl.DN(find(string(parTbl.PAD)=="cap sum with wires"))*1E-15;
    par.Rs_DP     = 11.886899e+03;
    par.Rs_DN     = 11.991468e+03;
                   
    % MEMS series resistors on CP and CN [Ohm]
    par.Rs_GP1    = 1.9869815e+03;
    par.Rs_GN1    = 2.1658862e+03;

    par.Rs_GP2    = 3.7589295e+03;
    par.Rs_GN2    = 3.7392400e+03;

    par.Rs_GP3    = 6.2928694e+03;
    par.Rs_GN3    = 6.3344795e+03;
 end
