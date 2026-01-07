function [ Hsd_model_diff, Hsd_model_dig, ABCD_model_SD1, ABCD_model_SD2 ] = rate_sd(obj)

    % Generate the State space model and ABCD matrix for the Rate Sigma Delta
    coeff = obj.imu.as.gyr.config.inf.sense.rate.coeff;
    % SD parameters
    
    z = ss(zpk('z',1/obj.imu.as.gyr.config.inf.fs));
    integ = 1/(1-(1/z));    %Integrator block
    
    % Main 2nd order SD:
    H_a0 = ss(coeff.a0);
    H_a0.InputName = 'u';
    H_a0.OutputName = 'u1';
    
    H_fb1 = ss(-coeff.f1);
    H_fb1.InputName = 'v';
    H_fb1.OutputName = 'f1';
    
    Sum1 = sumblk('e1 = u1+f1');
    
    H_int1 = integ*(1/z);
    H_int1.InputName = 'e1';
    H_int1.OutputName = 'o1';
    
    H_a1 = ss(coeff.a1);
    H_a1.InputName = 'o1';
    H_a1.OutputName = 'u2';
    
    H_a12 = ss((-coeff.a1/2)*(1-(1/z)));
    H_a12.InputName = 'o1';
    H_a12.OutputName = 'ff1';
    
    H_fb2 = ss(-coeff.b1*(1/z));
    H_fb2.InputName = 'v';
    H_fb2.OutputName = 'f2';
    
    Sum2 = sumblk('e2 = u2+f2+ff1');
    
    H_int2 = integ;
    H_int2.InputName = 'e2';
    H_int2.OutputName = 'o2';
    
    H_sum = ss(coeff.gainQ);
    H_sum.InputName = 'o2';
    H_sum.OutputName = 'y';
    
    Hsd_model_SD1 = minreal(...
        connect(H_a0,H_fb1,Sum1,H_int1,H_a1,H_a12,H_fb2,Sum2,H_int2,H_sum,...
        {'u', 'v'},{'y'}));
    
    % Create ABCD matrix
    ABCD_model_SD1 = [Hsd_model_SD1.A Hsd_model_SD1.B; ...
        Hsd_model_SD1.C Hsd_model_SD1.D];
    
    
    % Difference between input and output of comparator of first SD
    
    
    H_c1 = ss(-coeff.c1);
    H_c1.InputName = 'uv';
    H_c1.OutputName = 'uc1';
    
    H_c2 = ss(coeff.c2/coeff.gainQ);
    H_c2.InputName = 'uy';
    H_c2.OutputName = 'uc2';
    
    Sum_diff = sumblk('u2 = uc1+uc2');
    
    Hsd_model_diff = minreal(connect(H_c1,H_c2,Sum_diff,{'uv', 'uy'},{'u2'}));
    
    % Qnoise 2nd order SD:
    
    H_fb3 = ss(-coeff.f2);
    H_fb3.InputName = 'v2';
    H_fb3.OutputName = 'f3';
    
    Sum3 = sumblk('e3 = u2+f3');
    
    H_int3 = integ*(1/z);
    H_int3.InputName = 'e3';
    H_int3.OutputName = 'o3';
    
    H_a2 = ss(coeff.a2);
    H_a2.InputName = 'o3';
    H_a2.OutputName = 'u4';
    
    H_fb4 = ss(-coeff.b2*(1/z));
    H_fb4.InputName = 'v2';
    H_fb4.OutputName = 'f4';
    
    Sum4 = sumblk('e4 = u4+f4');
    
    H_int4 = integ;
    H_int4.InputName = 'e4';
    H_int4.OutputName = 'o4';
    
    H_sum2 = ss(coeff.gainQ2);
    H_sum2.InputName = 'o4';
    H_sum2.OutputName = 'y2';
    
    Hsd_model_SD2 = minreal(connect(H_fb3,Sum3,H_int3,H_a2,H_fb4,Sum4,H_int4,H_sum2,...
        {'u2', 'v2'},{'y2'}));
    
    
    % Create ABCD matrix
    ABCD_model_SD2 = [Hsd_model_SD2.A Hsd_model_SD2.B; ...
        Hsd_model_SD2.C Hsd_model_SD2.D];
    
    % Digital part
    
    H_d1 = ss(1/z);
    H_d1.InputName = 'din1';
    H_d1.OutputName = 'd1';
    
    H_d2 = ss(coeff.d2);
    H_d2.InputName = 'din2';
    H_d2.OutputName = 'dd2';
    
    H_s1 = ss(coeff.s1);
    H_s1.InputName = 'd1';
    H_s1.OutputName = 'ds1';
    
    Sum6 = sumblk('dd3 = dd2+ds1');
    
    H_d3 = ss((1-(1/z))^2);
    H_d3.InputName = 'dd3';
    H_d3.OutputName = 'd2';
    
    Sum7 = sumblk('ytot = d2+d1');
    
    Hsd_model_dig = minreal(connect(H_d1,H_d2,H_s1,Sum6,H_d3,Sum7,...
        {'din1', 'din2'},{'ytot'}));

end
