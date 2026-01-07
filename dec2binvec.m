function  out = dec2binvec(input, nbit, inv)
    qout = input;
    rem = zeros(nbit+1,length(qout));
    
    for jj = 1:length(qout)
        ii = 1+1;
        while (qout(jj) ~= 0)
            qout(jj) = floor(input(jj)/2);
            rem(ii,jj) = (input(jj)/2 - qout(jj))*2;
            input(jj) = qout(jj);
            ii = ii + 1;
        end
    end
    if inv == 1
        rem = not(rem);
    else
        rem = rem;
    end
    rem(1) = 1;
    
    out = reshape(rem',[],1);
end
