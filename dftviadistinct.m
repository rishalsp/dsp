
x = [0,7,0,69,0,7,0,-7,0,-69,0,-7];
dftdistinct(x);


function Xdis = dftdistinct(x)
    Xdis = [];
    Xdis(1,1) = 0;
    Xdis(1,2) = 0;
    i = 1;
    r = length(x);
    %tau = ndistinct(r);
    %count = 0;
    for D=2:r
        %if count<tau
            if mod(r,D)==0 %ie. D divides r
                val = 0;
                for d=2:r
                    if mod(r,d)==0 && mod(d,4)==0
                        val = val + x((r/d)+1)*cd(d,((r/D)+(d/4)));
                    end
                end
                %if ismember(val,Xdis(:,2))==0
                    i = i+1;
                    Xdis(i,1) = r/D;
                    Xdis(i,2) = 1j*val;
                    %count = count+1;
                %end
            end
        %end
    end
end


function tau = ndistinct(r)
%This function outputs the number of distinct nonzero values that xr(n)
%can take
%ndistinct = (m1+1)(m2+1)...(mk+1) where mi's are the powers of the
%prime factors of r/4
    primearr = factor(r/4)';
    marr = groupcounts(primearr);
    tau = marr(1)+1;
    if length(marr)>1
        for i=2:length(marr)
            tau = tau*(marr(i)+1);
        end
    end
end


function sum = cd(d, x) 
%Outputs cd(c+(d/4)) - the odd ramanujan sum
    sum = 0;
    for U=0:d-1
        if gcd(U,d)==1
            Wd = W(d);
            sum = sum + Wd^(-x*U);
        end
    end
end


function out = W(d)
%Outputs Wd
    out = exp(2*pi*1j/d);
end



    