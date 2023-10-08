
% Digital Signal Processing - Monsoon 2022 - Course Project
% Odd Ramanujan Sums of Complex Roots of Unity
% Project done by:
% Rishal - 2010110514
% Bhavani Namrata - 2010110189
% Akshat Sabavat - 2010110063

%-----------------------------------------

% This program computes the discrete fourier transform of 
% odd periodic signals as defined in the paper by using 
% odd ramanujan sums of complex roots of unity

x = [0,7,0,69,0,7,0,-7,0,-69,0,-7]; 
%sample input example as given in the paper
ourDft = dftall(x)' %taking dft using approach discussed in the paper
inbuiltDft = fft(x)' %taking dft using inbuilt fft algorithm
ourIdft = idftall(ourDft)' %taking idft using our approach
inbuiltIdft = ifft(inbuiltDft) %taking idft using inbuilt algorithm

%plotting the signal, fourier transforms, and the inverse transforms
subplot(5,1,1)
plot(x)
title('Main period of Odd Signal mod 12')
subplot(5,1,2)
plot(abs(ourDft))
title('DFT Using Odd Ramanujan Sums')
subplot(5,1,3)
plot(abs(inbuiltDft))
title('DFT Using Typical FFT Algo - Inbuilt')
subplot(5,1,4)
plot(abs(ourIdft))
title('IDFT Using Odd Ramanujan Sums')
subplot(5,1,5)
plot(abs(inbuiltIdft))
title('IDFT Using Typical IFFT Algo - Inbuilt')

%------------------------------------
% The functions which are used to implement the algorithm:

function x = idftall(X)
%Computes the idft of the fourier transform of an odd signal mod r
%using an approach making use of odd ramanujan sums
%This implements equation (14) given in the paper
    x = [];
    r = length(X);
    for n=0:r-1
        val = 0;
        for d=2:r
            if mod(r,d)==0 && mod(d,4)==0
                val = val + X((r/d)+1)*cd(d,n+(d/4));
            end
        end
        x = [x,-1j*val/r];
    end
end


function X = dftall(x)
%Computes the dft of the fourier transform of an odd signal mod r
%using an approach making use of odd ramanujan sums
%This implements equation (10) given in the paper
    X = [];
    r = length(x);
    for n=0:r-1
        val = 0;
        for d=2:r
            if mod(r,d)==0 && mod(d,4)==0
                val = val + x((r/d)+1)*cd(d,n+(d/4));
            end
        end
        X = [X,1j*val];
    end
end


function tau = ndistinct(r)
%This function outputs the number of distinct nonzero values 
%that xr(n) can take.
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
%this computes cd(n+(d/4)) - the odd ramanujan sum
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



    