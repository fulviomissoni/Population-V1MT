function y = Weibull(t,b,x)
%y = Weibull(p,x)
%
%Parameters:  p.b slope
%             p.t threshold yeilding ~80% correct
%             x   intensity values.

    g = 0;          % chance performance
    %e = (.5)^(1/3);   % threshold performance ( ~80%)
    e = 0.5;   % threshold performance ( ~75%)

    %here it is.
    k = (-log( (1-e)/(1-g)))^(1/b);
    y = 1- (1-g)*exp(- (k*x/t).^b);
end