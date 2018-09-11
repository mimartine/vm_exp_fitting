function y = expFun(parms, n)

a = parms(1);
b = parms(2);
c = parms(3);


if isscalar(n)
    x = 1:n;
else
    x = n;
end

y = a + b * exp(x*c);

end % of function... 