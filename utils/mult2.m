function Z = mult2(X,Y)

if verLessThan('matlab','9.1') % for MATLAB 2016a and earlier
    Z = bsxfun(@times,X,Y);
else % for MATLAB 2016b and later (faster)
    Z = X.*Y;
end