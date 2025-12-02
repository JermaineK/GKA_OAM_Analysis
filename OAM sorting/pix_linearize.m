% Mapping from measured pixel values to true, linear values
function result = pix_linearize(c,pix)
    result = zeros(size(pix));
    % polynomial terms with coefficients in vector c
    for it = 1:length(c)
        result = result + c(it).*pix.^(it);
    end
end