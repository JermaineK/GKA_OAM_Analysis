function [overlap, norm_const, integral] = calc_overlap(mode1, mode2)
    % Calculates the overlap between the two given modes
    % The function gives a normalization constant that can be used to
    % normalize a mode if mode1 and mode2 are the same
    
    integral = sum(conj(mode1).*mode2, 'all'); % This one is faster
    % Take real value since any imaginary part of the integral should be
    % a numerical error:
    norm_const = sqrt(complex(real(integral)));
    overlap = abs(integral).^2;
end