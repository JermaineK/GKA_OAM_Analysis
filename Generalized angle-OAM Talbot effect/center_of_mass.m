% A function to find the center of mass of a complex-valued matrix

function [cx, cy] = center_of_mass(A,threshold)
    % Take abs value
    A = abs(A);
    
    % Normalize to 0-1
    A = A./max(max(A));
    
    % Threshold
    A(A<threshold) = 0;
    
    % Dimensions
    [ny,nx] = size(A);
    
    % Normalization factor
    norm = sum(sum(A));
    
    cx = 0;
    cy = 0;
    
    % y-center
    for it = 1:ny
        cy = cy + it*sum(A(it,:))/norm;
    end
    
    % x-center
    for it = 1:nx
        cx = cx + it*sum(A(:,it))/norm;
    end
    
    figure(52)
    imagesc(A)
    hold on
    scatter(cx,cy);
    title('Center of mass');
end

    