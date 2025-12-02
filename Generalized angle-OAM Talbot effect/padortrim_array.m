%% A function that adds zeros or removes elements from array
% First input is the array and the second and third inputs are the desired x and y
% dimensions


function Array = padortrim_array(Array, nx, ny)
    A_size = size(Array);
    Ax = A_size(1);
    Ay = A_size(2);
    
    % First x-dir
    if Ax < nx
        Array = padarray(Array,[floor((nx-Ax)/2),0], 0, 'both');
        % Extra 1 row of padding if needed
        if size(Array,1) < nx
            Array = padarray(Array,[1 0],0,'pre');
        end
    elseif Ax > nx
        DiffX = abs(Ax-nx);
        % Scale axes up to equal and resize down
        Array = padarray(Array,[0, floor(DiffX/2)], 0, 'both');
        % Extra 1 row of padding if needed
        if size(Array,2) < Ay+DiffX
            Array = padarray(Array,[0 1],0,'pre');
        end
        % resizing down
        Array = imresize(Array, [nx,Ay], 'nearest');
    end
    
    % Second y-dir
    if Ay < ny
        Array = padarray(Array,[0, floor((ny-Ay)/2)], 0, 'both');
        % Extra 1 row of padding if needed
        if size(Array,2) < ny
            Array = padarray(Array,[0 1],0,'pre');
        end
    elseif Ay > ny
        DiffY = abs(Ay-ny);
        % Scale axes up to equal and resize down
        Array = padarray(Array,[floor(DiffY/2), 0], 0, 'both');
        % Extra 1 row of padding if needed
        if size(Array,2) < nx+DiffY
            Array = padarray(Array,[1 0],0,'pre');
        end
        % resizing down
        Array = imresize(Array, [nx,ny], 'nearest');
    end
    
end