%% intensity_phase_plot
% Displays the electric field phase with a colormap and the intensity overlayed as
% alpha (transparency) values.
function intensity_phase_plot(Mode_Efield)

    Int = abs(Mode_Efield).^2;
    Int = Int/max(max(Int));
    Ang = angle(Mode_Efield);
    
    % Colormap for angular phase:
%     colorscheme = othercolor('RdBu9',1000);
%     set(gcf,'Colormap',colorscheme);
    
%     grad=[colorGradient([0.3 0.3 0.3], [0.6980    0.0941    0.1686], round(256/9));...
%     colorGradient([0.6980    0.0941    0.1686], [0.8392    0.3765    0.3020], round(256/9));...
%     colorGradient([0.8392    0.3765    0.3020], [0.9569    0.6471    0.5098], round(256/9));...
%     colorGradient([0.9569    0.6471    0.5098], [0.9922    0.8588    0.7804], round(256/9));...
%     colorGradient([0.9922    0.8588    0.7804], [0.8196    0.8980    0.9412], round(256/9));...
%     colorGradient([0.8196    0.8980    0.9412], [0.5725    0.7725    0.8706], round(256/9));...
%     colorGradient([0.5725    0.7725    0.8706], [0.2627    0.5765    0.7647], round(256/9));...
%     colorGradient([0.2627    0.5765    0.7647], [0.1294    0.4000    0.6745], round(256/9));...
%     colorGradient([0.1294    0.4000    0.6745], [0.3 0.3 0.3], round(256/9))];
%     colormap(grad);
    colormap(hsv)
    imag = imagesc(Ang,[-pi pi]);
    
    %colormap(hot)
    %colormap(othercolor('OrRd3'))
    %imag = imagesc(Int,[0 1]);
    
    %colormap(hsv) % Optional colormap 
    set(imag, 'AlphaData', Int) %overlay with intensity
    set(gca,'color',[0 0 0]) %set background black
    %set(gcf, 'color', 'none');
    %set(gca, 'color', 'none');
    
    % Get rid of axes:
    set(gca,'xtick',[])
    set(gca,'ytick',[]) 

end