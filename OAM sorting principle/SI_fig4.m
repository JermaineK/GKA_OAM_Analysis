% Simulation of the OAM mode sorter, SI fig 4
% Jianqi Hu (jianqi.hu@lkb.ens.fr), Last updated (11/10/2024)

clear;
clc;
close all

% fiber length 
% p/q (p,q is coprime) of Talbot length (half of the Talbot length defined in Wiki) 
p= 1; 
q= 4;


% when q is odd, p is odd, the 0-th order starts from middle (\pi)
% when q is even, the 0-th order starts at 0;
% when q is odd, p is even, the 0-th order starts at 0;


fonts =24;
lines =3;

% oam order
for l0 =0:3

% phase
sampling_length =300; % total number of sampling point 2*q*sampling_length
phi = 0:1/(sampling_length*q)*pi: 2*pi-1/(sampling_length*q)*pi;

% input oam mode
e_in = exp(1i*l0*phi);

% find s/q 
if mod(q,2) == 1
   s = 2 * modinv(2*p,q);
else 
   s = modinv(p,2*q);
end

% talbot phase
for k =0:q-1
   theta(k+1) = -s*pi*k^2/q; 
end
 

% talbot phase mask
t_mask = reshape(repmat(theta,2*sampling_length,1),sampling_length*q*2,1)';

%figure,plot(1:length(t_mask),t_mask)

% output oam complex amplitude
e_out = e_in .* exp(1i*t_mask);

% oam fft
e_out_fft= fftshift(fft(e_out))/(2*q*sampling_length);

% oam x-axis
oam_range = -q*sampling_length:q*sampling_length-1;

% oam spectrum amplitdue plot
% figure,
% plot(oam_range, abs(e_out_fft))
% hold on
% plot(oam_range,abs(sinc((oam_range-l0)/q))/sqrt(q),'--','Linewidth',lines)
% xlim([-50,50])
% plot(oam_range,angle(e_out_fft))


ratio = p/q;
step = ratio;
step_num = length(step);

carpet_out = zeros(step_num,2*q*sampling_length);

for k =1:step_num
corr_pha = (oam_range).^2*pi*step(k);
e_out_fft_proc = e_out_fft.* exp(-1i*corr_pha);
e_out_fft_proc_ifft = ifft(ifftshift(e_out_fft_proc))*(2*q*sampling_length);

carpet_out(k,:) =abs(e_out_fft_proc_ifft);
end


%% output figure
figure2 = figure;
axes2= axes('Parent',figure2);
hold(axes2,'on');

%yyaxis left
phi_circle = [phi 2*pi]; % put 0 as 2*pi
e_out_fft_proc_ifft_circ = [e_out_fft_proc_ifft, e_out_fft_proc_ifft(1)];
plot(phi_circle,abs(e_out_fft_proc_ifft_circ).^2,'LineWidth',lines)
xlim([0 2*pi]), xlabel('Azimuthal angle (rad)'), ylabel('Intensity (a.u.)')

set(axes2,'FontSize',fonts,'XTick',...
    [0 pi/2 pi 3*pi/2 2*pi]);
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})


set(axes2,'FontSize',fonts,'YTick',...
    [0 1 2 3 4]);
yticklabels({'0','1','2','3','4'})
ylim([-0.25 4.25])

[~,phi_ind] = find(abs(e_out_fft_proc_ifft)>1e-4);
phi_new = phi(phi_ind);

pha_tai = angle(e_out_fft_proc_ifft(phi_ind)); 
% pha_tai = angle(e_out_fft_proc_ifft(phi_ind))+pi; 
% pha_tai = pha_tai - pha_tai(1); % remove the constant phase

% yyaxis right 
% plot(phi_new,pha_tai,'LineWidth',lines)
% ylim([-pi*1.125,pi*1.125]), ylabel('Phase (rad)')

box(axes2,'on');
hold(axes2,'off');

% set(axes2,'FontSize',fonts,'YTick',...
%     [-pi -pi/2 0 pi/2 pi]);
% yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})

ax = gca(); 
ax.XGrid = 'on'; 
yl = arrayfun(@(y)yline(ax, y,'LineStyle',ax.GridLineStyle,'Color',ax.GridColor,'Alpha',ax.GridAlpha),ax.YTick);

end 




