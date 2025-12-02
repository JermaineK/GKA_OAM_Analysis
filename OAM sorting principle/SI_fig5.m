% Simulation of the carpet for mode sorter, SI fig 5
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





% oam order
for l0 =0:4

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
% plot(oam_range,abs(sinc((oam_range-l0)/q))/sqrt(q),'--','Linewidth',2)
% xlim([-50,50])
% plot(oam_range,angle(e_out_fft))


ratio = p/q;
step = 0:ratio/3000 :ratio;
step_num = length(step);

carpet_out = zeros(step_num,2*q*sampling_length);

for k =1:step_num
corr_pha = (oam_range).^2*pi*step(k);
e_out_fft_proc = e_out_fft.* exp(-1i*corr_pha);
e_out_fft_proc_ifft = ifft(ifftshift(e_out_fft_proc))*(2*q*sampling_length);

carpet_out(k,:) =abs(e_out_fft_proc_ifft);
end

%% input fiugre
figure3 = figure;
axes3= axes('Parent',figure3);
hold(axes3,'on');

yyaxis left
phi_circle = [phi 2*pi]; % put 0 as 2*pi

plot(phi_circle,ones(size(phi_circle)),'LineWidth',2)
xlim([0 2*pi]), xlabel('Angle (rad)'), ylabel('Intensity (a.u.)')

set(axes3,'FontSize',20,'XTick',...
    [0 2.0943951023932 4.18879020478639 6.2831853071795]);
xticklabels({'0','2\pi/3','4\pi/3','2\pi'})

set(axes3,'FontSize',20,'YTick',...
    [0 1 2 3]);
yticklabels({'0','1','2','3'})
ylim([0 3*2.25/2])

yyaxis right 
plot(phi_circle,mod(l0*phi_circle,2*pi),'LineWidth',2)
ylim([0,2.25*pi]), ylabel('Phase (rad)')

hold on 
t_mask_circle = [t_mask t_mask(1)];
plot(phi_circle,mod(t_mask_circle,2*pi),'--','LineWidth',2,'color',[0.5 0.5 0.5])


box(axes3,'on');
hold(axes3,'off');

set(axes3,'FontSize',20,'YTick',...
    [0 2*pi/3 4*pi/3 2*pi]);
yticklabels({'0','2\pi/3','4\pi/3','2\pi'})

ax = gca(); 
ax.XGrid = 'on'; 
yl = arrayfun(@(y)yline(ax, y,'LineStyle',ax.GridLineStyle,'Color',ax.GridColor,'Alpha',ax.GridAlpha),ax.YTick);

%% output figure
figure2 = figure;
axes2= axes('Parent',figure2);
hold(axes2,'on');

yyaxis left
phi_circle = [phi 2*pi]; % put 0 as 2*pi
e_out_fft_proc_ifft_circ = [e_out_fft_proc_ifft, e_out_fft_proc_ifft(1)];
plot(phi_circle,abs(e_out_fft_proc_ifft_circ).^2,'LineWidth',2)
xlim([0 2*pi]), xlabel('Angle (rad)'), ylabel('Intensity (a.u.)')

set(axes2,'FontSize',20,'XTick',...
    [0 2.0943951023932 4.18879020478639 6.2831853071795]);
xticklabels({'0','2\pi/3','4\pi/3','2\pi'})

set(axes2,'FontSize',20,'YTick',...
    [0 1 2 3]);
yticklabels({'0','1','2','3'})
ylim([0 3*2.25/2])

[~,phi_ind] = find(abs(e_out_fft_proc_ifft)>1e-4);
phi_new = phi(phi_ind);

pha_tai = angle(e_out_fft_proc_ifft(phi_ind))+ pi;
% pha_tai = pha_tai - pha_tai(1); % remove the constant phase
yyaxis right 
plot(phi_new,pha_tai,'LineWidth',2)
ylim([0,2.25*pi]), ylabel('Phase (rad)')

box(axes2,'on');
hold(axes2,'off');

set(axes2,'FontSize',20,'YTick',...
    [0 2*pi/3 4*pi/3 2*pi]);
yticklabels({'0','2\pi/3','4\pi/3','2\pi'})

ax = gca(); 
ax.XGrid = 'on'; 
yl = arrayfun(@(y)yline(ax, y,'LineStyle',ax.GridLineStyle,'Color',ax.GridColor,'Alpha',ax.GridAlpha),ax.YTick);





figure1= figure,
axes1 = axes('Parent',figure1);
%plot(phi,abs(carpet_out(1,:)))
% mesh(phi,carpet_out)
imagesc(carpet_out')
%xlim([0,2*pi])
%hold on 
%plot(phi,angle(e_out_fft_proc_ifft))
figure1.Position = [600 600 1600 200]
% 
% set(axes1,'XAxisLocation','top','FontSize',18,'Layer','top','XTick',...
%      [1 1001 2001 3001 4001],'XTickLabel',...
%      {'0','1/4','1/2','3/4'},'YTick',0,...
%      'YTickLabel','');

% xlabel('p/q','FontAngle','italic','FontSize',18);
axis off
end 
