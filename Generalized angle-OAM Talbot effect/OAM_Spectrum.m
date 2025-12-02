function spectrum = OAM_Spectrum(field,Angle,Rad,OAM_limit)
% Given a 2D Complex field, calculate the OAM spectrum of the field

% OAM_limit = 20;

spectrum = zeros(1,2*OAM_limit+1);

OAMmode = zeros(size(Angle));

for l = -OAM_limit:OAM_limit
    OAMmode = exp(1i.*l.*Angle);
    OAMmode(Rad>0.5) = 0;
    [ov, ~, integral] = calc_overlap(field, OAMmode);
    spectrum(1,l+OAM_limit+1) = integral;
%     spectrum(1,l+OAM_limit+1) = ov.*exp(1i.*angle(integral));
end

% plot(-OAM_limit:OAM_limit,spectrum);