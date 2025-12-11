function W = weightingfunctions(xx,tt)

% sigma_t = 0.25;
sigma_t = 1;
sigmaGabor = sigma_t;

W = zeros([numel(xx),numel(xx),4]);

%cos weighting works only if I take into account the 16th orientations
% tt16 = tt(:).*xx(:);

W(:,:,1) = exp( ...
    -(xx(:).*cos(tt(:)-tt(:)') - xx(:)').^2./ ...
    (2*sigma_t.^2)); %Andre pesi originali (inviluppo exp(.))
%to be generalised to all numbers of cells
W(1:40,49:end,2) = -cos((tt(1:40)-tt(49:end)'));    %this changing allows for 
                                                    %take into account that 
                                                    %orientation channels in 
                                                    %my network spans only
                                                    %in the range between
                                                    %0 and pi and the
                                                    %others are encoded by
                                                    %negative velocities

W(49:end,49:end,2) = cos(tt(49:end)-tt(49:end)');

W(1:40,1:40,2) = (cos(tt(1:40)-tt(1:40)'));
W(49:end,1:40,2) = -(cos(tt(49:end)-tt(1:40)'));

W(:,:,3) = cos(xx(:).*cos(tt(:)'-tt(:)) - xx(:)'); %primo modello per introdurre l'opponenza (inviluppo coseno)
W(:,:,4) = exp(-(xx(:).*cos(tt(:)'-tt(:)) - xx(:)').^2./(2*sigmaGabor.^2)).* ...
    cos(2*pi*1./(4*xx(:)').*(xx(:).*cos(tt(:)'-tt(:)) - xx(:)')); %inviluppo Gabor per gestire la quantità e la posizione dei pesi negativi
% W2 = reshape(reshape(W2,8,11,8,11)./max(reshape(W2,8,11,8,11),[],4),88,88);
W(:,:,5) = (xx(:).*cos(tt(:)-tt(:)') -xx(:)'); %Andre pesi originali)


%
W(isnan(W)) = 0;
% W = W./repmat(sum(W,[1,2]),88,88,1);

