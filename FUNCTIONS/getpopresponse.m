function popresp = getpopresponse(W,e)

sz = size(W);
sze = size(e);

if ndims(W) > 2
    W = reshape(W,sz(1),sz(2),[],sz(end));
    e = repmat(reshape(e,prod(sze(1:2)),[]),1,1,sz(end));
    for ii = 1:prod(sz(3:end-1))
        for jj = 1:sz(end)
            pop_resp(:,ii,jj) = squeeze(W(:,:,ii,jj))*squeeze(e(:,ii,jj));
        end
    end
    popresp = reshape(pop_resp,[sze,sz(end)]);
end