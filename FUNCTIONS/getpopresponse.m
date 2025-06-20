function popresp = getpopresponse(W,e)

sz = size(W);
sze = size(e);

if ndims(W) > 2
    W = reshape(W,sz(1),sz(2),[]);
    e = reshape(e,prod(sze(1:2)),[]);
    for ii = 1:prod(sz(3:end))
        pop_resp(:,ii) = W(:,:,ii)*e(:,ii);
    end
    popresp = reshape(pop_resp,sze);
end