%% Ocean Potential Temperature first guess

%% use SST first guess (b0)
%  to make a first-guess reconstruction of
%  HMS Challenger data.

% reconstruction of HMS challenger
y0 = Gchall*b0(:);

% cost function value related to Challenger
% first term of 2-term cost function
% second term is actually zero because
% first-guess u = 0
J0 = challengercost(Gchall,y,b0(:),iW);

% the first-guess Atlantic and Pacific profiles
% DT = WOCE minus Challenger temperature difference
DT_pacz = (G_pacz_woce-G_pacz_chall)*b0(:);
DT_atlz = (G_atlz_woce-G_atlz_chall)*b0(:);
