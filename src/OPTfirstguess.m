%% Ocean Potential Temperature first guess

%% use b0 to make a first-guess reconstruction of
%  HMS Challenger data.
y0 = Gchall*b0(:);
J0 = challengercost(Gchall,y,b0(:),iW);

% the first-guess Atlantic and Pacific profiles
DT_pacz = (G_pacz_woce-G_pacz_chall)*b0(:);
DT_atlz = (G_atlz_woce-G_atlz_chall)*b0(:);
