
% need to clone TMI project?

if ~exist('../TMI') & ~exist('../../TMI')
    display('clone TMI project')
    !git clone https://github.com/ggebbie/TMI ../TMI
    TMIdir = '../TMI/';
elseif exist('../TMI')
    display('TMI project available')
    TMIdir = '../TMI/';
elseif exist('../../TMI')
    TMIdir = '../../TMI';
end
addpath(TMIdir)

cd ..
read_TMI_from_google_drive

load ~/mcode/transient/d_all_2deg.mat
load ~/mcode/TMI_v7/L_2deg_2012

  NZ = max(kt);
  NY = max(jt);
  NX = max(it);
  dx = LON(2)-LON(1);
  dy = LAT(2)-LAT(1);
  for ny = 1:length(LAT)
    distx(ny) = sw_dist([LAT(ny) LAT(ny)],[0 dx],'km').*1000;
  end
  disty = sw_dist([0 dy],[0 0],'km').*1000;
  distx = reshape(distx,length(distx),1);
  areaa = disty(ones(NY,NX)).*distx(:,ones(NX,1)); 
  zface= (DEPTH(1:end-1)+DEPTH(2:end))./2;

  % technically this dz is wrong: center depths not in center of cells.
  dz = ([zface(1) ; diff(zface); 500]);
  vl = permute(areaa(:,:,ones(NZ,1)),[3 1 2]) .*dz(:,ones(NY,1),ones(NX,1));
  vol = field_to_vector(vl,it,jt,kt);
