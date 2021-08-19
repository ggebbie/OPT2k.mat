%% instantiate TMI: download/upgrade TMI dependent package

% need to clone TMI project?
if ~exist('../TMI') && ~exist('../../TMI')
    display('clone TMI project')
    !git clone https://github.com/ggebbie/TMI ../TMI
    dir.TMI = '../TMI/src/';
else
    display('check for TMI project updates')
    if exist('../../TMI/')
        dir.TMI = '../../TMI/src/';
        % cd ../../TMI
    elseif exist('../TMI/')
        dir.TMI = '../TMI/src/';
        % cd ../TMI
    end
    
    % assumes origin is set to GitHub repository.
    !git fetch origin
    !git pull origin refs/heads/main
end
addpath(dir.TMI)
