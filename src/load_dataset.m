function [ dataset ] = load_dataset( pars )

t=tic;
fprintf('Loading %s dataset...', pars.dataset.name);

namefile = fullfile(pars.settings.rootFolder, 'data', ['dataset_' pars.dataset.name, '.mat']);
load(namefile);
fprintf('done in %.2f(s)\n', toc(t));

end

