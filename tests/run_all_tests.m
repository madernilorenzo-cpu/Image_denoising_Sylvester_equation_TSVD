%RUN_ALL_TESTS Runs all the tests of the project.
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, '..', 'src'));

test_build_sylvester_data();
test_bartels_stewart();
test_sylvester_tsvd();

fprintf('All tests passed.\n');
