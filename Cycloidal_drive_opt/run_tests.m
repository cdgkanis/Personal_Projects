function run_tests()
%RUN_TESTS Execute the project regression checks.
projectDir = fileparts(mfilename('fullpath'));
addpath(projectDir,fullfile(projectDir,'tests'));
run_all_tests();
end
