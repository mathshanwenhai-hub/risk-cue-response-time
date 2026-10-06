function github_smoke_test()
%GITHUB_SMOKE_TEST Short CI-only PDE run. Not publication data.
R=github_run_case(1,'smoke',24,.01,.20,'upwind',0);
assert(strcmp(R.status,'completed'));
assert(abs(R.time_final-.20)<1e-12);
assert(all(isfinite([R.P(:);R.N(:);R.W(:)])));
fprintf('GitHub MATLAB smoke test passed.\n');
end
