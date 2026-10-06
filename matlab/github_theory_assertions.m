function study = github_theory_assertions()
%GITHUB_THEORY_ASSERTIONS Reproduce benchmark thresholds and algebra checks.
study=run_siap_study('theory');
S=study.spectrum;
tin=S.entry.tau_minus;tout=S.exit.tau_plus;
assert(abs(tin-0.6087956390)<5e-7,'Unexpected lower response-time threshold.');
assert(abs(tout-6.6445349603)<5e-7,'Unexpected upper response-time threshold.');
entry=sort([S.entry.m,S.entry.n]);exitm=sort([S.exit.m,S.exit.n]);
assert(isequal(entry,[4 6]),'Unexpected entry mode family.');
assert(isequal(exitm,[4 8]),'Unexpected exit mode family.');
assert(all(S.tail_test_passed),'Spectral tail certification failed.');
assert(study.entry.eigen_residual<1e-7 && study.exit.eigen_residual<1e-7);
assert(study.entry.symmetric_error<1e-7 && study.exit.symmetric_error<1e-7);
fprintf('Theory assertions passed: tau_in=%.10f, tau_out=%.10f.\n',tin,tout);
end
