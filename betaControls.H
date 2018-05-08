const dictionary& betaControls = mesh.solverDict(beta.name());

label nBetaCorr(readLabel(betaControls.lookup("nBetaCorr")));

label nBetaSubCycles(readLabel(betaControls.lookup("nBetaSubCycles")));

bool MULESCorr(betaControls.lookupOrDefault<Switch>("MULESCorr", false));

bool betaOuterCorrectors
(
    betaControls.lookupOrDefault<Switch>("betaOuterCorrectors", false)
);

// Apply the compression correction from the previous iteration
// Improves efficiency for steady-simulations but can only be applied
// once the alpha field is reasonably steady, i.e. fully developed
bool betaApplyPrevCorr
(
    betaControls.lookupOrDefault<Switch>("betaApplyPrevCorr", false)
);

// Isotropic compression coefficient
scalar icBeta
(
    betaControls.lookupOrDefault<scalar>("icBeta", 0)
);
