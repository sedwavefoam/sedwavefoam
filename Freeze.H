// This version is modified so that Ua=0 in flow direction when alpha>alphaMax

    forAll (alpha, celll)
    {
       Ua[celll] = 0*Ua[celll];
    }
    Ua.correctBoundaryConditions();
// update phia
    phia = (fvc::interpolate(Ua) & mesh.Sf());
