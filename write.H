    if (runTime.outputTime())
    {
        volTensorField Tauf
        (
           IOobject
           (
                "Tauf",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
           ),
           nuEffb*(fvc::grad(Ub)+fvc::grad(Ub)().T())
        );
        volTensorField Taus
        (
           IOobject
           (
                "Taus",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
           ),
           (nuEffa)*(fvc::grad(Ua)+fvc::grad(Ua)().T()) // regarding nuEffa, check kineticTheoryModel.C
        );
    
    runTime.write();
    }
