
    Info<< "Reading field k\n" << endl;
    volScalarField k
    (
        IOobject
        (
            "k",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
       sgsModelb->k()
    );
    volScalarField epsilon
    (
        IOobject
        (
            "epsilon",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
       sgsModelb->epsilon()
    );
    volScalarField nutb
    (
        IOobject
        (
            "nutb",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
       sgsModelb->nuSgs()
    );

    Info<< "Calculating field nuEffa\n" << endl;
    volScalarField nuEffa
    (
        IOobject
        (
            "nuEffa",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sqr(Ct)*nutb + nua
    );

    Info<< "Calculating field nuEffb\n" << endl;
    volScalarField nuEffb
    (
        IOobject
        (
            "nuEffb",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        nutb + nub
    );

