    IOdictionary RASProperties
    (
        IOobject
        (
            "RASProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    Switch turbulence
    (
        RASProperties.lookup("turbulence")
    );

    Switch kEpsilon
    (
        RASProperties.lookup("kEpsilon")
    );

    Switch MixingLength
    (
        RASProperties.lookup("MixingLength")
    );

    dictionary kEpsilonDict
    (
        RASProperties.subDictPtr("kEpsilonCoeffs")
    );

    dimensionedScalar Cmu
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Cmu",
            kEpsilonDict,
            0.09
        )
    );

    dimensionedScalar C1
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "C1",
            kEpsilonDict,
            1.44
        )
    );

    dimensionedScalar C2
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "C2",
            kEpsilonDict,
            1.92
        )
    );

    dimensionedScalar alphak
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "alphak",
            kEpsilonDict,
            1.0
        )
    );

    dimensionedScalar alphaEps
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "alphaEps",
            kEpsilonDict,
            0.76923
        )
    );

    dictionary wallFunctionDict
    (
        RASProperties.subDictPtr("wallFunctionCoeffs")
    );

    dimensionedScalar kappa
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "kappa",
            wallFunctionDict,
            0.41
        )
    );

    dimensionedScalar E
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "E",
            wallFunctionDict,
            9.8
        )
    );

    dictionary MixingLengthDict
    (
        RASProperties.subDictPtr("MixingLengthCoeffs")
    );

    dimensionedScalar kappaLM
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "kappaLM",
            MixingLengthDict,
            0.41
        )
    );

dimensionedScalar alphaMaxLM
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "alphaMaxLM",
            MixingLengthDict,
            0.6
        )
    );


if(kEpsilon)Info<<"kEpsilon turbulence Model"<<nl<<endl;
  
if(MixingLength)Info<<"Mixing Length turbulence Model"<<nl<<endl;
  
if(kEpsilon && MixingLength) 
  {
    Info<<"Both turbulence model can not be turned on at the same time"<<nl
        <<"--> correct constant/RASProperties"<<nl<<endl; return 0;
  }
else if(turbulence && (!kEpsilon && !MixingLength)) 
  {
    Info<<"If turbulence model is turned on at least one model should be switched on"<<nl
        <<"--> correct constant/RASProperties"<<nl<<endl; return 0;
  }

    if (RASProperties.lookupOrDefault("printCoeffs", false))
    {
        Info<< "kEpsilonCoeffs" << kEpsilonDict << nl
	    << "MixingLengthCoeffs" << MixingLengthDict << nl
            << "wallFunctionCoeffs" << wallFunctionDict << nl<< endl;
    }


    nearWallDist y(mesh);


    Info<< "Reading field k\n" << endl;
    volScalarField k
    (
        IOobject
        (
            "k",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field epsilon\n" << endl;
    volScalarField epsilon
    (
        IOobject
        (
            "epsilon",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

     dimensionedScalar epsilonSmall0
        (
         "epsilonsmall0",
         dimensionSet(0, 2, -3, 0, 0, 0, 0),
         scalar(1e-10)
         );

    Info<< "Calculating field nutb\n" << endl;
    volScalarField nutb
    (
        IOobject
        (
            "nutb",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Cmu*sqr(k)/max(epsilon, epsilonSmall0)
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

    volScalarField nuFra
    (
        IOobject
        (
            "nuFra",
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
        nutb + nugl
    );

// added for modifed k-epsilon model
     dimensionedScalar SUS
    (
        RASProperties.lookup("SUS")
    );
      dimensionedScalar C3ep
    (
        RASProperties.lookup("C3ep")
    );
      dimensionedScalar C4ep
    (
        RASProperties.lookup("C4ep")
    );
            dimensionedScalar KE1
    (
        RASProperties.lookup("KE1")
    );
            dimensionedScalar KE2
    (
        RASProperties.lookup("KE2")
    );
            dimensionedScalar KE3
    (
        RASProperties.lookup("KE3")
    );
    dimensionedScalar KE4
    (
        RASProperties.lookup("KE4")
    );
    dimensionedScalar B
    (
        RASProperties.lookup("B")
    );

if(SUS.value()==0) Info<<"Turbulence suspension term is neglected"<<nl<<endl;
 else if(SUS.value()>0) Info<<"Turbulence suspension term is included"<<nl<<endl;
 else {Info<<"Turbulence suspension coeffcient SUS can't be negative"<<nl<<endl; exit(FatalError);}
