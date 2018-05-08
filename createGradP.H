Info<< "\nReading forceProperties for gradP\n" <<endl;

    IOdictionary forceProperties
    (
        IOobject
        (
            "forceProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar tilt
    (
        forceProperties.lookup("tilt")
    );

    dimensionedVector gradPAMP1
    (
        forceProperties.lookup("gradPAMP1")
    );
    dimensionedVector gradPAMP2
    (
        forceProperties.lookup("gradPAMP2")
    );

    dimensionedVector gradPMEAN
    (
        forceProperties.lookup("gradPMEAN")
    );

    dimensionedVector gradPOSC
    (
        forceProperties.lookup("gradPOSC")
    );

      dimensionedScalar oscpT
    (
        forceProperties.lookup("oscpT")
    );
      dimensionedScalar initTheta
    (
        forceProperties.lookup("initTheta")
    );
    dimensionedScalar InitFreeze
    (
        forceProperties.lookup("InitFreeze")
    );
    Info<< "Initializing with specified pressure gradient:"<<gradPOSC.value()<<endl
            << endl;
    // open files to write gradP info
    OFstream gradPOSCFile
    (
    "gradPOSC.txt"
    );
