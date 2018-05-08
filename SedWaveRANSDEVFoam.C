/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    twoPhaseEulerFoam

Description
    Solver for a system of 2 incompressible fluid phases with one phase
    dispersed, e.g. gas bubbles in a liquid or solid particles in a gas.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "MULES.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "singlePhaseTransportModel.H"
#include "nearWallDist.H"
#include "wallFvPatch.H"
#include "Switch.H"
#include "IFstream.H"
#include "OFstream.H"
#include "symmetryFvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "fixedFluxPressureFvPatchScalarField.H"
// added for wave2Foam
#include "relaxationZone.H"
#include "bound.H"
#include "dragModel.H"
#include "phaseModel.H"
#include "kineticTheoryModel.H"
#include "mathematicalConstants.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "fixedFluxPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    // added for wave2Foam
    #include "readWaveProperties.H"
    #include "createFields.H"
    #include "readPPProperties.H"
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
    #include "createGradP.H"
    pimpleControl pimple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    while (runTime.run())
    {
        #include "readTwoPhaseEulerFoamControls.H"
        #include "CourantNos.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
        // first solve the air-water interface
            beta = 1.0 - gamma; 
        // then solve for alpha, after that, get gamma = 1-beta
            #include "betaControls.H"
            if (pimple.firstIter() || betaOuterCorrectors)
            {
                #include "betaEqnSubCycle.H"
                //added for wave2Foam
                relaxing.correct(); // check createFields
                relaxing2.correct(); // check createFields
                interface.correct();
            }
        // then update gamma
            gamma = 1.0 - beta;
        // second: solve mixture-sediment 
            #include "alphaEqn.H"
        // gamma = air, alpha = sed, alpha1 = 1-gamma-alpha 
            alpha1 = 1.0 - gamma - alpha;
        // update the mixture properties
            rho == (gamma*rhoc + alpha1*rhob)/(scalar(1.0)-alpha);
            nugl=(alpha1*rhob*nub + gamma*rhoc*nuc)/(gamma*rhoc + alpha1*rhob);

            #include "liftDragCoeffs.H"
            if (pimple.turbCorr())
            {
                #include "kEpsilon.H"
            }
            #include "callKineticTheory.H"
            #include "UEqns.H"
            
            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
                if (correctAlpha && !pimple.finalIter())
                {
                    Info<< "Correct alpha after pstep" <<endl;
                    #include "alphaEqn.H"
                }
            }
            Info<< "Max Ua =" << max(Ua).value() <<",Min Ua =" << min(Ua).value()<<endl;
            Info<< "Max Ub =" << max(Ub).value() <<",Min Ub =" << min(Ub).value()<<endl;
            #include "DDtU.H"
        }
        #include "write.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
