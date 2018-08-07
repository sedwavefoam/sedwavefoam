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

\*---------------------------------------------------------------------------*/

#include "SyamlalOBrien.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SyamlalOBrien, 0);

    addToRunTimeSelectionTable
    (
        dragModel,
        SyamlalOBrien,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SyamlalOBrien::SyamlalOBrien
(
    const dictionary& interfaceDict,
    const volScalarField& alpha,
    const volScalarField& alpha1,
    const phaseModel& phasea,
    const phaseModel& phaseb,
    const phaseModel& phasec
)
:
    dragModel(interfaceDict, alpha, alpha1,phasea, phaseb,phasec)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SyamlalOBrien::~SyamlalOBrien()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::SyamlalOBrien::K
(
    const volScalarField& Ur
) const
{
    volScalarField voidfraction(max(scalar(1) - alpha_, scalar(1.0e-6)));
    volScalarField A(pow(voidfraction, 4.14));
    volScalarField B
    (
        neg(voidfraction - 0.85)*(0.8*pow(voidfraction, 1.28))
      + pos(voidfraction - 0.85)*(pow(voidfraction, 2.65))
    );

    // introduce a mixture visosity for air-fluid mixture
    volScalarField rhoMix = (alpha1_*phaseb_.rho()+(scalar(1)-alpha_-alpha1_)*phasec_.rho())/(scalar(1)-alpha_);
    volScalarField nuMix = (alpha1_*phaseb_.nu()*phaseb_.rho()+(scalar(1)-alpha_-alpha1_)*phasec_.nu()*phasec_.rho())/(alpha1_*phaseb_.rho()+(scalar(1)-alpha_-alpha1_)*phasec_.rho());
    volScalarField Re(max(Ur*phasea_.d()/nuMix, scalar(1.0e-3)));

    volScalarField Vr
    (
        0.5*
        (
            A - 0.06*Re + sqrt(sqr(0.06*Re) + 0.12*Re*(2.0*B - A) + sqr(A))
        )
    );

    volScalarField Cds(sqr(0.63 + 4.8*sqrt(Vr/Re)));

    return 0.75*Cds*rhoMix*Ur/(phasea_.d()*sqr(Vr));
}


// ************************************************************************* //
