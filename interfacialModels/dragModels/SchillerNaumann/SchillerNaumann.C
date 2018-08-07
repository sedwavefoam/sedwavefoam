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

#include "SchillerNaumann.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SchillerNaumann, 0);

    addToRunTimeSelectionTable
    (
        dragModel,
        SchillerNaumann,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SchillerNaumann::SchillerNaumann
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

Foam::SchillerNaumann::~SchillerNaumann()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::SchillerNaumann::K
(
    const volScalarField& Ur
) const
{
    // introduce a mixture visosity for air-fluid mixture
    volScalarField rhoMix = (alpha1_*phaseb_.rho()+(scalar(1)-alpha_-alpha1_)*phasec_.rho())/(scalar(1)-alpha_);
    volScalarField nuMix = (alpha1_*phaseb_.nu()*phaseb_.rho()+(scalar(1)-alpha_-alpha1_)*phasec_.nu()*phasec_.rho())/(alpha1_*phaseb_.rho()+(scalar(1)-alpha_-alpha1_)*phasec_.rho());
    volScalarField Re(max(Ur*phasea_.d()/nuMix, scalar(1.0e-3)));

    volScalarField Cds
    (
        neg(Re - 1000)*(24.0*(1.0 + 0.15*pow(Re, 0.687))/Re)
      + pos(Re - 1000)*0.44
    );

    return 0.75*Cds*rhoMix*Ur/phasea_.d();
}


// ************************************************************************* //
