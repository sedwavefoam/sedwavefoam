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

#include "GidaspowSchillerNaumann.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GidaspowSchillerNaumann, 0);

    addToRunTimeSelectionTable
    (
        dragModel,
        GidaspowSchillerNaumann,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GidaspowSchillerNaumann::GidaspowSchillerNaumann
(
    const dictionary& interfaceDict,
    const volScalarField& alpha,
    const volScalarField& alpha1,
    const phaseModel& phasea,
    const phaseModel& phaseb,
    const phaseModel& phasec
)
:
    dragModel(interfaceDict, alpha, alpha1, phasea, phaseb, phasec)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::GidaspowSchillerNaumann::~GidaspowSchillerNaumann()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::GidaspowSchillerNaumann::K
(
    const volScalarField& Ur
) const
{
    volScalarField voidfraction(max(scalar(1) - alpha_, scalar(1e-6)));
    //    volScalarField bp(pow(beta, -2.65));
    volScalarField bp(pow(voidfraction, -phasea_.hExp()));

    //volScalarField Re(max(beta*Ur*phasea_.d()/phaseb_.nu(), scalar(1.0e-3)));
    // introduce a mixture visosity for air-fluid mixture
    volScalarField rhoMix = (alpha1_*phaseb_.rho()+(scalar(1)-alpha_-alpha1_)*phasec_.rho())/(scalar(1)-alpha_);
    volScalarField nuMix = (alpha1_*phaseb_.nu()*phaseb_.rho()+(scalar(1)-alpha_-alpha1_)*phasec_.nu()*phasec_.rho())/(alpha1_*phaseb_.rho()+(scalar(1)-alpha_-alpha1_)*phasec_.rho());
    volScalarField Re(max(voidfraction*Ur*phasea_.d()*phasea_.sF()/nuMix, scalar(1.0e-9)));

    volScalarField Cds
    (
        neg(Re - 1000)*(24.0*(1.0 + 0.15*pow(Re, 0.687))/Re)
      + pos(Re - 1000)*0.44
    );

    return 0.75*Cds*rhoMix*Ur*bp/(phasea_.d()*phasea_.sF());
}


// ************************************************************************* //
