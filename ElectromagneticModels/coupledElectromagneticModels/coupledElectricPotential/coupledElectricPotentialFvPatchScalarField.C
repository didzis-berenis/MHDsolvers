/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "coupledElectricPotentialFvPatchScalarField.H"
#include "electromagneticModel.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"
#include "mappedPatchBase.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::coupledElectricPotentialFvPatchScalarField::getNbr
(
    tmp<scalarField>& sigmaNbr,
    tmp<scalarField>& sigmaEPotNbr
) const
{
    const electromagneticModel& em =
        patch().boundaryMesh().mesh()
       .lookupType<electromagneticModel>();

    sigmaNbr = em.sigma(patch().index());
    sigmaEPotNbr = sigmaNbr()*patchInternalField();
}

void Foam::coupledElectricPotentialFvPatchScalarField::add
(
    tmp<scalarField>& result,
    const tmp<scalarField>& field
) const
{
    if (result.valid())
    {
        result.ref() += field;
    }
    else
    {
        if (field.isTmp())
        {
            result = field;
        }
        else
        {
            result = field().clone();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupledElectricPotentialFvPatchScalarField::
coupledElectricPotentialFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    ePotName_(iF.name())
{
   mappedPatchBase::validateMapForField
    (
        *this,
        iF,
        dict,
        mappedPatchBase::from::differentPatch
    );

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        // default: valueFraction=1; (fixed value)
        valueFraction() = 1;
        refValue() = *this;
        // gradient condition not used
        refGrad() = 0;
    }
    //evaluate() is called after all solvers have been constructed
}


Foam::coupledElectricPotentialFvPatchScalarField::
coupledElectricPotentialFvPatchScalarField
(
    const coupledElectricPotentialFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    ePotName_(psf.ePotName_)
{}


Foam::coupledElectricPotentialFvPatchScalarField::
coupledElectricPotentialFvPatchScalarField
(
    const coupledElectricPotentialFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    ePotName_(psf.ePotName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coupledElectricPotentialFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp = mappedPatchBase::getMap(patch().patch());
    const label patchiNbr = mpp.nbrPolyPatch().index();
    const fvPatch& patchNbr =
        refCast<const fvMesh>(mpp.nbrMesh()).boundary()[patchiNbr];

    const fvPatchScalarField& ePotpNbr =
        patchNbr.lookupPatchField<volScalarField, scalar>(ePotName_);

    if (!isA<coupledElectricPotentialFvPatchScalarField>(ePotpNbr))
    {
        FatalErrorInFunction
            << "Patch field for " << internalField().name() << " on "
            << this->patch().name() << " is of type "
            << coupledElectricPotentialFvPatchScalarField::typeName
            << endl << "The neighbouring patch field "
            << internalField().name() << " on "
            << patchNbr.name() << " is required to be the same, but is "
            << "currently of type " << ePotpNbr.type() << exit(FatalError);
    }

    const coupledElectricPotentialFvPatchScalarField& coupledPotentialNbr =
        refCast<const coupledElectricPotentialFvPatchScalarField>(ePotpNbr);

    // Get neighbour contributions
    tmp<scalarField> sigma;
    tmp<scalarField> sigmaEPot;
    {
        tmp<scalarField> sigmaNbr;
        tmp<scalarField> sigmaEPotNbr;
        coupledPotentialNbr.getNbr(sigmaNbr, sigmaEPotNbr);

        add(sigmaEPot, mpp.fromNeighbour(sigmaEPotNbr));
        add(sigma, mpp.fromNeighbour(sigmaNbr));
    }

    // if sigma->0: ePot = ePotNbr; if sigmaNbr->0 => ePot = 0
    this->refValue() = sigmaEPot()/(sigma()+SMALL);

    mixedFvPatchScalarField::updateCoeffs();

    // Restore tag
    UPstream::msgType() = oldTag;
}

const Foam::scalarField& Foam::coupledElectricPotentialFvPatchScalarField::getNbr() const
{
    const electromagneticModel& em =
        patch().boundaryMesh().mesh()
       .lookupType<electromagneticModel>();

    // J_normal = sigma * ( grad(ePot) * n )
    // Return patch normal current.
    return em.sigma(patch().index())*snGrad();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        coupledElectricPotentialFvPatchScalarField
    );
}


// ************************************************************************* //
