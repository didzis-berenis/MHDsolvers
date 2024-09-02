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

#include "coupledCurrentDensityFvPatchVectorField.H"
#include "coupledElectricPotentialFvPatchScalarField.H"
#include "electromagneticModel.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"
#include "mappedPatchBase.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::coupledCurrentDensityFvPatchVectorField::getThis
(
    tmp<scalarField>& sigma
) const
{
    const electromagneticModel& em =
        patch().boundaryMesh().mesh()
       .lookupType<electromagneticModel>();

    sigma = em.sigma(patch().index());
}

void Foam::coupledCurrentDensityFvPatchVectorField::initCoupledPotential()
{
    //Allow for single initialization
    if (!coupled_)
    {
        const electromagneticModel& em =
            patch().boundaryMesh().mesh()
        .lookupType<electromagneticModel>();

        ePotName_ = em.getCoupledPotentialName(internalField().name());
        coupled_ = true;
    }
}

void Foam::coupledCurrentDensityFvPatchVectorField::add
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

Foam::coupledCurrentDensityFvPatchVectorField::
coupledCurrentDensityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF)
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
        refValue() = vectorField("refValue", dict, p.size());
        refGrad() = vectorField("refGradient", dict, p.size());
        valueFraction() = symmTensorField("valueFraction", dict, p.size());
    }
    else
    {
        //valueFraction() is symmTensorField:
        //"xx", "xy", "xz",
        //      "yy", "yz",
        //            "zz"
        valueFraction() = sqr(patch().nf());//fraction set to normal direction
        // Thus, refValue() is normalCurrent
        // Act as slip condition by default
        refValue() = Zero;
        //refGrad not used
        refGrad() = Zero;
    }
    //evaluate() is called after all solvers have been constructed
}


Foam::coupledCurrentDensityFvPatchVectorField::
coupledCurrentDensityFvPatchVectorField
(
    const coupledCurrentDensityFvPatchVectorField& psf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(psf, p, iF, mapper)
{}


Foam::coupledCurrentDensityFvPatchVectorField::
coupledCurrentDensityFvPatchVectorField
(
    const coupledCurrentDensityFvPatchVectorField& psf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(psf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coupledCurrentDensityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    if (!coupled_)
    {
        //Act as slip condition by default
        //Use default coefficients
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
            << coupledCurrentDensityFvPatchVectorField::typeName
            << endl << "The neighbouring patch field " << ePotName_
            << " on " << patchNbr.name() << " is required to be "
            << coupledElectricPotentialFvPatchScalarField::typeName
            << ", but is currently of type " << ePotpNbr.type()
            << exit(FatalError);
    }

    const coupledElectricPotentialFvPatchScalarField& coupledPotentialNbr =
        refCast<const coupledElectricPotentialFvPatchScalarField>(ePotpNbr);

    tmp<scalarField> sigma;
    // Get patch values
    getThis(sigma);
    tmp<scalarField> normalJ;
    // Add neighbour contributions
    {
        tmp<scalarField> normalJNbr;
        normalJNbr = coupledPotentialNbr.getNbr();
        add(normalJ, mpp.fromNeighbour(normalJNbr));
    }

    // default: J . n = - Jnbr . n = grad(ePotNbr)*sigmaNbr . n
    // if sigmaNbr->0 => Jnbr . n = 0
    // if sigma->0 => sigma()/(sigma() + SMALL) = 0 => J . n = 0
    this->refValue() = (sigma()/(sigma() + SMALL)) * ( - normalJ() ) * patch().nf();

    directionMixedFvPatchVectorField::updateCoeffs();

    // Restore tag
    UPstream::msgType() = oldTag;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        coupledCurrentDensityFvPatchVectorField
    );
}


// ************************************************************************* //
