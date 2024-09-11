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
#include "electromagneticModel.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"
#include "mappedPatchBase.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::coupledCurrentDensityFvPatchVectorField::getValues
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
    coupled_ = true;
}

Foam::word Foam::coupledCurrentDensityFvPatchVectorField::suffix() const
{
    const word Jname = internalField().name();
    if ( Jname != Jname_ &&
        Jname.size() >= Jname_.size())
    {
        return Jname.substr(Jname_.size());
    }
    return "";
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
    evaluate();// is called after all solvers have been constructed
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

    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;
    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    // Get patch values
    tmp<scalarField> sigma;
    getValues(sigma);

    // Get coupled electric potential patch field from this patch
    const fvPatchScalarField& ePotp =
        patch().lookupPatchField<volScalarField, scalar>(ePotName_ + suffix());

    tmp<scalarField> ePotPatch;
    tmp<scalarField> ePotInternal;
    ePotPatch = ePotp;//Returns patch value
    ePotInternal = ePotp.patchInternalField();//Returns nearest internal field value

    // Get velocity contribution field from this patch
    const fvPatchVectorField& deltaUxBp =
        patch().lookupPatchField<volVectorField, vector>(UxBname_+suffix());

    tmp<scalarField> deltaUxB;
    deltaUxB = deltaUxBp.patchInternalField() & patch().nf();

    this->refValue() =
        (
        -(ePotPatch() - ePotInternal()) * patch().deltaCoeffs()//Potential gradient
        +
        deltaUxB()
        ) * sigma() * patch().nf();

    //this->refValue() = patchInternalField();//assign nearest value

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
