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

/*void Foam::coupledCurrentDensityFvPatchVectorField::getValues
(
    tmp<scalarField>& sigma
) const
{
    const electromagneticModel& em =
        patch().boundaryMesh().mesh()
       .lookupType<electromagneticModel>();

    sigma = em.sigma(patch().index());
}*/

void Foam::coupledCurrentDensityFvPatchVectorField::initCoupledPotential()
{
    //Allow for single initialization
    if (!coupled_)
    {
        word Jname_ = internalField().name();
        const word Jname("deltaJ");
        if ( Jname_ != Jname &&
            Jname_.size() >= Jname.size())
        {
            ePotName_ = "PotE" + Jname_.substr(Jname_.size() - Jname.size());
        }
        else
        {
            ePotName_ = "PotE";
        }
            coupled_ = true;
    }
}

void Foam::coupledCurrentDensityFvPatchVectorField::getJfromPotential
(
    const coupledElectricPotentialFvPatchScalarField& psf,
    tmp<scalarField>& sigmaGradPotE
) const
{
    const electromagneticModel& em =
        psf.patch().boundaryMesh().mesh()
       .lookupType<electromagneticModel>();

    tmp<scalarField> gradPotE;
    psf.getGradientValues(gradPotE);
    sigmaGradPotE = em.sigma(psf.patch().index())* gradPotE();
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

    const fvPatchScalarField& ePotPatch =
        patch().lookupPatchField<volScalarField, scalar>(ePotName_);

    if (!isA<coupledElectricPotentialFvPatchScalarField>(ePotPatch))
    {
        FatalErrorInFunction
            << "Patch field for " << internalField().name() << " on "
            << this->patch().name() << " is of type "
            << coupledCurrentDensityFvPatchVectorField::typeName
            << endl << "The patch field " << ePotName_
            << " on " << this->patch().name() << " is required to be "
            << coupledElectricPotentialFvPatchScalarField::typeName
            << ", but is currently of type " << ePotPatch.type()
            << exit(FatalError);
    }

    const coupledElectricPotentialFvPatchScalarField& coupledPotential =
        refCast<const coupledElectricPotentialFvPatchScalarField>(ePotPatch);

    // Get patch values
    //tmp<scalarField> sigma;
    //getValues(sigma);
    tmp<scalarField> normalJ;
    getJfromPotential(coupledPotential,normalJ);

    this->refValue() = normalJ() * patch().nf();
    //(patchInternalField() & patch().nf())* patch().nf();//assign nearest normal component

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
